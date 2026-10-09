#!/usr/bin/env python3
"""Extract vertical sounding profiles from numerical weather model data.

The closest model grid point to the requested location and the model run /
forecast step closest to the requested time are used; no interpolation is done
in time or in the horizontal.  Two CSV files are written:

- ``<YYYYMMDD>_<HH>Z_<lat>_<lon>_orig.csv``: the original model levels
- ``<YYYYMMDD>_<HH>Z_<lat>_<lon>_200m.csv``: profile interpolated to constant
  altitude steps (200 m by default)

Dependencies: numpy and eccodes (``pip install numpy eccodes``).

New data sources are added by subclassing :class:`DataSource` and registering
the class in :data:`SOURCES`.  New data formats are added by implementing a
reader similar to :class:`GribReader`.
"""

import argparse
import datetime as dt
import re
import sys
import time
import urllib.error
import urllib.request
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

R_DRY = 287.05  # Specific gas constant of dry air, J/(kg K)
HEADER = ("# Altitude (m), Pressure (hPa), Temperature (deg C), "
          "Air density (kg/m3), Wind direction (deg), Wind speed (m/s)")


@dataclass
class Profile:
    """Vertical profile at a single grid point, sorted by increasing altitude.

    The first element is the surface, the rest are the model levels above it.
    """

    altitude: np.ndarray  # m above mean sea level
    pressure: np.ndarray  # Pa
    temperature: np.ndarray  # K
    u: np.ndarray  # m/s, eastward wind component
    v: np.ndarray  # m/s, northward wind component
    lat: float  # Grid point latitude
    lon: float  # Grid point longitude
    init_time: dt.datetime
    valid_time: dt.datetime

    @property
    def density(self):
        return self.pressure / (R_DRY * self.temperature)

    def interpolate(self, step):
        """Interpolate the profile to altitudes that are multiples of *step*.

        The lowest altitude is the surface altitude rounded to the nearest
        multiple of *step*.  If it is below the surface, surface values are used.
        """
        start = round(self.altitude[0] / step) * step
        altitude = np.arange(start, self.altitude[-1] + 1e-6, step)
        # Pressure decreases exponentially with height, so interpolate its logarithm
        pressure = np.exp(np.interp(altitude, self.altitude, np.log(self.pressure)))
        return Profile(
            altitude=altitude,
            pressure=pressure,
            temperature=np.interp(altitude, self.altitude, self.temperature),
            u=np.interp(altitude, self.altitude, self.u),
            v=np.interp(altitude, self.altitude, self.v),
            lat=self.lat, lon=self.lon,
            init_time=self.init_time, valid_time=self.valid_time,
        )

    def above_surface(self):
        """Return the profile without the surface level."""
        return Profile(self.altitude[1:], self.pressure[1:], self.temperature[1:],
                       self.u[1:], self.v[1:], self.lat, self.lon,
                       self.init_time, self.valid_time)

    def write_csv(self, path):
        speed = np.hypot(self.u, self.v)
        # Meteorological convention: the direction the wind blows from
        direction = np.round(np.degrees(np.arctan2(-self.u, -self.v))) % 360
        with open(path, "w") as fid:
            fid.write(HEADER + "\n")
            for row in zip(self.altitude, self.pressure, self.temperature,
                           self.density, direction, speed):
                z, p, t, rho, wdir, wspd = row
                # Use more decimals in the upper stratosphere and mesosphere
                p_fmt = ".1f" if p >= 100 else ".4f"
                rho_fmt = ".6f" if rho >= 1e-3 else ".4e"
                # Adding 0.0 turns a negative zero into a positive one
                t_c = round(t - 273.15, 1) + 0.0
                fid.write(f"{z:.0f},{p / 100:{p_fmt}},{t_c:.1f},{rho:{rho_fmt}},"
                          f"{wdir:.0f},{wspd:.1f}\n")


# ---------------------------------------------------------------------------
# Data format readers
# ---------------------------------------------------------------------------

class GribReader:
    """Read point values from GRIB (edition 1 or 2) messages using eccodes."""

    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon
        self.index = None
        self.grid_lat = None
        self.grid_lon = None

    def read(self, message):
        """Return the value of the closest grid point from a single GRIB message."""
        import eccodes

        gid = eccodes.codes_new_from_message(message)
        try:
            if self.index is None:
                nearest = eccodes.codes_grib_find_nearest(gid, self.lat, self.lon)[0]
                self.index = nearest.index
                self.grid_lat = nearest.lat
                self.grid_lon = (nearest.lon + 180) % 360 - 180
            return eccodes.codes_get_double_element(gid, "values", self.index)
        finally:
            eccodes.codes_release(gid)


# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

class DataSource(ABC):
    """Base class for model data sources.

    Subclasses define the model run cycle, the available forecast steps and
    how the data for a single run and step are read.
    """

    #: Hours between consecutive model runs
    cycle_hours = 6

    @abstractmethod
    def forecast_steps(self):
        """Return the available forecast steps in hours, in increasing order."""

    @abstractmethod
    def is_available(self, init_time, step):
        """Check whether data exist for the given run and forecast step."""

    @abstractmethod
    def get_profile(self, init_time, step, lat, lon):
        """Read a :class:`Profile` for the closest grid point."""

    def select_run(self, valid_time, max_lead):
        """Find the run and forecast step closest to *valid_time*.

        The candidates are ordered by the difference between the requested and
        the available time, and then by the forecast length, so that data
        closest to the initialization time are preferred.  The first one with
        available data is returned.
        """
        steps = self.forecast_steps()
        cycle = dt.timedelta(hours=self.cycle_hours)
        # Allow a run initialized a bit after the requested time when its
        # analysis is closer than any forecast of the previous run
        half_step = dt.timedelta(hours=steps[1] - steps[0]) / 2
        init = _floor_time(valid_time + half_step, cycle)
        earliest = valid_time - dt.timedelta(hours=max_lead)

        candidates = []
        while init >= earliest:
            lead = (valid_time - init).total_seconds() / 3600
            step = min(steps, key=lambda s: abs(s - lead))
            # Round to minutes to avoid floating point noise in the ordering
            diff = round(abs(step - lead) * 60)
            if step <= max_lead and diff <= (steps[-1] - steps[-2]) * 30:
                candidates.append((diff, step, init))
            init -= cycle
        candidates.sort(key=lambda c: (c[0], c[1]))

        for _, step, init in candidates:
            if self.is_available(init, step):
                return init, step
        raise RuntimeError(f"No data found for {valid_time:%Y-%m-%d %H:%M} UTC "
                           f"with forecast lengths up to {max_lead} hours")


class GFSSource(DataSource):
    """NOAA GFS from the AWS Open Data bucket (data available from 2021 onwards)."""

    base_url = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
    cycle_hours = 6
    level_variables = ("HGT", "TMP", "UGRD", "VGRD")
    surface_variables = {
        "HGT": "surface",
        "PRES": "surface",
        "TMP": "2 m above ground",
        "UGRD": "10 m above ground",
        "VGRD": "10 m above ground",
    }

    def __init__(self, resolution="0p25", workers=8):
        self.resolution = resolution
        self.workers = workers
        self._url_cache = {}

    def forecast_steps(self):
        if self.resolution == "0p25":
            return list(range(0, 120, 1)) + list(range(120, 385, 3))
        return list(range(0, 385, 3))

    def _candidate_urls(self, init_time, step):
        day = f"gfs.{init_time:%Y%m%d}/{init_time:%H}"
        fname = f"gfs.t{init_time:%H}z.pgrb2.{self.resolution}.f{step:03d}"
        # The "atmos" sub-directory was added in March 2021 (GFS v16)
        return [f"{self.base_url}/{day}/atmos/{fname}", f"{self.base_url}/{day}/{fname}"]

    def _data_url(self, init_time, step):
        key = (init_time, step)
        if key not in self._url_cache:
            self._url_cache[key] = None
            for url in self._candidate_urls(init_time, step):
                if _url_exists(url + ".idx"):
                    self._url_cache[key] = url
                    break
        return self._url_cache[key]

    def is_available(self, init_time, step):
        return self._data_url(init_time, step) is not None

    def get_profile(self, init_time, step, lat, lon):
        url = self._data_url(init_time, step)
        index = _parse_idx(_http_get(url + ".idx").decode())
        messages = self._download_messages(url, index)

        reader = GribReader(lat, lon % 360)
        values = {key: reader.read(msg) for key, msg in messages.items()}

        surface = {var: values[(var, level)] for var, level in self.surface_variables.items()}
        levels = sorted({level for (var, level) in values if level.endswith(" mb")},
                        key=lambda level: -float(level.split()[0]))
        rows = [(surface["HGT"], surface["PRES"], surface["TMP"],
                 surface["UGRD"], surface["VGRD"])]
        for level in levels:
            pressure = float(level.split()[0]) * 100
            height = values[("HGT", level)]
            # Skip levels that are extrapolated below the ground
            if pressure >= surface["PRES"] or height <= surface["HGT"]:
                continue
            rows.append((height, pressure, values[("TMP", level)],
                         values[("UGRD", level)], values[("VGRD", level)]))
        data = np.array(rows).T
        return Profile(*data, lat=reader.grid_lat, lon=reader.grid_lon,
                       init_time=init_time,
                       valid_time=init_time + dt.timedelta(hours=step))

    def _wanted(self, var, level):
        if self.surface_variables.get(var) == level:
            return True
        return var in self.level_variables and re.fullmatch(r"[\d.]+ mb", level) is not None

    def _download_messages(self, url, index):
        """Download the needed GRIB messages using HTTP range requests.

        Adjacent messages are merged into a single request.
        """
        # Keep only the first occurrence of each variable/level combination
        first = {}
        for i, entry in enumerate(index):
            if self._wanted(entry["var"], entry["level"]):
                first.setdefault((entry["var"], entry["level"]), i)
        wanted = sorted(first.values())
        groups = []
        for i in wanted:
            if groups and groups[-1][-1] == i - 1:
                groups[-1].append(i)
            else:
                groups.append([i])

        def fetch(group):
            start = index[group[0]]["offset"]
            end = index[group[-1]]["end"]
            data = _http_get(url, (start, end))
            return {(index[i]["var"], index[i]["level"]):
                    data[index[i]["offset"] - start:
                         None if index[i]["end"] is None else index[i]["end"] - start + 1]
                    for i in group}

        total = len(wanted)
        print(f"Downloading {total} GRIB messages from {url}", file=sys.stderr)
        messages = {}
        with ThreadPoolExecutor(self.workers) as executor:
            for result in executor.map(fetch, groups):
                messages.update(result)
        return messages


#: Available data sources, selectable with the --source command-line option
SOURCES = {
    "gfs": GFSSource,
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_idx(text):
    """Parse a wgrib2 style .idx inventory file."""
    entries = []
    for line in text.strip().splitlines():
        parts = line.split(":")
        entries.append({"offset": int(parts[1]), "var": parts[3], "level": parts[4]})
    for this, following in zip(entries, entries[1:] + [None]):
        this["end"] = None if following is None else following["offset"] - 1
    return entries


def _http_get(url, byte_range=None, retries=3):
    request = urllib.request.Request(url)
    if byte_range is not None:
        start, end = byte_range
        request.add_header("Range", f"bytes={start}-{'' if end is None else end}")
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError):
            if attempt == retries - 1:
                raise
            time.sleep(2 ** attempt)


def _url_exists(url):
    request = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(request, timeout=30):
            return True
    except urllib.error.HTTPError as err:
        if err.code in (403, 404):
            return False
        raise


def _floor_time(time_, interval):
    epoch = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)
    return epoch + (time_ - epoch) // interval * interval


def _format_coordinate(value, positive, negative):
    text = f"{abs(value):.2f}".rstrip("0")
    if text.endswith("."):
        text += "0"
    return text + (positive if value >= 0 else negative)


def _parse_time(text):
    time_ = dt.datetime.fromisoformat(text)
    if time_.tzinfo is None:
        return time_.replace(tzinfo=dt.timezone.utc)
    return time_.astimezone(dt.timezone.utc)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Save sounding profiles from weather model data as CSV files.")
    parser.add_argument("lat", type=float, help="Latitude in degrees, north positive")
    parser.add_argument("lon", type=float, help="Longitude in degrees, east positive")
    parser.add_argument("time", type=_parse_time,
                        help="Time in ISO format, e.g. 2024-01-01T06:00 (UTC if no time zone given)")
    parser.add_argument("--source", choices=SOURCES, default="gfs",
                        help="Model data source (default: %(default)s)")
    parser.add_argument("--resolution", choices=("0p25", "0p50", "1p00"), default="0p25",
                        help="GFS grid resolution in degrees (default: %(default)s)")
    parser.add_argument("--step", type=float, default=200,
                        help="Altitude step of the interpolated profile in meters "
                             "(default: %(default)s)")
    parser.add_argument("--max-lead", type=int, default=384,
                        help="Maximum forecast length in hours to accept (default: %(default)s)")
    parser.add_argument("--outdir", type=Path, default=Path("."),
                        help="Output directory (default: current directory)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.source == "gfs":
        source = GFSSource(resolution=args.resolution)
    else:
        source = SOURCES[args.source]()

    try:
        init_time, step = source.select_run(args.time, args.max_lead)
    except RuntimeError as err:
        sys.exit(str(err))
    print(f"Using {args.source.upper()} run {init_time:%Y-%m-%d %H} UTC, "
          f"forecast step +{step} h", file=sys.stderr)
    profile = source.get_profile(init_time, step, args.lat, args.lon)

    prefix = (f"{profile.valid_time:%Y%m%d_%H}Z_"
              f"{_format_coordinate(profile.lat, 'N', 'S')}_"
              f"{_format_coordinate(profile.lon, 'E', 'W')}")
    step_label = f"{args.step:g}m"
    args.outdir.mkdir(parents=True, exist_ok=True)
    orig_path = args.outdir / f"{prefix}_orig.csv"
    interp_path = args.outdir / f"{prefix}_{step_label}.csv"
    profile.above_surface().write_csv(orig_path)
    profile.interpolate(args.step).write_csv(interp_path)
    print(f"Grid point {profile.lat:.2f}, {profile.lon:.2f}; "
          f"valid time {profile.valid_time:%Y-%m-%d %H:%M} UTC", file=sys.stderr)
    print(f"Saved {orig_path} and {interp_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
