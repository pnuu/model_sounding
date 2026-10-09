# model_sounding.py

Download weather model data and save the vertical sounding profile at a
location as CSV files.

The script uses the model grid point closest to the given location, and the
model run and forecast step closest to the given time. Nothing is interpolated
in time or in the horizontal. If several model runs cover the requested time
equally well, the most recent one (shortest forecast length) is used. If the
newest run is not yet available, older runs are tried.

Currently supported source: **GFS** from the
[NOAA GFS AWS Open Data bucket](https://registry.opendata.aws/noaa-gfs-bdp-pds/).
The bucket holds data from 2021 onwards. Only the needed GRIB messages are
downloaded, using HTTP range requests.

## Installation

Requires Python 3.11 or newer, plus numpy and eccodes:

```sh
pip install numpy eccodes
```

## Usage

```sh
python model_sounding.py LAT LON TIME [options]
```

- `LAT`, `LON`: location in degrees. North and east are positive.
- `TIME`: ISO format time, such as `2024-01-01T06:00`. A time without a time
  zone is treated as UTC.

Options:

| Option | Default | Description |
|---|---|---|
| `--source` | `gfs` | Model data source |
| `--resolution` | `0p25` | GFS grid resolution: `0p25`, `0p50` or `1p00` degrees |
| `--step` | `200` | Altitude step of the interpolated profile, in meters |
| `--max-lead` | `384` | Longest forecast length to accept, in hours |
| `--outdir` | `.` | Output directory |

Example:

```sh
python model_sounding.py 40.5 0.0 2024-01-01T06:00
```

This writes these two files:

- `20240101_06Z_40.5N_0.0E_orig.csv` has the original model (pressure) levels
  above the ground.
- `20240101_06Z_40.5N_0.0E_200m.csv` has the profile interpolated to 200 m
  altitude steps.

The date, hour and coordinates in the file names are the valid time and the
grid point that were used. These can differ from the requested ones. The
script prints the selected run, forecast step and grid point.

With the 0.25 degree grid, one profile is about 100 MB of downloads. Use
`--resolution 1p00` for a much smaller and faster download. The 0.5 and 1.0
degree grids have 3-hourly forecast steps; the 0.25 degree grid is hourly
up to +120 h.

## Output format

```
# Altitude (m), Pressure (hPa), Temperature (deg C), Air density (kg/m3), Wind direction (deg), Wind speed (m/s)
1024,900.0,5.1,1.126923,314,5.1
```

- Altitudes are geopotential heights above mean sea level.
- Air density is computed from pressure and temperature, as for dry air.
- Wind direction is the direction the wind blows from (0 = north, 90 = east).
- Above 1 hPa, pressure and density get more decimals. GFS reaches about
  0.01 hPa (around 80 km).

The interpolated profile starts at the model surface altitude rounded to the
nearest step. Model surface values (surface pressure, 2 m temperature, 10 m
wind) form its lowest level. Temperature and wind components are interpolated
linearly in altitude, and pressure logarithmically.

## Adding data sources and formats

- **New source:** subclass `DataSource`. Implement `forecast_steps()`,
  `is_available()` and `get_profile()`, set `cycle_hours`, and add the class
  to the `SOURCES` dictionary. The base class handles selecting the run and
  forecast step.
- **New data format:** write a reader like `GribReader`, which returns the
  closest grid point value of a field. Use it from the source's
  `get_profile()`.
