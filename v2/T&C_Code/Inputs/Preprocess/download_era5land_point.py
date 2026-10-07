"""Download a multi-year ERA5-Land hourly point series and write it in the
shape Build_TC_Forcing_From_ERA5_Land_Timeseries expects.

Uses reanalysis-era5-land-timeseries, the ARCO point product, rather than the
gridded archive. That matters for two reasons:

  * It accepts a multi-decade date range in a single request, so this is about
    seven requests rather than the 2,520 that chunking the gridded product
    down to monthly field-count limits would require.
  * Its accumulated variables are served already de-accumulated, as hourly
    increments, which is what the MATLAB builder expects. The gridded product
    accumulates from 00 UTC and resets at midnight.

The request serialisation for location and date is not documented in a form
the client exposes, so the first request probes the plausible shapes and the
working one is reused for the rest.

CDS serves UTC and accepts whole dates only, while everything downstream works
in local standard time: Apply_CMIP6_Deltas_To_ERA5_Land groups the record by
year(utc + DeltaGMT) and rejects any partial year. A request framed as whole
UTC years therefore produces two broken local years -- a few hours belonging to
the year before the start, and the same number missing from the final year. So
the request brackets the target window in whole UTC days, and assemble() trims
it to the hour and checks every local year is complete before writing.

Requires: cdsapi, xarray, netcdf4, dask, numpy. A CDS account with ~/.cdsapirc
set up, and the dataset licence accepted on the CDS website.

Usage:  python -u download_era5land_point.py
"""

import os

import numpy as np
import xarray as xr

# ---- configuration ------------------------------------------------------
LAT, LON = 40.0625, -88.1963      # Energy Farm

DELTA_GMT = -6                 # must equal cfg.DeltaGMT in Delta_Forcing_Config.m
FIRST_LOCAL_YEAR, LAST_LOCAL_YEAR = 1995, 2024

_OFFSET = np.timedelta64(DELTA_GMT, "h")
UTC_FIRST = np.datetime64(f"{FIRST_LOCAL_YEAR}-01-01T00", "h") - _OFFSET
UTC_LAST = np.datetime64(f"{LAST_LOCAL_YEAR}-12-31T23", "h") - _OFFSET

START = str(UTC_FIRST.astype("datetime64[D]"))
END = str(UTC_LAST.astype("datetime64[D]"))

# The window is in the directory name, so a changed window can never be served
# silently from a previous download by the skip-if-present guard in download().
RAW_DIR = "era5land_raw_%s_%s" % (START.replace("-", ""), END.replace("-", ""))
OUT_FILE = f"energyfarm_{FIRST_LOCAL_YEAR}_{LAST_LOCAL_YEAR}.nc"
DATASET = "reanalysis-era5-land-timeseries"

VARIABLES = [
    "2m_temperature",
    "2m_dewpoint_temperature",
    "surface_pressure",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "total_precipitation",
    "surface_solar_radiation_downwards",
]

SHORT = {
    "2m_temperature": "t2m",
    "2m_dewpoint_temperature": "d2m",
    "surface_pressure": "sp",
    "10m_u_component_of_wind": "u10",
    "10m_v_component_of_wind": "v10",
    "total_precipitation": "tp",
    "surface_solar_radiation_downwards": "ssrd",
}


def request_shapes():
    """Plausible serialisations of location and date, most likely first."""
    locations = [
        {"latitude": LAT, "longitude": LON},
        {"n": LAT, "s": LAT, "e": LON, "w": LON},
    ]
    dates = [
        [f"{START}/{END}"],
        {"start": START, "end": END},
        [START, END],
    ]
    for loc in locations:
        for date in dates:
            yield loc, date


def unzip_in_place(path):
    """CDS may ignore download_format and return a zip regardless, so detect
    one by its magic bytes and replace it with the NetCDF it contains."""
    with open(path, "rb") as fh:
        if fh.read(2) != b"PK":
            return
    import shutil
    import tempfile
    import zipfile

    print("    returned a zip archive; extracting")
    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(path) as zf:
            members = [n for n in zf.namelist() if n.endswith(".nc")]
            if len(members) != 1:
                raise SystemExit(
                    f"Expected one NetCDF inside {path}, found: {zf.namelist()}"
                )
            zf.extract(members[0], tmp)
            extracted = os.path.join(tmp, members[0])
        os.remove(path)
        shutil.move(extracted, path)


def download():
    import cdsapi

    os.makedirs(RAW_DIR, exist_ok=True)
    client = cdsapi.Client()
    working = None
    for variable in VARIABLES:
        path = os.path.join(RAW_DIR, f"era5land_{SHORT[variable]}.nc")
        if os.path.exists(path) and os.path.getsize(path) > 0:
            print(f"  {variable}: already present, skipping")
            continue
        shapes = [working] if working else list(request_shapes())
        last_error = None
        for loc, date in shapes:
            request = {
                "variable": [variable],
                "location": loc,
                "date": date,
                "data_format": "netcdf",
                "download_format": "unarchived",
            }
            try:
                print(f"  {variable}: requesting {START} to {END}")
                client.retrieve(DATASET, request).download(path)
                unzip_in_place(path)
                working = (loc, date)
                break
            except Exception as err:
                if os.path.exists(path) and os.path.getsize(path) == 0:
                    os.remove(path)
                last_error = err
                if working:
                    raise
                print(f"    request shape rejected: {str(err)[:120]}")
        else:
            raise SystemExit(
                f"Every request shape was rejected for {variable}.\n"
                f"Last error: {last_error}\n"
                "Open the dataset's Download tab, click 'Show API request', "
                "and send me the generated request dict."
            )


def is_cumulative(values, hours):
    """Cumulative fields never decrease within a UTC day; increments do."""
    steps = np.diff(values)
    within_day = hours[1:] != 0
    inside = steps[within_day]
    inside = inside[np.isfinite(inside)]
    if inside.size == 0:
        return False
    return np.mean(inside < -1e-9) < 0.01


def deaccumulate(values, hours):
    out = np.diff(values, prepend=np.nan)
    start_of_day = hours == 0
    out[start_of_day] = values[start_of_day]
    out[0] = values[0]
    return np.maximum(out, 0.0)


def assemble():
    files = sorted(
        os.path.join(RAW_DIR, f)
        for f in os.listdir(RAW_DIR)
        if f.endswith(".nc")
    )
    if not files:
        raise SystemExit("No downloaded files found; run the download step first.")

    ds = xr.open_mfdataset(files, combine="by_coords").load()
    time_name = "valid_time" if "valid_time" in ds.dims else "time"

    lat_name = next((n for n in ("latitude", "lat") if n in ds.coords), None)
    lon_name = next((n for n in ("longitude", "lon") if n in ds.coords), None)
    if lat_name is None or lon_name is None:
        raise SystemExit(f"No coordinates found; variables present: {list(ds.coords)}")

    point = ds.squeeze()
    cell_lat = float(np.atleast_1d(point[lat_name].values)[0])
    cell_lon = float(np.atleast_1d(point[lon_name].values)[0])
    print(f"\nNearest ERA5-Land cell: {cell_lat:.4f} N, {cell_lon:.4f} E")
    print("Put these exact values into Delta_Forcing_Config.m")

    times = point[time_name].values
    hours = times.astype("datetime64[h]").astype(int) % 24

    gaps = np.unique(np.diff(times).astype("timedelta64[h]").astype(int))
    if not np.array_equal(gaps, [1]):
        raise SystemExit(f"Time axis is not strictly hourly; spacings found: {gaps}")

    out = xr.Dataset()
    out[time_name] = ((time_name,), times)
    for name in point.data_vars:
        values = np.asarray(point[name].values, dtype="float64").squeeze()
        if values.ndim != 1:
            continue
        if name in ("tp", "ssrd"):
            # This product should already be de-accumulated; check rather
            # than assume, in case the product convention ever changes.
            if is_cumulative(values, hours):
                print(f"  {name}: cumulative, de-accumulating")
                values = deaccumulate(values, hours)
            else:
                print(f"  {name}: hourly increments, as expected")
        out[name] = ((time_name,), values)

    units = {
        "t2m": "K", "d2m": "K", "sp": "Pa",
        "u10": "m s-1", "v10": "m s-1",
        "tp": "m", "ssrd": "J m-2",
    }
    missing = [v for v in units if v not in out]
    if missing:
        raise SystemExit(f"Download is missing required variables: {missing}")
    for name, unit in units.items():
        out[name].attrs["units"] = unit

    # ---- trim to whole local years --------------------------------------
    # Done after de-accumulation, so that logic still sees the UTC-day-aligned
    # record it expects. The 8760/8784 test below is the one MATLAB applies in
    # Apply_CMIP6_Deltas_To_ERA5_Land; running it here means a bad window fails
    # in seconds instead of surviving to the forcing build.
    keep = (times >= UTC_FIRST) & (times <= UTC_LAST)
    if not keep.any():
        raise SystemExit("Downloaded window %s..%s does not cover %s..%s"
                         % (times[0], times[-1], UTC_FIRST, UTC_LAST))
    out = out.isel({time_name: np.flatnonzero(keep)})
    times = times[keep]

    local_year = (times + _OFFSET).astype("datetime64[Y]").astype(int) + 1970
    for y in range(FIRST_LOCAL_YEAR, LAST_LOCAL_YEAR + 1):
        n = int((local_year == y).sum())
        want = 8784 if (y % 4 == 0 and y % 100 != 0) or y % 400 == 0 else 8760
        if n != want:
            raise SystemExit(
                "Local year %d has %d hours, expected %d; the requested UTC "
                "window does not cover whole local years." % (y, n, want)
            )
    print("  trimmed to %d hours = whole local years %d-%d"
          % (len(times), FIRST_LOCAL_YEAR, LAST_LOCAL_YEAR))

    out.attrs.update(
        source=DATASET,
        cds_request_utc="%s/%s" % (START, END),
        trimmed_window_utc="%s/%s" % (UTC_FIRST, UTC_LAST),
        local_years="%d-%d" % (FIRST_LOCAL_YEAR, LAST_LOCAL_YEAR),
        delta_gmt_hours=DELTA_GMT,
        requested_point="%.4f, %.4f" % (LAT, LON),
    )

    out["latitude"] = cell_lat
    out["latitude"].attrs["units"] = "degrees_north"
    out["longitude"] = cell_lon
    out["longitude"].attrs["units"] = "degrees_east"
    out[time_name].encoding.update(
        units="seconds since 1970-01-01 00:00:00", dtype="float64"
    )
    out.to_netcdf(OUT_FILE)
    print(f"\nWrote {OUT_FILE}: {len(times)} hourly records, "
          f"{str(times[0])[:10]} to {str(times[-1])[:10]}")


if __name__ == "__main__":
    print(f"Downloading {DATASET}, one request per variable:")
    download()
    print("\nAssembling point series:")
    assemble()
