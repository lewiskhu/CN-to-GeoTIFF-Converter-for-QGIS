#!/usr/bin/env python3
"""
tropomi_to_geotiff.py

Lightweight TROPOMI to GeoTIFF converter with safe dry-run and deferred imports.

Reads local gridded TROPOMI files (NetCDF or GeoTIFF) in a directory, computes
monthly means, annual mean, and seasonal metrics (amplitude & phase using
a simple harmonic fit), and writes GeoTIFF outputs ready for QGIS.

Default input dir:
    D:\Lewis\Global phenology maping\output_phen\sif_to_tiff

Default output dir:
    D:\Lewis\Global phenology maping\output_phen\ImageCollection

This file defers heavy GIS/scientific imports so you can run --dry-run to verify
file discovery even if xarray/rioxarray/rasterio are not installed in the current
Python environment. Full processing requires the packages in requirements.txt.
"""

import argparse
import os
from glob import glob
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Process gridded TROPOMI files to GeoTIFFs for QGIS")
    p.add_argument("--input-dir", required=False,
                   default=r"D:\\Lewis\\Global phenology maping\\output_phen\\sif_to_tiff",
                   help="Directory containing TROPOMI files (.nc, .tif). Default: D:\\Lewis\\Global phenology maping\\output_phen\\sif_to_tiff")
    p.add_argument("--out-dir", required=False,
                   default=r"D:\\Lewis\\Global phenology maping\\output_phen\\ImageCollection",
                   help="Directory to write output GeoTIFFs (default: D:\\Lewis\\Global phenology maping\\output_phen\\ImageCollection)")
    p.add_argument("--variable", default=None, help="Variable name to extract (for NetCDF). If not set, try to autodetect")
    p.add_argument("--time-coord", default=None, help="Time coordinate name in NetCDF (default: time or datetime)")
    p.add_argument("--skip-harmonic", action='store_true', help="Skip harmonic seasonal fit (amplitude/phase)")
    p.add_argument("--dry-run", action='store_true', help="Only list discovered input files and exit (no heavy imports)")
    return p.parse_args()


def find_files(input_dir):
    patterns = ["*.nc", "*.nc4", "*.tif", "*.tiff"]
    files = []
    for pat in patterns:
        files.extend(glob(os.path.join(input_dir, pat)))
    files.sort()
    return files


def safe_import(name, package=None):
    """Try to import a module and return it or print a helpful message and exit."""
    try:
        module = __import__(name) if package is None else __import__(package, fromlist=[name])
        return module
    except Exception as e:
        print(f"ERROR: required package '{name}' is not available: {e}")
        return None


def load_and_process(files, args):
    """Performs the heavy processing. Imports xarray/rioxarray and friends here so
    a --dry-run can work without them.
    """
    # Defer heavy imports
    xr = safe_import('xarray')
    np = safe_import('numpy')
    pd = safe_import('pandas')
    rioxarray = safe_import('rioxarray')
    scipy_opt = safe_import('scipy')
    if xr is None or np is None or pd is None or rioxarray is None or scipy_opt is None:
        print("One or more heavy dependencies are missing. Install requirements.txt (prefer conda-forge) and retry.")
        return 1

    from scipy.optimize import curve_fit

    # Now proceed with roughly the same logic as before, but with defensive checks
    datasets = []
    for f in files:
        if f.lower().endswith(('.nc', '.nc4')):
            ds = xr.open_dataset(f)
            var = args.variable if args.variable else list(ds.data_vars)[0]
            da = ds[var].squeeze()
            datasets.append(da)
        else:
            da = rioxarray.open_rasterio(f)
            if 'band' in da.dims and da.sizes['band'] > 1:
                da = da.isel(band=0)
            datasets.append(da)

    # time parsing
    times = []
    for i, da in enumerate(datasets):
        fname = os.path.basename(files[i])
        ts = None
        m = re.search(r'(20\d{2})[-_]?([01]?\d)(?:[-_]?([0-3]?\d))?', fname)
        if m:
            year = int(m.group(1)); month = int(m.group(2) or 1); day = int(m.group(3) or 1)
            try:
                ts = pd.Timestamp(year=year, month=month, day=day)
            except Exception:
                ts = None
        if ts is None:
            ts = pd.Timestamp(os.path.getmtime(files[i]), unit='s')
        times.append(ts)

    # align and concat
    base = datasets[0]
    aligned = []
    for da in datasets:
        try:
            da2 = da.rio.reproject_match(base)
        except Exception:
            try:
                da2 = da.interp_like(base)
            except Exception:
                da2 = da
        aligned.append(da2)

    stacked = xr.concat(aligned, dim='time')
    stacked = stacked.assign_coords(time=('time', times))

    # annual mean
    am = stacked.mean(dim='time')
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    def write_geotiff(arr2d, ref, out_path):
        arr = xr.DataArray(arr2d[np.newaxis, :, :], dims=('band','y','x'))
        try:
            arr = arr.assign_coords({'y': ref['y'], 'x': ref['x']})
            arr.rio.set_spatial_dims(x_dim='x', y_dim='y', inplace=True)
            if hasattr(ref.rio, 'crs') and ref.rio.crs is not None:
                arr.rio.write_crs(ref.rio.crs, inplace=True)
        except Exception:
            pass
        arr.rio.to_raster(out_path, driver='GTiff')

    am2 = am.values
    if am2.ndim == 3:
        am2 = am2[0,:,:]
    write_geotiff(am2, stacked.isel(time=0), os.path.join(out_dir, 'annual_mean.tif'))

    # monthly means
    monthly = stacked.groupby('time.month').mean(dim='time')
    for m in range(1,13):
        da_m = monthly.sel(month=m)
        arr = da_m.values
        if arr.ndim == 3:
            arr = arr[0,:,:]
        write_geotiff(arr, stacked.isel(time=0), os.path.join(out_dir, f'month_{m:02d}_mean.tif'))

    # harmonic metrics (optional)
    if not args.skip_harmonic:
        # very simple harmonic fit per-pixel
        months = np.arange(1,13)
        mdata = monthly.values
        if mdata.ndim == 4:
            mdata = mdata[:,0,:,:]
        ny, nx = mdata.shape[1], mdata.shape[2]
        m_flat = mdata.reshape(12, -1)
        means = np.nanmean(m_flat, axis=0)
        amps = np.full(means.shape, np.nan)
        phases = np.full(means.shape, np.nan)
        def fit_pixel(col):
            try:
                popt, _ = curve_fit(lambda t,a0,a1,b1: a0 + a1*np.cos(2*np.pi*t/12.) + b1*np.sin(2*np.pi*t/12.), months, col, p0=[np.nanmean(col),0,0], maxfev=2000)
                a0,a1,b1 = popt
                amp = np.sqrt(a1*a1 + b1*b1)
                phase = (np.arctan2(-b1,a1)/(2*np.pi)*12.0) % 12.0
                return amp, phase
            except Exception:
                return np.nan, np.nan
        for i in range(m_flat.shape[1]):
            col = m_flat[:,i]
            if np.all(np.isnan(col)):
                continue
            a,p = fit_pixel(col)
            amps[i] = a; phases[i] = p
        amp2 = amps.reshape(ny,nx)
        phase2 = phases.reshape(ny,nx)
        write_geotiff(amp2, stacked.isel(time=0), os.path.join(out_dir, 'harmonic_amplitude.tif'))
        write_geotiff(phase2, stacked.isel(time=0), os.path.join(out_dir, 'harmonic_phase_month.tif'))

    print('Processing complete. Outputs written to:', out_dir)
    return 0


def main():
    args = parse_args()
    files = find_files(args.input_dir)
    if args.dry_run:
        print('Dry-run: discovered files:')
        for f in files:
            print(' -', f)
        print(f'Total files: {len(files)}')
        return 0

    if len(files) == 0:
        print('No input files found in', args.input_dir)
        return 1

    # call heavy processing
    return load_and_process(files, args)


if __name__ == '__main__':
    rc = main()
    sys.exit(rc)