import netCDF4
import xarray as xr
from tqdm import tqdm
from pathlib import Path
import warnings
import sys

WMO_ID_MAPPING = {
    152:'goes-16',
    270:'goes-16',
    664:'goes-17',
    271:'goes-17',
    167:'himawari-8',
    173:'himawari-8',
    684:'meteosat-8',
    55:'meteosat-8',
    305:'meteosat-11',
    70:'meteosat-11'
}

def setup_dataset(orig, fname):
    nc = netCDF4.Dataset(fname, mode='w')
    nc.set_auto_maskandscale(False)
    # Define dimensions
    for d in orig.dimensions:
        nc.createDimension(d, orig.dimensions[d].size)
    nc.createDimension('layer', 3)
    # Define Variables
    for k,v in orig.variables.items():
        if hasattr(v, '_FillValue'):
            fill_value=v._FillValue
        else:
            fill_value = None
        dims = v.dimensions
        if dims == ('scan_lines_along_track_direction', 'pixel_elements_along_scan_direction'):
            dims = ('layer', *dims)
        nc.createVariable(k, v.datatype.str, dimensions=dims, zlib=True, fill_value=fill_value, complevel=1)
        for k2 in v.ncattrs():
            if k2 == '_FillValue':
                continue
            nc.variables[k].setncattr(k2, v.getncattr(k2))
    # Copy global attributes
    for k in orig.ncattrs():
        nc.setncattr(k, orig.getncattr(k))
    return nc


def prepare_inputs(l2_files):

    netcdfs = {}
    for f in l2_files:
        wmo_id = int(f.name.split('.')[0].split('_')[-1])
        name = WMO_ID_MAPPING[wmo_id]
        netcdfs[name] = netCDF4.Dataset(f)
        netcdfs[name].set_auto_maskandscale(False)
    orig = netcdfs[name]
    return netcdfs, orig


def make_index(wmo_file, sats):
    wmo_id_ds = xr.open_dataset(wmo_file)
    idx = {}
    l1g_mapping = {}

    for i in wmo_id_ds.wmo_id.attrs['satellite_names'].split(';'):
        id, name = i.split('=')
        l1g_mapping[name.lower()] = int(id)
        
    with tqdm(sats) as bar:
        for sat in bar:
            mask = (wmo_id_ds.wmo_id[0] == l1g_mapping[sat]).values
            idx[sat] = mask.nonzero()
    wmo_id_ds.close()
    return idx


def shuffle_data(nc, netcdfs, orig, idx):
    nc.set_auto_maskandscale(False)
    with tqdm(orig.variables.items()) as bar:
        # for each variable
        for k,v in bar:
            if v.dimensions == ('scan_lines_along_track_direction', 'pixel_elements_along_scan_direction'):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    data = nc.variables[k][:]
                # for each satellite
                for sat,l2 in netcdfs.items():
                    x,y,z = idx[sat]
                    bar.set_description(f'{k:30}{sat:20} read nc')
                    v = l2.variables[k][:]
                    bar.set_description(f'{k:30}{sat:20} shuffle')
                    data[x,y,z] = v[y,z]
                bar.set_description(f'{k:30}{"":20} write nc')
                nc.variables[k][:] = data
            else:
                nc.variables[k] = orig.variables[k][:]


def main(out_file, wmo_file, l2_files):
    # Open input files
    netcdfs, orig = prepare_inputs(l2_files)
    # Open output file and setup variables
    nc = setup_dataset(orig, out_file)
    # pre-compute shuffle
    print('Making index')
    idx = make_index(wmo_file, netcdfs.keys())
    # Do the shuffle
    print('Processing')
    shuffle_data(nc, netcdfs, orig, idx)
    nc.close()


def validate_args(args):
    out_file = Path(args.out_file)
    if out_file.exists():
        print(out_file,'already exists')
        sys.exit(1)
    wmo_file = Path(args.wmo_file)
    if not wmo_file.exists():
        print(wmo_file, 'does not exist')
        sys.exit(1)
    l2_files = []
    for f in args.l2_files:
        f = Path(f)
        if not f.exists():
            print(f, 'does not exist')
            sys.exit(1)
        l2_files.append(f)
    return out_file, wmo_file, l2_files
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('out_file', help='path to output netcdf')
    parser.add_argument('wmo_file', help='path to l1g wmo_id file')
    parser.add_argument('l2_files', nargs='+', help='paths to individual l2 files from each satellite')
    args = parser.parse_args()

    out_file, wmo_file, l2_files = validate_args(args)
    
    main(out_file, wmo_file, l2_files)

