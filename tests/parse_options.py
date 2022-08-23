import yaml


OLD_CONFIG = False
channel_map = None

VIIRS_CHANNELS = {f'M{i}':c for i,c in enumerate([8,9,3,4,1,15,2,5,26,6,7,20,22,29,31,32],1)}
VIIRS_CHANNELS.update({f'I{i}':i+38 for i in range(1,6)})
VIIRS_CHANNELS['DNB'] = 44

HIRS_CHANNELS = {'HIRS_4': 36, 'HIRS_5': 35, 'HIRS_6': 34, 'HIRS_7': 33, 'HIRS_9': 30, 'HIRS_11': 28, 'HIRS_12': 27, 'HIRS_14': 25, 'HIRS_15': 24, 'HIRS_18': 23}

AVHRR_CHANNELS = {
    'AVHRR_065':1,
    'AVHRR_086':2,
    'AVHRR_160':6,
}

FUSION_CHANNELS = {**HIRS_CHANNELS, **AVHRR_CHANNELS}
FUSION_CHANNELS['AVHRRHIRS_1100'] = 31
FUSION_CHANNELS['AVHRRHIRS_1200'] = 32
FUSION_CHANNELS['AVHRRHIRS_0375'] = 20

AHI_CHANNELS = {
'AHI_1':3, 'AHI_2':4, 'AHI_3':1, 'AHI_4':2, 'AHI_5':6, 'AHI_6':7,
'AHI_7':20, 'AHI_8':37, 'AHI_9':27, 'AHI_10':28, 'AHI_11':29, 'AHI_12':30,
'AHI_13':38, 'AHI_14':31, 'AHI_15':32, 'AHI_16':33
}


NEW_CHANNEL_MAP = {
    'viirs':VIIRS_CHANNELS,
    'avhrr_hirs':FUSION_CHANNELS,
    'ahi':AHI_CHANNELS
}


VERBOSITY = {
    'QUIET':0,
    'ERROR':1,
    'MINIMAL':2,
    'WARNING':4,
    'DEFAULT':5,
    'VERBOSE':9
}
CLOUD_MASK = {
    'Baseline':0,
    'ECM1':1,
    'ECM2':2
}

CCL_MODE = {
    'off':0,
    'top only':1,
    'top + base':2,
    'top + base + lower':3
}

CCL_TYPE = {
    'NOAT':0,
    'ISSCP':1,
    'NCEP':2
}

OUTPUT_FORMAT = {
    'hdf4':0,
    'netcdf4':1
}

NWP_MODEL = {
     'off':0,
     'gfs':1,
     'ncep':2,
     'cfsr':3,
     'gdas':4,
     'merra':5,
     'era':6,
     'gfs ait':7,
     'gfs fv3':8
}

NWP_MODE={
    'minimal':0,
    'all':1
}

RTM = {
    'crtm':0,
    'pfast':1,
    'rttov':2
}

NAV_OPTION = {
    'l1b':0,
    'external':1,
    'Reposnx':2
}

SNOW_MASK = {
    'no':0,
    'ims':1,
    'GlobSnow':2
}

REF_CAL = {
    'default':0,
    'recalibrate':1
}
THERM_CAL = {
    'default':0,
    'recalibrate':1
}

SFC_EMISS = {
    'umd':0,
    'rttov':1,
    'seebor':2,
    'crtm':3
}

MAPPING = {
    'verbosity':VERBOSITY,
    'cloud_mask':CLOUD_MASK,
    'ccl_mode':CCL_MODE,
    'ccl_type':CCL_TYPE,
    'output_format':OUTPUT_FORMAT,
    'nwp_model':NWP_MODEL,
    'nwp_mode':NWP_MODE,
    'rtm':RTM,
    'nav_option':NAV_OPTION,
    'snow_mask':SNOW_MASK,
    'therm_cal':THERM_CAL,
    'ref_cal':REF_CAL,
    'sfc_emiss':SFC_EMISS
}


def write_bounds(fp, options):
    enable = int(options['enable'])
    lat_north = options['lat_north']
    lat_south = options['lat_south']
    lon_west = options['lon_west']
    lon_east = options['lon_east']
    zen_min = options['zen_min']
    zen_max = options['zen_max']
    solzen_min = options['solzen_min']
    solzen_max = options['solzen_max']
    name = options['name']
    fp.write(f'{enable} {lat_south} {lat_north} {lon_west} {lon_east} {zen_min} {zen_max} {solzen_min} {solzen_max} {name}\n')

SAMPLE_MODES = {
    'sample':0,
    'average':1,
    'average+stats':2
}
def write_sample(fp, options):
    mode = SAMPLE_MODES[options['mode']]
    x_stride = options['x_stride']
    y_stride = options['y_stride']
    fp.write(f'{mode} {x_stride} {y_stride}\n')

def write_channels(fp, options):
    for r in [range(1,7), range(7,13), range(13,19), range(19,25), range(25,31), range(31,37)]:
        fp.write(' '.join([str(int(options['MODIS'][i])) for i in r])+'\n')
    fp.write(' '.join([str(int(i)) for i in [options['AxI'][37], options['AxI'][38],
                                             *[options['VIIRS'][f'I{j}'] for j in range(1,5)]
                                            ]
                      ])+'\n')
    fp.write(' '.join([str(int(i)) for i in [options['VIIRS']['I5'], options['VIIRS']['DNB'], *[0,0,0,0]]])+'\n')

def write_new_channels(fp, options):
    for r in [range(1,7), range(7,13), range(13,19), range(19,25), range(25,31), range(31,37), range(37,43)]:
        fp.write(' '.join([str(int(options[i])) for i in r])+'\n')
    fp.write(' '.join([str(int(i)) for i in [options[43], options[44], *[0,0,0,0]]])+'\n')

def write_options(fp, options):
    def convert(k, v):
        desc = ''
        if k in MAPPING:
            desc = ' ='+str(v)
            v = MAPPING[k][v]
        if isinstance(v, bool):
            return f'{str(int(v)):20} ! {k}{desc}\n'
        elif isinstance(v, int):
            return f'{str(v):20} ! {k}{desc}\n'
        else:
            print(f'No converter: {k}={v}')
            return f'{str(v):20} ! {k}{desc}\n'
    order1 = [
     'verbosity',
     'cloud_mask',
     'dcomp_alg',
     'acha_alg',
     'ccl_mode',
     'ccl_type',
     'enable_asos',
     'enable_nlcomp',
     'enable_aerosol',
     'enable_output',
     'output_format',
     'prc_cloud_flag',
     'scan_lines',
     'enable_sasrab',
     'nwp_model',
     'nwp_mode',
     'rtm',
     'nav_option',
     'compress_output',
     'aux_cloud_mask'
    ]
    order2 =  ['sfc_emiss',
     'see_emiss',
     'hires_sfc_type',
     'land_mask',
     'coast_mask',
     'elevation',
     'volcano',
     'snow_mask',
     'dark_composite',
    ]
    if OLD_CONFIG:
        order2[0] = 'seebor_emiss'
    order3 = ['lrc',
     'smooth_nwp',
     'process_undetected',
    ]
    fp.write(options['ancil_dir'].rstrip('/')+'/\n')
    fp.write(options['temp_dir'].rstrip('/')+'/\n')
    fp.write('9  ! Expert mode\n')
    for k in order1:
        fp.write(convert(k, options[k]))
    fp.write(options['lut']+'\n')
    for k in order2:
        fp.write(convert(k, options[k]))
    fp.write(convert('ref_cal', options['avhrr']['ref_cal']))
    fp.write(convert('therm_cal', options['avhrr']['therm_cal']))
    for k in order3:
        fp.write(convert(k, options[k]))
    write_bounds(fp, options['bounds'])
    write_sample(fp, options['sample'])
    if 'new_channels' in options:
        channels = new_channels(options['new_channels'])
        write_new_channels(fp, channels)
    else:
        write_channels(fp, options['channels'])


def new_channels(channels):
    # only look at this sensor
    new = {i:False for i in range(1,45)}
    sensor = channels['sensor']
    for c,b in channels[sensor].items():
        new[NEW_CHANNEL_MAP[sensor][c]] = b
    return new


def main():
    with open('clavrx_options.yml') as fp:
        x = yaml.load(fp, Loader=yaml.FullLoader)
    with open('clavrx_options', 'w') as fp:
        write_options(fp, x)

if __name__ == '__main__':
    main()


