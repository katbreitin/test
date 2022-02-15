from clavrx import run_clavrx, get_config, build_options_file, build_file_list
from pathlib import Path
from tempfile import mkdtemp
from shutil import rmtree
import subprocess
import pytest

CLAVRX = Path('../clavrx_bin/clavrxorb').absolute()
assert CLAVRX.exists(), str(CLAVRX)+' does not exist'

HERE = Path(__file__).absolute().parent


L2_LIST_CONTENT = "testing\nlatitude\nlongitude\ncloud_probability\n"

TEST_EXTRA = False

extra = pytest.mark.skipif(not TEST_EXTRA, reason="extra test")

def _run_it(main_l1b_file, aux_l1b_files=(), config_override=None):
    main_l1b_file = Path(main_l1b_file)
    tmpdir = Path(mkdtemp(dir=HERE))
    try:
        out_dir = tmpdir
        l1b_links = []
        for f in aux_l1b_files:
            f = Path(f)
            assert f.exists()
            newlink = tmpdir / f.name
            newlink.symlink_to(f)
            l1b_links.append(newlink)

        main_l1b_link = tmpdir / main_l1b_file.name
        main_l1b_link.symlink_to(main_l1b_file)

        print()
        args = ['ls','-l',str(tmpdir)]
        print(' '.join(args))
        subprocess.run(args)
        print()

        config = get_config()
        if config_override is not None:
            for k,v in config_override.items():
                if isinstance(v,dict):
                    for k2,v2 in v.items():
                        if isinstance(v2,dict):
                            config[k][k2].update(v2)
                        else:
                            config[k][k2] = v2
                else:
                    config[k] = v
        options_file_content = build_options_file(config)
        file_list_content = build_file_list(out_dir, [main_l1b_link])
        level2_list_content = L2_LIST_CONTENT
        run_clavrx(CLAVRX, file_list_content, options_file_content, level2_list_content, debug=False, log=False)
    finally:
        rmtree(tmpdir)


def test_viirs():
    ROOT = Path('/apollo/cloud/archive/Satellite_Input/VIIRS-N20/global/2018/098/')
    VIIRS_L1_FNAME = 'GMTCO_j01_d20180407_t2358242_e2359487_b01995_c20190226192728879159_noac_ops.h5'
    VIIRS_L1 = ROOT / VIIRS_L1_FNAME
    _run_it(VIIRS_L1)


def test_avhrr():
    AVHRR = Path('/arcdata/polar/noaa/noaa18/2020/2020_01_01_001/avhrr/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI')
    _run_it(AVHRR)


def test_fusion():
    FUSION = Path('/ships19/cloud/archive/Satellite_Input/HIRS-FUSION/NN/2020/001/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI.fusion.nc')
    AVHRR = Path('/arcdata/polar/noaa/noaa18/2020/2020_01_01_001/avhrr/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI')
    override = {'lut':'ecm2_lut_avhrr123_hirs_common_chs.nc'}
    _run_it(FUSION, [AVHRR], config_override=override)


def test_g16_fd():
    ROOT = Path('/arcdata/goes/grb/goes16/2020/2020_05_02_123/abi/L1b/RadF/')
    files = [next(ROOT.glob(f'OR_ABI-L1b-RadF-M6C{i:02d}_G16_s202012322201*.nc')) for i in range(1,17)]
    override = {'bounds':{
        'enable':True,
        'lat_north': 20.0,
        'lat_south': 10.0,
        'lon_west': -90.,
        'lon_east': -80.,
        'zen_min': 0.0,
        'zen_max': 90.0,
        'solzen_min': 0.0,
        'solzen_max': 180.0,
        'name': 'SUB'
        }}
    
    _run_it(files[0], files[1:], config_override=override)

def test_g16_conus():
    ROOT = Path('/arcdata/goes/grb/goes16/2020/2020_05_02_123/abi/L1b/RadC/')
    files = [next(ROOT.glob(f'OR_ABI-L1b-RadC-M6C{i:02d}_G16_s2020123000111*.nc')) for i in range(1,17)]
    override = {'bounds':{
        'enable':True,
        'lat_north': 20.0,
        'lat_south': 10.0,
        'lon_west': -90.,
        'lon_east': -80.,
        'zen_min': 0.0,
        'zen_max': 90.0,
        'solzen_min': 0.0,
        'solzen_max': 180.0,
        'name': 'SUB'
        }}
    
    _run_it(files[0], files[1:], config_override=override)


@extra
def test_process_cloud_off():
    # Use the avhrr template
    AVHRR = Path('/arcdata/polar/noaa/noaa18/2020/2020_01_01_001/avhrr/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI')
    override = { 'prc_cloud_flag':0 }
    _run_it(AVHRR, config_override=override)


def test_g17_fd():
    ROOT = Path('/arcdata/goes/grb/goes17/2020/2020_05_02_123/abi/L1b/RadF/')
    files = [next(ROOT.glob(f'OR_ABI-L1b-RadF-M6C{i:02d}_G17_s2020123222031*.nc')) for i in range(1,17)]
    override = {'bounds':{
        'enable':True,
        'lat_north': 20.0,
        'lat_south': 10.0,
        'lon_west': -130.,
        'lon_east': -120.,
        'zen_min': 0.0,
        'zen_max': 90.0,
        'solzen_min': 0.0,
        'solzen_max': 180.0,
        'name': 'SUB'
        }}
    
    _run_it(files[0], files[1:], config_override=override)


def test_g17_conus():
    ROOT = Path('/arcdata/goes/grb/goes17/2020/2020_05_02_123/abi/L1b/RadC/')
    files = [next(ROOT.glob(f'OR_ABI-L1b-RadC-M6C{i:02d}_G17_s2020123222117*.nc')) for i in range(1,17)]
    override = {'bounds':{
        'enable':True,
        'lat_north': 20.0,
        'lat_south': 10.0,
        'lon_west': -130.,
        'lon_east': -120.,
        'zen_min': 0.0,
        'zen_max': 90.0,
        'solzen_min': 0.0,
        'solzen_max': 180.0,
        'name': 'SUB'
        }}
    
    _run_it(files[0], files[1:], config_override=override)


def test_h8_fd_dat():
    ROOT = Path('/arcdata/nongoes/japan/himawari08/2021_01/2021_01_01_001/2300/')
    files = sorted(ROOT.glob('HS_H08_20210101_2300_B*_FLDK*'))
    override = {'bounds':{
        'enable':True,
        'lat_north': 20.0,
        'lat_south': 10.0,
        'lon_west': 130.,
        'lon_east': 140.,
        'zen_min': 0.0,
        'zen_max': 90.0,
        'solzen_min': 0.0,
        'solzen_max': 180.0,
        'name': 'SUB'
        }}
    
    _run_it(files[0], files[1:], config_override=override)


