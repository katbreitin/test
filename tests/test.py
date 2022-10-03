"""
test.py

Test cases can be added by creating a function with "test" in the name

A helper function "_run_it()" makes it easy to add a test.

@save decorator marks tests that we'll permanently save the output from (see save_record.py)

@extra decorator marks tests that won't be run unless we set TEST_EXTRA=True below)

"""

from clavrx import run_clavrx, get_config, build_options_file, build_file_list
from pathlib import Path
from tempfile import mkdtemp
from shutil import rmtree
import subprocess
import pytest
import signal
import os
import re

CLAVRX = Path('../run/bin/clavrxorb').absolute()
assert CLAVRX.exists(), str(CLAVRX)+' does not exist'

HERE = Path(__file__).absolute().parent

def get_all_l2_variables():
    l2_variables = set()
    pattern = re.compile(r'.*case\([\'"](\w+)[\'"]\).*')
    with open(HERE / '../main_src/level2_mod.f90') as fp:
        for line in fp:
            if 'case("' in line:
                match = pattern.match(line)
                if match is not None:
                    l2_variables.add(match.group(1))
    return sorted(l2_variables)


L2_LIST_CONTENT = '\n'.join(['testing',*get_all_l2_variables()])

TEST_EXTRA = False

DEFAULT_MODE=None

extra = pytest.mark.skipif(not TEST_EXTRA, reason="extra test")

save_funcs = []
def save(func):
    save_funcs.append(func)
    return func


def _run_it(main_l1b_file, aux_l1b_files=(), config_override=None, out_dir=None, mode=DEFAULT_MODE):
    """
    Run clavrx case

    main_l1b_file :: the primary l1b granule
    aux_l1b_files :: list of files that will get linked next to the primary file in a temp dir
    config_override :: dict of clavrx options to override the defaults (see options/default_options.yml)
    out_dir :: directory to persist output the l2 file (default is a temporary directory which is removed)
    
    returns a string containing the contents of the clavrx_options file
    """
    main_l1b_file = Path(main_l1b_file)
    tmpdir = Path(mkdtemp(dir=HERE))
    try:
        if out_dir is None:
            out_dir = tmpdir / 'out'
        else:
            out_dir = Path(out_dir)
        out_dir.mkdir(exist_ok=True, parents=True)
        l1b_links = []
        for f in aux_l1b_files:
            f = Path(f)
            assert f.exists()
            newlink = tmpdir / f.name
            newlink.symlink_to(f)
            l1b_links.append(newlink)

        main_l1b_link = tmpdir / main_l1b_file.name
        main_l1b_link.symlink_to(main_l1b_file)

        #print()
        #args = ['ls','-l',str(tmpdir)]
        #print(' '.join(args))
        #subprocess.run(args)
        #print()

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
        config['temp_dir'] = str(tmpdir) #str(out_dir / 'temp_dir')
        options_file_content = build_options_file(config)
        file_list_content = build_file_list(out_dir, [main_l1b_link])
        level2_list_content = L2_LIST_CONTENT
        options_file_path = tmpdir / 'clavrx_options'
        with open(options_file_path,'w') as fp:
            fp.write(options_file_content)
        level2_list_content = L2_LIST_CONTENT
        with open(tmpdir / 'level2_list', 'w') as fp:
            fp.write(level2_list_content)
        file_list_path = tmpdir / 'file_list'
        with open(file_list_path, 'w') as fp:
            fp.write(file_list_content)
        if mode=='perf':
            p = subprocess.run(['perf','record','-F99','-g',CLAVRX], cwd=tmpdir)
        elif mode=='gdb':
            p = subprocess.run(['gdb',CLAVRX], cwd=tmpdir)
        elif mode=='valgrind':
            p = subprocess.run(['valgrind','--suppressions='+str(Path('suppressions').absolute()), CLAVRX], cwd=tmpdir)
        else:
            p = subprocess.run([CLAVRX], cwd=tmpdir)
        assert p.returncode == 0
        return options_file_content
    finally:
        rmtree(tmpdir)

def test_pgroup_bug():
    main_l1b_file = Path('./dummy')
    tmpdir = Path(mkdtemp(dir=HERE))
    out_dir = tmpdir / 'out'
    out_dir.mkdir()
    try:
        config = get_config()
        config['temp_dir'] = str(out_dir / 'temp_dir')
        options_file_content = build_options_file(config)
        file_list_content = build_file_list(out_dir, [main_l1b_file])
        options_file_path = tmpdir / 'clavrx_options'
        with open(options_file_path,'w') as fp:
            fp.write(options_file_content)
        # Skip level2_list so the process ends early
        # This test is for a bug that causes fast failures
        #level2_list_content = L2_LIST_CONTENT
        #with open(tmpdir / 'level2_list', 'w') as fp:
            #fp.write(level2_list_content)
        file_list_path = tmpdir / 'file_list'
        with open(file_list_path, 'w') as fp:
            fp.write(file_list_content)
        # this returns immediately, and we need to keep the shell alive
        if os.fork() == 0:
            # child
            # ignore SIGINT
            signal.signal(2, signal.SIG_IGN)
            os.chdir(tmpdir)
            os.execv(str(CLAVRX),['clavrxorb'])
        else:
            pid,ws = os.wait()
            rc = ws >> 8
            assert rc != 125
    finally:
        rmtree(tmpdir)


@save
def test_noaa_viirs(out_dir=None):
    #ROOT = Path('/apollo/cloud/archive/Satellite_Input/VIIRS-N20/global/2018/098/')
    #VIIRS_L1_FNAME = 'GMTCO_j01_d20180407_t2358242_e2359487_b01995_c20190226192728879159_noac_ops.h5'
    ROOT = Path('/ships19/cloud/archive/clavrx_test_data/viirs/noaa')
    VIIRS_L1 = ROOT / 'GMTCO_npp_d20220221_t2356569_e0002373_b53481_c20220222005748538034_oebc_ops.h5'
    aux = set(ROOT.glob('*npp_d20220221_t2356569_e0002373*'))
    aux.remove(VIIRS_L1)
    return _run_it(VIIRS_L1, aux, out_dir=out_dir)

@save
def test_nasa_viirs(out_dir=None):
    ROOT = Path('/ships19/cloud/archive/clavrx_test_data/viirs/nasa')
    vnp03 = ROOT / 'VNP03MOD.A2019003.1700.002.2021102031552.nc'
    extra = set(ROOT.glob('*A2019003.1700*'))
    extra.remove(vnp03)
    return _run_it(vnp03, extra, out_dir=out_dir)

@extra
def test_nasa_viirs_rttov_slow(out_dir=None):
    ROOT = Path('/ships19/cloud/archive/clavrx_test_data/viirs/nasa')
    vnp03 = ROOT / 'VNP03MOD.A2019003.1700.002.2021102031552.nc'
    extra = set(ROOT.glob('*A2019003.1700*'))
    extra.remove(vnp03)
    config_override = {
        #'acha_alg':'off',
        #'dcomp_alg':0,
        #'prc_cloud_flag':0,
        'rtm':'rttov',
        'sfc_emiss':'rttov',
    }
    return _run_it(vnp03, extra, config_override=config_override, out_dir=out_dir)



@save
def test_avhrr(out_dir=None):
    AVHRR = Path('/arcdata/polar/noaa/noaa18/2020/2020_01_01_001/avhrr/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI')
    return _run_it(AVHRR, out_dir=out_dir)

def test_avhrr_get_goes_header_bug():
    # This file is empty
    AVHRR = Path('/arcdata/polar/noaa/noaa18/2022/2022_02_05_036/avhrr/NSS.GHRR.NN.D22036.S1847.E2013.B8614849.GC')
    try:
        return _run_it(AVHRR)
        assert False, 'This file should cause a nonzero returncode'
    except AssertionError as e:
        pass


@save
def test_fusion(out_dir=None):
    FUSION = Path('/ships19/cloud/archive/Satellite_Input/HIRS-FUSION/NN/2020/001/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI.fusion.nc')
    AVHRR = Path('/arcdata/polar/noaa/noaa18/2020/2020_01_01_001/avhrr/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI')
    override = {'lut':'ecm2_lut_avhrr123_hirs_common_chs.nc'}
    return _run_it(FUSION, [AVHRR], config_override=override, out_dir=out_dir)


@save
def test_g16_fd(out_dir=None):
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
    
    return _run_it(files[0], files[1:], config_override=override, out_dir=out_dir)

@save
def test_g16_conus(out_dir=None):
    ROOT = Path('/arcdata/goes/grb/goes16/2020/2020_05_02_123/abi/L1b/RadC/')
    files = [next(ROOT.glob(f'OR_ABI-L1b-RadC-M6C{i:02d}_G16_s2020123000111*.nc')) for i in range(1,17)]
    override = {'bounds':{
        'enable':True,
        'lat_north': 50.0,
        'lat_south': 40.0,
        'lon_west': -90.,
        'lon_east': -80.,
        'zen_min': 0.0,
        'zen_max': 90.0,
        'solzen_min': 0.0,
        'solzen_max': 180.0,
        'name': 'SUB'
        }}
    
    return _run_it(files[0], files[1:], config_override=override, out_dir=out_dir)


@extra
def test_process_cloud_off():
    # Use the avhrr template
    AVHRR = Path('/arcdata/polar/noaa/noaa18/2020/2020_01_01_001/avhrr/NSS.GHRR.NN.D20001.S0000.E0143.B7532324.WI')
    override = { 'prc_cloud_flag':0 }
    return _run_it(AVHRR, config_override=override)


@save
def test_g17_fd(out_dir=None):
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
    
    return _run_it(files[0], files[1:], config_override=override, out_dir=out_dir)


def _g17_conus(out_dir=None, sfc_emiss='seebor', rtm='pfast', see_emiss=True):
    ROOT = Path('/arcdata/goes/grb/goes17/2020/2020_05_02_123/abi/L1b/RadC/')
    files = [next(ROOT.glob(f'OR_ABI-L1b-RadC-M6C{i:02d}_G17_s2020123222117*.nc')) for i in range(1,17)]
    override = {'bounds':{
        'enable':True,
        'lat_north': 40.0,
        'lat_south': 30.0,
        'lon_west': -130.,
        'lon_east': -120.,
        'zen_min': 0.0,
        'zen_max': 90.0,
        'solzen_min': 0.0,
        'solzen_max': 180.0,
        'name': 'SUB'
        },
        'sfc_emiss':sfc_emiss,
        'see_emiss':see_emiss,
        'rtm':rtm
    }
    
    return _run_it(files[0], files[1:], config_override=override, out_dir=out_dir)

for rtm in ['pfast','rttov']:
    for sfc_emiss in ['umd','rttov','seebor']:
        for see_emiss in [False, True]:
            if see_emiss:
                see_emiss_str = 'IRseaemis'
            else:
                see_emiss_str = 'IRseaemis'
            exec(f"""
@save
def test_g17_conus_{rtm}_{sfc_emiss}_{see_emiss_str}(out_dir=None):
    return _g17_conus(out_dir=out_dir, rtm='{rtm}', sfc_emiss='{sfc_emiss}', see_emiss={see_emiss})
            """)

@save
def test_h8_fd_dat(out_dir=None):
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
    
    return _run_it(files[0], files[1:], config_override=override, out_dir=out_dir)


