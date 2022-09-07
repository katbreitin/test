from subprocess import run
from tempfile import TemporaryDirectory
import os
from pathlib import Path
import yaml
from parse_options import write_options
import io
from copy import deepcopy
import shutil
from uuid import uuid4
import re


print_orig = print

def print(*args, flush=True, **kwargs):
    return print_orig(*args, flush=True, **kwargs)

ROOT = Path(__file__).resolve().absolute().parent

CLAVRX_COMMIT = None

if False:
    DEFAULT_CLAVRX = (ROOT / 'clavrx_bin/clavrxorb.no_phase_f8').absolute()
else:
    DEFAULT_CLAVRX = (ROOT / 'clavrx_bin/clavrxorb.acha_tau_uncer').absolute()

DEFAULT_L2_LIST = ROOT / 'level2_lists/level2_list'


def run_clavrx(clavrx, file_list_content, options_file_content, level2_list_content, debug=False, log=True, scratch=None):
    print(f'Using {clavrx}')
    clavrx = Path(clavrx).absolute()
    assert clavrx.exists(), str(clavrx)
    clavrx = str(clavrx)
    if scratch is None:
        if Path('/scratch').is_dir():
            scratch = Path('/scratch/jobs') / os.environ['USER']
            scratch.mkdir(exist_ok=True, parents=True)
    with TemporaryDirectory(prefix='clavrx_tmp_', dir=scratch) as tmpdir:
        tmpdir = Path(tmpdir)
        print(f'Logs temporarily in {tmpdir}')
        def crash():
            print('Copying logs to ./ ')
            shutil.copytree(tmpdir, tmpdir.name)
        try:
            file_list = tmpdir / 'file_list'
            with open(file_list, 'w') as fp:
                fp.write(file_list_content)
            options_file = tmpdir / 'clavrx_options'
            with open(options_file,'w') as fp:
                fp.write(options_file_content)
            level2_list = tmpdir / 'level2_list'
            with open(level2_list,'w') as fp:
                fp.write(level2_list_content)
            if log:
                with open(tmpdir / 'clavrx.log','wb') as stdout:
                    with open(tmpdir / 'clavrx.err','wb') as stderr:
                        tmpdir.chmod(0o555)
                        p = run([clavrx], cwd=tmpdir, stdout=stdout, stderr=stderr)
                        tmpdir.chmod(0o755)
            else:
                tmpdir.chmod(0o555)
                p = run([clavrx], cwd=tmpdir)
                tmpdir.chmod(0o755)
            p.check_returncode()
        except Exception:
            crash()
            raise
        if debug:
            crash()
        

def build_file_list(out, files):
    fp = io.StringIO()
    out = Path(out)
    files = [Path(f) for f in files]
    for part in zip(*[f.absolute().parts[:-1] for f in files]):
        for i in part:
            assert i == part[0], 'Not all files in same directory'
    dir = files[0].absolute().parent
    print(dir)
    fp.write(str(dir)+'/\n')
    fp.write(str(out.absolute())+'/\n')
    for f in files:
        fp.write(f.name+'\n')
    fp.seek(0)
    file_list_content = fp.read()
    return file_list_content


def build_options_file(config):
    fp = io.StringIO()
    write_options(fp, x)
    fp.seek(0)
    options_file_content = fp.read()
    return options_file_content


def build_options_file(config):
    fp = io.StringIO()
    write_options(fp, config)
    fp.seek(0)
    options_file_content = fp.read()
    return options_file_content


def default_level2_list_content(override_list=None):
    if override_list is None:
        l = DEFAULT_L2_LIST
    else:
        l = Path(override_list)
    l = l.absolute()
    return get_level2_list_content(l)

def get_level2_list_content(file):
    print(f'Level2 List {file}')
    with open(file) as fp:
        level2_list_content = fp.read()
    return level2_list_content


def test_run_clavrx():
    clavrx = DEFAULT_CLAVRX
    config = get_config()
    options_file_content = build_options_file(config) 

    files = [
        '/arcdata/polar/noaa/noaa16/2001/2001_04_01_091/avhrr/NSS.GHRR.NL.D01091.S0019.E0212.B0270203.WI',
        '/arcdata/polar/noaa/noaa16/2001/2001_04_01_091/avhrr/NSS.GHRR.NL.D01091.S0207.E0401.B0270304.WI'
    ]
    file_list_content = build_file_list('/tmp/', files)
    level2_list_content = default_level2_list_content()
    run_clavrx(clavrx, file_list_content, options_file_content, level2_list_content)


def update_config(default, override):
    config = deepcopy(default)
    for k,v in override.items():
        if isinstance(v, dict):
            for k2,v2 in v.items():
                if isinstance(v2,dict):
                    config[k][k2].update(v2)
                else:
                    config[k][k2] = v2
        else:
            config[k] = v
    return config
    

def get_config(override_yml=None):
    with open(ROOT / 'options/default_options.yml') as fp:
        default = yaml.load(fp, Loader=yaml.FullLoader)
    if override_yml is not None:
        with open(override_yml) as fp:
            override = yaml.load(fp, Loader=yaml.FullLoader)
        config = update_config(default, override)
    else:
        config = default
    return config


def add_processing_attrs(f, clavrx, options_file_content=None, config=None):

    print('Adding processing attributes')
    try:
        clavrxorb_name = clavrx.resolve().name
        clavrx_commit = None
        if '.' in clavrxorb_name:
            _,suffix = clavrxorb_name.split('.')
            if re.match('^(dirty_)?[a-f0-9]{40}$', suffix):
                clavrx_commit = suffix

        print('import')
        import netCDF4
        import hashlib
        import yaml
        print('open')
        nc = netCDF4.Dataset(f,'a')
        nc.setncattr('clavrx_bin', clavrxorb_name)
        print('clavrx md5')
        # Get binary md5sum
        with open(clavrx, 'rb') as fp:
            clavrxorb_md5 = hashlib.md5(fp.read()).hexdigest()
        print('set clavrxorb_md5')
        nc.setncattr('clavrxorb_md5', clavrxorb_md5)
        if clavrx_commit is not None:
            print('set clavrxorb_md5')
            nc.setncattr('clavrx_commit', clavrx_commit)
        else:
            print('no clavrx_commit')
        if config is not None:
            print('set config_yaml')
            nc.setncattr('config_yaml', yaml.dump(config))
        else:
            print('no config')
        if options_file_content is not None:
            print('set options_file_content')
            nc.setncattr('clavrx_options', options_file_content)
        else:
            print('no options_file_content')
        print('close')
        nc.close()
        print('done')
    except Exception as e:
        print('!!!!!! Could not add processing attrs')
        print(e)
        print('------ Could not add processing attrs')


def main(l1b_files, out_dir, override_yml=None, override_list=None, clavrx=None, debug=False, log=True):

    if clavrx is None:
        clavrx = DEFAULT_CLAVRX
    elif isinstance(clavrx, str):
        p1 = Path(clavrx)
        p2 = DEFAULT_CLAVRX.parent / clavrx
        if p1.exists():
            clavrx = p1
        elif p2.exists():
            # just the name
            clavrx = p2
        
    for f in l1b_files:
        f = Path(f)
        if not f.is_file():
            print(f'{f} is not a file')
            print(f'Exiting')
            return 1

    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    config = get_config(override_yml)
    print(config)
    options_file_content = build_options_file(config) 

    print(f'Processing {len(l1b_files)} files')
    uuid = str(uuid4())
    tmp_out_dir = out_dir / f'processing_{uuid}'
    tmp_out_dir.mkdir(exist_ok=True)
    file_list_content = build_file_list(tmp_out_dir, l1b_files)
    level2_list_content = default_level2_list_content(override_list)
    run_clavrx(clavrx, file_list_content, options_file_content, level2_list_content, debug=debug, log=log)
    for f in tmp_out_dir.glob('*'):
        if f.name.endswith('.nc'):
            add_processing_attrs(f, clavrx, options_file_content=options_file_content, config=config)
        print('Move')
        shutil.move(str(f), str(out_dir))
    tmp_out_dir.rmdir()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-c',help='config', required=False, default=None)
    parser.add_argument('-l',help='l2_list', required=False, default=None)
    parser.add_argument('--debug',action='store_true')
    parser.add_argument('--nolog',action='store_true')
    parser.add_argument('--clavrx',required=False, default=None)
    parser.add_argument('out_dir')
    parser.add_argument('l1b_files',nargs='+')
    args = parser.parse_args()
    main(args.l1b_files, args.out_dir, override_yml=args.c, override_list=args.l, debug=args.debug, log=not args.nolog, clavrx=args.clavrx)
