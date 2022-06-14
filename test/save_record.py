from pathlib import Path
import test

RECORD_DIR = Path('/ships19/cloud/archive/clavrx_test_data/version_granules')

def main(version):
    version_dir = RECORD_DIR / version
    for func in test.save_funcs:
        name = func.__name__
        output_dir = version_dir / name
        if output_dir.exists() and len(list(output_dir.glob('*')))>0:
            print(output_dir, 'already exists')
        else:
            output_dir.mkdir(exist_ok=True,parents=True)
            print(output_dir)
            func(out_dir=output_dir)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('version')
    args = parser.parse_args()
    main(args.version)

