from pathlib import Path

HERE = Path(__file__)
ROOT = HERE.parent.parent

OUTPUT_FILE = ROOT / "src/main/default_options_file.f90"
INPUT_FILE = ROOT / 'run/clavrx_options'

with open(INPUT_FILE) as fp:
    lines = list(fp)

def string_chunks(chunksize=100):
    for line in lines:
        for i in range(0, len(line), chunksize):
            yield line[i:i+chunksize]

with open(OUTPUT_FILE,'w') as fp:
    fp.write('module default_options_file\n')
    fp.write('implicit none\ncontains\n')
    fp.write('subroutine print_default_options()\n')
    for chunk in string_chunks():
        if chunk.endswith('\n'):
            fp.write('WRITE(*,"(A)") "')
        else:
            fp.write('WRITE(*,"(A)",advance="no") "')
        fp.write(chunk.replace('\n',''))
        fp.write('"\n')
    fp.write("end subroutine print_default_options\n")
    fp.write("end module default_options_file\n")
print(OUTPUT_FILE)

