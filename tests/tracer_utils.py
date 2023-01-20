import os
os.environ['OMP_NUM_THREADS'] = '1'
import signal
import numpy as np
import subprocess
from pathlib import Path
import asyncio
import tempfile
import shutil
from contextlib import contextmanager
from collections import defaultdict, namedtuple

WP_L1B  = 1
WP_BEFORE_ECM2 = 5
WP_ECM2 = 2
WP_TYPE = 4
WP_FINAL = 3
WP_SAVE = 6
WP_ALL_FINISH = 7


SYMBOL_KEYS = {}
for t in ['i1','i2','i4','i8','f4']:
    for d in range(4):
        k = f'__tracer_MOD_symbol_names_{t}_{d}d_ptr'
        SYMBOL_KEYS[k] = {'type':t,'dims':d,'array':'names'}
        k = f'__tracer_MOD_symbol_shapes_{t}_{d}d_ptr'
        SYMBOL_KEYS[k] = {'type':t,'dims':d,'array':'shapes'}
        k = f'__tracer_MOD_symbol_ptrs_{t}_{d}d_ptr'
        SYMBOL_KEYS[k] = {'type':t,'dims':d,'array':'ptrs'}
        k = f'__tracer_MOD_num_symbols_{t}_{d}d'
        SYMBOL_KEYS[k] = {'type':t,'dims':d,'array':'len'}


def parse_elf(exe):
    p = subprocess.run(['nm','-g',str(exe)], capture_output=True)

    for l in p.stdout.decode().splitlines():
        try:
            addr,_,name = l.split()
            if name in SYMBOL_KEYS:
                SYMBOL_KEYS[name]['ptr'] = int(addr, 16)
            elif name == '__tracer_MOD_rtm_shape_3d_addr':
                global RTM_SHAPE_3D_ADDR_ADDR
                RTM_SHAPE_3D_ADDR_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_num_symbols':
                global NUM_SYMBOLS_ADDR
                NUM_SYMBOLS_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_wait_number':
                global WAIT_NUMBER_ADDR
                WAIT_NUMBER_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_signal_int':
                global SIGNAL_INT_ADDR
                SIGNAL_INT_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_child_pid':
                global CHILD_PID_ADDR
                CHILD_PID_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_num_rtm_symbols':
                global NUM_RTM_SYMBOLS_ADDR
                NUM_RTM_SYMBOLS_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_rtm_numel':
                global RTM_NUMEL_ADDR
                RTM_NUMEL_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_rtm_symbol_table_addr':
                global RTM_SYMBOL_TABLE_ADDR_ADDR
                RTM_SYMBOL_TABLE_ADDR_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_rtm_symbol_names_addr':
                global RTM_SYMBOL_NAMES_ADDR_ADDR
                RTM_SYMBOL_NAMES_ADDR_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_rtm_channel_allocated_addr':
                global RTM_CHAN_ALLOC_ADDR_ADDR
                RTM_CHAN_ALLOC_ADDR_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_number_of_lines':
                global NUM_LINES_ADDR
                NUM_LINES_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_skip_processing':
                global SKIP_PROCESSING_ADDR
                SKIP_PROCESSING_ADDR = int(addr, 16)
            elif name == '__tracer_MOD_sibling_pids_addr':
                global SIBLING_PIDS_ADDR_ADDR
                SIBLING_PIDS_ADDR_ADDR = int(addr, 16)
            elif name == '__pixel_common_mod_MOD_num_scans_level2_hdf':
                global NUM_SCANS_LEVEL2_HDF_ADDR
                NUM_SCANS_LEVEL2_HDF_ADDR  = int(addr, 16)
        finally:
            pass

    global TABLE_POINTERS
    TABLE_POINTERS = defaultdict(dict)
    for d in SYMBOL_KEYS.values():
        if 'ptr' in d:
            TABLE_POINTERS[d['type'],d['dims']][d['array']] = d['ptr']
    TABLE_POINTERS = dict(TABLE_POINTERS.items())


def read_all_tables(pid, base_addr=None):
    if base_addr is None:
        base_addr = get_base_addr(pid)
    arrays = {}
    ArrayPtr = namedtuple('ArrayPtr',('ptr','shape','type'))
    for (t,d),table_ptrs in TABLE_POINTERS.items():
        for name,(array_ptr,shape) in read_table(pid, table_ptrs, d, base_addr=base_addr).items():
            arrays[name] = ArrayPtr(array_ptr, shape, t)
    return arrays


def read_table(pid, pointers, dims, base_addr=None):
    if base_addr is None:
        base_addr = get_base_addr(pid)
    len_ptr_offset = pointers['len']
    names_ptr_ptr_offset = pointers['names']
    num_symbols_ptr = len_ptr_offset+base_addr
    names_ptr = peek(pid, names_ptr_ptr_offset+base_addr)
    ptr_ptr = peek(pid, pointers['ptrs']+base_addr)
    symbols = read_symbols(pid, num_symbols_ptr, ptr_ptr, names_ptr)
    if ('shapes' in pointers) and (dims > 0):
        shapes_ptr = peek(pid, base_addr+pointers['shapes'])
        shapes = read_array(pid, shapes_ptr, (len(symbols),dims), np.int64)
        for k,shape in zip(symbols,shapes):
            symbols[k] = symbols[k], shape
    else:
        for k in symbols:
            symbols[k] = symbols[k], (1,)
    return symbols
    
    

def mem_filename(pid):
    mem_f = f'/proc/{pid}/mem'
    assert os.path.isfile(mem_f)
    return mem_f

def map_filename(pid):
    # unused for now
    map_f = f'/proc/{pid}/maps'
    assert os.path.isfile(map_f)
    return map_f

def get_base_addr(pid):
    # If clavrxorb is static exe
    return 0
    map_f = map_filename(pid)
    with open(map_f) as fp:
        for line in fp:
            start_addr = int(line.split('-')[0], 16)
            return start_addr

def read_array(pid, ptr, shape, dtype):
    size = np.prod(shape)
    mem_f = mem_filename(pid)
    with open(mem_f,'rb') as fp:
        fp.seek(ptr)
        tmp = np.fromfile(fp, dtype=dtype, count=size)
    #tmp.resize(shape, refcheck=False)
    tmp = tmp.reshape(shape)
    return tmp

def peek(pid, ptr, dtype=np.int64):
    mem_f = mem_filename(pid)
    with open(mem_f,'rb') as fp:
        fp.seek(ptr)
        return np.fromfile(fp, dtype=dtype, count=1).item()

def read_symbols(pid, num_symbols_addr, table_addr, names_addr):
    num_symbols = peek(pid, num_symbols_addr)
    symbols = {}
    for i in range(num_symbols):
        ptr = peek(pid, table_addr + i*8)
        name_words = [peek(pid, names_addr + i*256 + j*8) for j in range(32)]
        name = np.array(name_words, dtype=np.int64).tobytes()
        name = name.decode().strip()
        symbols[name] = ptr
    return symbols

def read_arrayptr(pid, arrayptr):
    if arrayptr.type == 'i4':
        return read_array(pid, arrayptr.ptr, arrayptr.shape, np.int32)
    elif arrayptr.type == 'i8':
        return read_array(pid, arrayptr.ptr, arrayptr.shape, np.int64)
    elif arrayptr.type == 'i1':
        return read_array(pid, arrayptr.ptr, arrayptr.shape, np.int8)
    elif arrayptr.type == 'f4':
        return read_array(pid, arrayptr.ptr, arrayptr.shape, np.float32)
    else:
        raise ValueError(f'Unknown type {arrayptr.type}')

def read_all_data(pid):
    arrayptrs = read_all_tables(pid)
    return read_data(pid, arrayptrs)

def read_subset_data(pid, names):
    arrayptrs = read_all_tables(pid)
    subset = {}
    for k in arrayptrs:
        if k in names:
            subset[k] = arrayptrs[k]
    return read_data(pid, subset)

def read_data(pid, arrayptrs):
    ds = {}

    for k,arrayptr in arrayptrs.items():
        v = read_arrayptr(pid, arrayptr)
        ds[k] = np.ma.masked_array(v)

    mask_all(ds)

    return ds


def mask_all(ds):
    for k,v in ds.items():
        if isinstance(v, np.ma.masked_array):
            if v.dtype in [np.float32, np.float64]:
                v[v==-999.] = np.ma.masked
                v[v==-9999.] = np.ma.masked
            elif v.dtype in [np.int8]:
                v[v==-128.] = np.ma.masked
            elif v.dtype in [np.int32, np.int16]:
                v[v==-999] = np.ma.masked
                v[v==-9999] = np.ma.masked


def get_scan_line_number(pid):
    arrays = read_all_tables(pid)
    v = read_arrayptr(pid, arrays['scan_line_number'])
    return np.ma.masked_equal(v, -999)

def start_clavrx(name=None, num_clones=0, redirect=True, skip_output=True, exe='./clavrxorb'):
    parse_elf(exe)
    pid = os.fork()
    if pid == 0:
        if name is None:
            name = str(os.getpid())
        # child
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        env = os.environ
        env['CLAVRX_ENABLE_TRACER'] = '1'
        if skip_output:
            env['CLAVRX_SKIP_OUTPUT'] = '1'
        env['CLAVRX_TRACER_CLONES'] = str(num_clones)
        if redirect:
            # Create empty log files first because it's easier than
            # setting raw flags with os.open
            stdout_f = f'clavrx_{name}.stdout'
            stderr_f = f'clavrx_{name}.stderr'
            with open(stdout_f,'w') as fp:
                pass
            with open(stderr_f,'w') as fp:
                pass
            # Redirect stderr and stdout to log files
            os.close(1)
            assert os.open(stdout_f, os.O_WRONLY) == 1
            os.set_inheritable(1, True)
            os.close(2)
            assert os.open(stderr_f, os.O_WRONLY) == 2
            os.set_inheritable(2, True)
        #os.execve(exe, (exe,'-acha_x','3120','-acha_y','1'), env)
        os.execve(exe, (exe,), env)
    else:
        return pid
            

async def wait_until(pid, waitnum):
    while True:
        p,stat = os.waitpid(pid, os.WUNTRACED|os.WNOHANG)
        if (p,stat) == (0,0):
            await asyncio.sleep(.1)
        elif os.WIFSTOPPED(stat):
            if os.WSTOPSIG(stat) == 19:
                # Sometimes there are signals that aren't the wait
                # Check the wait number to make sure we're at a wait point
                wn = peek(pid, WAIT_NUMBER_ADDR, np.uint64)
                if wn == waitnum or (waitnum is None and wn > 0):
                    return wn
            os.kill(pid, signal.SIGCONT)
        elif os.WIFEXITED(stat):
            rc = os.WEXITSTATUS(stat)
            if rc not in [0,1]:
                print('clavrx error: exit with code',rc)
                print('check clavrx log')
            return -1

def resume(pid):
    os.kill(pid, signal.SIGCONT)

async def kill(pid):
    os.kill(pid, signal.SIGTERM)
    os.kill(pid, signal.SIGCONT)
    await wait_until(pid, -1)

def set_skip_processing(pid):
    mem_f = mem_filename(pid)
    with open(mem_f,'rb+') as fp:
        fp.seek(SKIP_PROCESSING_ADDR)
        fp.write(np.int32(1))


@contextmanager
def tempswap(out):
    tmp = out.parent / ('tmp.'+out.name)
    if tmp.exists():
        tmp.unlink()
    try:
        yield tmp
    except:
        if tmp.exists():
            tmp.unlink()
        raise
    if out.exists():
        out.unlink()
    if tmp.exists():
        tmp.rename(out)

async def get_or_done(queue, done):
    while True:
        try:
            return queue.get_nowait()
        except asyncio.QueueEmpty:
            if done.is_set():
                return None
            else:
                await asyncio.sleep(.1)

    
@contextmanager
def tmpdir(debug=False, **kwargs):
    tmp = Path(tempfile.mkdtemp(**kwargs))
    if debug:
        print(tmp)
    try:
        yield tmp
    finally:
        if not debug:
            shutil.rmtree(tmp)
    

