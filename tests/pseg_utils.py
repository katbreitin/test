import os
import signal
import asyncio
import numpy as np
from tracer_utils import read_array, start_clavrx, peek, resume, kill, set_skip_processing, mem_filename, wait_until, get_or_done, get_scan_line_number
import tracer_utils
from pathlib import Path
import subprocess
import sys
assert sys.version_info[0] == 3, 'Python 3 required'
assert sys.version_info[1] >= 7, 'Python >= 3.7 required'


async def handle_segment(pid, file_lock):
    resume(pid)
    wn = await wait_until(pid, 3)
    if wn == -1:
        print(pid, 'UNEXPECTEDLY TERMINATED BEFORE PROCESSING COMPLETED')
        return
    async with file_lock:
        resume(pid)
        wn = await wait_until(pid, 6)
    if wn == -1:
        print(pid, 'UNEXPECTEDLY TERMINATED BEFORE WRITE COMPLETED')
        return
    print('writing done, killing', pid)
    await kill(pid)


def increment_hdf_start(pid, num_scans):
    mem_f = mem_filename(pid)
    with open(mem_f,'rb+') as fp:
        fp.seek(tracer_utils.NUM_SCANS_LEVEL2_HDF_ADDR)
        current = np.fromfile(fp, dtype=np.int32, count=1).item()
        new = current + num_scans
        print('Setting min scan of',pid, 'from', current, 'to', new)
        fp.seek(tracer_utils.NUM_SCANS_LEVEL2_HDF_ADDR)
        fp.write(np.int32(new))

        
async def spawner(parent_pid, processing_done, worker_queue):
    created_file = False
    try:
        while True:
            wn = await wait_until(parent_pid, None)
            if wn == 7:
                processing_done.set()
                await kill(parent_pid)
                return
            elif wn == 1:
                sibling_pids_addr = peek(parent_pid, tracer_utils.SIBLING_PIDS_ADDR_ADDR)
                pids = read_array(parent_pid, sibling_pids_addr, 1, np.int32)
                wn = await wait_until(pids[0], 1)
                assert wn != -1, 'UNEXPECTEDLY TERMINATED BEFORE PROCESSING STARTED'

                if not created_file:
                    await kill(pids[0])
                else:
                    num_scans = len(get_scan_line_number(pids[0]))
                    increment_hdf_start(parent_pid, num_scans)
                    await worker_queue.put(pids[0])
                    # Skip processing on parent
                    set_skip_processing(parent_pid)
                created_file = True
            elif wn == -1:
                print('PARENT TERMINATED UNEXPECTEDLY')
                processing_done.set()
                return
            # resume parent
            resume(parent_pid)
    except:
        os.kill(parent_pid, signal.SIGTERM)
        processing_done.set()
        raise


async def main(exe='./clavrxorb'):
    file_lock = asyncio.Lock()
    processing_done = asyncio.Event()
    worker_pids = asyncio.Queue(1)

    # start parent
    parent_pid = start_clavrx('parent', num_clones=1, skip_output=False, redirect=False, exe=exe)
    segment_tasks = []
    spawner_task = asyncio.create_task(spawner(parent_pid, processing_done, worker_pids))

    while True:
        result = await get_or_done(worker_pids, processing_done)
        if result is None:
            await spawner_task
            break
        else:
            pid = result
        segment_tasks.append(asyncio.create_task(handle_segment(pid, file_lock)))

    await asyncio.gather(*segment_tasks)
        

