import os
import signal
import asyncio
from asyncio import CancelledError
import numpy as np
from tracer_utils import read_array, start_clavrx, peek, resume, kill, set_skip_processing, mem_filename, wait_until, get_or_done, get_scan_line_number
import tracer_utils
from pathlib import Path
import subprocess

async def handle_segment(pid, file_lock):
    try:
        alive = True
        resume(pid)
        wn = await wait_until(pid, 3)
        if wn == -1:
            alive = False
            print(pid, 'UNEXPECTEDLY TERMINATED BEFORE PROCESSING COMPLETED')
            raise ValueError(f'{pid} terminated before processing completed')
        async with file_lock:
            resume(pid)
            wn = await wait_until(pid, 6)
        if wn == -1:
            alive = False
            print(pid, 'UNEXPECTEDLY TERMINATED BEFORE WRITE COMPLETED')
            raise ValueError(f'{pid} terminated before output write completed')
        print('writing done')
    finally:
        if alive:
            print('killing worker', pid)
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
    alive = True
    try:
        while True:
            wn = await wait_until(parent_pid, None)
            if wn == 7:
                print('Parent Done')
                return
            elif wn == 1:
                sibling_pids_addr = peek(parent_pid, tracer_utils.SIBLING_PIDS_ADDR_ADDR)
                pids = read_array(parent_pid, sibling_pids_addr, 1, np.int32)
                wn = await wait_until(pids[0], 1)
                if wn == -1:
                    alive = False
                    print('UNEXPECTEDLY TERMINATED BEFORE PROCESSING STARTED')
                    raise ValueError('UNEXPECTEDLY TERMINATED BEFORE PROCESSING STARTED')

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
                alive = False
                print('PARENT TERMINATED UNEXPECTEDLY')
                raise ValueError('parent terminated unexpectedly')
            # resume parent
            resume(parent_pid)
    finally:
        processing_done.set()
        if alive:
            print('Kill parent', parent_pid)
            await kill(parent_pid)


async def make_worker(pid_queue, file_lock, processing_done):
    while True:
        result = await get_or_done(pid_queue, processing_done)
        if result is None:
            return
        else:
            pid = result
            await handle_segment(pid, file_lock)
        

async def main(exe='./clavrxorb'):
    loop = asyncio.get_event_loop()
    max_workers = 8
    max_parent_lead = 2
    file_lock = asyncio.Lock()
    processing_done = asyncio.Event()
    worker_pids = asyncio.Queue(max_parent_lead)
    workers = []

    # start parent
    parent_pid = start_clavrx('parent', num_clones=1, skip_output=False, redirect=False, exe=exe)
    try:
        spawner_task = loop.create_task(spawner(parent_pid, processing_done, worker_pids))

        for _ in range(max_workers):
            workers.append(loop.create_task(make_worker(worker_pids, file_lock, processing_done)))
        await spawner_task
        await asyncio.gather(*workers)
        print('PSeg Finished')
        return 0
    except (Exception, CancelledError) as main_e:
        print(str(main_e))
        print('PSeg clean-up')
        processing_done.set()
        spawner_task.cancel()
        for worker in workers:
            worker.cancel()
        for worker in workers:
            try:
                await worker
            except CancelledError:
                pass
            except Exception as e:
                print(str(e))
        print('workers clean')
        try:
            await spawner_task
        except CancelledError:
            pass
        except Exception as e:
            print(str(e))
        print('parent clean')

        while True:
            try:
                await kill(worker_pids.get_nowait())
            except asyncio.QueueEmpty:
                break
        print('queue clean')
        raise main_e
        

