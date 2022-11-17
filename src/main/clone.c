#define _GNU_SOURCE
#include <sched.h>          /* Definition of CLONE_* constants */
#include <sys/syscall.h>    /* Definition of SYS_* constants */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>


int getpid_nocache(){
    return syscall(SYS_getpid);
}

void reopen_files(int parent_pid){
    char proc_self[1024];
    off_t offset = 0;
    for(int fd=3; fd<64; fd++){
        int ret = fcntl(fd, F_GETFD);
        if(ret != -1){
            // Can't use dup(), we need a new "open file descriptors" not just fd
            // otherwise the offset is sync'd between procs
            sprintf(proc_self, "/proc/%d/fd/%d", parent_pid, fd);
            offset = lseek(fd, 0, SEEK_CUR);
            // Get open mode
            int flags = fcntl(fd, F_GETFL);
            flags &= (O_RDWR | O_WRONLY | O_RDONLY);
            close(fd);
            open(proc_self, ret | flags);
            lseek(fd, offset, SEEK_SET);
        }
    }
}

int sibling_clone(){
    void *stack = 0;
    unsigned long flags = CLONE_PARENT;
    long pid;

    pid = syscall(SYS_clone, flags, stack);
    return (int) pid;
}

