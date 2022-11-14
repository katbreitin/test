#include <linux/sched.h>    /* Definition of struct clone_args */
#include <sched.h>          /* Definition of CLONE_* constants */
#include <sys/syscall.h>    /* Definition of SYS_* constants */
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <linux/memfd.h>
#include <sys/mman.h>
#include <fcntl.h>


int getpid_nocache(){
    return syscall(SYS_getpid);
}

void find_heap_bounds(void ** heap_start, void ** heap_end){
   char maps_f[255];
   char * line;
   size_t linelen = 0;
   size_t n = 0;
   FILE * stream;

   sprintf(maps_f, "/proc/%d/maps", getpid_nocache());
   stream = fopen(maps_f, "r");
   while(linelen = getline(&line, &n, stream) > 0){
       if(strstr(line, "[heap]")){
           printf(line);
           char * line2;
           (*heap_start) = (void *) strtol(line, &line2, 16);
           line2++;
           (*heap_end) = (void *) strtol(line2, NULL, 16);
           printf("%x %x\n", (*heap_start), (*heap_end));
           free(line);
           n=0;
           break;
       }
       free(line);
       n=0;
   }
}


int mmap_heap(){
    int fd, pid;
    void * heap_start, *heap_end;
    // read heap start & end
    find_heap_bounds(&heap_start, &heap_end);
    fd = syscall(SYS_memfd_create, "heap", MFD_CLOEXEC);
    if(fd == -1){
       return -1;
    }
    mmap(heap_start, heap_end-heap_start, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
    return 0;
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

