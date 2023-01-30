#ifdef LINUX
#define _GNU_SOURCE
#endif
#include <sched.h>          /* Definition of CLONE_* constants */
#include <sys/syscall.h>    /* Definition of SYS_* constants */
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <signal.h>
#include <errno.h>


int getpid_nocache(){
    return syscall(SYS_getpid);
}

void reopen_files(int pid){
    char proc_self[1024];
    int close_after[64];
    off_t offset = 0;
    char path[1024];
    for(int fd=3; fd<64; fd++){
        close_after[fd] = 0;
        int ret = fcntl(fd, F_GETFD);
        if(ret != -1){
            // Can't use dup(), we need a new "open file descriptors" not just fd
            // otherwise the offset is sync'd between procs
            sprintf(proc_self, "/proc/%d/fd/%d", pid, fd);
            // get offset
            offset = lseek(fd, 0, SEEK_CUR);
            // Get open mode
            int flags = fcntl(fd, F_GETFL);

            // Read the file path from /proc/pid/fd/fd
            ssize_t size = readlink(proc_self, path, 1024);
            if(size == -1){
                perror("error reading link");
                printf("%s fd=%d\n", proc_self, fd);
                exit(1);
            } else if (size == 1024){
                printf("filename too long: %s\n", proc_self);
                exit(1);
            }
            // null-terminate
            path[size] = '\0';

            // Close the file descriptor
            if(close(fd) == -1){
                perror("error closing file");
                printf("%s fd=%d\n", proc_self, fd);
                exit(1);
            }
            // Re-open, new_fd should be the lowest unused fd
            // I need to make sure all fd < new_fd are used
            flags &= (O_RDWR | O_WRONLY | O_RDONLY);
            int new_fd = open(path, flags);
            if(new_fd == -1){
                perror("error reopening file");
                printf("%s\n", path);
                exit(1);
            }
            // set the file offset
            if(lseek(fd, offset, SEEK_SET)==-1){
                perror("Error seeking");
                printf("old_fd=%d, new_fd=%d\n", fd, new_fd);
                exit(1);
            }
        } else {
            // I need to make sure all file descriptors < fd are used
            open("/dev/null", O_RDONLY);
            close_after[fd] = 1;
        }
    }
    for(int fd=3; fd<64; fd++){
        if(close_after[fd] == 1){
            if(close(fd) == -1){
                perror("error closing file");
                printf("%d\n", fd);
                exit(1);
            }
        }
    }
}

int sibling_clone(){
#ifdef LINUX
    void *stack = 0;
    unsigned long flags = CLONE_PARENT;
    long pid;

    pid = syscall(SYS_clone, flags, stack);
    return (int) pid;
#else
    printf("ERROR: `sibling_clone` has been only implemented for Linux systems.");
    exit(1);
#endif
}


/*```````````````````````````````````````````````````````````````````*/
int sigstop_to_pid(long int pid_F) {   /* Input */

/*
!+ Send the SIGSTOP signal to a given process, relaying the return value of
    the operation.
*/

  pid_t pid;

  pid = (pid_t) pid_F;
  return kill(pid, SIGSTOP);
}
