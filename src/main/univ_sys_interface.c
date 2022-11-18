/*
!+ Implements useful operations that are difficult (or less portable) to
!   perform within other programming languages (e.g., Fortran).
*/

#include <stdio.h>
#include <stdlib.h>
#include <glob.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <signal.h>
#include <unistd.h>
#include "build_env_appind.h"  /* Macros, etc. from autoconf configure */


/* Compiler name-mangling workarounds: */
#if defined(FC_FUNC_)
#define UFCF FC_FUNC_(univ_filematch_count_f,UNIV_FILEMATCH_COUNT_F)
#define UFF FC_FUNC_(univ_filematch_f,UNIV_FILEMATCH_F)
#define USCF FC_FUNC_(univ_system_cmd_f,UNIV_SYSTEM_CMD_F)
#define UMPF FC_FUNC_(univ_mkdir_p_f,UNIV_MKDIR_P_F)
#define URF FC_FUNC_(univ_remove_f,UNIV_REMOVE_F)
#define RSIH FC_FUNC_(univ_reg_sigint_handler,UNIV_REG_SIGINT_HANDLER)
#define RSTH FC_FUNC_(univ_reg_sigterm_handler,UNIV_REG_SIGTERM_HANDLER)
#define EXIMM FC_FUNC_(univ__exit,UNIV__EXIT)
#else
#define UFCF univ_filematch_count_f
#define UFF univ_filematch_f
#define USCF univ_system_cmd_f
#define UMPF univ_mkdir_p_f
#define URF univ_remove_f
#define RSIH univ_reg_sigint_handler
#define RSTH univ_reg_sigterm_handler
#define EXIMM univ__exit
#endif

typedef void (*sighandler_t)(int);


/*```````````````````````````````````````````````````````````````````*/
void UFCF(int *pattfn_len, char *pattfn,  /* Input */
          int *mfn_cnt, int *ierr) {    /* Output */

/*
!+ Search filesystem for file matching a simple wildcard pattern
!   ('pattfn'; uses only '*' wildcard) -- assumes a Fortran-type input string
!     * If successful, ierr = 0, and sets number of matches ('mfn_cnt')
!     * If UNsuccessful, returns with ierr != 0
*/

  int F_pattfn_len, i;
  char *pattfn_C;
  glob_t glob_stuff;

  F_pattfn_len = *pattfn_len;
  pattfn_C = malloc(F_pattfn_len+1);
  for (i=0; i<F_pattfn_len; i++) {
    pattfn_C[i] = pattfn[i];
  }
  pattfn_C[F_pattfn_len] = '\0';   /* terminate C-type string */

  *ierr = glob(pattfn_C, 0, NULL, &glob_stuff);  /* Matching filename search */

  if (*ierr == 0) {    /* Search succeeded */
    *mfn_cnt = glob_stuff.gl_pathc;
  }
  free(pattfn_C);

}


/*```````````````````````````````````````````````````````````````````*/
void UFF(int *pattfn_len, char *pattfn, int *tmpfn_len, char *tmpfn, /* Input */
         int *ierr) {                 /* Output */

/*
!+ Search filesystem for file matching a simple wildcard pattern
!   ('pattfn'; uses only '*' wildcard) -- assumes a Fortran-type input string
!   * If successful, ierr = 0, opens a file (name='tmpfn') and writes any
!      matching filenames to it
!   * If UNsuccessful, no file is opened, and returns with ierr != 0
*/

  int F_pattfn_len, F_tmpfn_len, i, cnt;
  char *pattfn_C, *tmpfn_C;
  glob_t glob_stuff;
  FILE *fn_list;

  F_pattfn_len = *pattfn_len;
  pattfn_C = malloc(F_pattfn_len+1);
  for (i=0; i<F_pattfn_len; i++) {
    pattfn_C[i] = pattfn[i];
  }
  pattfn_C[F_pattfn_len] = '\0';   /* terminate C-type string */

  *ierr = glob(pattfn_C, 0, NULL, &glob_stuff);  /* Matching filename search */

  if (*ierr == 0) {    /* Search succeeded */
    F_tmpfn_len = *tmpfn_len;
    tmpfn_C = malloc(F_tmpfn_len+1);
    for (i=0; i<F_tmpfn_len; i++) {
      tmpfn_C[i] = tmpfn[i];
    }
    tmpfn_C[F_tmpfn_len] = '\0';   /* terminate C-type string */

    cnt = glob_stuff.gl_pathc;
    fn_list = fopen(tmpfn_C, "w+");
    for (i=0; i<cnt; i++) {  /* Write filename matches */
      fprintf(fn_list, "%s\n", glob_stuff.gl_pathv[i]);
    }
    i = fclose(fn_list);
    free(tmpfn_C);
  }
  free(pattfn_C);

}


/*```````````````````````````````````````````````````````````````````*/
void USCF(int *syscmd_strlen, char *syscmd_string,  /* Input */
          int *ierr) {                       /* Output */

/*
!+ Execute a system shell operation specified by the command string
!   'syscmd_string' and its length ('syscmd_strlen').
!   * If successful, returns with ierr = 0
!   * If UNsuccessful, returns with ierr != 0
*/

  int i_syscmd_strlen, i;
  char *syscmd_string_C;

  i_syscmd_strlen = *syscmd_strlen;
  syscmd_string_C = malloc(i_syscmd_strlen+1);
  for (i=0; i<i_syscmd_strlen; i++) {
    syscmd_string_C[i] = syscmd_string[i];
  }
  syscmd_string_C[i_syscmd_strlen] = '\0';   /* terminate C-type string */

  *ierr = system(syscmd_string_C);

  free(syscmd_string_C);
}


/*```````````````````````````````````````````````````````````````````*/
void UMPF(int *path_strlen, char *path_string,  /* Input */
          int *ierr) {                /* Output */

/*
!+ Recursively creates the path (given the string and its length).
!   * If successful, returns with ierr = 0
!   * If UNsuccessful, returns with ierr != 0
*/

  int i_path_strlen, i;
  char *path_string_C, *p;

  i_path_strlen = *path_strlen;
  path_string_C = malloc(i_path_strlen+1);
  for (i=0; i<i_path_strlen; i++) {
    path_string_C[i] = path_string[i];
  }
  path_string_C[i_path_strlen] = '\0';   /* terminate C-type string */

  *ierr = 0;   /* No error, yet */
  p = path_string_C+1;   /* Start at 2nd char */

  while (*p) {
    if (*p == '/') {
      *p = '\0';   /* Temporarily terminate string at '/' */

      *ierr = mkdir(path_string_C, S_IRWXU|S_IRWXG);
      if (*ierr == -1) {
	if (errno == EEXIST) { *ierr = 0; }  /* No error if directory exists */
        else { free(path_string_C); return; }
      }

      *p = '/';   /* Restore '/' to string */
    }
    p++;
  }

  if (*(p-1) != '/') {   /* Process final path element  */
    *ierr = mkdir(path_string_C, S_IRWXU|S_IRWXG);
    if (*ierr == -1 && errno == EEXIST) *ierr = 0;  /* No error if dir exists */
  }
  
  free(path_string_C);
}


/*```````````````````````````````````````````````````````````````````*/
void URF(int *to_remove_strlen, char *to_remove_string,  /* Input */
         int *ierr) {                       /* Output */

/*
!+ Removes/deletes a file or empty directory specified by the string
!   'to_remove_string' and its length ('to_remove_strlen').
!   * If successful, returns with ierr = 0
!   * If UNsuccessful, returns with ierr != 0
*/

  int i_to_remove_strlen, i;
  char *to_remove_string_C;

  i_to_remove_strlen = *to_remove_strlen;
  to_remove_string_C = malloc(i_to_remove_strlen+1);
  for (i=0; i<i_to_remove_strlen; i++) {
    to_remove_string_C[i] = to_remove_string[i];
  }
  to_remove_string_C[i_to_remove_strlen] = '\0';   /* terminate C-type string */

  *ierr = remove(to_remove_string_C);

  free(to_remove_string_C);
}


/*```````````````````````````````````````````````````````````````````*/
void RSIH(sighandler_t handler) {   /* Input */

/*
!+ Registers a user-specified routine to act as the handler for SIGINT signals.
*/

  signal(SIGINT, handler);
}


/*```````````````````````````````````````````````````````````````````*/
void RSTH(sighandler_t handler) {   /* Input */

/*
!+ Registers a user-specified routine to act as the handler for SIGTERM signals.
*/

  signal(SIGTERM, handler);
}


void EXIMM(int* status) {   /* Input */

/*
!+ Exit immediately, skip exit handlers
*/
    _exit(*status);
    
}
