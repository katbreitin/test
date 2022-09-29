#!/usr/bin/env python
"""
Build environment configuration script.

This program requires either Python version 2.6+ or 3.2+, and is importable as a
python module.
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
import os
import io
import subprocess
import getopt
import shlex
import shutil

#------------------------------------------------------------------------------
def rm_f(fpath):
    """
    Emulates 'rm -f' functionality.
    """

    if os.path.isfile(fpath) or os.path.islink(fpath):
        os.remove(fpath)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
if "check_output" not in dir(subprocess):
 # subprocess.check_output does not exist (Python 2.6), so use a backported
 #  version from Python 2.7 instead, and "add" it to subprocess:
    def check_output(*popenargs, **kwargs):
        """Run command with arguments and return its output."""

        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            error = subprocess.CalledProcessError(retcode, cmd)
            error.output = output
            raise error
        return output

    subprocess.check_output = check_output
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def main(cbenv_path, anchor_path, build_env_cfg_fpath, supp_defs_fpath,
         top_app_label):
    """Build environment configuration script."""

    # Define some useful parameters:
    cbenv_src_path = cbenv_path
    cbenv_work_path = os.path.join(anchor_path, "work_cbenv")

    # Begin build environment configuration:

    beg_cwd = os.getcwd()  # Record the current working directory

      # If needed, create the work directory:
    if not os.path.isdir(cbenv_work_path):
        os.mkdir(cbenv_work_path)

    os.chdir(cbenv_work_path)   # Enter the work directory

      # Include file(s) resulting from configuration process (and templates):
    pcfg_h = "cfghdr_appind.h"
    pcfg_sh = "config_appind.sh.include"
    pcfg_mk = "config_appind.mk.include"
    pcfg_h_fpath = os.path.join(cbenv_work_path, pcfg_h)
    pcfg_sh_fpath = os.path.join(cbenv_work_path, pcfg_sh)
    pcfg_mk_fpath = os.path.join(cbenv_work_path, pcfg_mk)

      # Check for template include files' existence:
    for ifile in [pcfg_h, pcfg_sh, pcfg_mk]:
        t_fpath = os.path.join(cbenv_src_path, ifile+"-template")
        if not os.path.isfile(t_fpath):
            print(top_app_label+" ERROR: Template file "+t_fpath+ \
                  " does not exist.")
            sys.exit(1)

      # Remove any old include files:
    for i_fpath in [pcfg_h_fpath, pcfg_sh_fpath, pcfg_mk_fpath]:
        rm_f(i_fpath)

      # Check for the existence of build environment configuration file; remove
      #  any prior "local" copy; make a copy of the file in the working
      #  directory:
    if os.path.isfile(build_env_cfg_fpath):
        local_build_env_cfg_fpath = os.path.join(cbenv_work_path,
                                                 "build_env_config")
        rm_f(local_build_env_cfg_fpath)  # Remove file
        shutil.copy(build_env_cfg_fpath, local_build_env_cfg_fpath)
    else:
        print(top_app_label+" ERROR: "+build_env_cfg_fpath+" does not exist.")
        sys.exit(1)

      # Read in build environment configuration file (for 'configure'
      #  command-line arguments below):
    becd = {}
    with open(build_env_cfg_fpath, "rt") as input_f:
        input_f_raw = input_f.read()
    for pair in shlex.split(input_f_raw, True, True):
        split_pair = (pair.replace(';', ' ')).split('=', 1)
        becd[split_pair[0]] = split_pair[1]

      # Remove any prior `configure` debris; force to start from "scratch":
    for ifile in ["config.log", "config.status"]:
        rm_f(os.path.join(cbenv_work_path, ifile))

      # Create working copies of various files:
    src_fn_to_copy_list = ["install-sh", "config.guess", "config.sub",
                           "configure", "config_appind.mk.include-template",
                           "config_appind.sh.include-template",
                           "cfghdr_appind.h-template"]
    for ifile in src_fn_to_copy_list:
        i_fpath = os.path.join(cbenv_src_path, ifile)
        shutil.copy(i_fpath, cbenv_work_path)

      # Run 'configure' to check whether current build environment is "sane":
    try:
        subprocess.check_call([os.path.join(cbenv_work_path, "configure"),
                               "CC="+becd["C_compiler"].strip(),
                               "FC="+becd["Fortran_compiler"].strip()])
    except:
        print (top_app_label+" ERROR: "+ \
               "build environment is broken; inspect "+cbenv_work_path+ \
               "config.log, and then modify build/env_settings/user_change_me.cfg?")
        raise

      # Ingest some configure results:
    pcfgd = {}
    with open(pcfg_sh_fpath, "rt") as input_f:
        input_f_raw = input_f.read()
    for pair in shlex.split(input_f_raw, True, True):
        split_pair = (pair.replace(';', ' ')).split('=', 1)
        pcfgd[split_pair[0]] = split_pair[1]

    Fortran_compiler = pcfgd["Fortran_compiler"].strip()
    Fortran_compiler_bn = os.path.basename(Fortran_compiler)
    C_compiler = pcfgd["C_compiler"].strip()
    C_compiler_bn = os.path.basename(C_compiler)

      # Append "MKDIR_P" definition to C header file:
    MKDIR_P_len_str = str(len(pcfgd["MKDIR_P"]))
    with open(pcfg_h_fpath, "at") as pcfg_h_f:
        text_to_write = "\n"+ \
"/* Define a (hopefully) working \'mkdir -p\' or equivalent workaround.  Note\n"+ \
"    that \'mkdir_p_str\' is defined without a terminating NULL character. */\n"+ \
"#define MKDIR_P_DECL_STRLEN int mkdir_p_strlen = "+MKDIR_P_len_str+";\n"+ \
"#define MKDIR_P_DECLARATION char mkdir_p_string["+MKDIR_P_len_str+"] = \""+pcfgd['MKDIR_P']+"\";\n"
        pcfg_h_f.write(text_to_write)

      # Files specifying the final application-INDEPENDENT build environment:
    appind_bldenvf_sh_fpath = os.path.join(anchor_path,
                                           "build_env_appind.sh.include")
    appind_bldenvf_mk_fpath = os.path.join(anchor_path,
                                           "build_env_appind.mk.include")
    appind_bldenvf_h_fpath = os.path.join(anchor_path, "build_env_appind.h")

      # Remove any prior final specification file(s):
    for i_fpath in [appind_bldenvf_sh_fpath, appind_bldenvf_mk_fpath,
                    appind_bldenvf_h_fpath]:
        rm_f(i_fpath)

      # Copy any final application-independent specification file(s):
    shutil.copy(pcfg_sh_fpath, appind_bldenvf_sh_fpath)
    shutil.copy(pcfg_mk_fpath, appind_bldenvf_mk_fpath)
    shutil.copy(pcfg_h_fpath, appind_bldenvf_h_fpath)

    #=== Create configuration files for utility programs:

      # Create F90+ compiler version information file:

    v_FC_opt_dict = dict( \
                         gfortran = "-dumpversion",
                         ifort = "-V",
                         ifc = "-V",
                         pgf90 = "-V",
                         pgf95 = "-V" \
                        )

    util_Fortran_v_fpath = os.path.join(cbenv_work_path,
                                       "build_env_configure-util_f90.compv.inf")
    rm_f(util_Fortran_v_fpath)

    try:
        text_to_write = subprocess.check_output([Fortran_compiler,
                  v_FC_opt_dict[Fortran_compiler_bn]], stderr=subprocess.STDOUT,
                                                universal_newlines=True)
    except:
        text_to_write = "UNKNOWN"

    with open(util_Fortran_v_fpath, "wt") as version_f:
        version_f.write(text_to_write)

      # Create C compiler version information file:

    v_CC_opt_dict = dict( \
                         gcc = "-v",
                         cc = "-v" \
                        )

    util_C_v_fpath = os.path.join(cbenv_work_path,
                                  "build_env_configure-util_c.compv.inf")
    rm_f(util_C_v_fpath)

    try:
        text_to_write = subprocess.check_output([C_compiler,
                        v_CC_opt_dict[C_compiler_bn]], stderr=subprocess.STDOUT,
                                                universal_newlines=True)
    except:
        text_to_write = "UNKNOWN"

    with open(util_C_v_fpath, "wt") as version_f:
        version_f.write(text_to_write)

      # Create F90+ utility configuration file:

    util_Fortran_cfg_fpath = os.path.join(cbenv_work_path,
                                          "build_env_configure-util_f90.cfg")
    rm_f(util_Fortran_cfg_fpath)

    try:
        text_to_write = unicode(top_app_label+"\n", "ascii")  # v2.6+
    except:
        text_to_write = top_app_label+"\n"  # v3.2+

    with io.open(util_Fortran_cfg_fpath, "wt", encoding="ascii") as cfg_f:
        cfg_f.write(text_to_write)

      # Create C utility configuration file:

    util_C_cfg_fpath = os.path.join(cbenv_work_path,
                                    "build_env_configure-util_c.cfg")
    rm_f(util_C_cfg_fpath)

    text_to_write = \
"sys_canonical_cputype "+pcfgd["sys_canonical_cputype"]+"\n"+ \
"sys_canonical_vendor "+pcfgd["sys_canonical_vendor"]+"\n"+ \
"sys_canonical_os "+pcfgd["sys_canonical_os"]+"\n"+ \
"sys_canonical_endianness "+pcfgd["sys_canonical_endianness"]+"\n"+ \
"compiler_name_C "+C_compiler_bn+"\n"+ \
"compiler_name_Fortran "+Fortran_compiler_bn+"\n"+ \
"compiler_C "+C_compiler+"\n"+ \
"compiler_Fortran "+Fortran_compiler+"\n"+ \
"comp_C_versionf "+util_C_v_fpath+"\n"+ \
"comp_Fortran_versionf "+util_Fortran_v_fpath+"\n"

    with open(util_C_cfg_fpath, "wt") as cfg_f:
        cfg_f.write(text_to_write)

    #=== Build and execute utility programs:

    util_prog_base = "build_env_configure-util"

    for utilp_type in ["F90", 'c']:
          # Check for source code file:
        utilp_src_fpath = os.path.join(cbenv_src_path,
                                       util_prog_base+'.'+utilp_type)
        if not os.path.isfile(utilp_src_fpath):
            print(top_app_label+" ERROR: "+utilp_src_fpath+" does not exist.")
            sys.exit(1)

         # Build utility executable:

        utilp_fpath = os.path.join(cbenv_work_path, 
                                   util_prog_base+'_'+utilp_type+"_exe")
        rm_f(utilp_fpath)
    
        if utilp_type == "F90":
            cmd_args = [Fortran_compiler, "-o", utilp_fpath,
                        pcfgd["Fortran_compiler_flags"],
                        pcfgd["Fortran_compiler_freeform_flag"],
                        utilp_src_fpath]
            tmp_fpath_A = os.path.join(anchor_path, "univ_kind_defs.f90")
            tmp_fpath_B = os.path.join(anchor_path, top_app_label+"_kinds.f90")
            rm_f(tmp_fpath_A) # rm prior file
            rm_f(tmp_fpath_B) # rm prior file
        elif utilp_type == 'c':
            cmd_args = [C_compiler, "-o", utilp_fpath,
                        pcfgd["C_compiler_flags"], utilp_src_fpath]

        cmd_args_clean = []
        for item in cmd_args:
            if not item.isspace():    # Remove any "blank" elements
                for i in item.split():  # Separate any compound arguments
                    cmd_args_clean.append(i)

          # Build utility program:
        subprocess.check_call(cmd_args_clean)

          # Run utility executable:
        utilp_out_fpath = os.path.join(cbenv_work_path,
                                      util_prog_base+'_'+utilp_type+"_outp.txt")
        rm_f(utilp_out_fpath)

        text_to_write = subprocess.check_output([utilp_fpath],
                                                stderr=subprocess.STDOUT,
                                                universal_newlines=True)
        print(text_to_write)
        with open(utilp_out_fpath, "wt") as utilp_out_f:
            utilp_out_f.write(text_to_write)

        if utilp_type == "F90":
            shutil.copy(os.path.join(cbenv_work_path, "kinds.A.f90"),
                        tmp_fpath_A)
            shutil.copy(os.path.join(cbenv_work_path, "kinds.B.f90"),
                        tmp_fpath_B)

    #=== Continue build environment configuration:

      # Files specifying the final application-INDEPENDENT build environment:
    appdep_bldenvf_sh_fpath = os.path.join(anchor_path,
                                           "build_env_appdep.sh.include")
    appdep_bldenvf_mk_fpath = os.path.join(anchor_path,
                                           "build_env_appdep.mk.include")

      # Remove any prior final specification file(s):
    for i_fpath in [appdep_bldenvf_sh_fpath, appdep_bldenvf_mk_fpath]:
        rm_f(i_fpath)

      # Determine system defines ("-D" option arguments):

    def_Fortran_dict = dict( \
                            gfortran = "-DFC_GNU",
                            ifort = "-DFC_INTEL",
                            ifc = "-DFC_INTEL",
                            pgf90 = "-DFC_PGI",
                            pgf95 = "-DFC_PGI",
                            unknown_default = ' ' \
                           )

    try:
        Fortran_compiler_def = def_Fortran_dict[Fortran_compiler_bn]
    except:
        Fortran_compiler_def = def_Fortran_dict["unknown_default"]

    def_C_dict = dict( \
                        gcc = "-DCC_GNU",
                        cc = ' ',
                        icc = "-DCC_INTEL",
                        unknown_default = ' ' \
                        )

    try:
        C_compiler_def = def_C_dict[C_compiler_bn]
    except:
        C_compiler_def = def_C_dict["unknown_default"]
                      
    sys_os_def = ' '   # Might be useful someday...

      # supp_defs_fpath is a sh-syntax file that contains any supplementary
      #  definitions to be added into the build environment at this final stage.

    sdd = dict( \
              supp_C_defs = ' ',
              supp_Fortran_defs = ' ',
              supp_C_includes = ' ',
              supp_Fortran_includes = ' ',
              supp_Fortran_modpaths = ' ',
              supp_C_linker_flags = ' ',
              supp_C_linker_libs = ' ',
              supp_Fortran_linker_flags = ' ',
              supp_Fortran_linker_libs = ' ' \
              )

    if os.path.isfile(supp_defs_fpath):     # Add any provided values
        with open(supp_defs_fpath, "rt") as input_f:
            input_f_raw = input_f.read()
        for pair in shlex.split(input_f_raw, True, True):
            split_pair = (pair.replace(';', ' ')).split('=', 1)
            sdd[split_pair[0]] = split_pair[1]

      # Write final sh-syntax output file:

    text_to_write = \
"# sh-syntax build environment include file (application-dependent)\n"+ \
"\n"+ \
"#~ auto-generated\n"+ \
"\n"+ \
"# C-language:\n"+ \
"\n"+ \
"app_C_compiler_flags=\"\";\n"+ \
"\n"+ \
"app_C_preproc_miscflags=\"\";\n"+ \
"app_C_preproc_defines=\""+sdd["supp_C_defs"]+' '+C_compiler_def+' '+sys_os_def+' '+sdd["supp_C_includes"]+"\";\n"+ \
"app_C_preproc_includes=\""+sdd["supp_C_includes"]+"\";\n"+ \
"\n"+ \
"app_C_linker_libs=\""+sdd["supp_C_linker_libs"]+"\";\n"+ \
"app_C_linker_flags=\""+sdd["supp_C_linker_flags"]+"\";\n"+ \
"\n"+ \
"# C++-language:\n"+ \
"\n"+ \
"app_CXX_compiler_flags=\"\";\n"+ \
"app_CXX_preproc_miscflags=\"\";\n"+ \
"app_CXX_preproc_defines=\"\";\n"+ \
"app_CXX_preproc_includes=\"\";\n"+ \
"\n"+ \
"app_CXX_linker_libs=\"\";\n"+ \
"app_CXX_linker_flags=\"\";\n"+ \
"\n"+ \
"# Fortran-language:\n"+ \
"\n"+ \
"app_Fortran_compiler_flags=\"\";\n"+ \
"\n"+ \
"app_Fortran_preproc_miscflags=\"\";\n"+ \
"app_Fortran_preproc_defines=\""+sdd["supp_Fortran_defs"]+' '+Fortran_compiler_def+' '+sys_os_def+"\";\n"+ \
"app_Fortran_preproc_includes=\""+sdd["supp_Fortran_includes"]+"\";\n"+ \
"app_Fortran_modpaths=\""+sdd["supp_Fortran_modpaths"]+"\";\n"+ \
"\n"+ \
"app_Fortran_linker_libs=\""+sdd["supp_Fortran_linker_libs"]+"\";\n"+ \
"app_Fortran_linker_flags=\""+sdd["supp_Fortran_linker_flags"]+"\";\n"

    with open(appdep_bldenvf_sh_fpath, "wt") as final_f:
        final_f.write(text_to_write)

      # Write final make-syntax output file:

    text_to_write = \
"# make-syntax build environment include file (application-dependent)\n"+ \
"\n"+ \
"#~ auto-generated\n"+ \
"\n"+ \
"# C-language:\n"+ \
"\n"+ \
"app_C_compiler_flags =\n"+ \
"\n"+ \
"app_C_preproc_miscflags =\n"+ \
"app_C_preproc_defines = "+sdd["supp_C_defs"]+' '+C_compiler_def+' '+sys_os_def+' '+sdd["supp_C_includes"]+"\n"+ \
"app_C_preproc_includes = "+sdd["supp_C_includes"]+"\n"+ \
"\n"+ \
"app_C_linker_libs = "+sdd["supp_C_linker_libs"]+"\n"+ \
"app_C_linker_flags = "+sdd["supp_C_linker_flags"]+"\n"+ \
"\n"+ \
"# C++-language:\n"+ \
"\n"+ \
"app_CXX_compiler_flags =\n"+ \
"\n"+ \
"app_CXX_preproc_miscflags =\n"+ \
"app_CXX_preproc_defines =\n"+ \
"app_CXX_preproc_includes =\n"+ \
"\n"+ \
"app_CXX_linker_libs =\n"+ \
"app_CXX_linker_flags =\n"+ \
"\n"+ \
"# Fortran-language:\n"+ \
"\n"+ \
"app_Fortran_compiler_flags =\n"+ \
"\n"+ \
"app_Fortran_preproc_miscflags =\n"+ \
"app_Fortran_preproc_defines = "+sdd["supp_Fortran_defs"]+' '+Fortran_compiler_def+' '+sys_os_def+"\n"+ \
"app_Fortran_preproc_includes = "+sdd["supp_Fortran_includes"]+"\n"+ \
"app_Fortran_modpaths = "+sdd["supp_Fortran_modpaths"]+"\n"+ \
"\n"+ \
"app_Fortran_linker_libs = "+sdd["supp_Fortran_linker_libs"]+"\n"+ \
"app_Fortran_linker_flags = "+sdd["supp_Fortran_linker_flags"]+"\n"

    with open(appdep_bldenvf_mk_fpath, "wt") as final_f:
        final_f.write(text_to_write)

    #=== Return to the initial directory:
    os.chdir(beg_cwd)
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def usage():
    usage_msg = """
usage: build_env_configure.py [-h]  anchor_path  build_env_cfg_fpath
                                    supp_defs_fpath  app_label

Build environment configuration script.

mandatory arguments:
  anchor_path   Directory where the calling script is located
  build_env_cfg_fpath   Build environment configuration file specifying library
                         and include paths, compiler options, etc.
  supp_defs_fpath   sh-syntax file that contains any supplementary definitions
                     to be added into the build environment at the final stage
  app_label   Label/name of the application associated with the calling script

optional arguments:
  -h, --help   show this help message and exit
"""
    print(usage_msg)
#------------------------------------------------------------------------------


if __name__ == "__main__":
    cbenv_path = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),
                                              ".."))

    # Process arguments:
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'h', ["help"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(1)
      # Default option values:
    # NONE
      # Act according to the options given:
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(1)
        else:
            assert False, "unhandled option"
      # Act according to the arguments given:
    if (len(args) > 4 or len(args) < 4):
        print("An incorrect number of arguments were provided.")
        usage()
        sys.exit(1)
    anchor_path = args[0]
    build_env_cfg_fpath = args[1]
    supp_defs_fpath = args[2]
    app_label = args[3]

    # Run driver:
    main(cbenv_path, anchor_path, build_env_cfg_fpath, supp_defs_fpath,
         app_label)
