#!/usr/bin/env python
"""
Top-level admin script for the compilation and build of the CLAVR-x package.

This program requires either Python version 2.6+ or 3.2+, and is importable as
a python module.
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
import os
import subprocess
import getopt
import shlex
import shutil
import glob
import re
import errno

if sys.version_info[0] >= 3:
    import configparser as cp
    from io import StringIO
else:
    import ConfigParser as cp
    from StringIO import StringIO

app_label = "CLAVRx"

#------------------------------------------------------------------------------
def mkdir_p(path):
    """Emulates 'mkdir -p' functionality"""
    if not os.path.isdir(path):
        os.makedirs(path)


#------------------------------------------------------------------------------
def rm_f(fpath):
    """Emulates 'rm -f' functionality."""
    if os.path.isfile(fpath) or os.path.islink(fpath):
        os.remove(fpath)


#------------------------------------------------------------------------------
def parse_via_configparser(cfg, section, entry, boolean=False):
    """Parse a given section & entry via configparser (+ Python-2 kludges)."""
    if not boolean:
        parsed = cfg.get(section, entry)
        if parsed == "____":
            parsed = ''
    else:
        parsed = cfg.getboolean(section, entry)
    return parsed


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
def extend_composite_pnf(anchor_path, pnf_fpath, sdwe_fpath, pnf_pfx,
                         available_ext_libs):
    """Augment a composite pathname file (list of filepaths)."""

    root_src_path = ''  # Initialize
    with open(pnf_fpath, "at") as pnf_f:   # Open for appending
        with open(sdwe_fpath, "rt") as input_f:
            input_f_line_list = input_f.readlines()  # Read a subdirectory list
        for line in input_f_line_list:
            trimmed_line = line.strip()
            if trimmed_line:     # If not completely "whitespace"
                if not root_src_path:   # Set root source code path
                    root_src_path = os.path.abspath( \
                        os.path.join(anchor_path, trimmed_line))
                else:      # A subdirectory to query
                    query_sub_pnf = True   # Assume it will be queried

                    if query_sub_pnf:   # Add to composite pathname list
                        sub_pnf_fpath = os.path.join(root_src_path,
                                                  trimmed_line, "etc", pnf_pfx)
                        if os.path.isfile(sub_pnf_fpath):
                            with open(sub_pnf_fpath, "rt") as sub_pnf_f:
                                  # Read a list of source code files:
                                sub_pnf_f_line_list = sub_pnf_f.readlines()
                            for pnf_line in sub_pnf_f_line_list:
                                if "%EXE" in pnf_line: # Executable source code
                                    exe_part, src_part = \
                                                pnf_line.partition("%:")[::2]
                                    exe_name = re.split("[\(\)]", exe_part)[1]
                                      # Write to a separate file, with the
                                      #  executable name as the filename's
                                      #  suffix:
                                    exe_pnf_fpath = pnf_fpath+'.'+exe_name
                                    line_to_write = os.path.join(root_src_path,
                                                                  src_part)
                                    with open(exe_pnf_fpath, "wt") as exe_pnf_f:
                                        exe_pnf_f.write(line_to_write)
                                elif "%LIB" in pnf_line: # Needs external lib
                                    lib_part, src_part = \
                                                  pnf_line.partition("%:")[::2]
                                    lib_name = re.split("[\(\)]", lib_part)[1]
                                    if lib_name in available_ext_libs:
                                        line_to_append = \
                                                    os.path.join(root_src_path,
                                                                 src_part)
                                        pnf_f.write(line_to_append)
                                elif "%AUXFLAGS" in pnf_line: # Compiler opts
                                    aux_part, src_part = \
                                                  pnf_line.partition("%:")[::2]
                                    aux_ID = re.split("[\(\)]", aux_part)[1]
                                    # Add aux_ID to appended line:
                                    line_to_append = \
                                        os.path.join(root_src_path,
                                                     src_part.strip())+ \
                                                        "@auxflags"+aux_ID+'\n'
                                    pnf_f.write(line_to_append)
                                else:            # Normal source code
                                    trimmed_pnf_line = pnf_line.strip()
                                    if trimmed_pnf_line:  # Not all "whitespace"
                                        # Append root source code path:
                                        line_to_append = \
                                           os.path.join(root_src_path, pnf_line)
                                        pnf_f.write(line_to_append)


#------------------------------------------------------------------------------
def build(anchor_path, *action_args):
    """Build driver."""

    os.chdir(anchor_path)   # Start from a "known" location
    
    #=== Define some useful parameters:

      # Frequently referenced parameters:

    local_build_path = os.path.abspath(anchor_path)
    env_autogen_path = os.path.join(anchor_path,
                                    "env_settings/autogenerated_from_cfg")
    work_path = os.path.join(anchor_path, "work")
    makedir_path = work_path    # Where source code is compiled 
    cbenv_work_path = work_path+"_cbenv"
    local_top_path = os.path.abspath(os.path.join(anchor_path, "..", ".."))
    local_run_path = os.path.join(local_top_path, "run")
    local_src_path = os.path.join(local_top_path, "src")
    local_runexe_path = os.path.join(local_run_path, "exe")
    local_runetc_path = os.path.join(local_run_path, "etc")
    local_runegrs_path = os.path.join(local_run_path, "example-runscripts")
    local_runegnl_path = os.path.join(local_run_path, "example-namelists")
    local_nmltmp_base_path = os.path.join(local_src_path, "base_model",
                                          "namelist_templates")
    local_nmltmp_post_path = os.path.join(local_src_path, "postprocessing",
                                          "namelist_templates")

    bec_other_fpath = os.path.join(anchor_path,
                                   "build_env_config.other_packages")

    build_envf_pfx = "build_env"
    bec_path = os.path.join(anchor_path, build_envf_pfx+"_config")
    pnf_pfx = "path_names"
    univ_kind_defs_fn = "univ_kind_defs.f90"
    app_kinds_fn = app_label+"_kinds.f90"
    src_dirs_w_etc_fn = "src_dirs_with_etc.cfg"
    make_exe = "make"   # Name of default system 'make' program to use

    makefile_fpath = os.path.join(work_path, "Makefile")
    crmf_include_fpath = os.path.join(makedir_path,
                                      build_envf_pfx+".mk.include")
    local_src_dirs_w_etc_fpath = os.path.join(anchor_path, src_dirs_w_etc_fn)
    univ_kind_defs_fpath = os.path.join(local_build_path, "src",
                                        univ_kind_defs_fn)
    app_kinds_fpath = os.path.join(local_build_path, "src", app_kinds_fn)
    build_info_admin_fpath = os.path.join(anchor_path, "build_info")

    mf_template_type_potvalues = ["normal", "debug"]

      # Further argument processing:

    arg_check_action = ["clean", "build", "query_build"]
    arg_check_clean = ["ALL"]
    arg_check_build = ["serial", "DM_only", "SM_only", "hybrid", "debug"]
    arg_check_unary = ["query_build"]
    if action_args[0] in arg_check_action:   # Is the specified action valid?
        this_action = action_args[0]
        len_action_args_list = len(action_args)
        if this_action in arg_check_unary:   # Should not have extraneous args
            if len_action_args_list > 1:   
                usage()   # Disallow extraneous arguments
                sys.exit(1)
        elif this_action == "clean":   # May have 0 or 1 extra argument
            clean_ALL = False   # Default
            for arg_item in action_args[1:]:
                if arg_item in arg_check_clean:  # Is the specified token valid?
                    clean_ALL = True
                else:                      # The specified token is invalid
                    usage()
                    sys.exit(1)
        elif this_action == "build":   # May have 1, 2, or 3 extra arguments
            build_debug = False              # Default build settings
            mf_template_type = "normal"      #
            build_parallel = False           #
            parallel_method = "serial"       #

            if len_action_args_list < 2:
                usage()   # Require at least 1 extra argument
                sys.exit(1)
            else:
                if action_args[1] in arg_check_build:
                    usage()   # Must not be a valid optional arg
                    sys.exit(1)
                else:
                    libexe_fpath = action_args[1]

            if len_action_args_list > 2:
                for arg_item in action_args[2:]:
                    if arg_item in arg_check_build:  # Is the given token valid?
                        if arg_item == "debug":
                            build_debug = True
                            mf_template_type = "debug"
                        elif arg_item == "serial":
                            build_parallel = False
                        else:
                            build_parallel = True
                            parallel_method = arg_item
                    else:                      # The specified token is invalid
                        usage()
                        sys.exit(1)
    else:                         # The specified action is invalid
        usage()
        sys.exit(1)

    #=== Perform specified action(s):

      # Parameters based on "external" paths:
    cbenv_path = os.path.join(anchor_path, "setup_autoconf")
    crmf_fpath = os.path.join(anchor_path, "produce_makefile.py")

    # Add useful Python module search path:
    sys.path.append(cbenv_path)

    if this_action == "clean":
        if os.path.isdir(makedir_path):  # If "make" work directory exists
            if os.path.isfile(makefile_fpath):  # If Makefile exists
                beg_cwd = os.getcwd()  # Record the current working directory
                os.chdir(makedir_path)   # Enter the build "work" directory
                subprocess.check_call([make_exe, "-f", makefile_fpath, "clean"])
                os.chdir(beg_cwd)   # Return to the original directory
            shutil.rmtree(makedir_path)   # Delete "make" work directory
        rm_f(build_info_admin_fpath)    # Remove file generated by admin script

        if clean_ALL:
            rm_f(univ_kind_defs_fpath)    # Remove old source file
            rm_f(app_kinds_fpath)    # Remove old source file
            dir_list_to_remove = [cbenv_work_path, env_autogen_path]
            for dir_path in dir_list_to_remove:
                if os.path.isdir(dir_path):
                    shutil.rmtree(dir_path)  # Delete directory+contents
            for template_type in mf_template_type_potvalues:
                tmp_fpath = os.path.join(anchor_path,
                                       build_envf_pfx+"_config."+template_type)
                rm_f(tmp_fpath)

    elif this_action == "build":
          # If needed, create the "make" work directory:
        if not os.path.isdir(makedir_path):
            os.mkdir(makedir_path)

          # Read the build_info_admin file, parse its contents, and store the
          #  result(s):
#@        biad = {}
#@        with open(build_info_admin_fpath, "rt") as input_f:
#@            input_f_raw = input_f.read()
#@        for pair in shlex.split(input_f_raw, True, True):
#@            split_pair = (pair.replace(';', ' ')).split('=', 1)
#@            biad[split_pair[0]] = split_pair[1]

        # Configure the build environment:
        build_env_cfg_fpath = os.path.join(anchor_path,
                                    build_envf_pfx+"_config."+mf_template_type)
        if not os.path.isfile(crmf_include_fpath):
            supp_defs_fpath = os.path.join(anchor_path,
                                           "supplementary_defs.sh.include")

              # Dictionary of additional needed parallel-related packages:
#@            pad = dict( \
#@                      serial = [' '],
#@                      DM_only = ["MPI_Fortran"],
#@                      SM_only = [' '],
#@                      hybrid = ["MPI_Fortran"] \
#@                      )

              # Write supplemental definitions file:
            text_to_write = \
"supp_C_defs=;\n"+ \
"supp_Fortran_defs=;\n"+ \
"supp_C_includes=\"-I"+makedir_path+"\";\n"+ \
"supp_Fortran_includes=\"-I"+makedir_path+"\";\n"+ \
"supp_Fortran_modpaths=\"-I"+makedir_path+"\";\n"+ \
"supp_C_linker_flags=;\n"+ \
"supp_C_linker_libs=;\n"+ \
"supp_Fortran_linker_flags=;\n"+ \
"supp_Fortran_linker_libs=;\n"
            with open(supp_defs_fpath, "wt") as text_f:
                text_f.write(text_to_write)

            rm_f(univ_kind_defs_fpath)   # Remove old "kind_defs" file
            rm_f(app_kinds_fpath)   # Remove old "kinds" file

              # Generate the "first-best-guess" build environment configuration
              #  file:

            ucm_cfg_fpath_parts = \
                  "env_settings/user_change_me.cfg".split('/')
            ucm_cfg_fpath = os.path.join(*ucm_cfg_fpath_parts)
    
            cfg = cp.RawConfigParser()

            if sys.version_info[0] < 3:
                # Workaround for Python-2 inflexibility:
                with open(ucm_cfg_fpath, "rt") as text_f:
                    ucm_cfg_str_l = []
                    for line in text_f:
                        ls_line = line.lstrip()
                        if len(line)-len(ls_line) < 10:
                            lineA = ls_line
                        else:
                            lineA = line
                        if len(lineA) >= 1:
                            if not (lineA[0] == '#' or lineA[0] == ';'):
                                parts = lineA.split(':')
                                if parts[-1].strip():
                                    if len(parts) > 1:
                                        rhs = parts[-1].split("#;")[0]
                                        lineB = parts[0]+': '+rhs
                                    else:
                                        lineB = parts[0]
                                else:
                                    lineB = lineA.strip()+' ____\n'
                                ucm_cfg_str_l.append(lineB)
                    ucm_cfg_flike = \
                                 StringIO('\n'.join(ucm_cfg_str_l).encode())
                    cfg.readfp(ucm_cfg_flike)
            else:
                cfg.read(ucm_cfg_fpath)

            # Read in some configuration parameters:

            cfg_aux_libs_mnk = []
            
            this_section = "Preferred Compilers"
            cfg_F_compiler = parse_via_configparser(cfg, this_section,
                                                    "Fortran compiler")
            cfg_C_compiler = parse_via_configparser(cfg, this_section,
                                                    "C compiler")

            if build_debug:
                this_section = "Compiler/Linker Options For Debugging"
            else:
                this_section = "Compiler/Linker Options For Production Use"
            cfg_F_preproc_opts = parse_via_configparser(cfg, this_section,
                                     "F preprocessor options (e.g. -D and -I)")
            cfg_F_preproc_opts += " -DISCLAVRX"
            cfg_F_compiler_opts = {}
            cfg_F_compiler_opts["nominal"] = parse_via_configparser(cfg,
                                                                  this_section,
                                   "F compiler options (no -I, -L, -l, or -D)")
            for i in range(10):
                key = "auxflags{:d}".format(i+1)
                token = "aux{:d}".format(i+1)
                try:
                    cfg_F_compiler_opts[key] = parse_via_configparser(cfg,
                                                                  this_section,
                            token+" F compiler options (no -I, -L, -l, or -D)")
                except:
                    cfg_F_compiler_opts[key] = None
            cfg_F_linker_opts = parse_via_configparser(cfg, this_section,
                                              "F linker options (e.g. -L, -l)")
            cfg_F_linker_opts_nol = ' '.join([x for x in
                                              cfg_F_linker_opts.split() if
                                              x[0:1] != "-l"])
            cfg_F_linker_opts_lib = ' '.join([x for x in
                                              cfg_F_linker_opts.split() if
                                              x[0:1] == "-l"])
            cfg_C_preproc_opts = parse_via_configparser(cfg, this_section,
                                     "C preprocessor options (e.g. -D and -I)")
            cfg_C_compiler_opts = parse_via_configparser(cfg, this_section,
                                   "C compiler options (no -I, -L, -l, or -D)")
            cfg_C_linker_opts = parse_via_configparser(cfg, this_section,
                                              "C linker options (e.g. -L, -l)")
            cfg_C_linker_opts_nol = ' '.join([x for x in
                                              cfg_C_linker_opts.split() if
                                              x[0:1] != "-l"])
            cfg_C_linker_opts_lib = ' '.join([x for x in
                                              cfg_C_linker_opts.split() if
                                              x[0:1] == "-l"])

            this_section = \
                        "Capabilities Which Require Special External Libraries"
            cfg_use_GRIB = parse_via_configparser(cfg, this_section,
                                            "NWP input can be in GRIB format?",
                                                  boolean=True)
            cfg_use_libHim = parse_via_configparser(cfg, this_section,
                         "Use libHimawari (needed for HSD-format files only)?",
                                                    boolean=True)
            cfg_use_RTTOV = parse_via_configparser(cfg, this_section,
                     "Use RTTOV (radiative transfer code/data) when possible?",
                                                   boolean=True)
            cfg_use_CRTM = parse_via_configparser(cfg, this_section,
                      "Use CRTM (radiative transfer code/data) when possible?",
                                                  boolean=True)

            cfg_deGRIB_txt = ''
            if cfg_use_GRIB:
                this_section = "Specs, De-GRIB NWP"
                cfg_aux_libs_mnk.append("deGRIB_Fortran")
                cfg_deGRIB_F_incdirs = parse_via_configparser(cfg, this_section,
                                                              "F include dirs")
                cfg_deGRIB_F_libdirs = parse_via_configparser(cfg, this_section,
                                                              "F library dirs")
                cfg_deGRIB_F_preproc_opts = parse_via_configparser(cfg,
                                                                   this_section,
                          "F preprocessor options (e.g. -D, but not -I or -L)")
                cfg_deGRIB_F_linker_opts = parse_via_configparser(cfg,
                                                                  this_section,
                                           "F linker options (e.g. libraries)")
                if cfg_deGRIB_F_preproc_opts:
                    cfg_F_preproc_opts += ' '+cfg_deGRIB_F_preproc_opts
                cfg_deGRIB_txt = \
"#--- de-GRIB software:\n"+ \
"\n"+ \
"deGRIB_Fortran_include_dirs=\""+cfg_deGRIB_F_incdirs+"\";\n"+ \
"deGRIB_Fortran_lib_dirs=\""+cfg_deGRIB_F_libdirs+"\";\n"+ \
"deGRIB_Fortran_libs=\""+cfg_deGRIB_F_linker_opts+"\";\n"

            cfg_libHim_txt = ''
            if cfg_use_libHim:
                this_section = "Specs, libHimawari"
                cfg_aux_libs_mnk.append("libHim_Fortran")
                cfg_libHim_F_incdirs = parse_via_configparser(cfg, this_section,
                                                              "F include dirs")
                cfg_libHim_F_libdirs = parse_via_configparser(cfg, this_section,
                                                              "F library dirs")
                cfg_libHim_F_preproc_opts = parse_via_configparser(cfg,
                                                                   this_section,
                          "F preprocessor options (e.g. -D, but not -I or -L)")
                cfg_libHim_F_linker_opts = parse_via_configparser(cfg,
                                                                  this_section,
                                           "F linker options (e.g. libraries)")
                if cfg_libHim_F_preproc_opts:
                    cfg_F_preproc_opts += ' '+cfg_libHim_F_preproc_opts
                cfg_libHim_txt = \
"#--- libHimawari library:\n"+ \
"\n"+ \
"libHim_Fortran_include_dirs=\""+cfg_libHim_F_incdirs+"\";\n"+ \
"libHim_Fortran_lib_dirs=\""+cfg_libHim_F_libdirs+"\";\n"+ \
"libHim_Fortran_libs=\""+cfg_libHim_F_linker_opts+"\";\n"

            cfg_RTTOV_txt = ''; RTTOV_def = ''
            if cfg_use_RTTOV:
                this_section = "Specs, RTTOV"
                cfg_aux_libs_mnk.append("RTTOV_Fortran")
                cfg_RTTOV_F_incdirs = parse_via_configparser(cfg, this_section,
                                                             "F include dirs")
                cfg_RTTOV_F_libdirs = parse_via_configparser(cfg, this_section,
                                                             "F library dirs")
                cfg_RTTOV_F_preproc_opts = parse_via_configparser(cfg,
                                                                  this_section,
                          "F preprocessor options (e.g. -D, but not -I or -L)")
                cfg_RTTOV_F_linker_opts = parse_via_configparser(cfg,
                                                                 this_section,
                                           "F linker options (e.g. libraries)")
                if cfg_RTTOV_F_preproc_opts:
                    cfg_F_preproc_opts += ' '+cfg_RTTOV_F_preproc_opts
                cfg_RTTOV_txt = \
"#--- RTTOV software:\n"+ \
"\n"+ \
"RTTOV_Fortran_include_dirs=\""+cfg_RTTOV_F_incdirs+"\";\n"+ \
"RTTOV_Fortran_lib_dirs=\""+cfg_RTTOV_F_libdirs+"\";\n"+ \
"RTTOV_Fortran_libs=\""+cfg_RTTOV_F_linker_opts+"\";\n"
                # Determine RTTOV major version number:
                tmp0_fpath = os.path.expandvars(cfg_RTTOV_F_libdirs.split()[0])
                try:
                    tmp_fpath = os.path.join(tmp0_fpath, "..", "src", "main",
                                             "rttov_const.F90")
                    with open(tmp_fpath, "r") as tmp_f:
                        for line in tmp_f:
                            if " version =" in line:
                                vp1 = line.split('=')[-1]
                            elif " release =" in line:
                                vp2 = line.split('=')[-1]
                            elif " minor_version =" in line:
                                vp3 = line.split('=')[-1]
                                break
                    RTTOV_mversion = vp1.strip()
                except:
                    tmp_fpath = os.path.join(tmp0_fpath, "f")
                    for fpath in \
                             glob.glob(tmp_fpath[:-1]+"librttov[0-9]*_main.a"):
                        fpath_bn = os.path.basename(fpath)
                        RTTOV_mversion = fpath_bn.split('_')[0][-2:]
                finally:
                    if RTTOV_mversion:
                        if int(RTTOV_mversion) <= 12:
                            RTTOV_def = "-DRTTOV_LE_V"+RTTOV_mversion
                        else:
                            RTTOV_def = "-DRTTOV_GE_V"+RTTOV_mversion
                    else:
                        print ("ERROR: Cannot determine the RTTOV version")
                        sys.exit(1)

            cfg_CRTM_txt = ''
            if cfg_use_CRTM:
                this_section = "Specs, CRTM"
                cfg_aux_libs_mnk.append("CRTM_Fortran")
                cfg_CRTM_F_incdirs = parse_via_configparser(cfg, this_section,
                                                            "F include dirs")
                cfg_CRTM_F_libdirs = parse_via_configparser(cfg, this_section,
                                                            "F library dirs")
                cfg_CRTM_F_preproc_opts = parse_via_configparser(cfg,
                                                                 this_section,
                          "F preprocessor options (e.g. -D, but not -I or -L)")
                cfg_CRTM_F_linker_opts = parse_via_configparser(cfg,
                                                                this_section,
                                           "F linker options (e.g. libraries)")
                if cfg_CRTM_F_preproc_opts:
                    cfg_F_preproc_opts += ' '+cfg_CRTM_F_preproc_opts
                cfg_CRTM_txt = \
"#--- CRTM software:\n"+ \
"\n"+ \
"CRTM_Fortran_include_dirs=\""+cfg_CRTM_F_incdirs+"\";\n"+ \
"CRTM_Fortran_lib_dirs=\""+cfg_CRTM_F_libdirs+"\";\n"+ \
"CRTM_Fortran_libs=\""+cfg_CRTM_F_linker_opts+"\";\n"

            this_section = "Specs, HDF4"
            cfg_aux_libs_mnk.append("HDF4_Fortran")
            cfg_HDF4_F_incdirs = parse_via_configparser(cfg, this_section,
                                                        "F include dirs")
            cfg_HDF4_F_libdirs = parse_via_configparser(cfg, this_section,
                                                        "F library dirs")
            cfg_HDF4_F_preproc_opts = parse_via_configparser(cfg, this_section,
                          "F preprocessor options (e.g. -D, but not -I or -L)")
            cfg_HDF4_F_linker_opts = parse_via_configparser(cfg, this_section,
                                           "F linker options (e.g. libraries)")
            if cfg_HDF4_F_preproc_opts:
                    cfg_F_preproc_opts += ' '+cfg_HDF4_F_preproc_opts
            cfg_HDF4_txt = \
"#--- HDF4 library:\n"+ \
"\n"+ \
"HDF4_Fortran_include_dirs=\""+cfg_HDF4_F_incdirs+"\";\n"+ \
"HDF4_Fortran_lib_dirs=\""+cfg_HDF4_F_libdirs+"\";\n"+ \
"HDF4_Fortran_libs=\""+cfg_HDF4_F_linker_opts+"\";\n"

            this_section = "Specs, HDF5"
            cfg_aux_libs_mnk.append("HDF5_Fortran")
            cfg_HDF5_F_incdirs = parse_via_configparser(cfg, this_section,
                                                        "F include dirs")
            cfg_HDF5_F_libdirs = parse_via_configparser(cfg, this_section,
                                                        "F library dirs")
            cfg_HDF5_F_preproc_opts = parse_via_configparser(cfg, this_section,
                          "F preprocessor options (e.g. -D, but not -I or -L)")
            cfg_HDF5_F_linker_opts = parse_via_configparser(cfg, this_section,
                                           "F linker options (e.g. libraries)")
            cfg_aux_libs_mnk.append("HDF5_C")
            cfg_HDF5_C_incdirs = parse_via_configparser(cfg, this_section,
                                                        "C include dirs")
            cfg_HDF5_C_libdirs = parse_via_configparser(cfg, this_section,
                                                        "C library dirs")
            cfg_HDF5_C_preproc_opts = parse_via_configparser(cfg, this_section,
                          "C preprocessor options (e.g. -D, but not -I or -L)")
            cfg_HDF5_C_linker_opts = parse_via_configparser(cfg, this_section,
                                           "C linker options (e.g. libraries)")
            if cfg_HDF5_F_preproc_opts:
                cfg_F_preproc_opts += ' '+cfg_HDF5_F_preproc_opts
            if cfg_HDF5_C_preproc_opts:
                cfg_C_preproc_opts += ' '+cfg_HDF5_C_preproc_opts
            cfg_HDF5_txt = \
"#--- HDF5 library:\n"+ \
"\n"+ \
"HDF5_Fortran_include_dirs=\""+cfg_HDF5_F_incdirs+"\";\n"+ \
"HDF5_Fortran_lib_dirs=\""+cfg_HDF5_F_libdirs+"\";\n"+ \
"HDF5_Fortran_libs=\""+cfg_HDF5_F_linker_opts+"\";\n"+ \
"HDF5_C_include_dirs=\""+cfg_HDF5_C_incdirs+"\";\n"+ \
"HDF5_C_lib_dirs=\""+cfg_HDF5_C_libdirs+"\";\n"+ \
"HDF5_C_libs=\""+cfg_HDF5_C_linker_opts+"\";\n"

            this_section = "Specs, NetCDF"
            cfg_aux_libs_mnk.append("NetCDF_Fortran")
            cfg_NetCDF_F_incdirs = parse_via_configparser(cfg, this_section,
                                                          "F include dirs")
            cfg_NetCDF_F_libdirs = parse_via_configparser(cfg, this_section,
                                                          "F library dirs")
            cfg_NetCDF_F_preproc_opts = parse_via_configparser(cfg,
                                                               this_section,
                          "F preprocessor options (e.g. -D, but not -I or -L)")
            cfg_NetCDF_F_linker_opts = parse_via_configparser(cfg, this_section,
                                           "F linker options (e.g. libraries)")
            cfg_aux_libs_mnk.append("NetCDF_C")
            cfg_NetCDF_C_incdirs = parse_via_configparser(cfg, this_section,
                                                          "C include dirs")
            cfg_NetCDF_C_libdirs = parse_via_configparser(cfg, this_section,
                                                          "C library dirs")
            cfg_NetCDF_C_preproc_opts = parse_via_configparser(cfg,
                                                               this_section,
                          "C preprocessor options (e.g. -D, but not -I or -L)")
            cfg_NetCDF_C_linker_opts = parse_via_configparser(cfg, this_section,
                                           "C linker options (e.g. libraries)")
            if cfg_NetCDF_F_preproc_opts:
                cfg_F_preproc_opts += ' '+cfg_NetCDF_F_preproc_opts
            if cfg_NetCDF_C_preproc_opts:
                cfg_C_preproc_opts += ' '+cfg_NetCDF_C_preproc_opts
            cfg_NetCDF_txt = \
"#--- NetCDF library:\n"+ \
"\n"+ \
"NetCDF_Fortran_include_dirs=\""+cfg_NetCDF_F_incdirs+"\";\n"+ \
"NetCDF_Fortran_lib_dirs=\""+cfg_NetCDF_F_libdirs+"\";\n"+ \
"NetCDF_Fortran_libs=\""+cfg_NetCDF_F_linker_opts+"\";\n"+ \
"NetCDF_C_include_dirs=\""+cfg_NetCDF_C_incdirs+"\";\n"+ \
"NetCDF_C_lib_dirs=\""+cfg_NetCDF_C_libdirs+"\";\n"+ \
"NetCDF_C_libs=\""+cfg_NetCDF_C_linker_opts+"\";\n"
            
            # Write sh-syntax output file (for use by 'configure'):
            F_compiler_opts_line_l = []
            for key in cfg_F_compiler_opts:
                if cfg_F_compiler_opts[key] is not None:
                    F_compiler_opts_line_l.append("Fortran_compiler_opts_"+ \
                                      key+"=\""+cfg_F_compiler_opts[key]+"\";")
            F_compiler_opts_lines = '\n'.join(F_compiler_opts_line_l)
            text_to_write = \
"#+ sh-syntax include file, specifying specific parameter values and tasks\n"+ \
"#   for 'configure'.\n"+ \
"\n"+ \
"#~ auto-generated\n"+ \
"\n"+ \
"# Compiler and general compiler/linker flags:\n"+ \
"\n"+ \
"#--- Fortran-language:\n"+ \
"\n"+ \
"Fortran_compiler="+cfg_F_compiler+";\n"+ \
F_compiler_opts_lines+'\n'+ \
"Fortran_preprocessor_opts=\""+cfg_F_preproc_opts+' '+RTTOV_def+"\";\n"+ \
"Fortran_linker_opts=\""+cfg_F_linker_opts_nol+"\";\n"+ \
"Fortran_libs=\""+cfg_F_linker_opts_lib+"\";\n"+ \
"\n"+ \
"#--- C-language:\n"+ \
"\n"+ \
"C_compiler="+cfg_C_compiler+";\n"+ \
"C_compiler_opts=\""+cfg_C_compiler_opts+"\";\n"+ \
"C_preprocessor_opts=\""+cfg_C_preproc_opts+"\";\n"+ \
"C_linker_opts=\""+cfg_C_linker_opts_nol+"\";\n"+ \
"C_libs=\""+cfg_C_linker_opts_lib+"\";\n"+ \
"\n"+ \
"# Information for external libraries/software:\n"+ \
"\n"+ \
cfg_deGRIB_txt+ \
"\n"+ \
cfg_libHim_txt+ \
"\n"+ \
cfg_RTTOV_txt+ \
"\n"+ \
cfg_CRTM_txt+ \
"\n"+ \
cfg_HDF4_txt+ \
"\n"+ \
cfg_HDF5_txt+ \
"\n"+ \
cfg_NetCDF_txt+ \
"\n"+ \
"# Other configuration information:\n"+ \
"\n"+ \
"# Programming languages needed:\n"+ \
"prog_languages_needed=\"Fortran C\";\n"+ \
"\n"+ \
"# Packages to check for:\n"+ \
"aux_libs_needed=\""+' '.join(cfg_aux_libs_mnk)+"\";\n"+ \
"\n"+ \
"# Fortran:\n"+ \
"\n"+ \
"LDFLAGS_Fortran_cfg_gen=\"${Fortran_linker_opts}\";\n"+ \
"LIBS_Fortran_cfg_gen=\"${Fortran_libs}\";\n"+ \
"FCFLAGS_cfg_gen=\"${Fortran_compiler_opts_nominal}\";\n"+ \
"CPPFLAGS_Fortran_cfg_gen=\"${Fortran_preprocessor_opts}\";\n"+ \
"\n"+ \
"# C:\n"+ \
"\n"+ \
"LDFLAGS_C_cfg_gen=\"${C_linker_opts}\";\n"+ \
"LIBS_C_cfg_gen=\"${C_libs}\";\n"+ \
"CFLAGS_cfg_gen=\"${C_compiler_opts}\";\n"+ \
"CPPFLAGS_C_cfg_gen=\"${C_preprocessor_opts}\";\n"

            with open(build_env_cfg_fpath, "wt") as text_f:
                text_f.write(text_to_write)

            # Determine the build environment:
            import build_env_configure
            build_env_configure.main(cbenv_path, anchor_path,
                                     build_env_cfg_fpath, supp_defs_fpath,
                                     app_label)

              # Move new "kind_defs" source code file to its proper location:
            fpath_to_move = os.path.join(anchor_path, univ_kind_defs_fn)
            shutil.move(fpath_to_move, univ_kind_defs_fpath)

              # Move new "kinds" source code file to its proper location:
            fpath_to_move = os.path.join(anchor_path, app_kinds_fn)
            shutil.move(fpath_to_move, app_kinds_fpath)

              # Move other needed files to their proper locations:
            fpath_globs_to_move = [ \
                      os.path.join(anchor_path, build_envf_pfx+"_app*.include"),
                      os.path.join(anchor_path, build_envf_pfx+"_app*.h")]
            for fpath_glob in fpath_globs_to_move:
                for fpath in glob.glob(fpath_glob):
                    shutil.copy(fpath, makedir_path)
                    rm_f(fpath)
            shutil.copy(supp_defs_fpath, makedir_path)
            rm_f(supp_defs_fpath)

        # Read build environment information files, parse their contents, and
        #  save the results in a dictionary for later reference:

        abd = {}

        with open(os.path.join(makedir_path,
                         build_envf_pfx+"_appind.sh.include"), "rt") as input_f:
            input_f_raw = input_f.read()
        for pair in shlex.split(input_f_raw, True, True):
            split_pair = (pair.replace(';', ' ')).split('=', 1)
            abd[split_pair[0]] = split_pair[1] 

        with open(os.path.join(makedir_path,
                         build_envf_pfx+"_appdep.sh.include"), "rt") as input_f:
            input_f_raw = input_f.read()
        for pair in shlex.split(input_f_raw, True, True):
            split_pair = (pair.replace(';', ' ')).split('=', 1)
            abd[split_pair[0]] = split_pair[1]

        becd = {}
        with open(build_env_cfg_fpath, "rt") as input_f:
            input_f_raw = input_f.read()
        for pair in shlex.split(input_f_raw, True, True):
            split_pair = (pair.replace(';', ' ')).split('=', 1)
            becd[split_pair[0]] = split_pair[1]

        # List of usable external library monikers:
        cfg_aux_libs_mnk = abd["usable_aux_libs"].split()

          # Dictionary of additional parallel-related flags:
          #   [ Fortran_defs, Fortran_flags ]
        pfd = dict( \
                  serial = [' ', ' '],
                  DM_only = ["-DUSE_MPI", ' '],
                  SM_only = [' ', abd["Fortran_OpenMP_flags"]],
                  hybrid = ["-DUSE_MPI", abd["Fortran_OpenMP_flags"]] \
                  )

    # Create composite pathname file (list); main code library:

        pnf_fpath = os.path.join(makedir_path, pnf_pfx)
        exe_pnf_fpath = pnf_fpath+'.'+os.path.basename(libexe_fpath)
        if not (os.path.isfile(pnf_fpath) and os.path.isfile(exe_pnf_fpath)):
            extend_composite_pnf(anchor_path, pnf_fpath,
                                 local_src_dirs_w_etc_fpath, pnf_pfx,
                                 cfg_aux_libs_mnk)

        # Create build environment configuration file for Makefile creation:
        F_compiler_opts_l = []
        for key in becd:
            if "Fortran_compiler_opts" in key:
                relkey = key.split('_')[-1]
                if relkey != "nominal":
                    F_compiler_opts_l.append("FFLAGS_"+relkey.upper()+" = "+ \
                                             becd[key]+' '+ \
                                           abd["Fortran_compiler_added_flags"])
        F_compiler_opts_lines = '\n'.join(F_compiler_opts_l)
        text_to_write = \
"# make-syntax build environment configuration file for Makefile creation\n"+ \
"\n"+ \
"#~ auto-generated\n"+ \
"\n"+ \
"FC = "+abd["Fortran_compiler"]+"\n"+ \
"FFLAGS = "+abd["Fortran_compiler_flags"]+abd["app_Fortran_compiler_flags"]+pfd[parallel_method][1]+"\n"+ \
F_compiler_opts_lines+'\n'+ \
"FPPDEFS = "+abd["Fortran_preproc_defines"]+abd["app_Fortran_preproc_defines"]+pfd[parallel_method][0]+"\n"+ \
"FPPFLAGS = "+abd["Fortran_preproc_miscflags"]+abd["app_Fortran_preproc_miscflags"]+abd["Fortran_preproc_includes"]+abd["app_Fortran_preproc_includes"]+pfd[parallel_method][0]+"\n"+ \
"FMODPATHS = "+abd["app_Fortran_modpaths"]+abd["Fortran_preproc_includes"]+abd["app_Fortran_preproc_includes"]+"\n"+ \
"\n"+ \
"CC = "+abd["C_compiler"]+"\n"+ \
"CFLAGS = "+abd["C_compiler_flags"]+abd["app_C_compiler_flags"]+"\n"+ \
"CPPDEFS = "+abd["C_preproc_defines"]+abd["app_C_preproc_defines"]+"\n"+ \
"CPPFLAGS = "+abd["C_preproc_miscflags"]+abd["app_C_preproc_miscflags"]+abd["C_preproc_includes"]+abd["app_C_preproc_includes"]+"\n"+ \
"\n"+ \
"LD = "+abd["Fortran_compiler"]+"\n"+ \
"LDFLAGS = "+abd["Fortran_linker_flags"]+abd["Fortran_linker_libs"]+abd["C_linker_flags"]+abd["C_linker_libs"]+abd["app_Fortran_linker_flags"]+abd["app_Fortran_linker_libs"]+abd["app_C_linker_flags"]+abd["app_C_linker_libs"]+"\n"+ \
"\n"
        with open(crmf_include_fpath, "wt") as text_f:
            text_f.write(text_to_write)

        # Generate a Makefile, then execute its instructions:

        import produce_makefile
        produce_makefile.crmf('', False, False, makefile_fpath, libexe_fpath,
                              crmf_include_fpath, '', '',
                              [pnf_fpath, exe_pnf_fpath])

        beg_cwd = os.getcwd()  # Record the current working directory
        os.chdir(makedir_path)   # Enter the build "work" directory
          # Use parallel compilation if available:
        try:
            subprocess.check_call([make_exe, "-f", makefile_fpath, "-j", '4'])
        except:
            subprocess.check_call([make_exe, "-f", makefile_fpath])

          # Special actions for 'prams':
        if "prams" in os.path.basename(libexe_fpath):
              # Create reference files needed to construct namelist template:
            fpath = os.path.join(makedir_path, "namelist_path-common")
            with open(fpath, "wt") as text_f:
                text_f.write(common_nmltmp_base_path)
            fpath = os.path.join(makedir_path, "namelist_path-bspec")
            with open(fpath, "wt") as text_f:
                text_f.write(local_nmltmp_base_path)

            c_p = []   # Initialize copy/processing list

              # Construct the namelist template file; copy it to the 'run' dir:
            dpf = local_nmltmp_base_path
            dpt = makedir_path
            c_p.append(([dpf, "PRAMS_IN.py"], dpt, []))
            c_p.append(([dpf, "PRAMS_IN.template"], dpt,
                        [univ_proc_tfile_fpath, "PRAMS_IN"]))
            c_p.append(([dpt, "PRAMS_IN"], local_runegnl_path, []))

              # Copy stock control scripts to the 'run/exe' directory:
            dpt = local_runexe_path
            dpf = os.path.join(common_src_path, "base_model", "control_scripts")
            c_p.append(([dpf, "reformat_PRAMS_namelist.py"], dpt, []))
            c_p.append(([dpf, "run_PRAMS.py"], dpt, []))
            dpf = os.path.join(local_src_path, "prepare_input_data",
                               "control_scripts")
            c_p.append(([dpf, "run_prep.py"], dpt, []))

              # Copy stock control scripts to the 'run' directory:
            dpt = local_run_path
            dpf = os.path.join(common_src_path, "base_model", "control_scripts")
            c_p.append(([dpf, "run_PRAMS"], local_runegrs_path, []))
            c_p.append(([dpf, "run_PRAMS-simulation.pbs"], local_runegrs_path,
                        []))
            c_p.append(([dpf, "run_PRAMS-gridcheck.pbs"], local_runegrs_path,
                        []))
            c_p.append(([dpf, "run_PRAMS-sfcfile.pbs"], local_runegrs_path, []))
            c_p.append(([dpf, "run_PRAMS-varfile.pbs"], local_runegrs_path, []))
            c_p.append(([dpf, "run_PRAMS-simulation.sbatch"],
                        local_runegrs_path, []))
            c_p.append(([dpf, "run_PRAMS-gridcheck.sbatch"], local_runegrs_path,
                        []))
            c_p.append(([dpf, "run_PRAMS-sfcfile.sbatch"], local_runegrs_path,
                        []))
            c_p.append(([dpf, "run_PRAMS-varfile.sbatch"], local_runegrs_path,
                        []))
            c_p.append(([dpf, "plot_PRAMS_vdt_from_log.py"], dpt, []))
            c_p.append(([dpf, "PRAMS_run_tool.py"], dpt, []))
            dpf = os.path.join(local_src_path, "base_model", "control_scripts")
            c_p.append(([dpf, "run_display_PRAMS_grids.py"], dpt, []))
            dpf = os.path.join(local_src_path, "prepare_input_data",
                               "control_scripts")
            c_p.append(([dpf, "run_prep"], local_runegrs_path, []))
            c_p.append(([dpf, "run_prep.pbs"], local_runegrs_path, []))
            c_p.append(([dpf, "run_prep.sbatch"], local_runegrs_path, []))
            dpf = os.path.join(local_src_path, "base_model",
                               "namelist_templates")
            c_p.append(([dpf, "run_display_PRAMS_grids.cfg"],
                        local_runegnl_path, []))

              # Copy stock info files to the 'run/etc' directory:
            dpt = local_runetc_path
            dpf = common_etc_path
            c_p.append(([dpf, "version_info"],
                        os.path.join(dpt, "version_info.common"), []))
            dpf = univ_etc_path
            c_p.append(([dpf, "version_info"],
                        os.path.join(dpt, "version_info.ulib"), []))
            dpf = local_etc_path
            c_p.append(([dpf, "version_info"],
                        os.path.join(dpt, "version_info.body"), []))
            dpf = os.path.join(local_src_path, "base_model")
            c_p.append(([dpf, "VTABLE"], dpt, []))

            # Perform list of copy/processing actions: 
            for item in c_p:
                copy_from_pl, copy_to, spcc_cmd = item
                shutil.copy(os.path.join(*copy_from_pl), copy_to)
                if spcc_cmd:
                    subprocess.check_call(spcc_cmd)

          # Special actions for 'postp':
        if "postp" in os.path.basename(libexe_fpath):
              # Create reference files needed to construct namelist template:
            fpath = os.path.join(makedir_path, "namelist_path-common")
            with open(fpath, "wt") as text_f:
                text_f.write(common_nmltmp_post_path)
            fpath = os.path.join(makedir_path, "namelist_path-bspec")
            with open(fpath, "wt") as text_f:
                text_f.write(local_nmltmp_post_path)

            c_p = []   # Initialize copy/processing list

              # Construct the namelist template file; copy it to the 'run' dir:
            dpf = local_nmltmp_post_path
            dpt = makedir_path
            c_p.append(([dpf, "POSTP_IN.py"], dpt, []))
            c_p.append(([dpf, "POSTP_IN.template"], dpt,
                        [univ_proc_tfile_fpath, "POSTP_IN"]))
            c_p.append(([dpt, "POSTP_IN"], local_runegnl_path, []))

              # Copy stock control scripts to the 'run/exe' directory:
            dpt = local_runexe_path
            dpf = os.path.join(common_src_path, "postprocessing",
                               "control_scripts")
            c_p.append(([dpf, "run_postp.py"], dpt, []))

              # Copy stock control scripts to the 'run' directory:
            dpt = local_run_path
            dpf = os.path.join(common_src_path, "postprocessing",
                               "control_scripts")
            c_p.append(([dpf, "run_postp"], local_runegrs_path, []))
            c_p.append(([dpf, "run_postp.pbs"], local_runegrs_path, []))
            c_p.append(([dpf, "run_postp.sbatch"], local_runegrs_path, []))

            # Perform list of copy/processing actions: 
            for item in c_p:
                copy_from_pl, copy_to, spcc_cmd = item
                shutil.copy(os.path.join(*copy_from_pl), copy_to)
                if spcc_cmd:
                    subprocess.check_call(spcc_cmd)

          # Return to the original directory:
        os.chdir(beg_cwd)

    if this_action == "query_build":
        if os.path.isfile(crmf_include_fpath):
            appind_sh_inc_fpath = os.path.join(makedir_path,
                                            build_envf_pfx+"_appind.sh.include")
            appind_mk_inc_fpath = os.path.join(makedir_path,
                                            build_envf_pfx+"_appind.mk.include")
            appind_c_h_fpath = os.path.join(makedir_path,
                                            build_envf_pfx+"_appind.h")
            extra_C_inc_fpath = os.path.join(anchor_path, "..", "..", "src",
                                             'C', "include")

              # Read build environment information files, parse their contents,
              #  and save the results in a dictionary for later reference:
            abd = {}
            with open(appind_sh_inc_fpath, "rt") as input_f:
                input_f_raw = input_f.read()
            for pair in shlex.split(input_f_raw, True, True):
                split_pair = (pair.replace(';', ' ')).split('=', 1)
                abd[split_pair[0]] = split_pair[1]

              # Construct and "deliver" query reply:
            text_to_print = \
"BUILT_LIBS_PATH=\""+makedir_path+"\";\n"+ \
"BUILT_MODS_PATH=\""+makedir_path+"\";\n"+ \
"\n"+ \
"APPIND_SH_INC=\""+appind_sh_inc_fpath+"\";\n"+ \
"APPIND_MK_INC=\""+appind_mk_inc_fpath+"\";\n"+ \
"APPIND_C_H=\""+appind_c_h_fpath+"\";\n"+ \
"\n"+ \
"C_INCLUDES=\"-I"+extra_C_inc_fpath+" "+abd["C_preproc_includes"]+"\";\n"+ \
"FORTRAN_INCLUDES=\""+abd["Fortran_preproc_includes"]+"\";\n"+ \
"\n"+ \
"C_LINKER=\""+abd["C_compiler"]+"\";\n"+ \
"FORTRAN_LINKER=\""+abd["Fortran_compiler"]+"\";\n"+ \
"\n"+ \
"C_LIBS_TO_LINK=\""+abd["C_linker_libs"]+"\";\n"+ \
"C_LINKER_FLAGS=\""+abd["C_linker_flags"]+"\";\n"+ \
"FORTRAN_LIBS_TO_LINK=\""+abd["Fortran_linker_libs"]+"\";\n"+ \
"FORTRAN_LINKER_FLAGS=\""+abd["Fortran_linker_flags"]+"\";\n"+ \
"\n"+ \
"FAILED_LIB_CHECKS=\""+abd["failed_aux_libs"]+"\";\n"+ \
"USABLE_LIB_CHECKS=\""+abd["usable_aux_libs"]+"\";\n"
            print(text_to_print)
        else:
            print("build_script.py ERROR: "+crmf_include_fpath+ \
                  " must exist for query_build mode to work.")
            sys.exit(1)


#------------------------------------------------------------------------------
def main(action):
    """Administrative task driver."""

    #=== Define some useful parameters:

    beg_cwd = os.getcwd()  # Record the current working directory

      # Determine fully-qualified filesystem location of this script:
    anchor_path = os.path.abspath(os.path.dirname(sys.argv[0]))
      # Start from a "known" location; in this case it should be 'build/'
    os.chdir(os.path.join(anchor_path, ".."))

      # Frequently referenced parameters:

    local_exe_path = os.path.join(anchor_path, "..", "run", "bin")
    local_runetc_path = os.path.join(anchor_path, "..", "run", "etc")
    local_runeg_l2_path = os.path.join(anchor_path, "..", "run",
                                       "examples-level2_list")
    local_runeg_co_path = os.path.join(anchor_path, "..", "run",
                                       "examples-clavrx_options")
    local_runeg_fl_path = os.path.join(anchor_path, "..", "run",
                                       "examples-file_list")
    local_build_path = anchor_path
#    version_fpath = os.path.join(local_build_path, "etc", "version_info")

      # Further argument processing:

    arg_check_action = ["version", "clean", "build", "create_msetup"]
    arg_check_clean = ["ALL"]
    arg_check_build = ["serial", "DM_only", "SM_only", "hybrid", "debug"]
    arg_check_unary = ["version", "create_msetup"]
    if action[0] in arg_check_action:   # Is the specified action valid?
        this_action = action[0]
        len_args_action_list = len(action)
        if this_action in arg_check_unary:   # Should not have extraneous args
            if len_args_action_list > 1:   
                usage()   # Disallow extraneous arguments
                sys.exit(1)
        elif this_action == "clean":   # May have 0 or 1 extra argument
            clean_ALL = False   # Default
            for arg_item in action[1:]:
                if arg_item in arg_check_clean:  # Is the specified token valid?
                    clean_ALL = True
                else:                      # The specified token is invalid
                    usage()
                    sys.exit(1)
        elif this_action == "build":   # May have 0, 1, or 2 extra arguments
            if len_args_action_list > 1:
                for arg_item in action[1:]:
                    if arg_item in arg_check_build:  # Is the token valid?
                        pass
                    else:                      # The specified token is invalid
                        usage()
                        sys.exit(1)
    else:                         # The specified action is invalid
        usage()
        sys.exit(1)

    #=== Perform specified action(s):

    if this_action == "version":
          # Read the application information file, parse its contents, and
          #  display the result(s):
        vid = {}
        with open(version_fpath, "rt") as version_f:
            version_f_raw = version_f.read()
        for pair in shlex.split(version_f_raw, True, True):
            split_pair = (pair.replace(';', ' ')).split('=', 1)
            vid[split_pair[0]] = split_pair[1]

        print(app_label+": "+vid["version"].strip())

    elif this_action == "create_msetup":
        # Create a machine_setup file from input configuration info:
        build(anchor_path, this_action)

    elif this_action == "build":

          # Invoke the admin script(s) for any other needed package(s), and
          #  create dictionary of 'query_build' results:
        qbd = {"C_INCLUDES": '', "FORTRAN_INCLUDES": '', "BUILT_MODS_PATH": '',
               "C_LINKER_FLAGS": '', "C_LIBS_TO_LINK": '',
               "BUILT_LIBS_PATH": '', "FORTRAN_LINKER_FLAGS": '',
               "FORTRAN_LIBS_TO_LINK": '', "BUILT_LIBS_PATH": '',
               "FAILED_LIB_CHECKS": '', "USABLE_LIB_CHECKS": ''}
        opn_list = []
        build_info_fpath = os.path.join(local_build_path, "build_info")

        text_to_write = \
"other_C_INCLUDES=\""+qbd["C_INCLUDES"]+"\"\n"+ \
"other_Fortran_INCLUDES=\""+qbd["FORTRAN_INCLUDES"]+"\"\n"+ \
"other_Fortran_MODPATHS=\"-I"+qbd["BUILT_MODS_PATH"]+"\"\n"+ \
"other_C_LDFLAGS=\""+qbd["C_LINKER_FLAGS"]+"\"\n"+ \
"other_C_LIBS=\""+qbd["C_LIBS_TO_LINK"]+" -L"+qbd["BUILT_LIBS_PATH"]+" ".join(opn_list)+"\"\n"+ \
"other_Fortran_LDFLAGS=\""+qbd["FORTRAN_LINKER_FLAGS"]+"\"\n"+ \
"other_Fortran_LIBS=\""+qbd["FORTRAN_LIBS_TO_LINK"]+" -L"+qbd["BUILT_LIBS_PATH"]+" ".join(opn_list)+"\"\n"+ \
"FAILED_EXT_LIBS=\""+qbd["FAILED_LIB_CHECKS"]+"\"\n"+ \
"USABLE_EXT_LIBS=\""+qbd["USABLE_LIB_CHECKS"]+"\"\n"
        with open(build_info_fpath, "wt") as text_f:
            text_f.write(text_to_write)

          # List of libraries and/or executables to build:
        libexe_fn_list = ["lib"+app_label+".a",
                          os.path.join(local_exe_path, "clavrxorb"),
                          os.path.join(local_exe_path, "comp_asc_des_level2b"),
                          os.path.join(local_exe_path, "drive_subset_level2b")]

          # Execute the build script as many times as necessary:
#@        mkdir_p(local_runetc_path)    #
#@        mkdir_p(local_runegrs_path)   # Ensure these directories exist
#@        mkdir_p(local_runegnl_path)   #
        for i in range(len(libexe_fn_list)):
            print("Building "+libexe_fn_list[i]+" ...")
              # Ensure the directory for the exe/lib exists:
            mkdir_p(os.path.abspath(os.path.dirname(libexe_fn_list[i])))
            args_list = [this_action, libexe_fn_list[i]]
            if len_args_action_list > 1:
                for arg in action[1:]:
                    args_list.append(arg)
            build(anchor_path, *args_list)

    elif this_action == "clean":
        # Execute the "cleaning" mode of build():
        print("Cleaning "+app_label+" debris ...")
        args_list = action
        build(anchor_path, *args_list)

#------------------------------------------------------------------------------
def usage():
    """Displays a usage message."""
    usage_msg = """
usage: py_admin.py [-h]  action

Top-level administrative script for the compilation and build of the CLAVR-x
package.  This program requires either Python version 2.6+ or 3.2+, and is
importable as a python module.

mandatory argument:
  action -- One of the following:
      version  :  output the version ID of this code package
      clean  :  Clean up local non-configure build process debris
      clean ALL  :  Clean up all known build process debris, including
                    that of configure
      build [debug]  :  Build package
        "   serial [debug]  :  (DEFAULT) all parallel features disabled
        "   DM_only [debug]  :  distributed-memory parallel (MPI) only
       " "  SM_only [debug]  :  shared-memory parallel (OpenMP) only
       " "  hybrid [debug]  :  hybrid parallel (MPI + OpenMP)
      create_msetup  :  Create a new machine_setup file from configuration

optional arguments:
  -h, --help   show this help message and exit
"""
    print(usage_msg)


#------------------------------------------------------------------------------
if __name__ == "__main__":
    # Retrieve command-line options/arguments:
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
    if (len(args) > 3 or len(args) < 1):
        print("An incorrect number of arguments were provided.")
        usage()
        sys.exit(1)
    action = args

    # Perform action(s):
    main(action)
