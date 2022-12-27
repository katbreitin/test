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
import copy
import distutils.dir_util

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
def parse_custom_srcdir(root_path, cfg, section, entry, default_evalue):
    """Return a tuple of absolute default and custom paths for an entry."""
    if cfg.has_option(section, entry):
        srcdir_tmp = parse_via_configparser(cfg, section, entry).strip()
        srcdir_default = os.path.abspath(os.path.join(root_path,
                                           *default_evalue.split(os.path.sep)))
        if len(srcdir_tmp) > 0:
            if os.path.isabs(srcdir_tmp):
                srcdir = os.path.abspath(srcdir_tmp)  # Remove any '.', ".."
            else:
                srcdir = os.path.abspath(os.path.join(root_path,
                                               *srcdir_tmp.split(os.path.sep)))
            if not os.path.isdir(srcdir):
                print ("ERROR: The custom source code directory ("+srcdir+\
                       ") specified in the build environment settings "+ \
                       "configuration file for the entry '"+entry+ \
                       "' does not exist.\n")
                sys.exit(1)
        else:
            srcdir = None  # No custom srcdir specified
        return (srcdir_default, srcdir) 
    else:
        return (None, None)  # Entry does not exist


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
                        sub_pnf_path = os.path.abspath(os.path.join(
                                                  root_src_path, trimmed_line))
                        sub_pnf_fpath = os.path.join(sub_pnf_path, "etc",
                                                     pnf_pfx)
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
                                    line_to_write = os.path.join(sub_pnf_path,
                                                         src_part.strip())+'\n'
                                    with open(exe_pnf_fpath, "wt") as exe_pnf_f:
                                        exe_pnf_f.write(line_to_write)
                                elif "%LIB" in pnf_line: # Needs external lib
                                    lib_part, src_part = \
                                                  pnf_line.partition("%:")[::2]
                                    lib_name = re.split("[\(\)]", lib_part)[1]
                                    if lib_name in available_ext_libs:
                                        line_to_append = \
                                                    os.path.join(sub_pnf_path,
                                                         src_part.strip())+'\n'
                                        pnf_f.write(line_to_append)
                                elif "%AUXFLAGS" in pnf_line: # Compiler opts
                                    aux_part, src_part = \
                                                  pnf_line.partition("%:")[::2]
                                    aux_ID = re.split("[\(\)]", aux_part)[1]
                                    # Add aux_ID to appended line:
                                    line_to_append = \
                                        os.path.join(sub_pnf_path,
                                                           src_part.strip())+ \
                                                        "@auxflags"+aux_ID+'\n'
                                    pnf_f.write(line_to_append)
                                else:            # Normal source code
                                    trimmed_pnf_line = pnf_line.strip()
                                    if trimmed_pnf_line:  # Not all "whitespace"
                                        # Append root source code path:
                                        line_to_append = \
                                             os.path.join(sub_pnf_path,
                                                         trimmed_pnf_line)+'\n'
                                        pnf_f.write(line_to_append)


#------------------------------------------------------------------------------
def build(abs_paths, abs_fpaths, *action_args):
    """Build driver."""

    os.chdir(abs_paths["anchor"])   # Start from a "known" location
    
    #=== Define some useful parameters:

    absl_paths = copy.deepcopy(abs_paths)    # Avoid side-effects
    absl_fpaths = copy.deepcopy(abs_fpaths)  #

      # Frequently referenced parameters:
    absl_paths["env_autogen"] = os.path.join(absl_paths["anchor"],
                                      "env_settings", "autogenerated_from_cfg")
    absl_paths["work"] = os.path.join(absl_paths["anchor"], "work")
    absl_paths["makedir"] = absl_paths["work"]  # Where source code is compiled
    absl_paths["cbenv_work"] = absl_paths["work"]+"_cbenv"

    absl_fpaths["bec_other"] = os.path.join(absl_paths["anchor"],
                                             "build_env_config.other_packages")

    build_envf_pfx = "build_env"
    pnf_pfx = "path_names"
    univ_kind_defs_fn = "univ_kind_defs.f90"
    app_kinds_fn = app_label+"_kinds.f90"
    src_dirs_w_etc_fn = "src_dirs_with_etc.cfg"
    make_exe = "make"  # Name of default system 'make' program to use

    absl_fpaths["makefile"] = os.path.join(absl_paths["work"], "Makefile")
    absl_fpaths["crmf_include"] = os.path.join(absl_paths["makedir"],
                                               build_envf_pfx+".mk.include")
    absl_fpaths["src_dirs_w_etc"] = os.path.join(absl_paths["anchor"],
                                                 src_dirs_w_etc_fn)
    absl_fpaths["autogen_src_dirs_w_etc"] = os.path.join(absl_paths["anchor"],
                                            "autogenerated_"+src_dirs_w_etc_fn)
    absl_fpaths["univ_kind_defs"] = os.path.join(absl_paths["build"], "src",
                                                 univ_kind_defs_fn)
    absl_fpaths["app_kinds"] = os.path.join(absl_paths["build"], "src",
                                            app_kinds_fn)

    mf_template_type_potvalues = ["normal", "debug"]

      # Further argument processing:

    arg_check_action = ["clean", "build"]
    arg_check_clean = ["ALL"]
    arg_check_build = ["serial", "DM_only", "SM_only", "hybrid", "debug"]
    arg_check_unary = []
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
                    absl_fpaths["libexe"] = action_args[1]

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
    absl_paths["cbenv"] = os.path.join(absl_paths["anchor"], "setup_autoconf")
    absl_fpaths["crmf"] = os.path.join(absl_paths["anchor"],
                                       "produce_makefile.py")

    # Add useful Python module search path:
    sys.path.append(absl_paths["cbenv"])

    if this_action == "clean":
        if os.path.isdir(absl_paths["makedir"]):  # If "make" work dir exists
            if os.path.isfile(absl_fpaths["makefile"]):  # If Makefile exists
                beg_cwd = os.getcwd()  # Record the current working directory
                os.chdir(absl_paths["makedir"])   # Enter the build "work" dir
                subprocess.check_call([make_exe, "-f", absl_fpaths["makefile"],
                                                                   "clean"])
                os.chdir(beg_cwd)   # Return to the original directory
            shutil.rmtree(absl_paths["makedir"])   # Delete "make" work dir
        rm_f(absl_fpaths["build_info"])  # Remove file generated by admin script

        if clean_ALL:
            rm_f(absl_fpaths["univ_kind_defs"])    # Remove old source file
            rm_f(absl_fpaths["app_kinds"])    # Remove old source file
            dir_list_to_remove = [absl_paths["cbenv_work"],
                                  absl_paths["env_autogen"]]
            for dir_path in dir_list_to_remove:
                if os.path.isdir(dir_path):
                    shutil.rmtree(dir_path)  # Delete directory+contents
            for template_type in mf_template_type_potvalues:
                tmp_fpath = os.path.join(absl_paths["anchor"],
                                       build_envf_pfx+"_config."+template_type)
                rm_f(tmp_fpath)

    elif this_action == "build":
          # If needed, create the "make" work directory:
        if not os.path.isdir(absl_paths["makedir"]):
            os.mkdir(absl_paths["makedir"])

          # Read the build_info_admin file, parse its contents, and store the
          #  result(s):
#@        biad = {}
#@        with open(build_info_admin_fpath, "rt") as input_f:
#@            input_f_raw = input_f.read()
#@        for pair in shlex.split(input_f_raw, True, True):
#@            split_pair = (pair.replace(';', ' ')).split('=', 1)
#@            biad[split_pair[0]] = split_pair[1]

        # Configure the build environment:
        absl_fpaths["build_env_cfg"] = os.path.join(absl_paths["anchor"],
                                    build_envf_pfx+"_config."+mf_template_type)
        if not os.path.isfile(absl_fpaths["crmf_include"]):
            absl_fpaths["supp_defs"] = os.path.join(absl_paths["anchor"],
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
"supp_C_includes=\"-I"+absl_paths["makedir"]+"\";\n"+ \
"supp_Fortran_includes=\"-I"+absl_paths["makedir"]+"\";\n"+ \
"supp_Fortran_modpaths=\"-I"+absl_paths["makedir"]+"\";\n"+ \
"supp_C_linker_flags=;\n"+ \
"supp_C_linker_libs=;\n"+ \
"supp_Fortran_linker_flags=;\n"+ \
"supp_Fortran_linker_libs=;\n"
            with open(absl_fpaths["supp_defs"], "wt") as text_f:
                text_f.write(text_to_write)

            rm_f(absl_fpaths["univ_kind_defs"])   # Remove old "kind_defs" file
            rm_f(absl_fpaths["app_kinds"])   # Remove old "kinds" file

              # Generate the "first-best-guess" build environment configuration
              #  file:

            ucm_cfg_fpath_parts = [absl_paths["anchor"], "env_settings",
                                   "user_change_me.cfg"]
            absl_fpaths["ucm_cfg"] = os.path.join(*ucm_cfg_fpath_parts)
    
            cfg = cp.RawConfigParser()

            if sys.version_info[0] < 3:
                # Workaround for Python-2 inflexibility:
                with open(absl_fpaths["ucm_cfg"], "rt") as text_f:
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
                    ucm_cfg_flike = StringIO('\n'.join(ucm_cfg_str_l).encode())
                    cfg.readfp(ucm_cfg_flike)
            else:
                cfg.read(absl_fpaths["ucm_cfg"])

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

            this_section = \
                  "Customizable Directories For Certain Source Code Components"
            customize_src_dirs = cfg.has_section(this_section)
            if customize_src_dirs:  # Build config file has this section
                csd_to_proc = []

                this_entry = "ACHA code directory to use"
                csd_to_proc.append(parse_custom_srcdir(absl_paths["pkg_top"],
                                                       cfg, this_section,
                                                       this_entry, "src/acha"))

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
                RTTOV_mversion = None
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
                    if RTTOV_mversion is not None:
                        if int(RTTOV_mversion) <= 12:
                            RTTOV_def = "-DRTTOV_LE_V"+RTTOV_mversion
                        else:
                            RTTOV_def = "-DRTTOV_GE_V"+RTTOV_mversion
                    else:
                        print ("ERROR: Cannot determine the RTTOV version")
                        print (tmp0_fpath)
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

            with open(absl_fpaths["build_env_cfg"], "wt") as text_f:
                text_f.write(text_to_write)

            # Determine the build environment:
            import build_env_configure
            build_env_configure.main(absl_paths["cbenv"],
                                     absl_paths["anchor"],
                                     absl_fpaths["build_env_cfg"],
                                     absl_fpaths["supp_defs"],
                                     app_label)

              # Move new "kind_defs" source code file to its proper location:
            fpath_to_move = os.path.join(absl_paths["anchor"],
                                         univ_kind_defs_fn)
            shutil.move(fpath_to_move, absl_fpaths["univ_kind_defs"])

              # Move new "kinds" source code file to its proper location:
            fpath_to_move = os.path.join(absl_paths["anchor"], app_kinds_fn)
            shutil.move(fpath_to_move, absl_fpaths["app_kinds"])

              # Move other needed files to their proper locations:
            fpath_globs_to_move = [ \
                      os.path.join(absl_paths["anchor"],
                                   build_envf_pfx+"_app*.include"),
                      os.path.join(absl_paths["anchor"],
                                   build_envf_pfx+"_app*.h")]
            for fpath_glob in fpath_globs_to_move:
                for fpath in glob.glob(fpath_glob):
                    shutil.copy(fpath, absl_paths["makedir"])
                    rm_f(fpath)
            shutil.copy(absl_fpaths["supp_defs"], absl_paths["makedir"])
            rm_f(absl_fpaths["supp_defs"])

            # Create customized source directory configuration file:
            if customize_src_dirs:
                with open(absl_fpaths["src_dirs_w_etc"], "rt") as f:
                    template = f.readlines()
                template_rootdir = template[0].strip()
                csd = {}
                try:
                    for default, custom in csd_to_proc:
                        if custom is not None:
                            rel_default = os.path.relpath(default,
                                                          template_rootdir)
                            rel_custom = os.path.relpath(custom,
                                                         template_rootdir)
                            csd[rel_default] = rel_custom
                except:
                    print("ERROR: Unable to determine relative path of a "+ \
                          "custom source directory.  Check "+ \
                          absl_fpaths["ucm_cfg"]+" and "+ \
                          absl_fpaths["src_dirs_w_etc"]+ \
                          " for errors. If you are working on a Microsoft "+ \
                          "Windows system, all source code must be located "+ \
                          "on the same drive.")
                if len(csd) > 0:  # Write a customized file (using a template)
                    with open(absl_fpaths["autogen_src_dirs_w_etc"], "wt") as f:
                        f.write(template[0])  # write root dir
                        for line in template[1:]:
                            line_key = line.strip()
                            if line_key in csd:
                                f.write(csd[line_key]+'\n')  # Customized
                            else:
                                f.write(line)  # Unchanged
                else:  # Simply copy the uncustomized file
                    shutil.copy(absl_fpaths["src_dirs_w_etc"],
                                absl_fpaths["autogen_src_dirs_w_etc"])
            else:
                shutil.copy(absl_fpaths["src_dirs_w_etc"],
                            absl_fpaths["autogen_src_dirs_w_etc"])

        # Read build environment information files, parse their contents, and
        #  save the results in a dictionary for later reference:

        abd = {}

        with open(os.path.join(absl_paths["makedir"],
                        build_envf_pfx+"_appind.sh.include"), "rt") as input_f:
            input_f_raw = input_f.read()
        for pair in shlex.split(input_f_raw, True, True):
            split_pair = (pair.replace(';', ' ')).split('=', 1)
            abd[split_pair[0]] = split_pair[1] 

        with open(os.path.join(absl_paths["makedir"],
                        build_envf_pfx+"_appdep.sh.include"), "rt") as input_f:
            input_f_raw = input_f.read()
        for pair in shlex.split(input_f_raw, True, True):
            split_pair = (pair.replace(';', ' ')).split('=', 1)
            abd[split_pair[0]] = split_pair[1]

        becd = {}
        with open(absl_fpaths["build_env_cfg"], "rt") as input_f:
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

    # Create composite pathname file (list); build main library:

        pnf_fpath = os.path.join(absl_paths["makedir"], pnf_pfx)
        exe_pnf_fpath = pnf_fpath+'.'+os.path.basename(absl_fpaths["libexe"])
        if not (os.path.isfile(pnf_fpath) and os.path.isfile(exe_pnf_fpath)):
            extend_composite_pnf(absl_paths["anchor"], pnf_fpath,
                                 absl_fpaths["autogen_src_dirs_w_etc"], pnf_pfx,
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
        with open(absl_fpaths["crmf_include"], "wt") as text_f:
            text_f.write(text_to_write)

        # Generate a Makefile, then execute its instructions:

        import produce_makefile
        produce_makefile.crmf('', False, False, absl_fpaths["makefile"],
                              absl_fpaths["libexe"],
                              absl_fpaths["crmf_include"], '', '',
                              [pnf_fpath, exe_pnf_fpath])

        beg_cwd = os.getcwd()  # Record the current working directory
        os.chdir(absl_paths["makedir"])   # Enter the build "work" directory
          # Use parallel compilation if available:
        try:
            subprocess.check_call([make_exe, "-f", absl_fpaths["makefile"],
                                   "-j", '4'])
        except:
            subprocess.check_call([make_exe, "-f", absl_fpaths["makefile"]])

          # Special actions for 'clavrxorb':
        if "clavrxorb" in os.path.basename(absl_fpaths["libexe"]):
            dirs_to_copy = ["examples-level2_list", "examples-clavrx_options",
                            "examples-file_list", "example_runscripts"]
            for item in dirs_to_copy:
                dir_from = os.path.abspath(os.path.join(abs_paths["src_ffrun"],
                                                        item))
                dir_to = os.path.abspath(os.path.join(abs_paths["run"],
                                                      item))
                if os.path.isdir(dir_to):
                    shutil.rmtree(dir_to)   # Remove any identical existing dir
                mkdir_p(dir_to)  # Create directory if needed
                distutils.dir_util.copy_tree(dir_from, dir_to, preserve_mode=1,
                                      preserve_times=1, preserve_symlinks=1)

          # Return to the original directory:
        os.chdir(beg_cwd)


#------------------------------------------------------------------------------
def main(action):
    """Administrative task driver."""

    #=== Define some useful parameters:

    beg_cwd = os.getcwd()  # Record the current working directory
    abs_paths = {}
    abs_fpaths = {}

      # Determine fully-qualified filesystem location of this script:
    abs_paths["anchor"] = os.path.abspath(os.path.dirname(sys.argv[0]))
      # Start from a "known" location; in this case it should be 'build/'
    os.chdir(os.path.join(abs_paths["anchor"], ".."))

      # Construct frequently-referenced paths; store in a dictionary:
    paths_cfg = [("build", ("anchor", ['.'])),
                 ("pkg_top", ("anchor", [".."])),
                 ("built_exe", ("anchor", ["..", "run", "bin"])),
                 ("run", ("anchor", ["..", "run"])),
                 ("run_etc", ("run", ["etc"])),
                 ("src_ffrun", ("anchor", ["..", "src", "files_for_runtime"]))]
    for pkey, psegments in paths_cfg:
        psegs = []
        for segment in psegments:
            if isinstance(segment, list):
                for seg in segment:
                    psegs.append(seg)
            else:
                psegs.append(abs_paths[segment])
        abs_paths[pkey] = os.path.abspath(os.path.join(*psegs))

    fpaths_cfg = [("build_info", ("build", ["build_info"]))]
    for pkey, psegments in fpaths_cfg:
        psegs = []
        for segment in psegments:
            if isinstance(segment, list):
                for seg in segment:
                    psegs.append(seg)
            else:
                psegs.append(abs_paths[segment])
        abs_fpaths[pkey] = os.path.abspath(os.path.join(*psegs))

#    version_fpath = os.path.join(local_build_path, "etc", "version_info")

      # Further argument processing:

    arg_check_action = ["version", "clean", "build"]
    arg_check_clean = ["ALL"]
    arg_check_build = ["serial", "DM_only", "SM_only", "hybrid", "debug"]
    arg_check_unary = ["version"]
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

    elif this_action == "build":

          # Invoke the admin script(s) for any other needed package(s), and
          #  create dictionary of any 'query_build' results:
        qbd = {"C_INCLUDES": '', "FORTRAN_INCLUDES": '', "BUILT_MODS_PATH": '',
               "C_LINKER_FLAGS": '', "C_LIBS_TO_LINK": '',
               "BUILT_LIBS_PATH": '', "FORTRAN_LINKER_FLAGS": '',
               "FORTRAN_LIBS_TO_LINK": '', "BUILT_LIBS_PATH": '',
               "FAILED_LIB_CHECKS": '', "USABLE_LIB_CHECKS": ''}
        opn_list = []

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
        with open(abs_fpaths["build_info"], "wt") as text_f:
            text_f.write(text_to_write)

          # List of libraries and/or executables to build:
        libexe_fn_list = ["lib"+app_label+".a",
                          os.path.join(abs_paths["built_exe"], "clavrxorb")]

          # Ensure that needed paths exist; if not, create them:
        tmp_l = [abs_paths["run_etc"]]
        for path in tmp_l:
            mkdir_p(path)  # Ensure that this path exists

          # Execute the build script as many times as necessary:
        for i in range(len(libexe_fn_list)):
            print("Building "+libexe_fn_list[i]+" ...")
              # Ensure the directory for the exe/lib exists:
            mkdir_p(os.path.abspath(os.path.dirname(libexe_fn_list[i])))
            args_list = [this_action, libexe_fn_list[i]]
            if len_args_action_list > 1:
                for arg in action[1:]:
                    args_list.append(arg)
            build(abs_paths, abs_fpaths, *args_list)

    elif this_action == "clean":
        # Execute the "cleaning" mode of build():
        print("Cleaning "+app_label+" debris ...")
        args_list = action
        build(abs_paths, abs_fpaths, *args_list)


#------------------------------------------------------------------------------
def usage():
    """Displays a usage message."""
    usage_msg = """
usage: admin OR py_admin.py  [-h] action

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
