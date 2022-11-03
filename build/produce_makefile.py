#!/usr/bin/env python
"""
Generate a Makefile after determining dependencies, given certain information.

Generate a Makefile after automatically determining source code dependencies,
given the source filepaths and other information (for Fortran 90+ and C source).

This program requires either Python version 2.6+ or 3.2+, and is importable as a
python module.  Only standard library types/methods are utilized for greater
portability.  May need to be invoked via 'produce_makefile.py ...' on some
systems.
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
import re
import os
import getopt
import operator

#--------------------------------------------------------------------------
def file_keyword_scan(obj, src, items_parsed_d, obj_with_module_d,
                      modules_used_d, files_included_d, arg_v, arg_d):
    """
    Scan source code file for 'include', 'use', and 'module' statements.
    """

    if src in items_parsed_d:  # File has already been processed; nothing to do
        return

    module_re = re.compile(r"^\s*module\s+(\w+)", re.IGNORECASE)
    use_re = re.compile(r"^\s*use\s+(\w+)", re.IGNORECASE)
    include_re = re.compile(r"^[#\s]*include\s+['\"<]([\w/\-\.]+)['\">]",
                            re.IGNORECASE)

    if arg_v:
        print("Processing file "+src+" of object "+obj+" ...")

    try:
        with open(src, "rt") as src_f:
            line_list = src_f.readlines()
    except UnicodeDecodeError:
        with open(src, "r", encoding='ISO-8859-1') as src_f:
            line_list = src_f.readlines()

    for line in line_list:
          # Search for Fortran module definitions:
        module_res = module_re.search(line)
        if module_res:
            mod_nm = (module_res.group(1)).lower()
            if mod_nm in obj_with_module_d and mod_nm != "procedure":
                print("\nModule "+mod_nm+" found in "+obj+ \
                      " source as well as in "+obj_with_module_d[mod_nm]+ \
                      " -- AMBIGUOUS.\n")
                sys.exit()
            obj_with_module_d[mod_nm] = obj
          # Search for Fortran 'use' directives:
        use_res = use_re.search(line)
        if use_res:
            modules_used_d[obj] += ' '+(use_res.group(1)).lower()
          # Search for include directives (both C- and Fortran-style):
        include_res = include_re.search(line)
        if include_res:
            if src not in files_included_d:
                files_included_d[src] = ''
            files_included_d[src] += ' '+include_res.group(1)
    items_parsed_d[src] = 1

    if arg_d:
        if modules_used_d[obj]:
            print("   requires modules: "+modules_used_d[obj])
        if src in files_included_d:
            print("   requires includes: "+files_included_d[src])
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
def delve_for_includes(obj, src, srci_sfx_tuple, det_dirs, arg_pathpfx,
                       obj_with_incf_d, mftgt_line_list, includes_d,
                       nonlocal_files_d, paths_to_include, items_parsed_d,
                       obj_with_module_d, modules_used_d, files_included_d,
                       arg_v, arg_d):
    """Recursively scan source code file for 'include' statements/info."""

    include_path_list = []

    if src in files_included_d:
        for incf in files_included_d[src].split():
            if arg_d:
                print("obj = "+obj+" src = "+src+" inc = "+incf)
            (i_head, i_tail) = os.path.split(incf)
            if i_tail.endswith(srci_sfx_tuple):
                if i_head:
                    if i_head == '.':
                        i_head = ''
                    if i_head[0] == os.sep:
                        include_path_list = i_head
                    else:
                        include_path_list = det_dirs
                else:
                    include_path_list = det_dirs
                for incp in include_path_list:
                    incp_tmp = os.path.join(incp, i_head)
                    if incp_tmp == '.'+os.sep:
                        incp_tmp = ''
                    incf_tmp = os.path.join(incp_tmp, i_tail)
                    if arg_pathpfx and incp_tmp[0] != os.sep:
                        incp_tmp = os.path.join("$(SRC_ROOT)", incp_tmp)
                    if arg_d:
                        print("Is "+incf_tmp+" in "+incp+'?')
                    if os.path.isfile(incf_tmp):
                        process_this = False
                        if incf_tmp in obj_with_incf_d:
                            if obj_with_incf_d[incf_tmp] != obj:
                                process_this = True
                        else:
                            process_this = True
                        if process_this:
                            if arg_v:
                                print("   Found "+incf_tmp)
                            mftgt_line_list.append(os.path.join(incp_tmp,
                                                                i_tail))
                            includes_d[incf_tmp] = 1
                            if incp_tmp:
                                nonlocal_files_d[incf_tmp] = 1
                            else:
                                incp_tmp = '.'
                            if incp_tmp not in paths_to_include:
                                paths_to_include.append(incp_tmp)
                            file_keyword_scan(obj, incf_tmp, items_parsed_d,
                                          obj_with_module_d, modules_used_d,
                                          files_included_d, arg_v, arg_d)
                            obj_with_incf_d[incf_tmp] = obj
                            delve_for_includes(obj, incf_tmp, srci_sfx_tuple,
                                   det_dirs, arg_pathpfx, obj_with_incf_d,
                                   mftgt_line_list, includes_d,
                                   nonlocal_files_d, paths_to_include,
                                   items_parsed_d, obj_with_module_d,
                                   modules_used_d, files_included_d, arg_v,
                                   arg_d)
                            break
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
def crmf(arg_pathpfx, arg_d, arg_v, arg_mffpath, arg_finaltgt,
         arg_mfpreamble, arg_incpaths_other, arg_compflags_other, arg_items):
    """ Main method."""

      # arg_d implies arg_v:
    if arg_d:
        arg_v = True

    if arg_v:
        print("Determining/creating the Makefile ("+arg_mffpath+')')

    arg_items_list = ['.']  # Always include current working dir first
    arg_items_list[1:] = arg_items

      # Tuple containing suffixes denoting "source" files:
    srcf_sfx_tuple = (".F", ".F90", ".F95", ".F03", ".F08",
                      ".f", ".f90", ".f95", ".f03", ".f08",
                      ".c")
      # Tuple containing suffixes denoting "include" files:
    srci_sfx_list = [".H", ".fh", ".inc", ".h90", ".h"]
    for i in range(100):
        srci_sfx_list.append(".inc{:02d}".format(i))  # Add the form '.incXX'
    srci_sfx_tuple = tuple(srci_sfx_list)
      # The final target is assumed to be an executable unless its name has one
      #  of the following suffixes:
    ftgt_sfx_tuple = (".a")

    compile_str_d = dict([ \
           (".F", "$(FC) $(FPPDEFS) $(FPPFLAGS) $(FMODPATHS) $(FFLAGS) $(OTHERFLAGS) -c"),
           (".F_auxflags1", "$(FC) $(FPPDEFS) $(FPPFLAGS) $(FMODPATHS) $(FFLAGS_AUXFLAGS1) $(OTHERFLAGS) -c"),
           (".f", "$(FC) $(FMODPATHS) $(FFLAGS) $(OTHERFLAGS) -c"),
           (".f_auxflags1", "$(FC) $(FMODPATHS) $(FFLAGS_AUXFLAGS1) $(OTHERFLAGS) -c"),
           (".F90", "$(FC) $(FPPDEFS) $(FPPFLAGS) $(FMODPATHS) $(FFLAGS) $(OTHERFLAGS) -c"),
           (".F90_auxflags1", "$(FC) $(FPPDEFS) $(FPPFLAGS) $(FMODPATHS) $(FFLAGS_AUXFLAGS1) $(OTHERFLAGS) -c"),
           (".f90", "$(FC) $(FMODPATHS) $(FFLAGS) $(OTHERFLAGS) -c"),
           (".f90_auxflags1", "$(FC) $(FMODPATHS) $(FFLAGS_AUXFLAGS1) $(OTHERFLAGS) -c"),
           (".c", "$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c") \
                       ] )
    str_delimiter_d = dict([ ("'", "'"), ('"', '"'), ('<', '>') ])

    # Begin writing the preamble of the Makefile:
    with open (arg_mffpath, 'w') as mf_f:
        mf_f.write("# Makefile (generated)\n\n")
        mf_f.write("SRC_ROOT = "+arg_pathpfx+"\n\n")
        mf_f.write("OTHERFLAGS = "+arg_compflags_other+"\n\n")
        if arg_mfpreamble:
            mf_f.write("include "+arg_mfpreamble+"\n\n")
        mf_f.write(".DEFAULT:\n\t-echo target $@ is unknown.\n")
        mf_f.write("all: "+arg_finaltgt+"\n")

    items_parsed_d = {}
    orig_src_fpath_d = {}
    src_fpath_d = {}
    paths_to_include = []
    paths_to_include_base = []
    order_of_proc = 0  
    
    for item in arg_items:
        if arg_pathpfx and (item[0] != '/'):    # Prepend path prefix
            item[:0] = arg_pathpfx
        if arg_v:
            print("item = "+item+"\n")

        if os.path.isdir(item) and item not in items_parsed_d:
          # Item is a directory; search for source files:
            if arg_v:
                print("Scanning directory "+item+"\n")
            files_in_dir = os.listdir(item)
            for fid in files_in_dir:
                (i_head, i_tail) = os.path.split(os.path.join(item, fid))
                if i_tail.endswith(srci_sfx_tuple):
                    paths_to_include_base.append(item)  # This is an include dir
                i_sfx_b = i_tail.endswith(srcf_sfx_tuple)
                (object_f, i_ext) = os.path.splitext(i_tail)
                object_f += ".o"
                if i_sfx_b and object_f not in orig_src_fpath_d:
                    if i_head == '.':
                        i_head = ''
                    if arg_pathpfx and (i_head[0] != '/'):
                        src_fpath_d[object_f] = "$(SRC_ROOT)"+os.sep+ \
                                        os.path.join(i_head, i_tail)
                        i_head = os.path.join(arg_pathpfx, i_head)
                    orig_src_fpath_d[object_f] = os.path.join(i_head,
                                                              i_tail)
                    if not object_f in src_fpath_d:
                        src_fpath_d[object_f] = orig_src_fpath_d[object_f]
            order_of_proc += 1
            items_parsed_d[item] = order_of_proc
        elif os.path.isfile(item):
          # Item is a regular file:
            (i_head, i_tail) = os.path.split(item)
            i_sfx_b = i_tail.endswith(srcf_sfx_tuple)
            (object_f, i_ext) = os.path.splitext(i_tail)
            object_f += ".o"
            if not object_f in orig_src_fpath_d:
                if i_sfx_b:
                  # Item is a source file:
                    if i_head == '.':
                        i_head = ''
                    if arg_pathpfx and (i_head[0] != '/'):
                        src_fpath_d[object_f] = "$(SRC_ROOT)"+os.sep+ \
                                        os.path.join(i_head, i_tail)
                        i_head = os.path.join(arg_pathpfx, i_head)
                    orig_src_fpath_d[object_f] = os.path.join(i_head, i_tail)
                    if not object_f in src_fpath_d:
                        src_fpath_d[object_f] = orig_src_fpath_d[object_f]
                else:
                  # Item is likely a file that contains a list of source files:
                    i_sfx_b = i_tail.endswith(srci_sfx_tuple)
                    if not i_sfx_b:   # Only if not an include file
                        if arg_v:
                            print("Reading source file list from "+item+"...")
                        with open (item, 'r') as pnf_f:
                            line_list = pnf_f.readlines()
                            for line in line_list:
                                words = line.split()
                                fif = words[-1]
                                (i_head, i_tail0) = os.path.split(fif)
                                i_tail = i_tail0.split('@')[0]
                                i_sfx_b = i_tail.endswith(srcf_sfx_tuple)
                                if arg_d:
                                    print("fif = "+fif)
                                (object_f, i_ext) = os.path.splitext(i_tail)
                                object_f += ".o"
                                if i_sfx_b and object_f not in orig_src_fpath_d:
                                    if i_head == '.':
                                        i_head = ''
                                    if arg_pathpfx and (i_head[0] != '/'):
                                        src_fpath_d[object_f] = "$(SRC_ROOT)"+ \
                                               os.sep+os.path.join(i_head,
                                                                   i_tail)
                                        i_head = os.path.join(arg_pathpfx,
                                                              i_head)
                                    orig_src_fpath_d[object_f] = \
                                              os.path.join(i_head, i_tail)
                                    if not object_f in src_fpath_d:
                                        src_fpath_d[object_f] = \
                                                     orig_src_fpath_d[object_f]
                                    if not i_head in items_parsed_d:
                                        order_of_proc += 1
                                        items_parsed_d[i_head] = order_of_proc
                                    try:
                                        i_cmd = i_tail0.split('@')[1].strip()
                                    except:
                                        i_cmd = None
                                    if i_cmd:
                                        tmp, i_ext = os.path.splitext(i_tail)
                                        compile_str_d[i_tail] = \
                                                 compile_str_d[i_ext+'_'+i_cmd]
                                        if arg_v:
                                            print("Custom compile command (for "+ \
                                                  i_tail+"): "+i_cmd)
                                if not i_sfx_b:
                                  # Is this an include file?
                                    i_sfx_b = i_tail.endswith(srci_sfx_tuple)
                                    if i_sfx_b and i_head not in items_parsed_d:
                                        order_of_proc += 1
                                        items_parsed_d[i_head] = order_of_proc

    sorted_keys = sorted(items_parsed_d.items(), key=operator.itemgetter(1)) 
    det_dirs = [i[0] for i in sorted_keys]
    det_srcf = src_fpath_d.values()
    det_objf = src_fpath_d.keys()
    if arg_d:
        print("det_dirs =")
        print(det_dirs)
        print("det_srcf =")
        print(det_srcf)
        print("det_objf =")
        print(det_objf)

    obj_with_module_d = {}
    modules_used_d = {}
    files_included_d = {}

    for obj in det_objf:
        modules_used_d[obj] = ''
        file_keyword_scan(obj, orig_src_fpath_d[obj], items_parsed_d,
                          obj_with_module_d, modules_used_d, files_included_d,
                          arg_v, arg_d)

    includes_d = {}
    nonlocal_files_d = {}
    obj_used_by_obj_d = {}

    with open (arg_mffpath, 'a') as mf_f:
        for obj in sorted(det_objf):
            obj_is_used_d = {}
            obj_with_incf_d = {}
            mftgt_line_list = []
            obj_is_used_d[obj] = 1   # Init to avoid recursion
            if arg_v:
                print("Performing dependency collection for "+obj+" ...")
            mftgt_line_list.append(obj+": "+src_fpath_d[obj])
            (i_head, i_tail) = os.path.split(orig_src_fpath_d[obj])
            if i_head != '.' or i_head != '':
                nonlocal_files_d[src_fpath_d[obj]] = 1

            paths_to_include = paths_to_include_base
            delve_for_includes(obj, orig_src_fpath_d[obj], srci_sfx_tuple,
                           det_dirs, arg_pathpfx, obj_with_incf_d,
                           mftgt_line_list, includes_d, nonlocal_files_d,
                           paths_to_include, items_parsed_d, obj_with_module_d,
                           modules_used_d, files_included_d, arg_v, arg_d)
 
            for mod_nm in modules_used_d[obj].split():
                if mod_nm in obj_with_module_d:
                    obj_tmp = obj_with_module_d[mod_nm]
                    if obj_tmp:
                        if obj_tmp not in obj_is_used_d:
                            obj_is_used_d[obj_tmp] = 1
                            mftgt_line_list.append(obj_tmp)
                            obj_used_by_obj_d[obj_tmp] = 1
                            if arg_v:
                                print("   Found module ("+mod_nm+ \
                                      ") in object ("+obj_tmp+")")

              # Write the "target: dependencies" line:
            mf_f.write(' '.join(mftgt_line_list)+"\n")

              # Write the command line(s):
            if i_tail in compile_str_d:
                mf_f.write("\t"+compile_str_d[i_tail])
            else:
                tmp, i_ext = os.path.splitext(i_tail)
                mf_f.write("\t"+compile_str_d[i_ext])
            for incp in paths_to_include:
                mf_f.write(" -I"+incp)
            if arg_incpaths_other:
                for incp in arg_incpaths_other.split():
                    mf_f.write(" -I"+incp)
            mf_f.write("\t"+src_fpath_d[obj]+"\n")

          # Write command(s) to create local copies of source code 
        for key in nonlocal_files_d.keys():
            ftmp = os.path.basename(key)
            ftmp = ftmp.replace("$(SRC_ROOT)", '')
            mf_f.write(os.path.join(".", ftmp)+": "+key+"\n\tcp "+key+" .\n")

          # Create list of objects which are not USEd by code in other objects:
        obj_unused_by_obj_list = []
        for obj in det_objf:
            if obj not in obj_used_by_obj_d:
                obj_unused_by_obj_list.append(obj)

          # Write lists of file/objects to Makefile:
        tmpl = includes_d.keys()
        mf_f.write("SRC = "+' '.join(det_srcf)+' '.join(tmpl)+"\n")
        mf_f.write("OBJ = "+' '.join(det_objf)+"\n")
        if len(nonlocal_files_d) > 0:
            tmpl = nonlocal_files_d.keys()
            mf_f.write("NONLOC = "+' '.join(tmpl)+"\n")

          # Write Makefile postamble targets:
        mf_f.write("clean: neat\n\t-rm -f $(OBJ) "+arg_finaltgt+"\n")
        mf_f.write("neat:\n\t-rm -f $(TMPFILES)\n")
        if len(nonlocal_files_d) > 0:
            mf_f.write("localize: $(NONLOC)\n\tcp $(NONLOC) .\n")
        if arg_finaltgt.endswith(ftgt_sfx_tuple):
            mf_f.write(arg_finaltgt+": $(OBJ)\n\t$(AR) $(ARFLAGS) "+ \
                       arg_finaltgt+" $(OBJ)\n")
        else:
            mf_f.write(arg_finaltgt+": $(OBJ)\n\t$(LD) $(OBJ) -o "+ \
                       arg_finaltgt+" $(LDFLAGS)\n")

    print("Finished determining/creating "+arg_mffpath)
#--------------------------------------------------------------------------


#------------------------------------------------------------------------------
def usage():
    usage_msg = """
usage: produce_makefile.py [-h][-d][-v] [--pathpfx Pathpfx] [--mffpath Mffpath]
                    [--finaltgt Finaltgt] [--mfpreamble Mfpreamble]
                    [--otherincp Incpaths_other] [--otherflags Compflags_other]
                    arg_item1 arg_item2 ... arg_itemN

Generate a Makefile after automatically determining source code dependencies,
given the source filepaths and other information (for Fortran 90+ and C source).
This program requires either Python version 2.6+ or 3.2+, and is importable as
a python module.

mandatory argument(s):
  arg_item(s)   List of source files (i.e., regular files having supported
                suffixes/file_extensions) -AND/OR- directories that contain
                source files -AND/OR- files that contain a list of source files
                (source filepaths are the last contiguous, w.r.t. whitespace,
                strings on each line, and any preceding data is assumed to be
                the compile command to apply to that source file.

optional arguments:
  -h, --help    Show this help message and exit
  --pathpfx     Path appended to the front of all relative sourcefile paths
  -d            Output copious debugging information
  -v            Be verbose about the operations taking place
  --mffpath     Name (full path) of Makefile to create
  --finaltgt    Name of final target (executable or library)
  --mfpreamble  File containing a preamble for the Makefile
  --otherincp   Other include paths
  --otherflags  Other compiler flags
"""
    print(usage_msg)
#------------------------------------------------------------------------------


if __name__ == "__main__":
    # Process arguments:
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hdv",
                                   ["help", "pathpfx=", "mffpath=",
                                    "finaltgt=", "mfpreamble=", "otherincp=",
                                    "otherflags="])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(1)
      # Default option values:
    arg_d = False
    arg_v = False
    pathpfx = ''
    mffpath = "Makefile"
    finaltgt = "a.out"
    mfpreamble = ''
    incpaths_other = ''
    compflags_other = ''
      # Act according to the options given:
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(1)
        elif o == "-d":
            arg_d = True
        elif o == "-v":
            arg_v = True
        elif o == "--pathpfx":
            pathpfx = a
        elif o == "--mffpath":
            mffpath = a
        elif o == "--finaltgt":
            finaltgt = a
        elif o == "--mfpreamble":
            mfpreamble = a
        elif o == "--otherincp":
            incpaths_other = a
        elif o == "--otherflags":
            compflags_other = a
        else:
            assert False, "unhandled option"
      # Act according to the arguments given:
    if (len(args) < 1):
        print("An incorrect number of arguments were provided.")
        usage()
        sys.exit(1)
    arg_items = args

    crmf(pathpfx, arg_d, arg_v, mffpath, finaltgt, mfpreamble, incpaths_other,
         compflags_other, arg_items)
