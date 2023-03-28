#!/usr/bin/env python

#------------------------------------------------------------------------------
# CLAVR-x (CLouds from AVHRR - eXtended) PROCESSING SOFTWARE
#
# COPYRIGHT
# THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
# DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
# FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
# INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
# EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
# DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
# THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
# SUPPORT TO USERS.
#------------------------------------------------------------------------------

"""
Write files based on a template file and a substitution keys/values file.

This is a Python script intended to simplify/automate the creation of one
or more files containing one or more text blocks that are sufficiently similar
("templatable") -- all from two minimalistic input files.  This code requires a
single command-line argument, the prefix (including path, if needed) of the
two input files {input_filename_prefix}.

{input_filename_prefix}.py  :  A Python-syntax file that contains substitution
                               keys/values that will be iteratively applied to
                               the templated text.
    Three named objects are defined by this file:
        'defines_dict'  :  A dictionary of keys with only scalar values that
                           will be used to substitute into the given values of
                           'sdl_string' and 'output_fn_template'.  Special
                           placeholders (e.g., '%SUBFILE%:filepath') will
                           substitute the (string) contents of 'filepath'.
        'sdl_string'  :  A quoted list of dictionaries which have keys with
                         lists as their values.  It must remain quoted
                         (rendering it as a long string) to permit substitution
                         of any tokens within it using the 'defines_dict'
                         dictionary.  Each of the value elements of each
                         dictionary key implies an iteration.  Multiple
                         dictionaries can be used to produce even
                         more complex output.  A single dictionary with key
                         values having 3 elements will generate 3 iterations
                         of template substitution.  Two dictionaries with
                         3 and 4 elements, respectively, will generate 12
                         iterations of template substitution. Special
                         placeholders (e.g., '%SUBFILE%:filepath') will
                         substitute the (string) contents of 'filepath'.
        'output_fn_template'  :  A quoted string denoting the output filename
                                 (with path, if desired).  This will undergo
                                 substitution of any tokens within it using 
                                 both the 'defines_dict' dictionary and a single
                                 iteration of the 'sdl_string' dictionary. Any
                                 missing parent directories will be created, if
                                 possible.

{input_filename_prefix}.template  :  A free-syntax file that contains templated
                                     text with substitution tokens indicated by
                                     $token or ${token} ('$$' results in a
                                     literal '$' character).

This program requires either Python version 2.6+ or 3.2+, and is importable as a
Python module.  Only standard library types/methods are utilized for greater
portability.  May need to be invoked via 'python process_template_file.py ...'
on some systems.
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
import os
import string
import shutil
import io
import getopt

#--------------------------------------------------------------------------
def create_any_needed_dirs(new_dir):
    """ Create a directory, ensuring all parent folders exist.

        * This method may execute itself recursively
        * Any missing parent directories will be created
        * If final directory already exists, return to caller
        * If there is another filesystem object with the same name,
           raise an exception
    """

    if os.path.isdir(new_dir):
        pass
    elif os.path.isfile(new_dir):
        raise OSError("Cannot create directory '{0}'; another filesystem object"
                      " with the same name already exists.".format(new_dir))
    else:
        head, tail = os.path.split(new_dir)
        if head and not os.path.isdir(head):
            create_any_needed_dirs(head)
        if tail:
            os.mkdir(new_dir)
#--------------------------------------------------------------------------


#--------------------------------------------------------------------------
def process_template(input_fn_prefix_raw):
    """Main method."""

    # Expand and parse input filename prefix string:
    input_fn_prefix_full = os.path.abspath(input_fn_prefix_raw)
    input_fn_prefix_dir = os.path.dirname(input_fn_prefix_full)
    input_fn_prefix_base = os.path.basename(input_fn_prefix_full)

    # Add location of input file(s) to the beginning of the module search path;
    # import the relevant file (*.py), which defines the objects 'defines_dict',
    #  'sdl_string', and 'output_fn_template':
    sys.path.insert(0, input_fn_prefix_dir)
    mod_to_import = input_fn_prefix_base
    need_to_delete_tmp_mod = False
    if '.' in mod_to_import:  # Work around a '.' in the input prefix string
        mod_to_import = "tmp_3838299"
        tmp_mod_fn_base = os.path.join(input_fn_prefix_dir, mod_to_import)
        need_to_delete_tmp_mod = True
        shutil.copy2(input_fn_prefix_full+".py", tmp_mod_fn_base+".py")
    inpmod = __import__(mod_to_import,
                        ["output_f_codec", "defines_dict", "sdl_string",
                         "output_fn_template"], -1)

    # For inpmod.defines_dict, substitute file contents (string) for special
    #  placeholders (e.g., "%SUBFILE%:filepath"):
    for (key, value) in inpmod.defines_dict.items():
        if value.startswith("%SUBFILE%:"):
            sub_fn = (value.split(':'))[1]
            if (sub_fn.strip()).startswith(os.sep):
                abs_path = sub_fn.strip()
            else:
                rel_path = os.path.join(input_fn_prefix_dir, sub_fn.strip())
                abs_path = os.path.abspath(rel_path)
            with open (abs_path, "rt") as sub_f:
                inpmod.defines_dict[key] = (sub_f.read()).rstrip()

    # Perform template value substitution on 'inpmod.sdl_string' using the
    #  dictionary 'inpmod.defines_dict'; evaluate resulting string to create a
    #  full-fledged dictionary:
    s = string.Template(inpmod.sdl_string)
    subbed_sdl_string = s.substitute(inpmod.defines_dict)
    sdl = eval(subbed_sdl_string)

    n_dicts = len(sdl)   # Number of dictionaries in 'sdl'

    # For sdl, substitute file contents (string) for special
    #  placeholders (e.g., "%SUBFILE%:filepath"):
    for d_id in range(0,n_dicts):
        for (key, value) in sdl[d_id].items():
            for i_id in range(0,len(value)):
                if value[i_id].startswith("%SUBFILE%:"):
                    sub_fn = (value[i_id].split(':'))[1]
                    if (sub_fn.strip()).startswith(os.sep):
                        abs_path = sub_fn.strip()
                    else:
                        rel_path = os.path.join(input_fn_prefix_dir,
                                                sub_fn.strip())
                        abs_path = os.path.abspath(rel_path)
                    with open(abs_path, "rt") as sub_f:
                        s = string.Template((sub_f.read()).rstrip())
                        subbed_string = s.substitute(inpmod.defines_dict)
                        sdl[d_id][key][i_id] = subbed_string

    # Determine number of iterations needed by each dictionary in 'sdl':
    dict_n_iters = []
    for idx_dict in range(n_dicts):
        tmp_str = list(sdl[idx_dict].keys())[0]   # Sample key in dictionary
        dict_n_iters.append(len(sdl[idx_dict][tmp_str]))

    # Construct an iteration control list:
    iter_control_list = []
    str_to_exec = ''
    icl_str = ''
    for i in range(n_dicts):
        i_str = 'i'+str(i)
        icl_str = icl_str+i_str
        if i < n_dicts-1:
            icl_str = icl_str+", "
        indent = ' '*(i*4)
        str_to_exec = str_to_exec+indent+"for "+i_str+ \
            " in range(dict_n_iters["+str(i)+"]):\n"
    indent = ' '*(n_dicts*4)
    str_to_exec = str_to_exec+indent+"iter_control_list.append(["+icl_str+"])\n"
    exec(str_to_exec)  # Run dynamically-created code to fill iter_control_list

    # Open input template file; read its contents; close the file:
    with open(input_fn_prefix_full+".template", "rt") as template_file:
        s = string.Template(template_file.read())

    # Appropriately substitute tokens in template text; write to resulting text
    #  to filename(s) as configured:
    ofn_list = []
    iter_dict = {}
    iter_dict.update(inpmod.defines_dict)   # Add all entries in 'defines_dict'
    this_is_first_iter = True
    for iter_control in iter_control_list:
        # Finish constructing this iteration's template substitution dictionary:
        for idx_dict in range(n_dicts):
            d = sdl[idx_dict]
            for key in d.keys():
                iter_dict[key] = d[key][iter_control[idx_dict]]

        # Make use of the current template substitution dictionary; write
        #  results to file (starting files and appending to them as
        #  necessary/appropriate):
        if this_is_first_iter:
            s_fn = string.Template(inpmod.output_fn_template)
            tmp_str = s_fn.substitute(iter_dict)
            output_filename = tmp_str
            ofn_list.append(output_filename)
            create_any_needed_dirs(os.path.dirname(output_filename))
            # New file:
            if inpmod.output_f_codec:
                output_file = io.open(output_filename, "wt",
                                      encoding=inpmod.output_f_codec)
            else:   # Default encoding for locale
                output_file = io.open(output_filename, "wt")
            this_is_first_iter = False
        else:
            tmp_str = s_fn.substitute(iter_dict)
            if tmp_str in ofn_list:   # Already wrote a file of this name
                if output_filename != tmp_str: # Must reopen file for appending
                    output_file.close()
                    output_filename = tmp_str
                    if inpmod.output_f_codec:
                        output_file = io.open(output_filename, "a+t",
                                              encoding=inpmod.output_f_codec)
                    else:   # Default encoding for locale
                        output_file = io.open(output_filename, "a+t")
            else:   # Write as a new file
                output_file.close()
                output_filename = tmp_str
                ofn_list.append(output_filename)   # Add filename to checklist
                create_any_needed_dirs(os.path.dirname(output_filename))
                if inpmod.output_f_codec:
                    output_file = io.open(output_filename, "wt",
                                          encoding=inpmod.output_f_codec)
                else:   # Default encoding for locale
                    output_file = io.open(output_filename, "wt")
        print(s.substitute(iter_dict)+"\n", file=output_file)

    output_file.close()   # Close final output file

    # Remove temporary files if necessary:
    if need_to_delete_tmp_mod:  # Remove temporary .py and bytecode files
        try:
            os.remove(tmp_mod_fn_base+".py")
        except: pass   # Ignore exceptions
        try:
            os.remove(tmp_mod_fn_base+".pyc")
        except: pass   # Ignore exceptions
    else:               # Only remove the Python module bytecode file
        try:
            os.remove(input_fn_prefix_full+".pyc")
        except: pass   # Ignore exceptions
    try:   # Remove any bytecode cache
        shutil.rmtree(os.path.join(input_fn_prefix_dir, "__pycache__"))
    except: pass   # Ignore exceptions

    # Print list of files that were generated:
    for ofn in ofn_list:
        print(ofn)

    # "Unload" custom module (otherwise, a module with the same name will not
    #  load):
    del sys.modules[mod_to_import]
#--------------------------------------------------------------------------


#------------------------------------------------------------------------------
def usage():
    usage_msg = """
usage: process_template_file.py [-h]  input_filename_prefix

Write files based on a template file and a substitution keys/values file.  This
is a Python script intended to simplify/automate the creation of one or more
files containing one or more text blocks that are sufficiently similar
("templatable") -- all from two minimalistic input files.  This program requires
either Python version 2.6+ or 3.2+, and is importable as a Python module.

mandatory argument:
  input_filename_prefix   The prefix (including path, if needed) of the two
                           input files ({}.py and {}.template)

optional arguments:
  -h, --help   show this help message and exit
"""
    print(usage_msg)
#------------------------------------------------------------------------------


if __name__ == "__main__":
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
    if (len(args) > 1 or len(args) < 1):
        print("An incorrect number of arguments were provided.")
        usage()
        sys.exit(1)
    input_filename_prefix = args[0]

    # Run driver:
    process_template(input_filename_prefix)
