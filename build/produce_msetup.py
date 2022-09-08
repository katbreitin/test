#!/usr/bin/env python
"""
Produces a machine setup file from information in a configuration file.

This program requires either Python version 2.6+ or 3.2+, and is importable as
a python module.
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
import os
import getopt

if sys.version_info[0] >= 3:
    import configparser as cp
else:
    import ConfigParser as cp


#------------------------------------------------------------------------------
def main():
    """Driver."""

    machine_setup_fpath_parts = \
             "env_settings/autogenerated_from_cfg/bec-machine_setup".split('/')
    machine_setup_fpath = os.path.join(*machine_setup_fpath_parts)

    msetup_cfg_fpath_parts = \
             "env_settings/user_change_me.cfg".split('/')
    msetup_cfg_fpath = os.path.join(*msetup_cfg_fpath_parts)
    
    cfg = cp.RawConfigParser()
    cfg.read(msetup_cfg_fpath)

    this_section = "Modules To Load"
    precursor_script = cfg.get(this_section,
                            "Shell script to load BEFORE (un)loading modules")
    purge_modules_first = cfg.getboolean(this_section,
                                         "Purge any currently-loaded modules?")
    module_commands = cfg.get(this_section,
                              "Modules to (un)load (in order)")
    postcursor_script = cfg.get(this_section,
                              "Shell script to load AFTER (un)loading modules")

    text_to_write = "##!/bin/sh\n"+"#\n"+"# auto-generated from "+ \
                    msetup_cfg_fpath+"\n"
    text_to_write +="""#
# Intended to contain instructions to properly set up the computing environment
# (especially on large managed clusters/supercomputers).  These will be
# executed during code compilation and/or at runtime.
#
# This code has been written using "heirloom" sh-syntax, which should
#  allow the greatest (easily) possible portability among "modern" computing
#  systems.  All modifications to this code should therefore employ only
#  "heirloom" sh-syntax (e.g., *NO* syntax, assumptions, or structures
#  specific to bash, ksh, tcsh, zsh, or any nonstandard "heirloom" sh).\n
"""

    tmp_str = precursor_script.strip()
    if tmp_str:
        text_to_write += \
              "# First, load any precursor script:\n. "+tmp_str+"\n\n"

    if purge_modules_first:
        text_to_write += \
              "# Remove all currently loaded modules:\n"+"module purge\n\n"

    text_to_write += \
          "# Load and/or unload modules:\n"+module_commands.strip()+"\n\n"

    tmp_str = postcursor_script.strip()
    if tmp_str:
        text_to_write += \
              "# Last, load any postcursor script:\n ."+tmp_str+"\n\n"

    with open(machine_setup_fpath, "wt") as text_f:
        text_f.write(text_to_write)
            
    sys.exit(0)
    

#------------------------------------------------------------------------------
def usage():
    """Displays a usage message."""
    usage_msg = """
usage: produce_msetup.py [-h]

Produces a machine setup file from information in a configuration file.  This
program requires either Python version 2.6+ or 3.2+, and is importable as a
Python module.

optional argument:
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
    if (len(args) > 0):
        print("An incorrect number of arguments were provided.")
        usage()
        sys.exit(1)

    # Perform action(s):
    main()
