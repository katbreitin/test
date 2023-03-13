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
This is a Python-syntax configuration file for template substitutions.

The Python-syntax blocks below define substitution keywords and their values
(over multiple iterations) -- for use with the associated template (*.template)
and 'process_template_file.py' (which uses the standard string.Template class),
producing the final associated file(s).
"""

# Output file text encoding -- set to
#  any valid Python codec name -OR- "" (uses default encoding for locale)
#   NOTE: for Fortran source code, set to "ascii"
output_f_codec = "ascii"

# A dictionary of keys with only scalar values that will be used to substitute
# into the given values of 'sdl_string' and 'output_fn_template'.

defines_dict = dict( \
                   use_kd_mod = 'use univ_kind_defs_mod, only:'
                   )

# A quoted list of dictionaries which have keys with lists as their values.  It
# must remain quoted (rendering it as a long string) to permit substitution of
# any tokens within it using the 'defines_dict' dictionary.  Each of the value
# elements of each dictionary key implies an iteration.  Multiple (up to 3)
# dictionaries can be used to produce even more complex output.  A single
# dictionary with key values having 3 elements will generate 3 iterations of
# template substitution.  Two dictionaries with 3 and 4 elements, respectively,
# will generate 12 iterations of template substitution.
#
# Every key list below (i.e., within the '[]' after each '=') must have the same
# number of elements.  Only strings allowed (e.g., no unquoted numeric values).

sdl_string = '''[ \
    dict(
        kind_idx_n = ['f4', 'f8'],
        kind_idx_use = ['f4', 'f8'],
        kind_idx_dec = ['f4', 'f8'],
        kind_desc = ['4-byte floating-point', '8-byte floating-point'],
        kind_mod = ['${use_kd_mod}', '${use_kd_mod}'],
        ), \
    dict(
        dim_n = ['1d', '2d', '3d', '4d', '5d', '6d', '7d'],
        dim_desc = ['1-dimensional', '2-dimensional', '3-dimensional',
                    '4-dimensional', '5-dimensional', '6-dimensional',
                    '7-dimensional'],
        dim_alloc = ['(size(A, 1))', '(size(A, 1), size(A, 2))',
                       '(size(A, 1), size(A, 2), size(A, 3))',
                       '(size(A, 1), size(A, 2), size(A, 3), size(A, 4))',
                       '(size(A, 1), size(A, 2), size(A, 3), size(A, 4), size(a, 5))',
                       '(size(A, 1), size(A, 2), size(A, 3), size(A, 4), size(a, 5), size(A, 6))',
                       '(size(A, 1), size(A, 2), size(A, 3), size(A, 4), size(a, 5), size(A, 6), size(A, 7))'],
        dim_rank_dec = ['(:)', '(:,:)', '(:,:,:)', '(:,:,:,:)', '(:,:,:,:,:)',
                        '(:,:,:,:,:,:)', '(:,:,:,:,:,:,:)']
        ) ]'''

# A quoted string denoting the output filename (with path, if desired).  This
# will undergo substitution of any tokens within it using both the
# 'defines_dict' dictionary and a single iteration of the 'sdl_string'
# dictionary.

output_fn_template = \
        '''univ_fp_comparison.f90.inc02'''
