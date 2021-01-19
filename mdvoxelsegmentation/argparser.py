#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:38:35 2020

@author: bart
This should contain all functionality regarding input. Both an input file 
and the terminal argparsing is included. If any functionality would require
user input, that should be handled here.

"""
import argparse
import psutil
import os

def write_default_input(path):
    """
    Writes a default selections input file for most standard Martini 
    molecules. Common lipids in the lipidome as well as basic protein and DNA 
    are included.

    Parameters
    ----------
    path : str
        The path to the intended location.

    Returns
    -------
    None.
    """
    message = (
"""# The selections should be in the MDAnalysis selection format. They need to be one line, directly
#  after their header, anything else is ignored.

[martini_heads]
(name PO4 NC3 CNO NH3 TAP GL1 GL2 AM1 AM2 GM1 GM2 GM3 GM4 GM5 GM6 GM7 GM8 GM9 GM10 GM11 GM12 GM13 GM14 GM15 GM16 GM17) or (resname CHOL and name ROH) or (resname PAPI PIPI POP1 POP2 POP3 POPI PUPI and name C1 C2 C3 P1 P2 P3)

[martini_tails]
(name D6A D6B C6A C6B D5A D5B C5A C5B D4A D4B C4A C4B D3A D3B C3A C3B) or (resname CHOL XHOL and name C1)

[martini_linkers]
(name PO4 GL1 GL2 D1A D1B C1A C1B AM1 AM2 GM1 GM2 COO COOH) or (resname CHOL and name ROH) or (resname PAPI PIPI POP1 POP2 POP3 POPI PUPI and name C1 C2 C3 P1 P2 P3) or (name BB)

[martini_proteins]
(name BB SC1 SC2 SC3 SC4 SC5)

[martini_dna]
(name BB1 BB2 SC1 SC2 SC3 SC4)

[none]
False
"""
)
    with open(path, 'w') as f:
        f.write(message)
    return

def read_arguments():
    """
    Parses the input arguments from the command line.

    Returns
    -------
    args = NameSpace

    """
    auto_threads = psutil.cpu_count()
    
    # Generating the argparser object.
    parser = argparse.ArgumentParser(add_help=False)
    required_grp = parser.add_argument_group(title='required arguments')
    optional_grp = parser.add_argument_group(title='optional arguments')
    # REQUIRED
    required_grp.add_argument(
        '-f', '--reference', nargs='?', required=True, type=str,
        help='an MDAnalysis coordinate file including atom names, resids and positions (e.g. GRO, TPR or PDB)',
        )
    required_grp.add_argument(
        '-x', '--trajectory', nargs='?', required=True, type=str,
        help='an MDAnalysis compatible trajectory file (e.g. XTC)',
        )
    
    # OPIONAL
    #TODO Make it so a selection input is read.
    optional_grp.add_argument(
        '-si', '--selections_input', nargs='?', default='selections.inp', 
        type = str,
        help='the selection input file (default=selections.inp)',
        )    
    optional_grp.add_argument(
        '-hg', '--headgroups', nargs='?', default='martini_heads', type=str,
        help='the selection name in [selections.inp] used for the headgroups (default=martini_heads)',
        )
    optional_grp.add_argument(
        '-tg', '--tailgroups', nargs='?', default='martini_tails', type=str,
        help='the selection name in [selections.inp] used for the tails (default=martini_tails)',
        )
    optional_grp.add_argument(
        '-lg', '--linkergroups', nargs='?', default='martini_linkers', 
        type=str,
        help='the selection name in [selections.inp] used for the linkers (default=martini_linkers)',
        )
    optional_grp.add_argument(
        '-eg', '--exclusiongroups', nargs='?', default='martini_proteins',
        type=str,
        help='the selection name in the [selections.inp] used for the exclusions (default=martini_proteins)',
        )
    
    optional_grp.add_argument(
        '-res', '--resolution', nargs='?', default=0.5, type=float,
        help='the binning resolution in the same units as the reference file (default=0.5)',
        )
    optional_grp.add_argument(
        '-hres', '--hyper_resolution', nargs='?', default=0.5, type=float,
        help='blurs the coordinates by a fraction of the bin dimension (default=0.5)',
        )
    # old default is 10
    optional_grp.add_argument(
        '-rd', '--recursion_depth', nargs='?', default=0, type=int,
        help='amount of iterations for forced segmentation (default=10; 0 is off)',
        )
    # old default is 20
    optional_grp.add_argument(
        '-fs', '--force_segmentation', nargs='?', default=0, type=float,
        help='forces segmentation within set radius, the units are the same as in the reference file (default=20)',
        )
    optional_grp.add_argument(
        '-fi', '--force_info', nargs='?', default=False, type=bool,
        help='set force segmentation information printing True/False (default=False)',
        )
    optional_grp.add_argument(
        '-min', '--minimum_size', nargs='?', default=50, type=int,
        help='the minimum size of a segment in the amount of beads (default=50)',
        )
    optional_grp.add_argument(
        '-b', '--begin', nargs='?', default=0, type=int,
        help='set starting frame (default=0)',
        )
    optional_grp.add_argument(
        '-e', '--end', nargs='?', default=None, type=int,
        help='set end frame (default=None)',
        )
    optional_grp.add_argument(
        '-s', '--stride', nargs='?', default=1, type=int,
        help='set stride for reading trajectory (default=1)',
        )
    optional_grp.add_argument(
        '-nt', '--threads', nargs='?', default=auto_threads, type=int,
        help='the maximum number of threads available (detected={})'.format(auto_threads),
        )
    optional_grp.add_argument(
        '-bs', '--bit_size', nargs='?', default='uint32', type=str,
        help='set cluster array bit size (default=uint32)',
        )
    # optional_grp.add_argument(
    #     '-p', '--plotting', nargs='?', default=False, 
    #     help='turn on plots for testing (default=False)',
    #     )
    # optional_grp.add_argument(
    #     '-rp', '--reduce_points', nargs='?', default=1, 
    #     help='set the stride in points for plotting (default=1)',
    #     )
    #TODO Change this so its called segmentations
    optional_grp.add_argument(
        '-o', '--output', nargs='?', default='clusters', type=str,
        help='change segmentation array name (default=clusters), changing this setting breaks the program',
        )
    optional_grp.add_argument(
        '-v', '--verbose', nargs='?', default=False, type=bool,
        help='set verbose True/False (default=False)',
        )
    optional_grp.add_argument(
        '-h', '--help', action="help",
        help='show this help message and exit',
        )
    # parse the arguments into the name space
    args = parser.parse_args()
    return args

def read_selections(args):
    """
    Adds selection queries to the args based on the selection groups.

    Parameters
    ----------
    args : NameSpace
        The terminal arguments.
    
    Returns
    -------
    args : NameSpace
        The MDAnalysis selection syntax for the given selections is added to the NameSpace.

    """
    # Use and write default input if the input file does not exist.
    if not os.path.isfile(args.selections_input):
        write_default_input(args.selections_input)
        
    selection_headers = (args.headgroups, args.tailgroups, 
                         args.linkergroups, args.exclusiongroups)
    selection_header_patterns = ['[{}]'.format(x) for x in selection_headers]
    selection_strings = []
    for selection_header_pattern in selection_header_patterns:
        with open(args.selections_input, 'r') as f:
            for line in f:
                if line.strip() == selection_header_pattern:
                    selection_strings.append(f.readline().strip())
                    break
            else:
                raise IndexError('{} could not be found in {}'.format(selection_header_pattern[1:-1], args.selections_input))
            
    args.headgroups_selection_query = selection_strings[0]
    args.tails_selection_query = selection_strings[1]
    args.linkers_selection_query = selection_strings[2]
    args.exclusions_selection_query = selection_strings[3]
    return args


def main():
    """
    Returns the parsed arguments and write them to a log file.

    Returns
    -------
    args : NameSpace
        The complete input arguments.

    """
    args = read_arguments()
    args = read_selections(args)
    with open('input_arguments.log', 'w') as f:
        f.write(str(args))
    return args

if __name__=='__main__':
    main()
