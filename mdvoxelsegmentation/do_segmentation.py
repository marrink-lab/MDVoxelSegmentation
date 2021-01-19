#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:17:19 2020

@author: bart
Running this file is all that is needed to perform a basic leaflet
segmentation for Martini lipids.
"""
import time
import os
from shutil import copy
import numpy as np
import multiprocessing as mp
import MDAnalysis as mda
from functools import partial
from mdvoxelsegmentation import argparser
from mdvoxelsegmentation import leaflets
from mdvoxelsegmentation import settests
from mdvoxelsegmentation import plotting
     

def main():
    # Reading in the terminal commands/input files
    args = argparser.main()
    
    # Reading trajectory
    print('\nReading trajectory...')
    universe = mda.Universe(args.reference, args.trajectory)
    if args.end == None:
        args.end = len(universe.trajectory)
    
    # Finding the total amount of frames
    args.frames = int((args.end - args.begin) / args.stride)
    
    # Setting the amount of threads to the lowest value [detected, trajlen]
    if args.frames < args.threads:
        args.threads = args.frames
            
    
    # Starting the multithreaded leaflet segmentation
    #print('Actual segmentation..')
    start = time.time()
    
    pool = mp.Pool(args.threads)
    segments = pool.map(partial(leaflets.mf_leaflets_threaded, args=args), range(args.threads))
    segments = np.asarray(segments)
    segments = np.concatenate(segments, axis=0)
    pool.close()
    pool.join()
    
    #print(segments, len(segments))
    
    print()
    print('Segmentation took: {:.2f} seconds.\n'.format(time.time() - start))
    print('The segmentation array was written to [{}.npy].'.format(args.output))
    
    #TODO Writing the ouput at once, this should become a per frame append!
    np.save(args.output, segments.astype(args.bit_size))

    #TODO Segment the output over time (handle the user input in the argparser)
    settests.main()

    # Copies the default VMD related files for easy visualization
    file_path = "{}/example_inputs".format('/'.join(__file__.split('/')[:-2]))
    input_PY = 'vmd_clusters_visualization.py' 
    input_VMD = 'vmd_clusters_visualization.vmd'
    cwd = os.getcwd()
    for path in [input_PY, input_VMD]:
        copy('{}/{}'.format(file_path, path), cwd)
    
    # Replace the default filenames with the active filenames
    args_dict = vars(args)
    with open('vmd_clusters_visualization.py', 'r') as f:
        sample1 = f.read().replace('your_ref.gro', args_dict['reference'])
        sample2 = sample1.replace('your_ref.xtc', args_dict['trajectory'])
    with open('vmd_clusters_visualization.py', 'w') as f:
        f.write(sample2)

    # Copy the plotting script into the folder where the analysis is performed
    #  and runs the default plotting procedure.
    file_path = "{}".format('/'.join(__file__.split('/')[:-1]))
    input_PY = 'plotting.py'
    copy('{}/{}'.format(file_path, input_PY), cwd)
    plotting.main()
    
 
if __name__=='__main__':
    main()
