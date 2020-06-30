#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:17:19 2020

@author: bart
Running this file is all that is needed to perform a basic leaflet
segmentation for Martini lipids.
"""
import time
import numpy as np
import multiprocessing as mp
import MDAnalysis as mda
from functools import partial
from mdvoxelsegmentation import argparser
from mdvoxelsegmentation import leaflets
from mdvoxelsegmentation import settests

     

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
    
if __name__=='__main__':
    main()
