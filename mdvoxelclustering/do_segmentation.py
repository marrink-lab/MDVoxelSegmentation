#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 16:17:19 2020

@author: bart
Running this file is all that is needed to perform a basic leaflet
segmentation for Martini lipids.
"""
import time
import multiprocessing as mp
import MDAnalysis as mda
from mdvoxelclustering import argparser
from mdvoxelclustering import leaflets

def main():
    # Reading in the terminal commands/input files
    args = argparser.main()
    
    # Processing input
    print('Reading trajectory...')
    universe = mda.Universe(args.reference, args.trajectory)
    
    # Setting the amount of threads to the lowest value [detected, trajlen]
    if len(universe.trajectory) < args.threads:
        args.threads = len(universe.trajectory)
    
    # Starting the multithreaded leaflet segmentation
    print('Actual segmentations..')
    start = time.time()
    
    pool = mp.Pool(args.threads)
    segments = pool.map(leaflets.mf_leaflets_threaded, range(args.threads))
    segments = np.asarray(segments)
    segments = np.concatenate(segments, axis=0)
    pool.close()
    pool.join()
    
    print(segments, len(segments))
    
    print('Segmentation took: {}.'.format(time.time() - start))
    
    #TODO Writing the ouput at once, this should become a per frame append!
    np.save(args.output, segments.astype(args.bit_size))
    
    if args.plotting:
        plot_clusters(universe, segments, args.stride, args.reduce_points, min_size = 1,
                      start_frame = args.begin, stop_frame = args.end)
    
if __name__=='__main__':
    main()