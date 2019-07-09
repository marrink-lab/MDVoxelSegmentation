#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:03 2019

@author: bart
"""
from collections import deque
import numpy as np
import matplotlib.pyplot as plt
import pickle

identity_threshold = 0.0000001
cluster_mutations = {}


def getlargestset(keylist,setdict):
    largestkey = 0
    setsize = 0
    for key in keylist:
        if len(setdict[key]) > setsize:
            largestkey = key
            setsize = len(setdict[key])
    return largestkey

def add_cluster_change(time,change):
   # print(cluster_mutations)
    if time in cluster_mutations:
        cluster_mutations[time].append(change)
    else:
        cluster_mutations[time] = [change]


def compare_clusters (clustersold,clustersnew,deprecated_clusters,clusternumber,frame_number):
    output = {}
    mergelist = {k: [] for k in range(len(clustersnew))}
    for key,value in clustersold.items():
        if key == 0:
            mergelist[0] = [key]
            continue
            
        likeness_score=index=0
        for idx, clusternew in enumerate(clustersnew):
            temp_score = len(value.intersection(clusternew))/max(len(value),len(clusternew),1) 
            if temp_score > likeness_score:
                likeness_score = temp_score
                index = idx
        mergelist[index].append(key)


        
    for key,value in mergelist.items():
        if len(value) == 1:
            output[value[0]] = clustersnew[key]
        if len(value) > 1:
            largestset = getlargestset(value,clustersold)
            output[largestset] = clustersnew[key]
            value.remove(largestset)
            for item in value:
                 deprecated_clusters[item] = clustersold[item]
                 print("cluster " + str(item) + " was merged into cluster " + str(largestset) + " at frame " + str(frame_number) + "  1")   
                 add_cluster_change(frame_number,(item,largestset,"m"))
        if len(value) == 0:
            likeness_score=index=0
            for key2,value2 in clustersold.items():
                temp_score = len(value2.intersection(clustersnew[key]))/max(len(value2),len(clustersnew[key]),1) 
                if temp_score > likeness_score:
                    likeness_score = temp_score
                    index = key2
            likeness_scoredep=indexdep=0
            for key2,value2 in deprecated_clusters.items():
                temp_score = len(value2.intersection(clustersnew[key]))/max(len(value2),len(clustersnew[key]),1) 
                if temp_score > likeness_scoredep  and temp_score > identity_threshold :
                    likeness_scoredep = temp_score
                    indexdep = key2
            if indexdep != 0:
                del deprecated_clusters[indexdep]
                output[indexdep] = clustersnew[key]
                print ("cluster " + str(indexdep) + " was recreated from cluster " + str(index) + " at frame " + str(frame_number))
                add_cluster_change(frame_number,(indexdep,index,"rc"))
            else: 
                newkey = clusternumber + 1
                clusternumber = newkey
                output[newkey] = clustersnew[key]
                print ("cluster " + str(newkey) + " was created from cluster " + str(index) + " at frame " + str(frame_number))
                add_cluster_change(frame_number,(newkey,index,"c"))
               
    return output,clusternumber, deprecated_clusters  

def sort_clusters(data):
    clustercount = 0
    olddict = {0 : set()}
    deprecated_clusters = {}
    plotinfo = {0 : []}
    for idx, frame in enumerate(data):
        newcluster = deque()
        plotvariables = []
        for cluster in set(frame):     
            newcluster.append(set(np.where(frame == cluster)[0]))  
        olddict, clustercount, deprecated_clusters = compare_clusters(olddict, newcluster, deprecated_clusters,clustercount,idx)
        for key, value in olddict.items():
            for index in value:
                frame[index] = key
            plotvariables.append((key,len(value)))
        plotinfo[idx] = plotvariables
    return data,plotinfo


def main():
    data = np.load('clusters.npy')
    
    data, plotinfo = sort_clusters(data)
    
    np.save('clusters_ordered', data)
    
    with open('visualization_data.pickle', 'wb') as handle:
        pickle.dump(plotinfo, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    with open('cluster_mutations.pickle', 'wb') as handle:
        pickle.dump(cluster_mutations, handle, protocol=pickle.HIGHEST_PROTOCOL)

    time = range(len(plotinfo))
    ordered_clusters = plotinfo.values()
    print(time)
    print(ordered_clusters) 

if __name__=='__main__':
    main()