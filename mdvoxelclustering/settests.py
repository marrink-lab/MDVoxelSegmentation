#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:03 2019

@author: albert
"""
from collections import deque
import numpy as np
import matplotlib.pyplot as plt
import pickle

identity_threshold = 0.5
cluster_mutations = {}


def getlargestset(keylist,setdict):
    # Finds the largest set from setdict with the keys in keylist
    largestkey = 0
    setsize = 0
    for key in keylist:
        if len(setdict[key]) > setsize:
            largestkey = key
            setsize = len(setdict[key])
    return largestkey

def add_cluster_change(time,change):
    """ Adds a cluster change to the timeline of changes in cluster_mutations. If this time point already exists, it adds the change to the list of changes in the time point."""
    if time in cluster_mutations:
        cluster_mutations[time].append(change)
    else:
        cluster_mutations[time] = [change]


def get_match(target_key, matchlist):
    """
    Returns the key of the cluster stored in matchlist and removes that value from the dictionary
    """
    for key, matched_clusters in matchlist.items():
        if target_key in matched_clusters:
            matched_clusters.remove(target_key)
    return key

def compare_clusters (clustersold,clustersnew,deprecated_clusters,clusternumber,frame_number):
    """
    Compares the sorterd clusters in clusters old to the new clusters in clustersnew. Clusters are segmentend by assiging each new cluster to an old cluster, based on the highest number of common elements.
    New clusters are created if there is no fit and merges between clusters are detected. The cluster with key == 0 is seen as the void and is not considered in clustering.

    Returns
    ----------
    This function returns the sorted clusters, the history of previous existing clusters and the index number of clusters.



    """
    output = {}
    matchlist = {k: [] for k in range(len(clustersnew))}
    for key,value in clustersold.items():
        # The cluster 0 is autmaticly matched to the new 0 as it is the void
        if key == 0:
            matchlist[0] = [key]
            continue
        # The likeness score range from 0 to 1 and determines how much 1 cluster is like the other
        likeness_score=index=0
        #assert 0 in matchlist[0], "Zero is not zero..."
        for idx, clusternew in enumerate(clustersnew):
            # Compare each cluster to the new clusters, assign tempscore
            # if the likeness is higher than the previous tempscore
            temp_score = (len(value.intersection(clusternew))
                /max(len(value),len(clusternew),1) )
            #if the new caculated score is higher, this becomes the score to beat
            if (temp_score > likeness_score) and (temp_score > identity_threshold):
                likeness_score = temp_score
                index = idx
        # Add the most likely combination to the list of matches
        matchlist[index].append(key)

    # Check for newly created clusters
    for key,value in matchlist.items():
        #print(str(matchlist) + "   " + str(output))
        # if the key has been removed, do not consider it
        # if key != 0 :
        #     assert 0 not in matchlist[key], "Zero is not being a hero..."
        if  "Blocked" in value:
            continue
        # if exactly 1 old cluster has matched to new cluster, write it to the output file    
        if len(value) == 1:
            output[value[0]] = clustersnew[key]
        # if more then one old clusters matched to this new cluster, there is a merge
        if len(value) > 1:
        # check which set is largest, this becomes the new cluster
            if 0 in value:
                largestset = 0
            else:
                largestset = getlargestset(value,clustersold)
            output[largestset] = clustersnew[key]
            value.remove(largestset)
        # all other clusters are written to depricated clusters 
            for item in value:
                 deprecated_clusters[item] = clustersold[item]
                 print("cluster " + str(item) + " was merged into cluster " + str(largestset) + " at frame " + str(frame_number) + "  1")  
                 # write the merge to the log 
                 add_cluster_change(frame_number,(item,largestset,"m"))
            value.append(largestset)
        # If no old cluster is matched to the new cluster, check for a split
        if len(value) == 0:
        # find the  old cluster from which this new cluster has split
            likeness_score_bestmatch=index_bestmatch=0
            for key2,value2 in clustersold.items():
                temp_score = (len(value2.intersection(clustersnew[key]))
                    /max(len(value2),len(clustersnew[key]),1))
                if temp_score > likeness_score_bestmatch:
                    likeness_score_bestmatch = temp_score
                    index_bestmatch = key2
            likeness_scoredep=indexdep=0
        # Look through the list of deprecated clusters and see if there is a match
            for key2,value2 in deprecated_clusters.items():
                temp_score = (len(value2.intersection(clustersnew[key]))
                    /max(len(value2),len(clustersnew[key]),1) )
                if temp_score > likeness_scoredep  and\
                temp_score > identity_threshold :
                    likeness_scoredep = temp_score
                    indexdep = key2
        # If there is a match, recreate the deprecated cluster and fill it with the new cluster
            if indexdep != 0:
                del deprecated_clusters[indexdep]
                output[indexdep] = clustersnew[key]
                print ("cluster " + str(indexdep) + " was recreated from cluster " + str(index_bestmatch) + " at frame " + str(frame_number))
                add_cluster_change(frame_number,(indexdep,index_bestmatch,"rc"))
            else: 
        # If there is no match, check if the new cluster occupying the parent old cluster has a match
                if index_bestmatch != 0:
                    indexdep_bestmatch = likeness_scoredep_bestmatch = 0
                    fight_cluster = None
                    for key3, matched_clusters in matchlist.items():
                        if index_bestmatch in matched_clusters:
                            fight_cluster = key3

                    if fight_cluster is None:
                        newkey = clusternumber + 1
                        clusternumber = newkey
                        output[newkey] = clustersnew[key]
                        print ("cluster " + str(newkey) + " was created from cluster " + str(index_bestmatch) + " at frame " + str(frame_number))
                        add_cluster_change(frame_number,(newkey,index_bestmatch,"c"))
                    else:
                        for key2,value2 in deprecated_clusters.items():
                            temp_score = (len(value2.intersection(clustersnew[fight_cluster]))
                                /max(len(value2),len(clustersnew[fight_cluster]),1) )
                            if temp_score > likeness_scoredep  and\
                            temp_score > identity_threshold :
                                likeness_scoredep_bestmatch = temp_score
                                indexdep_bestmatch = key2
                                print(key, value, indexdep_bestmatch, likeness_scoredep_bestmatch,index_bestmatch)
                        # if there is a match, match the new cluster occupying the parent cluster to this new match
                        if indexdep_bestmatch != 0:
                            # delete the new cluster in history
                            del deprecated_clusters[indexdep_bestmatch]

                                # print("key " + str(key3) + " matched_clusters " + str(matched_clusters) + " index_bestmatch " + str(index_bestmatch))
                            matchlist[fight_cluster].remove(index_bestmatch)
                                # and remove this cluster from the matching list
                            if matchlist[fight_cluster] == []:
                                matchlist[fight_cluster] = ["Blocked"]
                                    #print("key " + str(key3) + " matched_clusters " + str(matched_clusters) + " index_bestmatch " + str(index_bestmatch))
                                output[indexdep_bestmatch] = clustersnew[fight_cluster]
                                print ("cluster " + str(indexdep_bestmatch) + " was recreated specially for you from cluster " + str(index_bestmatch) + " at frame " + str(frame_number))
                                add_cluster_change(frame_number,(indexdep_bestmatch,index_bestmatch,"rc"))
                                # add the current viewed new cluster as match to it's best matching old cluster
                                output[index_bestmatch] = clustersnew[key]
                        else:

                            # If occupying new cluster does not have a match, create a new cluster
                            newkey = clusternumber + 1
                            clusternumber = newkey
                            output[newkey] = clustersnew[key]
                            print ("cluster " + str(newkey) + " was created from cluster " + str(index_bestmatch) + " at frame " + str(frame_number))
                            add_cluster_change(frame_number,(newkey,index_bestmatch,"c"))
                else:

                    # If occupying new cluster does not have a match, create a new cluster
                    newkey = clusternumber + 1
                    clusternumber = newkey
                    output[newkey] = clustersnew[key]
                    print ("cluster " + str(newkey) + " was created from cluster " + str(index_bestmatch) + " at frame " + str(frame_number))
                    add_cluster_change(frame_number,(newkey,index_bestmatch,"c"))
               
    return output, clusternumber, deprecated_clusters  

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
    return data,plotinfo, clustercount


def compare_2masks(mask1, mask2, selection):
    condition1 = mask1[mask1 == selection]
    condition2 = mask2[mask2 == selection]
    comparision = np.all(condition1 == condition2)
    return comparision


def main():
    data = np.load('clusters.npy')
    data_old = np.copy(data)
    data, plotinfo, n_clusters = sort_clusters(data)
    for index, _ in enumerate(data):
        if not compare_2masks(data[index],data_old[index], 0):
            print(index)



    np.save('clusters_ordered', data)
    
    with open('visualization_data.pickle', 'wb') as handle:
        pickle.dump(plotinfo, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    with open('cluster_mutations.pickle', 'wb') as handle:
        pickle.dump(cluster_mutations, handle, protocol=pickle.HIGHEST_PROTOCOL)
    with open('n_clusters.pickle', 'wb') as handle:
        pickle.dump(n_clusters, handle, protocol=pickle.HIGHEST_PROTOCOL)

    time = range(len(plotinfo))
    ordered_clusters = plotinfo.values()
    print(time)
    #print(ordered_clusters) 

if __name__=='__main__':
    main()