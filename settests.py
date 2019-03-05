#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:17:03 2019

@author: bart
"""
from collections import deque


cluster1 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
cluster2 = {'b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}
cluster3 = {'c1','c2','c3','c4','c5','c6','c7','c8','c9'}
cluster4 = {'d1','d2','d3'}
cluster5 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
cluster6 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
cluster7 = {'g1','g2','g3','g4'}
cluster8 = {'d7','d8','d9','d10','d4','d5'}

olddict1 = {
        1 : cluster1,
        2 : cluster2,
        3 : cluster3,
        4 : cluster4,
        5 : cluster5,
        6 : cluster6,
        7 : cluster7}


olddict2 = {
        1 : cluster1,
        2 : cluster2,
        3 : cluster3,
        4 : cluster4,
        5 : cluster5,
        6 : cluster6,
        7 : cluster7,
        8 : cluster8   
        }


olddict3 = {
        1 : cluster1,
        2 : cluster2,
        3 : cluster3,
        4 : cluster8,
        5 : cluster5,
        6 : cluster6,
        7 : cluster7,
        8 : cluster4   
        }

cluster11 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
cluster12 = {'b1','b2','b3','b4','b5'}
cluster13 = {'c1','c2','c3','c4','c5','c6'}
cluster14 = {'d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'}
cluster15 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
cluster16 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
cluster17 = {'g1','g2','g3'}
cluster18 = {'b7','b8','b9','b10'}


frame21 = deque([cluster11,cluster12,cluster13,cluster14,cluster15,cluster16,cluster17,cluster18])
frame22 = deque([cluster11,cluster12,cluster13,cluster14,cluster15,cluster16,cluster17])
frame23 = deque([cluster11,cluster12,cluster13,cluster14,cluster15,cluster16,cluster17])

test1 = {1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 4: {'d9', 'd5', 'd6', 'd8', 'd10', 'd2', 'd7', 'd1', 'd3', 'd4'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'}, 6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'}, 7: {'g3', 'g1', 'g2'}, 8: {'b8', 'b9', 'b7', 'b10'}}
test2 = {1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'}, 6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'},  7: {'g3', 'g1', 'g2'}, 8: {'d9', 'd5', 'd8', 'd10', 'd7', 'd4'}}
test3 = {1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 4: {'d9', 'd5', 'd6', 'd8', 'd10', 'd2', 'd7', 'd1', 'd3', 'd4'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'},  6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'}, 7: {'g3', 'g1', 'g2'}}

def compare_clusters (clustersold,clustersnew):
    output = {}
    
    
    for key,value in clustersold.items():
        likeness_score=index=0
        if not clustersnew:
             for key_out,value_out in output.items():
                 temp_score = len(value.intersection(value_out))/max(len(value),len(value_out)) 
                 if temp_score > likeness_score:
                        likeness_score = temp_score
                        index = key_out
             if len(value) > len(value_out):
                 output[key] = value 
                 del output[index]
                 print("cluster " + str(index) + " was merged into cluster " + str(key))
             else:
                 print("cluster " + str(key) + " was merged into cluster " + str(index))   
        else:
            for idx, clusternew in enumerate(clustersnew):
                temp_score = len(value.intersection(clusternew))/max(len(value),len(clusternew)) 
                if temp_score > likeness_score:
                    likeness_score = temp_score
                    index = idx
    
            clustersnew.rotate(-index)
            output[key] = clustersnew.popleft()
    while clustersnew:
        likeness_score=index=0
        for key,value in clustersold.items():
            temp_score = len(value.intersection(clusternew))/max(len(value),len(clusternew)) 
            if temp_score > likeness_score:
                likeness_score = temp_score
                index = key
        newkey = max(clustersold,key=int) + 1
        output[newkey ] = clustersnew.popleft()
        print ("cluster " + str(newkey) + " was created from cluster " + str(index))
    return output

    
                
                  
                
            
            
if  compare_clusters(olddict1,frame21)== test1 :
    print("test1 ok :)")

if  compare_clusters(olddict2,frame22)== test2 :
    print("test2 ok :)")
    
if  compare_clusters(olddict3,frame23) == test3:
   print("test3 ok :)")
        
        
    
    

