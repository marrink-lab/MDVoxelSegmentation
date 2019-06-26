import os
import pickle
import numpy

def assignsets(clusters):
    frame_sets = []
    for idx in range(len(clusters)):
        setassignments = {}
        for idx in range(len(clusters[0])):
            print(idx)
            setassignments[idx] = set(clusters[0][idx].atoms.ids)
        frame_sets.append(setassignments)
    with open('frame_sets.pkl', 'wb') as f:
        pickle.dump(frame_sets, f)
    with open('frame_sets.pkl', 'rb') as f:
        testlist = pickle.load(f)
    print(testlist)
    

def ordersets():
    with open('../../frame_sets.pkl', 'rb') as f:
        testlist = pickle.load(f)
    print(type(testlist))
    print(len(testlist))

ordersets()