import mdvoxelclustering as mdv
from mdvoxelclustering import settests as sets
from collections import deque

def test_standard(olddict1, frame20):
    result,_,_ = sets.compare_clusters(olddict1,frame20,{},8,1)
    test0 = {0: set(), 1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 4: {'d9', 'd5', 'd6', 'd8', 'd10', 'd2', 'd7', 'd1', 'd3', 'd4'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'}, 6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'}, 7: {'g3', 'g1', 'g2'}}
    assert(result == test0), "Standard merging failed"    


  
def test_newcreate(olddict1,frame21):    
    result,resultnumber,_ = sets.compare_clusters(olddict1,frame21,{},8,1)
    test1 = {0: set(), 1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 4: {'d9', 'd5', 'd6', 'd8', 'd10', 'd2', 'd7', 'd1', 'd3', 'd4'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'}, 6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'}, 7: {'g3', 'g1', 'g2'}, 9: {'b8', 'b9', 'b7', 'b10'}}
    assert (result == test1 and resultnumber == 9), "Create new cluster failed" 
  
def test_recreate(olddict1,frame21,deprecated_dict):
    
    result,resultnumber,resultdep =sets.compare_clusters(olddict1,frame21,deprecated_dict,8,1)

    test4 = {0: set(), 1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 4: {'d9', 'd5', 'd6', 'd8', 'd10', 'd2', 'd7', 'd1', 'd3', 'd4'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'}, 6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'}, 7: {'g3', 'g1', 'g2'}, 8: {'b8', 'b9', 'b7', 'b10'}}

    assert (result == test4 and resultnumber == 8 and resultdep == {}), "Recreate old cluster failed"
   
def test_merge(olddict2,frame22):
    result,_,_ = sets.compare_clusters(olddict2,frame22,{},8,1)
    test2 = {0: set(), 1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'}, 6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'},  7: {'g3', 'g1', 'g2'}, 8: {'d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'}}
    assert(result == test2), "Merging clusters failed" 
     
    
def test_deprecation(olddict3,frame23,test3deprecated):
    result,_,resultdep = sets.compare_clusters(olddict3,frame23,{},8,1)
    test3 = {0: set(), 1: {'a8', 'a5', 'a6', 'a4', 'a9', 'a2', 'a10', 'a7', 'a3', 'a1'}, 2: {'b1', 'b4', 'b3', 'b5', 'b2'}, 3: {'c4', 'c5', 'c6', 'c1', 'c2', 'c3'}, 4: {'d9', 'd5', 'd6', 'd8', 'd10', 'd2', 'd7', 'd1', 'd3', 'd4'}, 5: {'e7', 'e5', 'e6', 'e1', 'e2', 'e3', 'e8', 'e4'},  6: {'f2', 'f1', 'f3', 'f5', 'f6', 'f7', 'f9', 'f10', 'f8', 'f4'}, 7: {'g3', 'g1', 'g2'}}
    assert(result == test3 and resultdep == test3deprecated), "Deprecating cluster failed"

def test_bilayer(bilayerdict,frame_bilayer,deprecated_bilayer_dict):
    result, resultnumber, resultdep = sets.compare_clusters(bilayerdict,frame_bilayer,deprecated_bilayer_dict,2,1)
    test_bilayer ={0: set(), 1: {'a1','a2','a3','a4','a5','a6','a7','a8','a9'}, 2: {'b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','a10'}}
    assert(result == test_bilayer and resultdep == {} and resultnumber == 2 ), " Bilayer recreationg failed"


def test_zero(zerodict,frame_zero):
    result,_,_ = sets.compare_clusters(zerodict,frame_zero,{},1,1)
    test_zero = {0 : {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}}
    assert result == test_zero, " Zero is no longer zero"
