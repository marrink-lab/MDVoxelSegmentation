from collections import deque
import pytest

@pytest.fixture
def olddict1():
    cluster1 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
    cluster2 = {'b1','b2','b3','b4','b5','b6','b7'}
    cluster3 = {'c1','c2','c3','c4','c5','c6','c7','c8','c9'}
    cluster4 = {'d9', 'd5', 'd6', 'd8', 'd10', 'd2', 'd7'}
    cluster5 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster6 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster7 = {'g1','g2','g3','g4'}


    olddict1 = {
            0 : set(),
            1 : cluster1,
            2 : cluster2,
            3 : cluster3,
            4 : cluster4,
            5 : cluster5,
            6 : cluster6,
            7 : cluster7}
    return olddict1

@pytest.fixture    
def olddict2():
    cluster1 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
    cluster2 = {'b1','b2','b3','b4','b5','b6','b7'}
    cluster3 = {'c1','c2','c3','c4','c5','c6','c7','c8','c9'}
    cluster4 = {'d1','d2','d3'}
    cluster5 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster6 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster7 = {'g1','g2','g3','g4'}
    cluster8 = {'d7','d8','d9','d10','d4','d5'}

    olddict2 = {
            0 : set(),
            1 : cluster1,
            2 : cluster2,
            3 : cluster3,
            4 : cluster4,
            5 : cluster5,
            6 : cluster6,
            7 : cluster7,
            8 : cluster8   
            }
    return olddict2

@pytest.fixture  
def olddict3():
    cluster1 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
    cluster2 = {'b1','b2','b3','b4','b5','b6','b7'}
    cluster3 = {'c1','c2','c3','c4','c5','c6','c7','c8','c9'}
    cluster4 = {'d1','d2','d3'}
    cluster5 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster6 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster7 = {'g1','g2','g3','g4'}
    cluster8 = {'d7','d8','d9','d10','d4','d5'}
    olddict3 = {
            0 : set(),
            1 : cluster1,
            2 : cluster2,
            3 : cluster3,
            4 : cluster8,
            5 : cluster5,
            6 : cluster6,
            7 : cluster7,
            8 : cluster4   
            }
    return olddict3

@pytest.fixture  
def frame20():
    cluster11 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
    cluster12 = {'b1','b2','b3','b4','b5'}
    cluster13 = {'c1','c2','c3','c4','c5','c6'}
    cluster14 = {'d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'}
    cluster15 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster16 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster17 = {'g1','g2','g3'}
    cluster0 = set()
    frame20 = deque([cluster0,cluster11,cluster12,cluster13,cluster14,cluster15,cluster16,cluster17])

    return frame20

@pytest.fixture  
def frame21():
    cluster11 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
    cluster12 = {'b1','b2','b3','b4','b5'}
    cluster13 = {'c1','c2','c3','c4','c5','c6'}
    cluster14 = {'d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'}
    cluster15 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster16 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster17 = {'g1','g2','g3'}
    cluster18 = {'b7','b8','b9','b10'}
    cluster0 = set()
    frame21 = deque([cluster0,cluster11,cluster12,cluster13,cluster14,cluster15,cluster16,cluster17,cluster18])

    return frame21

@pytest.fixture  
def frame22():
    cluster11 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
    cluster12 = {'b1','b2','b3','b4','b5'}
    cluster13 = {'c1','c2','c3','c4','c5','c6'}
    cluster14 = {'d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'}
    cluster15 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster16 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster17 = {'g1','g2','g3'}
    cluster0 = set()
    frame22 = deque([cluster0,cluster11,cluster12,cluster13,cluster14,cluster15,cluster16,cluster17])

    return frame22

@pytest.fixture  
def frame23():
    cluster11 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10'}
    cluster12 = {'b1','b2','b3','b4','b5'}
    cluster13 = {'c1','c2','c3','c4','c5','c6'}
    cluster14 = {'d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'}
    cluster15 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster16 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster17 = {'g1','g2','g3'}
    cluster0 = set()
    frame23 = deque([cluster0,cluster11,cluster12,cluster13,cluster14,cluster15,cluster16,cluster17])

    return frame23

@pytest.fixture
def multimergeframe():
    cluster11 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5','c6','d1','d2','d3','d4','d5','d6','d7','d8','d9','d10'}
    cluster15 = {'f1','f2','f3','f4','f5','f6','f7','f8','f9','f10'}
    cluster16 = {'e1','e2','e3','e4','e5','e6','e7','e8'}
    cluster17 = {'g1','g2','g3'}
    cluster0 = set()
    multimergeframe = deque([cluster0,cluster11,cluster15,cluster16,cluster17])

    return multimergeframe

@pytest.fixture
def testingclusternumber():
    testingclusternumber = 8
    return testingclusternumber


@pytest.fixture
def test3deprecated():
    test3deprecated = {8: {'d1','d2','d3'}}
    return test3deprecated

@pytest.fixture
def deprecated_dict():
    
    deprecated_dict = {
            8 : {'b7','b8','b9'}
    }
    return deprecated_dict


@pytest.fixture  
def bilayerdict():
    cluster1 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}
  
    bilayerdict = {
            0 : set(),
            1 : cluster1,}
    return bilayerdict


@pytest.fixture  
def deprecated_bilayer_dict():   
    deprecated_bilayer_dict = {
            2 : {'b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}
    }
    return deprecated_bilayer_dict

@pytest.fixture 

def frame_bilayer():
    cluster0 = set()
    cluster11 ={'a1','a2','a3','a4','a5','a6','a7','a8','a9'}
    cluster12 = {'b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','a10'}
    bilayer_frame = deque([cluster0,cluster11,cluster12])
    return bilayer_frame

@pytest.fixture 

def zerodict():
    cluster1 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7'}
  
    zerodict = {
            0 : {'b8','b9','b10'},
            1 : cluster1,}
    return zerodict

@pytest.fixture  

def frame_zero():
    cluster0 = {'a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}
    zero_frame = deque([cluster0])
    return zero_frame    