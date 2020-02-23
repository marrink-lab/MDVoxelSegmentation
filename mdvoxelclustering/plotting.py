import numpy as np
import matplotlib.pyplot as plt
import pickle
from matplotlib.collections import LineCollection

my_dpi = 300
scale = 0.01
with open('cluster_mutations.pickle', 'rb') as handle:
    cluster_mutations = pickle.load(handle)
with open('n_clusters.pickle', 'rb') as handle:
    n_clusters = pickle.load(handle)

n_clusters = n_clusters

with open('visualization_data.pickle', 'rb') as handle:
    visualization_data = pickle.load(handle)

frames = len(visualization_data)
threshold = 50
plottingdata = {0: [[[0, 0], [0, 200]]]}

for key, value in cluster_mutations.items():
    for change in value:
        if change[2] == 'c':
            plottingdata[change[0]] =[[[change[0],change[1]],[key,key]],[[change[0],change[0]],[key,frames]]]
            
        if change[2] == 'm':
            plottingdata[change[0]][-1][1][1] = key
            plottingdata[change[0]].append([[change[0],change[1]],[key,key]])
            
        if change[2] == 'rc':
            plottingdata[change[0]].extend([[[change[0],change[1]],[key,key]],[[change[0],change[0]],[key,frames]]])


fig, ax = plt.subplots()
ax.set_xlim(0,6)
ax.set_ylim(0,len(visualization_data))
ax.set_xticks(range(1,6))
ax.set_aspect
thicknessratio = 0
for value in visualization_data[0]:
    thicknessratio += value[1]
    
thicknessratio = thicknessratio * scale

thickness = np.zeros((len(visualization_data),len(plottingdata)))
for key, value in visualization_data.items():
    for point in value:
        thickness[key][point[0]] = point[1]



for key, value in plottingdata.items():
    if key == 0:
        continue
    for linedata in value:
        y = np.linspace(linedata[1][0],linedata[1][1],linedata[1][1] - linedata[1][0] ,False, dtype = int)
        x = np.linspace(linedata[0][0],linedata[0][1],abs(linedata[0][1] - linedata[0][0]) ,False, dtype = int)
        points = []
        linethickness = []
        if len(y) > len(x):
            for i in y:
                points.append([(linedata[0][0],i),(linedata[0][0],i+1)])
                linethickness.append(thickness[i][linedata[0][0]]/thicknessratio)
            line_segments = LineCollection(points,linethickness)
            ax.add_collection(line_segments)
        else:
            if thickness[linedata[1][0]][(linedata[0][0]-1)] > threshold:
                print(thickness[linedata[1][0]][(linedata[0][0]-1)])
                print(linedata[0],linedata[1])
                print(linedata[1][0])
                plt.plot(linedata[0],linedata[1],'-k')
            else:
                if key < 10:
                    print(thickness[linedata[1][0]][(linedata[0][0]-2)])
                    print(linedata[0],linedata[1])
                    print(linedata[1][0])

    #            for i in x:
    #                points.append([(i,linedata[1][0]),(i+1,linedata[1][0])])
    #                if thickness[linedata[1][0]][linedata[0][0]]/thicknessratio == 0:
    #                    linethickness.append(0.000000001)
    #                else: 
    #                    linethickness.append(thickness[linedata[1][0]][linedata[0][0]]/thicknessratio)


            
          


plt.savefig('bababoem.png')

#plt.bar(time, ordered_clusters, align='center')
#plt.xticks(range(len(plotinfo)), list(plotinfo.keys()))
#        



# -*- coding: utf-8 -*-

