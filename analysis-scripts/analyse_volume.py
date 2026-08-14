"""
* * * * * * * * * * * * * * * * * * * * *
*                           o---------x *
*     G e c K o         x-o |  #  w w | *
*                     o-|-x |       w | *
*       Grid          x-|---|---o     | *
*       Kirby     x-----o   |   |   # | *
*                 |     x---|---|-----o *
*             x---o     |   |   | x-o   *
*             |     o---|---x o-|-|-x   *
*         o-x |     |   |     | x-o     *
*       x-|-|-o x---|---o     |         *
*       o-x |   |   |     o---x         *
* o---------|---|---x     |             *
* |         |   |   o-----x             *
* |     o---|---|---|-x                 *
* |     |   |   | x-|-o        By       *
* |   x-|---o   | o-x         Lucy      *
* | o-|-x       |             Tobin     *
* x-|-o         |                       *
*   x-----------o                       *
* * * * * * * * * * * * * * * * * * * * *
"""

# Script for plotting hyperbolic volume and counting likely hyperbolic 3-manifolds, from
# both decomposition data files (sample3, exhaustive3) and volume data files (sample3-volume).

# Includes one data file for each size up to the specified max_size where one exists, preferencing
# volume data files.

# Note that when decomposition files are used this will run significantly slower than with
# volume files (where the necessary calculations have already been done), and involves some
# randomness in retriangulating so will not give precisely the same result when re-run

import functools
import os.path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from regina import *
import gecko_utils as gecko

volumes = {}

max_size = 12 # set manually
knot = False # whether to read link or knot data, set manually
sample_size = 10000 # set manually

knot_ex1 = ""
knot_ex2 = ""
if (knot):
    knot_ex1 = "knot_"
    knot_ex2 = "knots_"

for n in range(2,max_size+1):
    sample_volumes = []
    if os.path.isfile("../sample3_data/volumes_"+knot_ex1+"n"+str(n)+"_"+str(sample_size)+"_0.txt"):
        samples = gecko.read_volume_data("../sample3_data/volumes_"+knot_ex1+"n"+str(n)+"_"+str(sample_size)+"_0.txt")
        sample_volumes = [pair[0] for pair in samples[0]]
    elif os.path.isfile("../sample3_data/"+knot_ex2+"n"+str(n)+"_"+str(sample_size)+"_0.txt"):
        samples = gecko.read_sample_data("../sample3_data/"+knot_ex2+"n"+str(n)+"_"+str(sample_size)+"_0.txt")
        sample_volumes = []
        for (sigs,grids) in samples:
            hyp = True
            if len(sigs) == 1:
                trig = Triangulation3(sigs[0])
                sptrig = SnapPeaTriangulation(trig)
                volume_ests = []
                for i in range(100):
                    sptrig.randomise()
                    vol = sptrig.volume()
                    if vol < 0.5:
                        hyp = False
                        break
                    else:
                        volume_ests.append(vol)
                if (hyp):
                    volume = np.mean(volume_ests)
                    if volume > 0.9:
                        for i in range(len(grids)):
                            sample_volumes.append(volume)
    if len(sample_volumes) > 0:
        print("n = ",n," | likely hyperbolic count = ",len(sample_volumes)," | mean = ",np.mean(sample_volumes)," | median = ",np.median(sample_volumes))
        volumes[n] = sample_volumes

fig0 = plt.figure(0,figsize=(8,5))
plt.errorbar(volumes.keys(),[np.mean(volumes[k]) for k in volumes.keys()],yerr=[np.std(volumes[k]) for k in volumes.keys()],marker='o',capsize=8,label="Mean +/- Standard Deviation")
plt.title("Average Volume of Likely Hyperbolic Manifolds")
plt.xticks(range(5,max_size+5,5))
plt.xlim(0,max_size+5)
plt.ylim(bottom=0)
plt.xlabel("Grid Size")
plt.ylabel("Volume")
plt.legend()

fig1 = plt.figure(1,figsize=(8,5))
plt.plot(volumes.keys(),[max(volumes[k]) for k in volumes.keys()],marker='o',label="Max")
plt.plot(volumes.keys(),[np.quantile(volumes[k],0.75) for k in volumes.keys()],marker='o',label="Third Quartile")
plt.plot(volumes.keys(),[np.median(volumes[k]) for k in volumes.keys()],marker='o',label="Median")
plt.plot(volumes.keys(),[np.quantile(volumes[k],0.25) for k in volumes.keys()],marker='o',label="First Quartile")
plt.plot(volumes.keys(),[min(volumes[k]) for k in volumes.keys()],marker='o',label="Min")
plt.title("Maximum Hyperbolic Volume")
plt.xticks(range(5,max_size+5,5))
plt.xlim(0,max_size+5)
plt.ylim(bottom=0)
plt.xlabel("Grid Size")
plt.ylabel("Volume")
plt.legend()

fig2 = plt.figure(2,figsize=(8,5))
plt.plot(volumes.keys(),[len(volumes[k]) for k in volumes.keys()],marker='o')
plt.plot([0,max_size+5],[sample_size,sample_size],color='black',ls='--')
plt.title("Number of Likely Hyperbolic Manifolds")
plt.xticks(range(5,max_size+5,5))
plt.xlim(0,max_size+5)
plt.ylim(bottom=0)
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))

fig3 = plt.figure(3)
x = []
y = []
for n in [k for k in volumes.keys() if k <= 30]:
    x = x + [n for i in range(len(volumes[n]))]
    y = y + volumes[n]

# NOTE: values on this line should be set manually
# this is for size up to 30 (the largest for which I had data for every intermediate grid size in the initial run)
# vmax was set to the largest value for knots in my data, so that both knots and links used the same scale
_,_,_,im = plt.hist2d(x,y,bins=[[n-0.5 for n in range(2,30+1)]+[30+0.5],np.linspace(0,80,400)],vmax=178) 
plt.title("Distribution of Hyperbolic Volume")
plt.xlabel("Grid Size")
plt.ylabel("Volume")
plt.colorbar(im)

plt.show()