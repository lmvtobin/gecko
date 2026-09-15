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

# Script for counting and plotting occurence of Lens spaces in a data file from sample3 or exhaustive3

import numpy as np
import matplotlib.pyplot as plt

from regina import *
import gecko_utils as gecko

max_size = 20 # set manually
knot = False # whether to read link or knot data, set manually
sample_size = 10000 # set manually

names = []
counts = {}

for n in range(2,max_size+1):
    samples = []
    if (knot):
        samples = gecko.read_sample_data("../data/sample3/knots_n"+str(n)+"_"+str(sample_size)+"_0.txt")
    else:
        samples = gecko.read_sample_data("../data/sample3/n"+str(n)+"_"+str(sample_size)+"_0.txt")
    sample_counts = {}
    for (sigs,grids) in samples:
        if len(sigs) > 1:
            continue
        trig = Triangulation3(sigs[0])
        name = ""

        recog = StandardTriangulation.recognise(trig)
        if recog:
            name = recog.manifold().name()
            if name == "RP3":
                name = "L(2,1)"
            if name.startswith("L("):
                if name in sample_counts:
                    sample_counts[name] += len(grids)
                else:
                    sample_counts[name] = len(grids)
                if not name in names:
                    names.append(name)
            continue

        lookup = Census.lookup(trig)
        if lookup:
            name = lookup[0].name().split(" : #")[0]
            if name == "RP3":
                name = "L(2,1)"
            if name.startswith("L("):
                if name in sample_counts:
                    sample_counts[name] += len(grids)
                else:
                    sample_counts[name] = len(grids)
                if not name in names:
                    names.append(name)
            continue

    print("n = ",n," | lens spaces: ",len(sample_counts.keys()))
    counts[n] = sample_counts

print("\ntotal lens spaces:",len(names))
name_counts = {}
for name in names:
    name_counts[name] = sum([counts[n][name] if name in counts[n].keys() else 0 for n in counts.keys()])

for name in sorted(names,key=lambda x:name_counts[x],reverse=True):
    print(name,":",name_counts[name])

fig0 = plt.figure(0,figsize=(8,5))
for name in sorted(names,key=lambda x:name_counts[x],reverse=True)[2:12]: # next 10 after S1xS2 by total number
# for name in sorted(names,key=lambda x: 0 if not x in counts[max_size].keys() else counts[max_size][x],reverse=True)[2:12]: # next 10 after S1xS2 by number in largest grid size
# for name in sorted([name for name in names if name.endswith(",1)")],key=lambda x:int(x[2:].split(',')[0])): # L(n,1) only
    plt.plot(range(2,max_size+1), [counts[n][name] if name in counts[n] else None for n in range(2,max_size+1)], marker='o', label=name)
plt.title("Identified Lens Spaces")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.xlim(1,max_size+1)
plt.xticks(list(range(2,max_size+1)))
plt.legend()

fig1 = plt.figure(1,figsize=(8,5))
for name in sorted(names,key=lambda x:name_counts[x],reverse=True)[2:12]: # next 10 after S1xS2 by total number
# for name in sorted(names,key=lambda x: 0 if not x in counts[max_size].keys() else counts[max_size][x],reverse=True)[2:12]: # next 10 after S1xS2 by number in largest grid size
# for name in sorted([name for name in names if name.endswith(",1)")],key=lambda x:int(x[2:].split(',')[0])): # L(n,1) only
    plt.plot(range(2,max_size+1), [counts[n][name] if name in counts[n] else None for n in range(2,max_size+1)], marker='o', label=name)
plt.title("Identified Lens Spaces")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.xlim(1,max_size+1)
plt.xticks(list(range(2,max_size+1)))
plt.semilogy()
plt.legend()

plt.show()