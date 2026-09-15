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

# Script for counting and plotting occurence of connected sums of S1xS2
# in a data file from sample3 or exhaustive3

import numpy as np
import matplotlib.pyplot as plt

from regina import *
import gecko_utils as gecko

max_size = 20 # set manually
knot = False # whether to read link or knot data, set manually
sample_size = 10000 # set manually

names = ["S3", "S1 x S2"]
counts = {}

for n in range(2,max_size+1):
    samples = []
    if knot:
        samples = gecko.read_sample_data("../data/sample3/knots_n"+str(n)+"_"+str(sample_size)+"_0.txt")
    else:
        samples = gecko.read_sample_data("../data/sample3/n"+str(n)+"_"+str(sample_size)+"_0.txt")
    sample_counts = {}
    for (sigs,grids) in samples:
        if len(sigs) == 1:
            trig = Triangulation3(sigs[0])

            recog = StandardTriangulation.recognise(trig)
            if recog:
                name = recog.manifold().name()
                if (name == "S3"):
                    if "S3" in sample_counts:
                        sample_counts["S3"] += len(grids)
                    else:
                        sample_counts["S3"] = len(grids)
                if (name == "S2 x S1"):
                    if "S1 x S2" in sample_counts:
                        sample_counts["S1 x S2"] += len(grids)
                    else:
                        sample_counts["S1 x S2"] = len(grids)
                continue

            lookup = Census.lookup(trig)
            if lookup:
                name = lookup[0].name().split(" : #")[0]
                if (name == "S3"):
                    if "S3" in sample_counts:
                        sample_counts["S3"] += len(grids)
                    else:
                        sample_counts["S3"] = len(grids)
                if (name == "S2 x S1"):
                    if "S1 x S2" in sample_counts:
                        sample_counts["S1 x S2"] += len(grids)
                    else:
                        sample_counts["S1 x S2"] = len(grids)
                continue

            if not trig.homology(1).isTrivial() and not trig.homology(1).isZ():
                continue

            groupname = trig.group().recogniseGroup
            if (groupname != "" and groupname != "0" and groupname != "Z"):
                continue

            if trig.isSphere():
                if "S3" in sample_counts:
                    sample_counts["S3"] += len(grids)
                else:
                    sample_counts["S3"] = len(grids)
                continue
            elif not trig.isIrreducible():
                if "S1 x S2" in sample_counts:
                    sample_counts["S1 x S2"] += len(grids)
                else:
                    sample_counts["S1 x S2"] = len(grids)
                continue

        else:
            is_handlebody = True
            for sig in sigs:
                trig = Triangulation3(sigs[0])

                name = ""

                recog = StandardTriangulation.recognise(trig)
                if recog:
                    name = recog.manifold().name()
                    if (name == "S2 x S1"):
                        continue
                    else:
                        is_handlebody = False
                        break

                lookup = Census.lookup(trig)
                if lookup:
                    name = lookup[0].name().split(" : #")[0]
                    if (name == "S2 x S1"):
                        continue
                    else:
                        is_handlebody = False
                        break

                if not trig.homology(1).isZ():
                    is_handlebody = False
                    break

                groupname = trig.group().recogniseGroup
                if (groupname == "Z"):
                    continue
                elif (groupname != ""):
                    is_handlebody = False
                    break

                if trig.isIrreducible():
                    is_handlebody = False
                    break

            if is_handlebody:
                name = "#^"+str(len(sigs))+" (S1 x S2)"
                if name in sample_counts:
                    sample_counts[name] += len(grids)
                else:
                    sample_counts[name] = len(grids)
                if not name in names:
                    names.append(name)


    print("n =",n,"| handlebody boundaries: ",len(sample_counts.keys()))
    print(sample_counts)
    counts[n] = sample_counts

fig0 = plt.figure(0,figsize=(8,5))
plt.plot(counts.keys(), [sum([counts[n][key] for key in counts[n]]) for n in counts.keys()], marker='o', label="Total")
for name in names:
    plt.plot(counts.keys(), [counts[n][name] if name in counts[n] else None for n in counts.keys()], marker='o', label=name)
plt.title("Identified Connected Sums of S1 x S2")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.xticks(list(counts.keys()))
plt.legend()

fig1 = plt.figure(1,figsize=(8,5))
plt.plot(counts.keys(), [sum([counts[n][key] for key in counts[n]]) for n in counts.keys()], marker='o', label="Total")
for name in names:
    plt.plot(counts.keys(), [counts[n][name] if name in counts[n] else None for n in counts.keys()], marker='o', label=name)
plt.title("Identified Connected Sums of S1 x S2")
plt.xlabel("Grid Size") 
plt.ylabel("Count / "+str(sample_size))
plt.xticks(list(counts.keys()))
plt.semilogy()
plt.legend()

plt.show()