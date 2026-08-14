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

# Script for calculating and plotting occurence of first and second homology groups
# in a data file from sample3 or exhaustive3

import numpy as np
import matplotlib.pyplot as plt

from regina import *
import gecko_utils as gecko

max_size = 20 # set manually
prime_only = True # whether to only count prime manifolds, set manually
knot = False # whether to read link or knot data, set manually
sample_size = 10000 # set manually

H1s = set()
H2s = set()
H1s_triv_H2 = set()

H1_counts = {}
H2_counts = {}
H1_counts_triv_H2 = {}

for n in range(2,max_size+1):
    samples = []
    if (knot):
        samples = gecko.read_sample_data("../sample3_data/knots_n"+str(n)+"_"+str(sample_size)+"_0.txt")
    else:
        samples = gecko.read_sample_data("../sample3_data/n"+str(n)+"_"+str(sample_size)+"_0.txt")
    sample_H1_counts = {}
    sample_H2_counts = {}
    sample_H1_counts_triv_H2 = {}
    for (sigs,grids) in samples:
        if prime_only and len(sigs) > 1:
            continue
        H1 = AbelianGroup()
        H2 = AbelianGroup()
        for sig in sigs:
            trig = Triangulation3(sig)
            H1.addGroup(trig.homology(1))
            H2.addGroup(trig.homology(2))
        if H1.str() in sample_H1_counts:
            sample_H1_counts[H1.str()] += len(grids)
        else:
            sample_H1_counts[H1.str()] = len(grids)
        if H2.str() in sample_H2_counts:
            sample_H2_counts[H2.str()] += len(grids)
        else:
            sample_H2_counts[H2.str()] = len(grids)
        H1s.add(H1.str())
        H2s.add(H2.str())
        if (H2.isTrivial()):
            if H1.str() in sample_H1_counts_triv_H2:
                sample_H1_counts_triv_H2[H1.str()] += len(grids)
            else:
                sample_H1_counts_triv_H2[H1.str()] = len(grids)
            H1s_triv_H2.add(H1.str())
    print("n = ",n," | H1s: ",len(sample_H1_counts.keys())," | H2s: ",len(sample_H2_counts.keys()))
    H1_counts[n] = sample_H1_counts
    H2_counts[n] = sample_H2_counts
    H1_counts_triv_H2[n] = sample_H1_counts_triv_H2

H1s_sorted = sorted(H1s,key=lambda s: sum([H1_counts[n][s] if s in H1_counts[n] else 0 for n in range(2,max_size+1)]),reverse=True)
H2s_sorted = sorted(H2s,key=lambda s: sum([H2_counts[n][s] if s in H2_counts[n] else 0 for n in range(2,max_size+1)]),reverse=True)
H1s_triv_H2_sorted = sorted(H1s_triv_H2,key=lambda s: sum([H1_counts_triv_H2[n][s] if s in H1_counts_triv_H2[n] else 0 for n in range(2,max_size+1)]),reverse=True)

fig0 = plt.figure(0)
for H1_string in H1s_sorted:
    plt.plot(range(2,max_size+1), [H1_counts[n][H1_string] if H1_string in H1_counts[n] else None for n in range(2,max_size+1)], marker='o', label=H1_string)
plt.title("First Homology")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.semilogy()
plt.legend()

fig1 = plt.figure(1)
for H2_string in H2s_sorted:
    plt.plot(range(2,max_size+1), [H2_counts[n][H2_string] if H2_string in H2_counts[n] else None for n in range(2,max_size+1)], marker='o', label=H2_string)
plt.title("Second Homology")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.semilogy()
plt.legend()

fig2 = plt.figure(2)
for H1_string in list(H1s_sorted)[:10]:
    plt.plot(range(2,max_size+1), [H1_counts[n][H1_string] if H1_string in H1_counts[n] else None for n in range(2,max_size+1)], marker='o', label=H1_string)
plt.title("First Homology\n(showing 10 most common only)")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.semilogy()
plt.legend()

fig3 = plt.figure(3)
plt.stackplot(range(2,max_size+1), [[H1_counts[n][H1_string] if H1_string in H1_counts[n] else 0 for n in range(2,max_size+1)] for H1_string in H1s_sorted], labels=H1s_sorted)
plt.title("First Homology")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.legend()

fig4 = plt.figure(4)
for H1_string in H1s_triv_H2_sorted:
    plt.plot(range(2,max_size+1), [H1_counts_triv_H2[n][H1_string] if H1_string in H1_counts_triv_H2[n] else None for n in range(2,max_size+1)], marker='o', label=H1_string)
plt.title("First Homology for H2=0")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.semilogy()
plt.legend()

fig5 = plt.figure(5)
for H1_string in list(H1s_triv_H2_sorted)[:10]:
    plt.plot(range(2,max_size+1), [H1_counts_triv_H2[n][H1_string] if H1_string in H1_counts_triv_H2[n] else None for n in range(2,max_size+1)], marker='o', label=H1_string)
plt.title("First Homology for H2=0\n(showing 10 most common only)")
plt.xlabel("Grid Size")
plt.ylabel("Count / "+str(sample_size))
plt.semilogy()
plt.legend()

plt.show()