from __future__ import print_function
import mdtraj as md
import numpy as np
import os
import math

# This script will first calculate centroids, then RMSF and B-factor from MD simulation output
for f in os.listdir('./run_0/'):
    if f.endswith("_1.dcd"):
        # load the file
        print(f)
        # setup traj and topology
        traj = md.load('./run_0/' + f, top="../build_system/MM/" + f.split("_")[0]+"_"+f.split("_")[1]+"_"+f.split("_")[2] + ".pdb")
        topology = traj.topology
        # Calculate centroid
        atom_indices = [a.index for a in traj.topology.atoms if a.element.symbol != 'H']
        distances = np.empty((traj.n_frames, traj.n_frames))
        for i in range(traj.n_frames):
            distances[i] = md.rmsd(traj, traj, i, atom_indices=atom_indices)
        beta = 1
        index = np.exp(-beta*distances / distances.std()).sum(axis=1).argmax()
        centroid = traj[index]
        # Calculate RMSD
        rmsf = md.rmsf(traj, centroid, 0)
        # Calculate B_factor
        B_factor = 8 / 3 * math.pi * math.pi * ((rmsf)**2)
        # Zip atom and B_factor to tables
        atom = [str(atom).split("-")[1] for atom in topology.atoms]
        B_factor_table  = list(zip(atom, B_factor))
        print('B factor table: ', B_factor_table)
