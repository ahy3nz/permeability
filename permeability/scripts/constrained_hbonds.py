import os
import sys
import numpy as np
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import mdtraj as mdt
import MDAnalysis as mda
import MDAnalysis.analysis.hbonds as hbonds
import pandas as pd


""" Compute hydrogen bonds for a Z-constrained MD trajectory



"""
# In the simulation folder
curr_dir = os.getcwd()
directory_prefix = 'sweep'
sweep_dirs = [f for f in os.listdir() if directory_prefix in f and os.path.isdir(f)]
# In the sweep directory
for sweep_folder in sweep_dirs:
    os.chdir(os.path.join(os.getcwd(), sweep_folder))

    tracers = np.loadtxt('tracers.out')
    n_tracer = np.shape(tracers)[0]

    z_windows = np.loadtxt('z_windows.out')
    n_windows = np.shape(z_windows)[0]

    sim_folders = [sim for sim in os.listdir() if 'Sim' in sim and os.path.isdir(sim)]
    n_sims = len(sim_folders)

    # Need to identify the tracer
    # In the Sim folder
    for sim_index, sim in enumerate(sim_folders):
        os.chdir(os.path.join(os.getcwd(), sim))
        pdbfile = [f for f in os.listdir() if '.pdb' in f and 'Stage5' in f]
        grofile = [f for f in os.listdir() if '.gro' in f and 'Stage5' in f]
        tprfile = [f for f in os.listdir() if '.tpr' in f and 'Stage5' in f]

        # If no pdbfile, make one
        if len(pdbfile) == 0:
            with open('gmx_pdb.log', 'w') as outfile:
                p = subProcess.Popen('gmx trjconv -f {0} -s {1} -o {2} -conect -pbc mol'.format(grofile, tprfile, pdbfile), shell=True, stdout=outfile, stderr=outfile)
                p.wait()

        # Compute hbonds for every tracer and output to a different file
        for tracer_index, tracer in enumerate(tracers):
            u = mda.Universe(pdbfile)
            hbond_analysis = hbonds.HydrogenBondAnalysis(u, 'resid {}'.format(tracer), 'resname HOH or resname SOL', distance=2.5, angle=150)
            hbond_analysis.run()
            hbond_analysis.generate_table()
            hbond_info = pd.DataFrame.from_records(hbond_analysis.table)
            hbond_dict={}
            for donor, acceptor in zip(hbond_info['donor_resnm'], hbond_info['acceptor_resnm']):
                key = donor + '-' + acceptor
                if key in hbond_dict:
                    hbond_dict[key] += 1
                else:
                    hbond_dict[key] = 1

            file_index = (tracer_index * n_sims) + sim_index
            with open('hbonds_mda_{0}.dat'.format(file_index), 'w') as f:
               for key,value in hbond_dict.iteritems():
                   f.write('{:<20s}: {} \n'. format(key, value))
    os.chdir(curr_dir)

