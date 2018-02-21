import numpy as np
import pdb
from multiprocessing import Pool
import os
import glob
import bilayer_analysis_functions
import mdtraj
from collections import OrderedDict
import itertools

def _get_z_variation(sweep, sim):
    os.chdir(os.path.join(curr_dir,sweep,sim))
    xtc_file = glob.glob("Stage4_Eq*.xtc")[0]
    gro_file = glob.glob("Stage4_Eq*.gro")[0]
    traj = mdtraj.load(xtc_file, top=gro_file)
    topol = traj.topology
    lipid_dict, headgroup_dict = bilayer_analysis_functions.get_lipids(topol)
    lipid_tails,lipid_heads = bilayer_analysis_functions.get_lipid_tails(topol, lipid_dict)

    bot_leaf, top_leaf = bilayer_analysis_functions.identify_leaflets(traj, topol, lipid_dict)

    n_lipid = len(lipid_dict.keys())
    n_lipid_tails = len(lipid_tails.keys())
    n_tails_per_lipid = n_lipid_tails/n_lipid

    # Iterate through the bot leaflet or top leaflet
    # Get an atom, its residue, and its residue index
    # If the residue index hasn't been tabulated, add the average z coordinate
    already_examined_resids = set()
    top_avgs = []
    top_stds = []
    bot_avgs = []
    bot_stds = []
    for i, leaflet in enumerate((bot_leaf, top_leaf)):
        all_avgs = []
        all_stds = []
        for atomid in leaflet:
            resindex = topol.atom(atomid).residue.index
            if str(resindex) in lipid_heads.keys():
                if resindex not in already_examined_resids:
                    z_stats = np.array([[np.mean(xyz), np.std(xyz)] for xyz in \
                         traj.xyz[:,lipid_heads[str(resindex)],2]])
                    already_examined_resids.add(resindex)
                    all_avgs.append(z_stats[:,0])
                    all_stds.append(z_stats[:,1])
        if i == 0:
            bot_avgs = np.mean(np.asarray(all_avgs),axis=0)
            bot_stds = np.std(np.asarray(all_avgs),axis=0) 
        else:
            top_avgs = np.mean(np.asarray(all_avgs),axis=0)
            top_stds = np.std(np.asarray(all_avgs),axis=0) 
    return bot_stds, top_stds


all_sweeps = [thing for thing in os.listdir() if os.path.isdir(thing) and 'sweep' in thing[0:6]]
curr_dir = os.getcwd()
bot_to_print = []
top_to_print = []
for sweep in all_sweeps:
    print("----{}----".format(sweep))
    os.chdir(os.path.join(curr_dir, sweep))
    all_sims = [thing for thing in os.listdir() if os.path.isdir(thing) and 'Sim' in thing]
    #sim_bot_stds = []
    #sim_top_stds = []
    with Pool() as p:
        all_stds = p.starmap(_get_z_variation, zip(itertools.repeat(sweep), all_sims))
    all_stds = np.asarray(all_stds)
    sim_bot_stds = all_stds[:,0]
    sim_top_stds = all_stds[:,1]
        
        #sim_bot_stds.append(np.mean(bot_stds))
        #sim_top_stds.append(np.mean(top_stds))
    
    bot_to_print.append("{} {}\n".format(sweep[-7:], np.mean(sim_bot_stds)))
    top_to_print.append("{} {}\n".format(sweep[-7:], np.mean(sim_top_stds)))
    print(np.mean(sim_bot_stds))
    print(np.mean(sim_top_stds))

