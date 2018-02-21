import numpy as np
import time
import pdb
import os
import glob
import json
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import mdtraj
from multiprocessing import Pool
import itertools
plt.style.use('bmh')

def _find_nn(traj, tracerid, number_neighbors=6):
    """ Find nearest neighbors of a given tracer id
    Utilize mdTraj topology selection to make ignore an entire molecule
    if its atom is already one of the nearest neighbors
    Hardcoded to find 6 neighbors

    """
    atoms = [a for a in traj.topology.atoms]
    frame_by_frame_nn = []
    for frame in range(traj.n_frames):
        nearest_neighbors=[]
        new_selection_string=""
        #base_selection_string = "not resname HOH and not resid {} and not resid ".format(tracerid)
        #group1 = traj.topology.select('not resname HOH and not resid {}'.format(tracerid))
        base_selection_string = "not resid {} and not resid ".format(tracerid)
        group1 = traj.topology.select('not resid {}'.format(tracerid))

        group2 = traj.topology.select('resid {}'.format(tracerid))
        ignore_resids = []

        for i in range(number_neighbors):
            a1, a2, dist = mdtraj.find_closest_contact(traj, group1, group2, frame=frame)
            # Find the residue of atom1
            res1_id = atoms[a1].residue.index
            res2_id = atoms[a2].residue.index
            ignore_resids.append(res1_id)
            ignore_resids.append(res2_id)
            new_selection_string = base_selection_string + " and not resid ".join([str(a) for a in ignore_resids])
            group1 = traj.topology.select(new_selection_string)
            
            #nearest_neighbors.append(atoms[a1].residue)
            #nearest_neighbors.append(atoms[a2].residue)
            nearest_neighbors.append(atoms[a1].residue.name)
            nearest_neighbors.append(atoms[a2].residue.name)

        nearest_neighors=set(nearest_neighbors)
        frame_by_frame_nn.append(nearest_neighbors)
    return frame_by_frame_nn

def _extract_resname_timeseries(nearest_neighbors, resname=None):
    """ Given a timeseries that lists the nearest neighbors,
    pluck out the molecules with the particular resname"""
    #nearest_neighbors = open(tracer_nn, 'r').readlines()
    res_neighbor_timeseries = []
    for frame_neighbor_list in nearest_neighbors:
        #items = line.split(",")
        extracted_items = [item for item in frame_neighbor_list if resname in item]
        #extracted_items = [item for item in frame_neighbor_list if resname in item.name]
        res_neighbor_timeseries.append(len(extracted_items))
    return res_neighbor_timeseries


curr_dir = os.getcwd()
all_sweeps = [thing for thing in os.listdir() if os.path.isdir(thing) and 'sweep' in thing[0:5]]
sample_z_windows = np.loadtxt("{}/z_windows.out".format(all_sweeps[0]))
all_neighbor_dicts = []
start=time.time()
for sweep in all_sweeps:
    print("Nearest neighbors for {}".format(sweep))
    os.chdir(os.path.join(curr_dir, sweep))
#if __name__ == "__main__":
    #curr_dir = os.getcwd()
    z_windows = np.loadtxt('z_windows.out')
    # Overarching datatype will be dictionaries within dictionaries
    # 1st key is zwindow
    # 2nd key is resname
    # The value is the timeseries of n_neighbors
    # 2nd key is 'all_neighbors'
    # The value is the list of all neighbors over all frames
    neighbor_dict = {key : {} for key in  z_windows}

    # Create a list of all possible residue names, which is later turned to set
    all_res_names = set()

    sims = [thing for thing in os.listdir() if 'Sim' in thing and os.path.isdir(thing)]
    n_sims = len(sims)
    # Iterate through sim folders, relating window to all nearest neighbors
    #for sim in range(1):
    for sim in range(n_sims):
        sim_folder = "Sim{}".format(sim)
        os.chdir(os.path.join(curr_dir, sweep,sim_folder))
        tracers = np.loadtxt('tracers.out', dtype=int)
        z_windows_i = z_windows[sim::n_sims]

        # load trajectory
        xtcfile = glob.glob("Stage4_Eq*.xtc")[0]
        grofile = glob.glob("Stage4_Eq*.gro")[0]
        traj = mdtraj.load(xtcfile, top=grofile)
        tracer_data = np.loadtxt('tracers.out',dtype=int)
        residues = [res for res in traj.topology.residues]
        all_res_names.update(res.name for res in residues)

        #for r in residues:
            #all_res_names.update(r.name)

        atoms = [a for a in traj.topology.atoms]

        # Parallel Version
        with Pool() as p:
            all_nn = p.starmap(_find_nn, zip(itertools.repeat(traj),tracer_data-1))
        for (window, nn) in zip(z_windows_i, all_nn):
            neighbor_dict[window].update({'neighbors':nn})

        # Serial Version
        #for i, (tracer,window) in enumerate(zip(tracer_data, z_windows_i)):
        #    print("*"*20)
        #    print("Analyzing tracer {0} ({1}/{2})".format(tracer, str(i+1), str(len(tracer_data))))
        #    # Find nearest neighbors
        #    nn = _find_nn(traj, tracer-1)
        #    neighbor_dict[window].update({'neighbors':nn})

    # Iterate through neighbor_dict
    for key, val  in neighbor_dict.items():
        # Extract the number of `resname` neighbors for each frame
        for resname in all_res_names:
            if 'neighbors' in val.keys():
                val.update({resname: _extract_resname_timeseries(val['neighbors'],
                    resname=resname)})
            else:
                print("Issue with window: {:.2f}".format(key))

    
    # For a given resname, plot the neighbor profiles
    for resname in all_res_names:
        fig, ax = plt.subplots(1,1, figsize=(8,6))

        ax.errorbar(z_windows, 
                [np.mean(neighbor_dict[window][resname]) for window in neighbor_dict.keys()],
                yerr=[np.std(neighbor_dict[window][resname]) for window in neighbor_dict.keys()], 
                fmt='-o')
        ax.set_ylabel("Number of {} neighbors".format(resname), fontsize=20)
        ax.set_xlabel("Z (nm)", fontsize=20)
        fig.savefig("neighbors_{}.svg".format(resname))
        plt.close()

    with open("nearest_neighbors.dat", 'w') as f:
        json.dump(neighbor_dict, f)
    all_neighbor_dicts.append(neighbor_dict)

os.chdir(curr_dir)
summary_neighbors = {}
for window in z_windows:
    if window not in summary_neighbors.keys():
        summary_neighbors.update({window: {}})
    for resname in all_res_names:
        if resname not in summary_neighbors[window].keys():
            summary_neighbors[window].update({resname: []})

        for neighbor_dict in all_neighbor_dicts:
            avg_neighbor = np.mean(neighbor_dict[window][resname])
            summary_neighbors[window][resname].append(avg_neighbor)
for resname in all_res_names:
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    ax.errorbar(z_windows, [np.mean(summary_neighbors[window][resname]) for window in summary_neighbors.keys()],
        yerr=[np.std(summary_neighbors[window][resname]) for window in summary_neighbors.keys()], fmt='-o')
    ax.set_ylabel("Numer of {} neighbors".format(resname),fontsize=20)
    ax.set_xlabel("Z (nm)", fontsize=20)
    fig.savefig("neighbors_{}.svg".format(resname))
    plt.close()
end = time.time()
print("NN analysis took {}".format(end-start))
