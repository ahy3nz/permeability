import os
import permeability as prm
import matplotlib.pyplot as plt
import pickle as pickle
from os import system 
import pdb
import numpy as np
import subprocess

data_dir = os.getcwd()
to_print = []
sweeps = [thing for thing in os.listdir() if os.path.isdir(thing) and 'sweep' in thing[0:5]]

#truncation_pts = {133336: "400ps", 166669:"500ps", 200002:"600ps",
        #233336:"700ps", 266669:"800ps", 300002:"900ps", None:"1000ps"}
truncation_pts = {None:"1000ps"}

for trunc_end, label in truncation_pts.items():
    print("*"*20)
    print("Analyzing {}...".format(label))
    # Make the appropraite directory
    preamble = "{}_sampling3".format(label)
    p = subprocess.Popen("mkdir {}".format(preamble), shell=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    # analyze sweeps, passing truncation
    prm.analyze_sweeps(data_dir, timestep=1.0, verbosity=0, directory_prefix='sweep', n_sweeps=26, correlation_length=250, end=trunc_end)

    # Perform integrals
    output = prm.analyze_force_acf_data(data_dir, 305.0, timestep=1, verbosity=0, directory_prefix='sweep',n_sweeps=26, kB=1.987e-3)
    pickle.dump(output, open('{}/output.p'.format(preamble), 'wb'))

    # Save permeability
    to_print.append("{}\t{}\t{}".format(label, output['permeability'], output['perm_err']))

    # Plotting
    prm.plot_forces(output['z'], output['forces'], fig_filename='{}/forces.pdf'.format(preamble),sweep_alpha=0.2)
    prm.plot_free_energy_z(output['z'], output['dG'], fig_filename='{}/delta_G.pdf'.format(preamble))
    prm.plot_force_acfs_time(output['time'], output['facf_windows'], fig_filename="{}/force_acf.png".format(preamble), normalize=True)
    prm.plot_int_acfs_time(output['time'], output['int_facf_windows'], fig_filename="{}/int-acf.png".format(preamble))
    prm.plot_symmetrized_free_energy(output['z'], output['dG_sym'], 
            output['dG_sym_err'],savefig=True, fig_filename="{}/delG-sym.pdf".format(preamble))
    prm.plot_sym_diffusion_coefficient_z(output['z'], output['d_z_sym'], 
            output['d_z_sym_err'],savefig=True, fig_filename="{}/d-sym_z.pdf".format(preamble))
    prm.plot_resistance_z(output['z'], output['R_z'], output['R_z_err'], savefig=True, fig_filename="{}/R_z.pdf".format(preamble))
    prm.plot_sym_exp_free_energy(output['z'], output['dG_sym'], output['dG_sym_err'], output['d_z_sym'], 
            output['d_z_sym_err'], output['R_z'], output['R_z_err'], 305,
        fig_filename="{}/expdelG-sym.pdf".format(preamble))





    # Look at the range of permeabilities
    all_perm = []
    for sweep in sweeps:
        print("Individually looking at {}".format(sweep))
        output = prm.analyze_force_acf_data(data_dir, 305.0, timestep=1, verbosity=1, directory_prefix=sweep,n_sweeps=1, kB=1.987e-3, parallel=False)
        all_perm.append(output['permeability'])
    np.savetxt("{}/separate_perm.dat".format(preamble), all_perm)


with open("1000ps.dat", 'w') as f:
    for line in to_print:
        f.write("{}\n".format(line))
