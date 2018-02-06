import os
import permeability as prm
import matplotlib.pyplot as plt
import pickle as pickle
from os import system 
import pdb
import numpy as np
import mdtraj as mdt

#data_dir = '/Users/rmhartkamp/Dropbox/PostDoc_2_Vanderbilt/Simulation/Permeability/DSPC_C12OH_3_1'
#data_dir = '/raid6/homes/ahy3nz/Trajectories/Data/11_DSPC_C18OH/DSPC-50_alc18-50_5-4a'
data_dir = os.getcwd()
n_sweeps = 26
#n_sweeps = 30
preamble = "folder"

prm.analyze_sweeps(data_dir, timestep=1.0, verbosity=2, directory_prefix='sweep', n_sweeps=n_sweeps, correlation_length=300)
#prm.analyze_sweeps(data_dir, timestep=1000.0, verbosity=2, directory_prefix='sweep', n_sweeps=n_sweeps)


#forcetime = prm.force_timeseries(data_dir, timestep=2.0, n_windows=40, start_window=15, directory_prefix='sweep')
#prm.plot_timeseries(forcetime['time'], forcetime['forces']) 

 
output = prm.analyze_force_acf_data(data_dir, 305.0, timestep=1, verbosity=2, directory_prefix='sweep',n_sweeps=n_sweeps, kB=1.987e-3)
pickle.dump(output, open('{}/output.p'.format(preamble), 'wb'))

#output = pickle.load(open('output.p', 'rb'))
#pdb.set_trace()

#system('rm *.pdf *.png')

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
print('Permeability (cm/sec): {} ({})'.format(output['permeability'], output['perm_err']))

#prm.plot_forces(output['z'], output['forces'], fig_filename='forces.pdf',sweep_alpha=0.2)
#prm.plot_free_energy_z(output['z'], output['dG'], fig_filename='delta_G.pdf')
#prm.plot_force_acfs_time(output['time'], output['facf_windows'], normalize=True)
#prm.plot_int_acfs_time(output['time'], output['int_facf_windows'])
#prm.plot_symmetrized_free_energy(output['z'], output['dG_sym'], 
#        output['dG_sym_err'],savefig=True)
#prm.plot_sym_diffusion_coefficient_z(output['z'], output['d_z_sym'], 
#        output['d_z_sym_err'],savefig=True)
#prm.plot_resistance_z(output['z'], output['R_z'], output['R_z_err'], savefig=True)
#prm.plot_sym_exp_free_energy(output['z'], output['dG_sym'], output['dG_sym_err'], output['d_z_sym'], 
#        output['d_z_sym_err'], output['R_z'], output['R_z_err'], 305)
#print('Permeability (cm/sec): {} ({})'.format(output['permeability'], output['perm_err']))
##system('open -a preview *.pdf *.png')
