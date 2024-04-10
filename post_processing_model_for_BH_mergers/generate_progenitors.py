import sys
sys.path.append('/home/aklantbhowmick/anaconda_prev/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda_prev/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda_prev/envs/nbodykit-env/lib/python3.6/site-packages/')
import mdot_to_Lbol
import arepo_package
import scipy.interpolate
radiative_efficiency=0.2
total_conv=mdot_to_Lbol.get_conversion_factor_arepo(radiative_efficiency)
import h5py
import illustris_python as il
import numpy 
import pandas as pd
import scipy.integrate
from scipy.optimize import curve_fit
import os
from Config import *

title_fontsize=30
import warnings



'''
path_to_output='/orange/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS///' # this is the folder containing the simulation run
run='/L3p125n512/AREPO/' # name of the simulation runs
output='/output_ratio1000_SFMFGM5_seed3.19_bFOF/' 
basePath=path_to_output+run+output
desired_redshift=7
'''

BH_Progs,o=arepo_package.get_particle_property(basePath,'BH_Progs',5,desired_redshift)

indices=numpy.arange(0,len(BH_Progs))


def create_directory_if_not_exists(directory_name):
    # Check if the directory exists
    if not os.path.exists(basePath + directory_name):
        # Create the directory if it doesn't exist
        os.makedirs(basePath + '/' + directory_name)
        print(f"Directory '{directory_name}' created.")
    else:
        print(f"Directory '{directory_name}' already exists.")
    
create_directory_if_not_exists('merger_progenitors')

for N in indices:
    rootblackhole,N_mergers,indices_of_included_events=arepo_package.get_blackhole_progenitors(basePath,N,desired_redshift,0.01)
    numpy.save(basePath+'/merger_progenitors/'+'progenitor_merger_indices_z%d_%d.npy'%(desired_redshift,N),indices_of_included_events)
