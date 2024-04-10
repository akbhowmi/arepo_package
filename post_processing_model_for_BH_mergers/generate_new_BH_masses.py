import sys
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
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
warnings.filterwarnings("ignore")

label_fontsize=40

'''
path_to_output='/orange/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS///' # this is the folder containing the simulation run
run='/L3p125n512/AREPO/' # name of the simulation runs
output='/output_ratio1000_SFMFGM5_seed3.19_bFOF/' 
basePath=path_to_output+run+output


desired_redshift=7
'''


SubhaloBHMass,o=arepo_package.get_group_property(basePath,'GroupBHMass',desired_redshift,file_format='fof_sub_subfind')
SubhaloMassType,o=arepo_package.get_group_property(basePath,'GroupMassType',desired_redshift,file_format='fof_sub_subfind')
SubhaloLenType,o=arepo_package.get_group_property(basePath,'GroupLenType',desired_redshift,file_format='fof_sub_subfind')
SubhaloStellarMass=SubhaloMassType[:,4]

mask=SubhaloBHMass>0

SubhaloIndices=numpy.arange(0,len(SubhaloBHMass))

SubhaloBHIndices=SubhaloIndices[SubhaloBHMass>0]
SubhaloBHStellarMass=SubhaloStellarMass[SubhaloBHMass>0]
SubhaloBHBHMass=SubhaloBHMass[SubhaloBHMass>0]

SubhaloIndex=SubhaloIndices[0]
BH_Mass,o=arepo_package.get_particle_property(basePath,'BH_Mass',5,desired_redshift)
BH_Progs,o=arepo_package.get_particle_property(basePath,'BH_Progs',5,desired_redshift)
BH_ID,o=arepo_package.get_particle_property(basePath,'ParticleIDs',5,desired_redshift)

SubhaloBHLen=SubhaloLenType[:,5]
SubhaloBHBHLen=SubhaloBHLen[SubhaloBHMass>0]

SM_space=[]
BH_ID_space=[]
BH_Progs_space=[]
BH_Mass_space=[]

start=0
for length,SM in list(zip(SubhaloBHLen,SubhaloStellarMass)):
    BH_Mass_Subhalo=BH_Mass[start:start+length]
    BH_Progs_Subhalo=BH_Progs[start:start+length]
    BH_ID_Subhalo=BH_ID[start:start+length]
    start=start+length
    if(length==0):
        continue
    elif(length==1):
        SM_space.append(SM)
        #print(BH_ID_Subhalo)
        BH_ID_space.append(BH_ID_Subhalo[0])
        BH_Mass_space.append(BH_Mass_Subhalo[0])
        BH_Progs_space.append(BH_Progs_Subhalo[0])
    elif(length>1):
        SM_space.append(SM)
        BH_ID_space.append(BH_ID_Subhalo[BH_Mass_Subhalo==numpy.amax(BH_Mass_Subhalo)][0])
        BH_Progs_space.append(BH_Progs_Subhalo[BH_Mass_Subhalo==numpy.amax(BH_Mass_Subhalo)][0])
        BH_Mass_space.append(BH_Mass_Subhalo[BH_Mass_Subhalo==numpy.amax(BH_Mass_Subhalo)][0])
        #BH_Mass_space.append(sum(BH_Mass_Subhalo))
    
SM_space=numpy.array(SM_space) 
BH_ID_space=numpy.array(BH_ID_space) 
BH_Mass_space=numpy.array(BH_Mass_space) 
BH_Progs_space=numpy.array(BH_Progs_space) 
    




'''
new_directory='postprocessing_mergers'
fixed_delay_Myr=125
'''


#BH_Mass_approx=[]
BH_Progs_approx=[]
BH_Progs_approx_delayed=[]
#BH_Mass_approx2=[]
merger_redshifts,new_merger_redshifts,primary_id_sorted2,secondary_id_sorted2,primary_mass_sorted2,secondary_mass_sorted2,event_multiplicities,event_indices=numpy.load(basePath + '/' + new_directory + '/merger_events_upto_redshift_%d_time_delay_%d_Myr.npy'%(desired_redshift,fixed_delay_Myr))
event_indices=event_indices.astype(int)
BH_indices=numpy.arange(0,len(BH_ID))

#seedmass=1.56e-7

for cur_ID in BH_ID_space[:]:
    cur_index=BH_indices[cur_ID==BH_ID][0]
    #print(cur_index)
    
    indices_of_included_events=numpy.load(basePath+'/merger_progenitors/'+'progenitor_merger_indices_z%d_%d.npy'%(desired_redshift,cur_index))
    indices_of_included_events_r=indices_of_included_events[1:]
    #print(indices_of_included_events)
    
    all_positions=numpy.arange(0,len(event_indices))
    positions =  numpy.array([(all_positions[event_indices==ii])[0] for ii in indices_of_included_events_r])

 
    m2_sorted=numpy.amin([primary_mass_sorted2,secondary_mass_sorted2],axis=0)
    m1_sorted=numpy.amax([primary_mass_sorted2,secondary_mass_sorted2],axis=0)
    
    try:
        m2_final=m2_sorted[positions]
        m1_final=m1_sorted[positions]
        merger_redshifts_final=merger_redshifts[positions]
        new_merger_redshifts_final=new_merger_redshifts[positions]
        mask3=merger_redshifts_final>=o
        mask4=new_merger_redshifts_final>=o
        BH_Progs_approx.append(len(m2_final[mask3]))
        
        BH_Progs_approx_delayed.append(len(m2_final[mask4]))
    except IndexError:
        BH_Progs_approx.append(0)
        BH_Progs_approx_delayed.append(0)

    
    
    
    
    #BH_Mass_approx2.append((len(m2_sorted[mask1 & mask3]) + len(m2_sorted[mask2 & mask3]))*seedmass)
    
#BH_Mass_approx=numpy.array(BH_Mass_approx)
BH_Progs_approx=numpy.array(BH_Progs_approx)+1
BH_Progs_approx_delayed=numpy.array(BH_Progs_approx_delayed)+1


#BH_Mass_approx=BH_Progs_approx*seedmass
    
    
numpy.save(basePath+'/' + new_directory + '/new_BH_populations_redshift%d_timedelayMyr%d.npy'%(desired_redshift,fixed_delay_Myr),[SM_space,BH_ID_space,BH_Mass_space,BH_Progs_space,BH_Progs_approx,BH_Progs_approx_delayed])


