#-------------------------------------------------Loading the necessary packages----------------------------------
import sys
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import mdot_to_Lbol
import arepo_package
import scipy.interpolate
import illustris_python
import h5py
import illustris_python.sublink
import numpy
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#plt.use('Agg')
#-----------------------------------------------------------------------------------------------------------------
path_to_simulation='/blue/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS/'
run='density_and_metallicity_based_criterion'
basePath=path_to_simulation+run+'/AREPO/output/'

subhalo_index=0
desired_redshift=10
load_output_path='/home/aklantbhowmick/Aklant/arepo_code_development/descendant_outputs/'
SubfindID_descendants,SnapNum_descendants=numpy.load(load_output_path+'/%s_index_%d_redshift_%.3f.npy'%(run,subhalo_index,desired_redshift))

f,axx=plt.subplots(6,3,figsize=(30,50))
snap_list,redshift_list=arepo_package.get_snapshot_redshift_correspondence(basePath)
j=0
i=0
for snap,index in zip(SnapNum_descendants,SubfindID_descendants):
    ax=axx[i,j]
    desired_redshift=redshift_list[snap_list==snap][0]
    ax.set_title('$z=%.1f$'%desired_redshift,fontsize=30)    
    DM_particle_positions_subhalo,DM_particle_positions_group,output_redshift=arepo_package.get_particle_property_within_groups(basePath,'Coordinates',1,desired_redshift,index,group_type='subhalo')
    #DM_particle_positions_subhalo/=1e3
    #DM_particle_positions_group/=1e3
    nostars=0
    try:
        nostars=1
        star_particle_positions_subhalo,star_particle_positions_group,output_redshift=arepo_package.get_particle_property_within_groups(basePath,'Coordinates',4,desired_redshift,index,group_type='subhalo')
        #star_particle_positions_subhalo/=1e3
        #star_particle_positions_group/=1e3
    except:
        nostars=0
        print("No stars")
    NBINS=200

    
    ax.hist2d(DM_particle_positions_group[:,1],DM_particle_positions_group[:,2], bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap='Greys_r');
    if (nostars):
        ax.hist2d(star_particle_positions_subhalo[:,1],star_particle_positions_subhalo[:,2], bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap='Reds_r');
        ax.hist2d(DM_particle_positions_group[:,1],DM_particle_positions_group[:,2], bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap='Reds_r',alpha=0);
    ax.tick_params(labelsize=30)
    j+=1   
    print(i,j)
    if (j==3):
        j=0
        i+=1

plt.savefig(load_output_path+'/%s_index_%d_redshift_%.3f_image.png'%(run,subhalo_index,desired_redshift),bbox_inches='tight')

