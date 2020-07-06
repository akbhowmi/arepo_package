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

#-------------------------------------Key inputs-------------------------------------------------------------------
path_to_simulation='/blue/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS/'
run='density_and_metallicity_based_criterion'
basePath=path_to_simulation+run+'/AREPO/output/'
subhalo_index=100
desired_redshift=5
output_redshift,output_snapshot=arepo_package.desired_redshift_to_output_redshift(basePath,desired_redshift)
save_output_path='/home/aklantbhowmick/Aklant/arepo_code_development/descendant_outputs/'
#------------------------------------------------------------------------------------------------------------------
#-------------------------------------Loading the merger tree------------------------------------------------------
tree=h5py.File(basePath+'/postprocessing/tree_extended.hdf5','r')
SubfindID=tree.get('SubfindID')[:]
SubhaloID=tree.get('SubhaloID')[:]
DescendantID=tree.get('DescendantID')[:]
TreeID=tree.get('TreeID')[:]
SnapNum=tree.get('SnapNum')[:]
#------------------------------------------------------------------------------------------------------------------

#-----------Find the selected subhalo in the tree data, and extracting the corresponding tree to which it belongs-------
find_original_subhalo=(SubfindID==subhalo_index)&(SnapNum==output_snapshot)
if (len(SubfindID[find_original_subhalo])!=1):
    print("Warning: The number of selected subhaloes must be 1")
Target_TreeID=TreeID[find_original_subhalo]
extract_tree=Target_TreeID==TreeID
SubhaloID_Tree=SubhaloID[extract_tree]
SnapNum_Tree=SnapNum[extract_tree]
DescendantID_Tree=DescendantID[extract_tree]
SubfindID_Tree=SubfindID[extract_tree]
find_subhalo_on_tree=(SubfindID_Tree==subhalo_index)&(SnapNum_Tree==output_snapshot)
#----------------------------------------------------------------------------------------------------------------------

#----------Initializing the arrays that'll save the descendant information, starting with the originally selected subhalo------
SubfindID_descendants=[SubfindID[find_original_subhalo][0]]
SnapNum_descendants=[SnapNum[find_original_subhalo][0]]
#------------------------------------------------------------------------------------------------------------------------------

#------------------------Auxilliary variables to navigate the while loop, change if needed-------------------------------------
MAX_ITERATIONS=100
i=0
Last_Snap_Num=200
#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------Loop that computes the all the descendants------------------------------------------------------------
while (i<MAX_ITERATIONS):
    Subhalo_ID_tracking=DescendantID_Tree[find_subhalo_on_tree]
    if (Subhalo_ID_tracking[0]==[-1]):
        break
    fetch_subhalo=SubhaloID_Tree==Subhalo_ID_tracking 
    SubfindID_descendants.append(SubfindID_Tree[fetch_subhalo][0])
    SnapNum_descendants.append(SnapNum_Tree[fetch_subhalo][0])
    find_subhalo_on_tree=fetch_subhalo   
    print(i)
    i+=1
#-----------------------------------------------------------------------------------------------------------------------------
#----------------------------Saving the output--------------------------------------------------------------------------------
numpy.save(save_output_path+'/%s_index_%d_redshift_%.3f.npy'%(run,subhalo_index,desired_redshift),[SubfindID_descendants,SnapNum_descendants])
#-----------------------------------------------------------------------------------------------------------------------------

