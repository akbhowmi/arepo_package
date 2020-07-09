#-------------------------------------------------------Load all necessary packages-------------------------------------
import sys
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import mdot_to_Lbol
import arepo_package
import scipy.interpolate
import illustris_python
#%pylab inline
import h5py
import illustris_python.sublink
import numpy
import pickle
#-----------------------------------------------------------------------------------------------------------------------
#------------------Class definition for the binary tree used to store the progenitor info------------------------------ 
class Subhalo:
    def __init__(self):
        self.Index=-1
        self.MostMassiveProgenitor = -1
        self.NextMostMassiveProgenitor = -1
        self.Snap=-1 
#-----------------------------------------------------------------------------------------------------------------------
#----------------------------------This function fills up the progenitor tree------------------------------------------               
def function_fill_progenitor_tree(subhalo_index,currentsubhalo,current_subhalo_ID):
    print(current_subhalo_ID)
    fetch_subhalo=current_subhalo_ID==SubhaloID_Tree
    currentsubhalo.Index=SubfindID_Tree[fetch_subhalo][0]
    currentsubhalo.Snap=SnapNum_Tree[fetch_subhalo][0]
    MostMassiveProgenitorID=FirstProgenitorID_Tree[fetch_subhalo][0]
    NextMostMassiveProgenitorID=NextProgenitorID_Tree[fetch_subhalo][0]
     
    MostMassiveProgenitorIndex=SubfindID_Tree[MostMassiveProgenitorID==SubhaloID_Tree]
    NextMostMassiveProgenitorIndex=SubfindID_Tree[NextMostMassiveProgenitorID==SubhaloID_Tree]
    
    
    currentsubhalo.MostMassiveProgenitor=Subhalo()
    currentsubhalo.NextMostMassiveProgenitor=Subhalo()

    
    if (MostMassiveProgenitorID!=-1):
        print("Following the most massive progenitor branch")
        function_fill_progenitor_tree(MostMassiveProgenitorIndex,currentsubhalo.MostMassiveProgenitor,MostMassiveProgenitorID)
    if (NextMostMassiveProgenitorID!=-1):
        print("Following the next most massive progenitor branch")
        function_fill_progenitor_tree(NextMostMassiveProgenitorIndex,currentsubhalo.NextMostMassiveProgenitor,NextMostMassiveProgenitorID)    
    print("Found no more progenitors: returning to previous recursion")
    return 
#-------------------------------------------------------------Simulation path -------------------------------------------------------
path_to_simulation='/blue/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS/'
run='density_and_metallicity_based_criterion'
basePath=path_to_simulation+run+'/AREPO/output/'
tree=h5py.File(basePath+'/postprocessing/tree_extended.hdf5','r')
save_output_path='/home/aklantbhowmick/Aklant/arepo_code_development/progenitor_outputs/'
#------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------Details of the root subhalo-------------------------------------------------
subhalo_index=1
desired_redshift=0
#------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------Extracting all the necessary data from the sublink catalogs and snapshots----
output_redshift,output_snapshot=arepo_package.desired_redshift_to_output_redshift(basePath,desired_redshift)
SubfindID=tree.get('SubfindID')[:]
SubhaloID=tree.get('SubhaloID')[:]
FirstProgenitorID=tree.get('FirstProgenitorID')[:]
NextProgenitorID=tree.get('NextProgenitorID')[:]
TreeID=tree.get('TreeID')[:]
SnapNum=tree.get('SnapNum')[:]
#------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------Extracting the tree to which the selected subhalo belongs-----------------------------------
find_original_subhalo=(SubfindID==subhalo_index)&(SnapNum==output_snapshot)
if (len(SubfindID[find_original_subhalo])!=1):
    print("Warning: The number of selected subhaloes must be 1")
Target_TreeID=TreeID[find_original_subhalo]
extract_tree=Target_TreeID==TreeID
SubhaloID_Tree=SubhaloID[extract_tree]
SnapNum_Tree=SnapNum[extract_tree]
FirstProgenitorID_Tree=FirstProgenitorID[extract_tree]
NextProgenitorID_Tree=NextProgenitorID[extract_tree]
SubfindID_Tree=SubfindID[extract_tree]
#-------------------------------------Finding the corresponding unique identifier of the selected subhalo-----------------------------------
find_subhalo_on_tree=(SubfindID_Tree==subhalo_index)&(SnapNum_Tree==output_snapshot)
RootSubhaloID=SubhaloID_Tree[find_subhalo_on_tree][0]
#-----Setting recursion limit: If the recursion limit is reached, try increasing this value, but that might also be a hint of a bug--------------------------------------------------------------------
sys.setrecursionlimit(1000)
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------Initializing the root subhalo of the progenitor tree and calling the function to extract all progenitors ----------------------------------------
rootsubhalo = Subhalo()
function_fill_progenitor_tree(subhalo_index,rootsubhalo,RootSubhaloID)
#----------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------Saving the output--------------------------------------------------------------------------------
#pickle.dump(open(save_output_path+'/%s_index_%d_redshift_%.3f.pickle'%(run,subhalo_index,desired_redshift),'w'),rootsubhalo)
numpy.save(save_output_path+'/%s_index_%d_redshift_%.3f.pickle'%(run,subhalo_index,desired_redshift),[rootsubhalo])

#-----------------------------------------------------------------------------------------------------------------------------



    
