#-----------------------This code takes a set of objects and their positions and groups them based on Friends of Friends linking length criterion
#-----------------------------------import required packages--------------------------------------------------------------------
import sys
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
#sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import numpy
import time
import h5py
import pickle
import sys
import arepo_package
sys.setrecursionlimit(20000)
PROJECTION_MODE=True
SPECTROSCOPIC_MODE=True

run='L75n1820TNG'
basePath='/n/hernquistfs3/IllustrisTNG/Runs/'+run+'/output/'


#basePath='/ufrc/lblecha/aklantbhowmick/arepo_runs_aklant/L25_n256/output/'
#--------------------------------------------------------------------------------------------

BOXSIZE=arepo_package.get_box_size(basePath)/1e3    # Simulation Box size: Default value set for Massive Black II. Note:- Periodic boundary conditions apply 

global FOF_tag  # Array of friends of freinds tag assigned as a global variable and can be modified by the various  functions called within the code 

def min_dis(x1,x2,BOXSIZE=BOXSIZE): 
        xdis_option1=(x1-x2)*(x1-x2)
        #print xdis_option1
        xdis_option2=(x1-x2+BOXSIZE)*(x1-x2+BOXSIZE)
        xdis_option3=(x1-x2-BOXSIZE)*(x1-x2-BOXSIZE)
        return numpy.amin(numpy.array([xdis_option1,xdis_option2,xdis_option3]),axis=0)

def calc_dis_square(pos1,pos2): # Calculates the distance between the two bjects in the box using periodic boundary condition   

        xdis2 = min_dis(pos1[0],pos2[0])
        ydis2 = min_dis(pos1[1],pos2[1])
        zdis2 = min_dis(pos1[2],pos2[2])
        calc_dis = xdis2+ydis2+zdis2
        return calc_dis
    
def calc_dis_square_vectorized(pos1,pos2):# Calculates the distance between the two objects in the box using periodic boundary condition, the second argument can be a vector   

        
        xdis2 = min_dis(pos1[0],pos2[:,0])
        ydis2 = min_dis(pos1[1],pos2[:,1])
        zdis2 = min_dis(pos1[2],pos2[:,2])
        
        #print xdis2
        if (PROJECTION=='z'):
            par_dis=xdis2+ydis2
            los_dis=zdis2
        if (PROJECTION=='x'):
            par_dis=ydis2+zdis2
            los_dis=xdis2
        if (PROJECTION=='y'):
            par_dis=xdis2+zdis2
            los_dis=ydis2
        if (PROJECTION_MODE):    
            return par_dis,los_dis
        else:
            return par_dis+los_dis
    
def make_cuts(quantities,cut): #selects an array of quantities (argument 1) and makes cuts (argument 2)
    cutted_quantities=[quantity[cut] for quantity in quantities]
    return cutted_quantities

def get_indices_to_place_tag(primary_position,other_positions_untagged): # takes an object's position and identofies the neighbours
            #print other_positions_untagged
            if (PROJECTION_MODE):
                par_distances_to_other_positions,los_distances_to_other_positions=numpy.sqrt(calc_dis_square_vectorized(primary_position,other_positions_untagged))
                los_LINKING_LENGTH=10**los_log_LINKING_LENGTH
            else:
                distances_to_other_positions=numpy.sqrt(calc_dis_square_vectorized(primary_position,other_positions_untagged))

            LINKING_LENGTH=10**log_LINKING_LENGTH
            
            if (TEST_MODE):
                LINKING_LENGTH=TEST_LENGTH+0.01
                los_LINKING_LENGTH=TEST_LENGTH_los+0.01
            
            if (PROJECTION_MODE): 
                extract_neighbors=(par_distances_to_other_positions<=LINKING_LENGTH)& (los_distances_to_other_positions<=los_LINKING_LENGTH)             
                extract_neighbors_but_self=extract_neighbors&(par_distances_to_other_positions+los_distances_to_other_positions>0.)
            else:
                extract_neighbors=distances_to_other_positions<=LINKING_LENGTH
                extract_neighbors_but_self=extract_neighbors&(distances_to_other_positions>0.)
            return extract_neighbors,extract_neighbors_but_self # returns a boolean array which can used to extract the array indices of the neighbors

def generate_tags(log_bhmass_cut,redshift,log_LINKING_LENGTH,START): 
# This is the main function of the code, it selects objects based on luminosities and redshifts and places the tags on the various friends of freinds groups
    global FOF_tag
    def tag_neighbors(primary_object_position): # This is function which is recursively called upon to tag the neighbors 
        other_positions_untagged,indices_untagged=(make_cuts([other_positions,indices],FOF_tag==-1)) # Find the indices of untagged objects
        
        extract_neighbors,extract_neighbors_but_self=get_indices_to_place_tag(primary_object_position,other_positions_untagged)  # Determine the indices the place the tag
        indices_to_place_tag=indices_untagged[extract_neighbors]
        neighbor_positions=other_positions_untagged[extract_neighbors_but_self]
        FOF_tag[indices_to_place_tag]=current_tag # place the current tag on the index positions
        if (FOF):  # IF FOF mode is activated, find the neighbors of neighbors too  
            if (len(neighbor_positions)>0): #recursively calling the function to iteratively finding the neighbors at each step
                #print "There are",len(neighbor_positions) ,"neighbors"
                for neighbor_position_primary in neighbor_positions:
                    tag_neighbors(neighbor_position_primary)  
            else:
                return #return back if no- more neighbours are left to be tagged in this branch
            
#-------------------------This whole section reads the properties of objects----------------------------------          
    #g=h5py.File('/nfs/nas-0-1/akbhowmi/quasar_properties/bh_lum_host_halo_with_id_and_mass_z%.2f_and_sub_halo_id_added_velocity'%redshift)
    #g2=h5py.File('/nfs/nas-0-1/akbhowmi/quasar_properties/bh_lum_host_halo_with_id_and_mass_z%.2f_and_halo_id_added_velocity'%redshift)


#-------------------------This section performs cuts based on what the user desires----------------------------------     
#    mask_luminosity=(bolometric_luminosity>10**log_luminosity_cuts)
    mask_bhmass=(bhmass>10**log_bhmass_cut)


    subsampled_quantities=make_cuts([position_vector,velocity_vector_gadget_units,blackhole_id,subhalo_id,hosthalo_id],mask_bhmass)
    position_vector_cut,velocity_vector_gadget_units_cut,blackhole_id_cut,subhalo_id_cut,hosthalo_id_cut=subsampled_quantities
    if (TEST_MODE): #--------Stack of points for testing the codes-----------------------------------------------
        position_vector_cut1=numpy.array([numpy.array([0.,i*TEST_LENGTH,0.]) for i in numpy.arange(1,10).astype('float')]) 
        position_vector_cut2=numpy.array([numpy.array([0.,i*TEST_LENGTH,1.]) for i in numpy.arange(1,10).astype('float')]) 
        position_vector_cut3=numpy.array([numpy.array([0.,i*TEST_LENGTH,200.]) for i in numpy.arange(1,10).astype('float')]) 
        position_vector_cut4=numpy.array([numpy.array([0.,i*TEST_LENGTH,201.]) for i in numpy.arange(1,10).astype('float')]) 

        position_vector_cut=numpy.append(position_vector_cut1,position_vector_cut2,axis=0)
        position_vector_cut=numpy.append(position_vector_cut,position_vector_cut3,axis=0)
        position_vector_cut=numpy.append(position_vector_cut,position_vector_cut4,axis=0)
#-------------------------------------------------------------------------------------------------------
    
    FOF_tag=numpy.array([-1]*len(position_vector_cut)) # Initializing the FOF_tag array
    
    start=time.time()

    original_number=len(position_vector_cut)
    print("Total number of blackholes to start with: ",original_number)

    MAX_ITERATIONS=50000.

    
    other_positions=position_vector_cut  # Assigning 'other_positions' as all the available positions 
    other_positions_untagged=position_vector_cut  # Initializing 'other_positions_untagged' as all the available positions 
    indices=numpy.arange(0,len(other_positions))
    indices_untagged=indices
    
    i=0
    while i<MAX_ITERATIONS: # Continue until the maximum number number of iterations, if that happens, something is probably wrong

        current_tag=i       # Current tag to placed. Each iteration identifies an FOF group and places a unique tag
        other_positions_untagged,indices_untagged=(make_cuts([other_positions,indices],FOF_tag==-1)) # Find the indices of untagged objects

        if (len(indices_untagged)==0):
            break       # If all tags have been placed, exit the loop. the calculation is done
        if(START=='LAST'):# This decides which object in the array do we start with. Not too important
            IP=-1
        else:
            IP=0
        
        primary_object_position=other_positions_untagged[IP] #Amongst the untagged objects, choose any primary object
        tag_neighbors(primary_object_position)# place the current tag on the index positions

                       
        i+=1
        if (i==MAX_ITERATIONS):  #Exit if maximum number of iterations have been reached
            print("Warning: Maximum number of iterations reached!!")
            break
            
        if (i%1000==0):
            print("%d groups identified"%(i))
        
    end=time.time()

    print("Time taken in seconds:", end-start)
    return FOF_tag,subsampled_quantities  # returns the FOF tags of the objects along with all the subsampled quantities
    


omk=0
oml=0.725
omm=0.275
cc=3*10**5
H0=100.
vmax=600.

PROJECTION_MODE=0
SPECTROSCOPIC_MODE=0
line_of_sight_squashing_factor=1.
if ((PROJECTION_MODE==1)&(SPECTROSCOPIC_MODE==1)):
    print("Both modes cannot be on at the same time")
    exit()
    

PROJECTION='y'

TEST_MODE=0#--------Activate test mode if ON.  
TEST_LENGTH=0.3 
TEST_LENGTH_los=10000.
#--------Testing linking length if ON. See construction of points for details


FOF=1   # Activate if you want Friends of Friends
redshift_space=[4,3,2,1,0]
log_LINKING_LENGTH_space=numpy.linspace(-2,2,40)
log_luminosity_cuts_space=[43.0,43.5,44.0,44.5,45.0,45.5]
log_bhmass_cuts_space=[8,7,6]
#log_luminosity_cuts_space=[42.5,42.0]

output_path_for_FOF_tags='/n/home00/aklantbhowmick/illustris_processing/FOF_tag_outputs/'
output_path_for_RICHNESS='/n/home00/aklantbhowmick/illustris_processing/richness_outputs/'

#run='L75n1820TNG'
#basePath='/n/hernquistfs3/IllustrisTNG/Runs/'+run+'/output/'













for redshift in redshift_space:
    p_type=5
    desired_redshift=redshift    
    output_redshift,output_snap=arepo_package.desired_redshift_to_output_redshift(basePath,desired_redshift)
 
    print('Reading positions')   
    Coordinates,output_redshift=arepo_package.get_particle_property(basePath, 'Coordinates', p_type,desired_redshift, list_all=False)
    pos_x=Coordinates[:,0]/1e3
    pos_y=Coordinates[:,1]/1e3
    pos_z=Coordinates[:,2]/1e3
    #print(pos_x)
    print('Reading velocities')
    Velocities,output_redshift=arepo_package.get_particle_property(basePath, 'Velocities', p_type,desired_redshift, list_all=False)

    vel_x=Velocities[:,0]
    vel_y=Velocities[:,1]
    vel_z=Velocities[:,2]
    if (PROJECTION_MODE | SPECTROSCOPIC_MODE):
        print("Warning: Performing redshift space distortions!!!")
        if (PROJECTION=='x'):
                pos_x=(pos_x+vel_x*numpy.sqrt(scale)/(scale*H))*line_of_sight_squashing_factor
        if (PROJECTION=='y'):
                pos_y=(pos_y+vel_y*numpy.sqrt(scale)/(scale*H))*line_of_sight_squashing_factor
        if (PROJECTION=='z'):
                pos_z=(pos_z+vel_z*numpy.sqrt(scale)/(scale*H))*line_of_sight_squashing_factor

    print('Reading particle ids')
    blackhole_id,output_redshift=arepo_package.get_particle_property(basePath, 'ParticleIDs', p_type,desired_redshift, list_all=False)
    print('Reading bh masses')
    bhmass,output_redshift=arepo_package.get_particle_property(basePath, 'BH_Mass', p_type,desired_redshift, list_all=False)
    bhmass=bhmass*1e10
    print('Reading/generating group ids')
    hosthalo_id=arepo_package.generate_group_ids(basePath,desired_redshift,5,save_output_path='./'+run+'_group_ids/')
    subhalo_id=hosthalo_id
    print('linking blackholes')
    position_vector=numpy.transpose([pos_x,pos_y,pos_z])
    velocity_vector_gadget_units=numpy.transpose([vel_x,vel_y,vel_z])
    
    H=H0*numpy.sqrt((1.+redshift)**3*omm+oml)
    scale=1./(1.+redshift)
    window=2.*vmax/H/scale
    los_log_LINKING_LENGTH=numpy.log10(window)
    print("Window width in Mpc",window)
    
    for log_bhmass_cut in log_bhmass_cuts_space:
        for log_LINKING_LENGTH in log_LINKING_LENGTH_space:
            print("-----------------------------------------")
            print("Outputing for redshift ", redshift,"; log_bhmass_cut ",log_bhmass_cut,"; log_LINKING_LENGTH ",log_LINKING_LENGTH)
            FOF=1
            final_FOF_tag,subsampled_quantities=generate_tags(log_bhmass_cut,redshift,log_LINKING_LENGTH,'FIRST')
            #if (~TEST_MODE):
            #    pickle.dump([final_FOF_tag,subsampled_quantities],open('./FOF_tag_outputs/log_luminosity_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f.pickle'%(log_luminosity_cuts,redshift,log_LINKING_LENGTH),'w'))
           
            if (PROJECTION_MODE):
                pickle.dump([final_FOF_tag,subsampled_quantities],open('/nfs/nas-0-1/akbhowmi/FOF_tag_outputs_bhmass_cuts/log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f_projected_%s_window_%.1f_new_run.pickle'%(log_bhmass_cut,redshift,log_LINKING_LENGTH,PROJECTION,vmax),'w'))
            elif (SPECTROSCOPIC_MODE):
                pickle.dump([final_FOF_tag,subsampled_quantities],open('/nfs/nas-0-1/akbhowmi/FOF_tag_outputs_bhmass_cuts/log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f_projected_%s_window_%.1f_spectroscopic_squashing_%.2f.pickle'%(log_bhmass_cut,redshift,log_LINKING_LENGTH,PROJECTION,vmax,line_of_sight_squashing_factor),'w'))
            else:
                numpy.save(output_path_for_FOF_tags+'%s_log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f.npy'%(run,log_bhmass_cut,redshift,log_LINKING_LENGTH),[final_FOF_tag,subsampled_quantities])
                #print(final_FOF_tag)
                #print("------------")
                
                
                
                final_FOF_tag,subsampled_quantities=numpy.load(output_path_for_FOF_tags+'%s_log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f.npy'%(run,log_bhmass_cut,redshift,log_LINKING_LENGTH))
                position_vector_cut,velocity_vector_gadget_units_cut,blackhole_id_cut,subhalo_id_cut,hosthalo_id_cut=subsampled_quantities
                unique_FOF_tag=numpy.unique(final_FOF_tag)
                RICHNESS=[]
                for tag in unique_FOF_tag:
                    extract_FOFs=(final_FOF_tag==tag)
                    
                    #if (len(final_FOF_tag[extract_FOFs])==5)&(log_luminosity_cuts==44.0):
                    #    print 10**log_LINKING_LENGTH
                    #    print "---------------"
                     #   print hosthalo_id_cut[extract_FOFs]
                     #   print "---"
                     #   print subhalo_id_cut[extract_FOFs]
                        
                        
                    RICHNESS.append(len(final_FOF_tag[extract_FOFs]))
                    
                    #if (len(final_FOF_tag[extract_FOFs])>10):
                    #    print len(final_FOF_tag[extract_FOFs]),tag
                        
                    
                RICHNESS=numpy.array(RICHNESS)
                #print(RICHNESS) 
                numpy.save(output_path_for_RICHNESS+'%s_richness_log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f.npy'%(run,log_bhmass_cut,redshift,log_LINKING_LENGTH),RICHNESS)

                
                #RICHNESS=numpy.load(output_path_for_RICHNESS+'%s_richness_log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f.npy'%(log_bhmass_cut,redshift,log_LINKING_LENGTH))
                #if(log_LINKING_LENGTH==-1.794871794871795):
                #    print(RICHNESS)    
                #    break
                
                


