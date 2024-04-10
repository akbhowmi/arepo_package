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
import numpy as np
import pandas as pd
import scipy.integrate
from scipy.optimize import curve_fit
import os
from Config import *


import warnings
warnings.filterwarnings("ignore")
from scipy.integrate import quad
'''
omk=0
oml=0.6911
omm=0.3089
c=3*10**5
h=0.6771
'''

def DC(z0, z1):
    # Comoving distance in Mpc                                                                                          $
    def integrand(z):
        return 1.0/E(z)
    return c/(100.0)*quad(integrand,z0,z1)[0]

def dx_dz(z):
    # Comoving distance in Mpc                                                                                          $
    def integrand(z):
        return 1.0/E(z)
    return c/(100.0)*integrand(z)


def E(z):
        #normalized hubble parameter                                                                                    $
    return np.sqrt(omm*(1.0+z)**3 + oml + omk*(1.0+z)**2)

def DL(z):
    #Luminosity distance in Mpcs                                                                                        $
    return (1.0+z)*DC(0,z)

def DM(z):
        #Distance Modulus                                                                                               $
        return 5*np.log10((DL(z)*1e6)/10.0)

def mtoM(m,z):
    return m - DM(z) - K(z)

def Mtom(M, z):
    return M + DM(z) + K(z)

def LtoM(L):
    M = -2.5 * np.log10(L*1e-7)+34.1
    return M

def MtoL(M):
    return 10**((-1/2.5*(M-34.1)))*1e7

def K(z):
    #K correction, assuming same filter size in emitted and observed frame                                              $
    return 2.5*np.log10(1.0/(1+z))

def T(z0,z1):
    sectoyrs=3.171e-08
    Mpctokm=3.0857E+19
    def integrand(z):
        return 1.0/E(z)/(1+z)
    return 1./(100.0)*quad(integrand,z0,z1)[0]*Mpctokm*sectoyrs

def update_events(events,previous_id1, previous_id2, previous_mass1, previous_mass2, M,actual_merger_time_of_previous_binary):
        global event_multiplicities, new_merger_times, event_processed,recursion_depth
        global triple_id1s, triple_id2s, triple_id3s, triple_mass1s, triple_mass2s, triple_mass3s, triple_formation_time, triple_event_sorted_index
        if(recursion_depth==0):
            id_1, id_2, mass_1, mass_2, time, event_index = events[i_event]
        else:
            id_1, id_2, mass_1, mass_2, time, event_index = events[0]


        if(M>2):   #This means that previous iteration(recursion) identified this event as a triple or higher multiple
            triple_id1s.append(previous_id1)
            triple_id2s.append(previous_id2)
            triple_id3s.append(unretained_id)

            triple_formation_time.append(time)

            triple_event_sorted_index.append(event_index)

            triple_mass1s.append(previous_mass1)
            triple_mass2s.append(previous_mass2)  
            triple_mass3s.append(unretained_mass)
            Ms.append(M)
        #Let us first apply the usual isolated binary delay time.     
        if (FIXED_DELAY==1):    
            new_time = time + fixed_delay 

        elif(RANDOM_DELAY==1):
            if(recursion_depth==0):
                new_time = time + random_delay[i_event]
            else:
                new_time = time + random_delay[0]
        #Update the merger time if the binary would encounter a third blackhole from the next 'simulated' mergere  
        if(new_time > next_merger_time_of_retained_id):
            new_time = next_merger_time_of_retained_id + triple_induced_delay
            #if this event was already a multiple, let's reduce the multiplicity if the two of the BHs from that triple would merge before
            if(M>2):
              if(new_time > actual_merger_time_of_previous_binary):
                #print("We have subtracted multiplicity")
                M-=1 
            M+=1

        new_merger_times[event_index]=new_time
        # Find subsequent events with the same IDs to determine the retained ID
        subsequent_events_id1 = [(i1, i2, m1, m2, t, index) for i1, i2, m1, m2,  t, index in all_events if (i1 == id_1 or i2 == id_1) and t > time]
        subsequent_events_id2 = [(i1, i2, m1, m2, t, index) for i1, i2, m2, m2,  t, index in all_events if (i1 == id_2 or i2 == id_2) and t > time]
        # Determine the retained ID based on subsequent events
        if subsequent_events_id1:
            retained_id = id_1  # If no subsequent events found, retain id_1
            subsequent_events = subsequent_events_id1
            unretained_id = id_2
            unretained_mass = mass_2
        else:
            retained_id = id_2  # If no subsequent events found, retain id_1
            subsequent_events = subsequent_events_id2
            unretained_id = id_1
            unretained_mass = mass_1

        event_processed[event_index]=1
        #print(subsequent_events)
        if not subsequent_events:
            recursion_depth-=1
            return
        
        next_merger_time_of_retained_id = subsequent_events[0][4]
        next_event_index_of_retained_id = subsequent_events[0][5]

        next_merger_id1_of_retained_id = subsequent_events[0][0]
        next_merger_id2_of_retained_id = subsequent_events[0][1]


        #print(recursion_depth,new_time,actual_merger_time_of_previous_binary)
        # Calculate new time by adding the fixed delay 'dt'

        recursion_depth+=1    
        update_events(subsequent_events,id_1, id_2, mass_1, mass_2, M,new_time) 
        


global triple_id1s, triple_id2s, triple_id3s, triple_mass1s, triple_mass2s, triple_mass3s, triple_formation_time, triple_event_sorted_index
triple_id1s=[]
triple_id2s=[]
triple_id3s=[]

triple_mass1s=[]
triple_mass2s=[]
triple_mass3s=[]
Ms=[]

triple_formation_time=[]
triple_event_sorted_index=[]

redshift_grid=numpy.linspace(0,25,50)
time_grid=numpy.array([(T(0,110000)-T(0,redshift))/h/1e9 for redshift in redshift_grid])
gen_log_time = scipy.interpolate.interp1d(redshift_grid,time_grid,fill_value='extrapolate')
gen_redshift = scipy.interpolate.interp1d(time_grid,redshift_grid,fill_value='extrapolate')


#path_to_output='/orange/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS///' # this is the folder containing the simulation run
#run='/L3p125n512/AREPO/' # name of the simulation runs
#output='/output_ratio1000_SFMFGM5_seed3.19_bFOF/' 
#basePath=path_to_output+run+output



upto_redshift=desired_redshift
scale_fac_complete_sorted,primary_mass_sorted,secondary_mass_sorted,primary_id_sorted,secondary_id_sorted,file_id_complete_sorted,N_empty=arepo_package.get_merger_events_from_snapshot(basePath,upto_redshift,HOSTS=0)


#scale_fac_complete_sorted,primary_mass_sorted,secondary_mass_sorted,primary_id_sorted,secondary_id_sorted,file_id_complete_sorted,N_empty=arepo_package.get_merger_events(basePath)

event_index=numpy.arange(0,len(scale_fac_complete_sorted))

event_index_sorted2 = arepo_package.sort_X_based_on_Y(event_index,scale_fac_complete_sorted)

scale_fac_complete_sorted2=scale_fac_complete_sorted[event_index_sorted2]
primary_id_sorted2=primary_id_sorted[event_index_sorted2]
secondary_id_sorted2=secondary_id_sorted[event_index_sorted2]
primary_mass_sorted2=primary_mass_sorted[event_index_sorted2]
secondary_mass_sorted2=secondary_mass_sorted[event_index_sorted2]


#FIXED_DELAY=1
#RANDOM_DELAY=0

global event_multiplicities, new_merger_times, event_processed,recursion_depth
merger_redshifts=1./scale_fac_complete_sorted2-1.
merger_times=gen_log_time(merger_redshifts)
new_merger_times = merger_times
new_merger_times_triple_induced = merger_times
event_multiplicities = numpy.array([2]*len(merger_times))

event_processed = numpy.array([0]*len(merger_times))
blackhole_id_1 = primary_id_sorted2
blackhole_id_2 = secondary_id_sorted2

blackhole_mass_1 = primary_mass_sorted2
blackhole_mass_2 = secondary_mass_sorted2

event_indices=numpy.arange(0,len(merger_times))

#fixed_delay_Myr=numpy.float(sys.argv[1])

fixed_delay = fixed_delay_Myr/1e3 # Example fixed delay time

#min_delay=0
#max_delay=0
random_delay = numpy.array([numpy.random.uniform(min_delay, max_delay) for _ in range(len(merger_redshifts))])

all_events = list(zip(blackhole_id_1, blackhole_id_2, blackhole_mass_1, blackhole_mass_2 ,merger_times,event_indices))

M=2
max_iterations=100000
actual_merger_time_of_previous_binary=0
for i_event in range(0,len(all_events)): 
    if(event_processed[i_event]==0):
        recursion_depth=0
        update_events(all_events,0,0,0,0,M,actual_merger_time_of_previous_binary)
        event_processed[i_event]=1
    if (len(merger_times[event_processed==0])==0):
        print("finished_processing_all_mergers")
        break
    i_event+=1
    if(i_event>max_iterations):
        print("max_iterations reached")
        break
    





def create_directory_if_not_exists(directory_name):
    # Check if the directory exists
    if not os.path.exists(basePath + directory_name):
        # Create the directory if it doesn't exist
        os.makedirs(basePath + '/' + directory_name)
        print(f"Directory '{directory_name}' created.")
    else:
        print(f"Directory '{directory_name}' already exists.")

# Specify the directory name you want to create
#new_directory = "postprocessing_mergers"  # Replace with your desired directory name

# Call the function to create the directory if it doesn't exist
create_directory_if_not_exists(new_directory)
new_merger_redshifts=gen_redshift(new_merger_times)



triple_id1s = numpy.array(triple_id1s)
triple_id2s = numpy.array(triple_id2s)
triple_id3s = numpy.array(triple_id3s)

triple_mass1s = numpy.array(triple_mass1s)
triple_mass2s = numpy.array(triple_mass2s)
triple_mass3s = numpy.array(triple_mass3s)
Ms = numpy.array(Ms)

triple_formation_time = numpy.array(triple_formation_time)
triple_event_sorted_index = numpy.array(triple_event_sorted_index)






numpy.save(basePath + '/' + new_directory + '/all_events_upto_redshift_%d_time_delay_%d_Myr_triple_%d_Myr.npy'%(desired_redshift,fixed_delay_Myr,triple_induced_delay),[merger_redshifts,new_merger_redshifts,primary_id_sorted2,secondary_id_sorted2,primary_mass_sorted2,secondary_mass_sorted2,event_multiplicities,event_index_sorted2])


numpy.save(basePath + '/' + new_directory + '/' + 'triple_data_upto_redshift_%d_time_delay_%d_Myr_triple_%d_Myr.npy'%(desired_redshift,fixed_delay_Myr,triple_induced_delay), {
    'triple_id1s': triple_id1s,
    'triple_id2s': triple_id2s,
    'triple_id3s': triple_id3s,
    'triple_mass1s': triple_mass1s,
    'triple_mass2s': triple_mass2s,
    'triple_mass3s': triple_mass3s,
    'Ms': Ms,
    'triple_formation_time': triple_formation_time,
    'triple_event_sorted_index': triple_event_sorted_index
})

