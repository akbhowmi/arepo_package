import numpy
import sys
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import mdot_to_Lbol
import arepo_package
import scipy.interpolate
radiative_efficiency=0.1
total_conv=mdot_to_Lbol.get_conversion_factor_arepo(radiative_efficiency)
import h5py
import os
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

path_to_output='/ufrc/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_ZOOM_RUNS2/'
run='density_and_metallicity_based_criterion_zoom_levelmin7_levelmax9_haloindex4_redshift5.00_logbhseedmass5.90_NSC/AREPO' # name of the simulation runs
desired_redshift_of_selected_halo=1 #your desired redshift. Note that the code will choose the closest redshift available in the data 
basePath=path_to_output+run+'/output/'
DESIRED_NO_OF_COLUMNS=6

output_file_names=numpy.array(os.listdir(basePath+'blackhole_details/'))

mask=numpy.array(['blackhole_details' in name for name in output_file_names])
details_files=output_file_names[mask]
details_files=numpy.sort(details_files)

number_of_files_for_each_rank=int(len(details_files)/(size-1))

number_of_files_for_last_rank=len(details_files)%(size-1)

if (rank==size-1):
    details_files_per_rank=details_files[rank * number_of_files_for_each_rank:rank * number_of_files_for_each_rank + number_of_files_for_last_rank]
else:
    details_files_per_rank=details_files[rank * number_of_files_for_each_rank:(rank + 1) * number_of_files_for_each_rank]

print(rank,details_files_per_rank)



path_to_cleaned_up_details=basePath+'blackhole_details_cleaned/'
if not os.path.exists(path_to_cleaned_up_details):
    print("Making directory for storing groupids")
    os.makedirs(path_to_cleaned_up_details)

for current_file in details_files_per_rank:
    #current_file='blackhole_details_190.txt' 
    g=open(path_to_cleaned_up_details + current_file,'w')
    with open(basePath+'blackhole_details/' + current_file) as f:
        for line in f:
            splitted_line=numpy.array(line.split(' ')) 
            remove_empty=splitted_line != ''
            splitted_line=splitted_line[remove_empty]
            if(len(splitted_line) == DESIRED_NO_OF_COLUMNS):
                if((splitted_line[0])[0:2] == 'BH'):
                   if ('nan' not in line):
#                       if ('-nan' not in splitted_line):
                           g.write(line) 
    g.close()
