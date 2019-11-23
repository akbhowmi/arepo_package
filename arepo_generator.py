import sys
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import arepo_package
import scipy.interpolate
import h5py
import os
import numpy
#%pylab inline

path_to_generated_files='/ufrc/lblecha/aklantbhowmick/NEW_AREPO_RUNS/'
BOXSIZE_IN_MPC=25
BASEGRID_PARAMETER=7
N=2**BASEGRID_PARAMETER
TYPE='UNIFORM'
MUSIC_OUTPUT_FILENAME='IC.hdf5'


FOLDERNAME='L%dn%d_%s'%(BOXSIZE_IN_MPC,N,TYPE)
if not os.path.exists(path_to_generated_files+FOLDERNAME):
        print("Making new directory")
        os.makedirs(path_to_generated_files+FOLDERNAME)
        
if not os.path.exists(path_to_generated_files+FOLDERNAME+'/'+'MUSIC'):
        print("Making new directory")
        os.makedirs(path_to_generated_files+FOLDERNAME+'/'+'MUSIC')
        
if not os.path.exists(path_to_generated_files+FOLDERNAME+'/'+'arepo'):
        print("Making new directory")
        os.makedirs(path_to_generated_files+FOLDERNAME+'/'+'arepo')
        
variables_to_be_edited=numpy.array(['levelmin','Omega_m','Omega_L','Omega_b','H0','sigma_8','filename'])
parameter_values=numpy.array(['8','0.3','0.7','0.04','70','0.82',MUSIC_OUTPUT_FILENAME])
#parameter_values=parameter_values.astype(str)

if(len(variables_to_be_edited)!=len(parameter_values)):
    print("ERROR: must have equal number of parameters and their values")
    exit

prototypefilepath = '/ufrc/lblecha/aklantbhowmick/MUSIC/final.conf'
fp=open(prototypefilepath)
generated_music_config_filename = 'MUSIC_CONFIG.txt'
line = fp.readline()
music_config=open(path_to_generated_files+FOLDERNAME+'/'+'MUSIC/'+generated_music_config_filename,'w')
while line:
    #print(line)
    line = fp.readline()
    
    line_to_write=line
    
    if ('=' in line):
        line_splitted=line.split('=')
        parameter_name,parameter_value=line_splitted[0],line_splitted[1]
        #print(parameter_name)
        for name,value in list(zip(variables_to_be_edited,parameter_values)):
            #print(name, parameter_name)
            #print(name in parameter_name)
            
            if (parameter_name in name)|(name in parameter_name):
                #print("Match")
                line_to_write='%s\t\t=\t%s\n'%(name,value)
                
                continue
                #   
    
    #if ()
    
    music_config.write(line_to_write) 
    print(line_to_write)
   
    
    
    
music_config.close()
