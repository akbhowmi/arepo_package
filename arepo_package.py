import sys
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')

sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')

#%pylab inline
import h5py
import numpy
import illustris_python as il
import os
from kdcount import correlate
import scipy

def get_snapshot_redshift_correspondence(output_path):
    output_file_names=os.listdir(output_path)
    snapshot_space=[]
    redshift_space=[]
    for name in output_file_names:
        if ('groups' in name):
            snapshot_number=int(name[7:])
            snapshot_space.append(snapshot_number)
    snapshot_space=numpy.sort(numpy.array(snapshot_space))
    for snapshot_number in snapshot_space:
            header=il.groupcat.loadHeader(output_path,snapshot_number)
            redshift=header.get('Redshift')   
            redshift_space.append(redshift)   
    return numpy.array(snapshot_space),numpy.array(redshift_space)


def get_box_size(output_path):
    output_file_names=os.listdir(output_path)
    snapshot_space=[]
    redshift_space=[]
    for name in output_file_names:
        if ('groups' in name):
            snapshot_number=int(name[7:])
            header=il.groupcat.loadHeader(output_path,snapshot_number)
            box_size=header.get('BoxSize')   
            return box_size/1000.
        
        
        
def get_cosmology(output_path):
    output_file_names=os.listdir(output_path)
    snapshot_space=[]
    redshift_space=[]
    for name in output_file_names:
        if ('groups' in name):
            snapshot_number=int(name[7:])
            header=il.groupcat.loadHeader(output_path,snapshot_number)
            om0=header.get('Omega0')
            oml=header.get('OmegaLambda')
            h=header.get('HubbleParam')            
            return om0,oml,h
        
        
def load_snapshot_header(output_path,desired_redshift):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)
    with h5py.File(il.snapshot.snapPath(output_path, output_snapshot)) as f:
        header = dict(f['Header'].attrs.items())
        print(header)
        return header

def desired_redshift_to_output_redshift(output_path,desired_redshift,list_all=True): 
    snapshot_space,redshift_space=get_snapshot_redshift_correspondence(output_path)
    redshift_difference=numpy.abs(redshift_space-desired_redshift)
    min_redshift_difference=numpy.amin(redshift_difference)
    output_snapshot=snapshot_space[redshift_difference==min_redshift_difference][0]
    output_redshift=redshift_space[redshift_difference==min_redshift_difference][0]
    if (list_all):
        print("Desired redshift: ",desired_redshift)
        print("Output redshift: ",output_redshift)
        print("Output snapshot: ",output_snapshot)            
    return output_redshift,output_snapshot      
        
def get_group_property(output_path,group_property,desired_redshift,list_all=True):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)
    halos = il.groupcat.loadHalos(output_path,output_snapshot,fields=None)
    if (list_all):
        print('Below are the list of properties')
        print(halos.keys())
    return halos.get(group_property),output_redshift

def get_subhalo_property(output_path,subhalo_property,desired_redshift,list_all=True):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)      
    subhalos = il.groupcat.loadSubhalos(output_path,output_snapshot,fields=None)
    if (list_all):
        print('Below are the list of properties')
        print(subhalos.keys())
    return subhalos.get(subhalo_property),output_redshift

def get_particle_property(output_path,particle_property,p_type,desired_redshift,list_all=True):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift,list_all)
    if (list_all):
        print('Below are the list of properties for ptype ',p_type)
        print(il.snapshot.loadSubset(output_path,output_snapshot,p_type).keys())
    return il.snapshot.loadSubset(output_path,output_snapshot,p_type)[particle_property],output_redshift




def mass_function(HM,box_size,Nbins,log_HM_min,log_HM_max):
    #print(HM)
    
    def extract(HM_min,HM_max):
        mask=(HM>HM_min/1e10)&(HM<HM_max/1e10)
        return (HM_min+HM_max)/2,len(HM[mask])

    HM_bin=numpy.logspace(log_HM_min,log_HM_max,Nbins,endpoint=True)
    out=[extract(HM_bin[i],HM_bin[i+1]) for i in range(0,len(HM_bin)-1)]
    #return out
    centers=numpy.array(list(zip(*out))[0])
    counts=numpy.array(list(zip(*out))[1])
    dM=centers*numpy.diff(numpy.log(centers))[0]
    HMF=counts/dM/box_size**3
    dHMF=numpy.sqrt(counts)/dM/box_size**3
    return centers,HMF,dHMF


def get_mass_function(category,object_type,desired_redshift,output_path,Nbins,log_mass_min,log_mass_max,list_all=True):
    box_size=get_box_size(output_path)
      #print (box_)
    if (object_type=='group'):
        if (category=='total'):
            mass,output_redshift=get_group_property(output_path,'GroupMass',desired_redshift,list_all=list_all)
            centers,HMF,dHMF=mass_function(mass,box_size,Nbins,log_mass_min,log_mass_max)
            return centers,HMF,dHMF,output_redshift
        
    if (object_type=='subhalo'):
        if (category=='total'):
            mass,output_redshift=get_subhalo_property(output_path,'SubhaloMass',desired_redshift,list_all=list_all)
            centers,HMF,dHMF=mass_function(mass,box_size,Nbins,log_mass_min,log_mass_max)
            return centers,HMF,dHMF,output_redshift
        else:
            if (category=='stellar'):
                p_type=4
            if (category=='bh'):
                p_type=5               
            if (category=='dark'):
                p_type=1          
            if (category=='gas'):
                p_type=0                
                
            mass,output_redshift=get_subhalo_property(output_path,'SubhaloMassType',desired_redshift,list_all=list_all)
            #print(mass)
            
            centers,HMF,dHMF=mass_function(mass[:,p_type],box_size,Nbins,log_mass_min,log_mass_max)
            return centers,HMF,dHMF,output_redshift        

def get_particle_history(z_latest,z_earliest,z_no_of_bins,p_type,p_id_to_be_tracked,desired_property,output_path):
    z_space=numpy.linspace(z_latest,z_earliest,z_no_of_bins)
    prperty_history=[]
    z_history=[]
    for z in z_space: 
        p_id_current_z,output_redshift=get_particle_property(output_path,'ParticleIDs',p_type,z,list_all=False)
        extract_target_id=p_id_to_be_tracked==p_id_current_z
        prperty_space,output_redshift=get_particle_property(output_path,desired_property,p_type,z,list_all=False)
        temp=prperty_space[extract_target_id]
        if (len(temp)==1):
            prperty_history.append(temp[0])
            z_history.append(output_redshift)
    return numpy.array(prperty_history),numpy.array(z_history)

def poiss(rmin,rmax,BOXSIZE):
    p=4./3*scipy.pi*(rmax**3-rmin**3)/BOXSIZE**3
    return p 

    
def correlate_info(data, NBINS, RMIN, RMAX, BOXSIZE, WRAP):
    if data is not None:
        if RMAX is None:
            RMAX = BOXSIZE
        
        if WRAP:
            wrap_length = BOXSIZE
        else:
            wrap_length = None
        
        dataset = correlate.points(data, boxsize = wrap_length)  
        
        binning = correlate.RBinning(numpy.logspace(numpy.log10(RMIN),numpy.log10(RMAX),NBINS+1))
        
#	RR=N**2*numpy.asarray([poiss(rbin[i],rbin[i+1]) for i in range(0,nbins)])
        DD = correlate.paircount(dataset, dataset, binning, np=16)
        DD = DD.sum1
        
#        print 'Done correlating'
        r = binning.centers
        rbin=binning.edges
        N=len(data)
        
        RR=(N**2-N)*numpy.asarray([poiss(rbin[i],rbin[i+1],BOXSIZE) for i in range(0,NBINS)])
    
    
        return r, DD,RR
    else:
        return None, None

def get_dark_matter_correlation_function(output_path,input_redshift,NBINS, RMIN, RMAX, WRAP,subsample_factor):
    BOXSIZE=get_box_size(output_path)
    positions,output_redshift=get_particle_property(output_path,'Coordinates',1,input_redshift)
    positions=positions/1e3
    r,DD,RR=correlate_info(positions[::subsample_factor],NBINS, RMIN, RMAX, BOXSIZE, WRAP)
    xi=DD/RR-1
    return r,DD,RR,xi,output_redshift





