import sys
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')

#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
#sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')

#%pylab inline
import h5py
import numpy
import illustris_python as il
import os
from kdcount import correlate
import scipy
import matplotlib as mpl

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
            return box_size
        
        
        
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
        
def make_cuts(quantities,cut): #selects an array of quantities (argument 1) and makes cuts (argument 2)
    cutted_quantities=[quantity[cut] for quantity in quantities]
    return cutted_quantities

def get_group_property(output_path,group_property,desired_redshift,list_all=True):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift,list_all=False)
    halos = il.groupcat.loadHalos(output_path,output_snapshot,fields=None)
    if (list_all):
        print('Below are the list of properties')
        print(halos.keys())
    return halos.get(group_property),output_redshift

def get_subhalo_property(output_path,subhalo_property,desired_redshift,list_all=True):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift,list_all=False)      
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
    box_size_Mpc=box_size/1000.
    #print(HM)
    
    def extract(HM_min,HM_max):
        mask=(HM>HM_min)&(HM<HM_max)
        return (HM_min+HM_max)/2,len(HM[mask])

    HM_bin=numpy.logspace(log_HM_min,log_HM_max,Nbins,endpoint=True)
    out=[extract(HM_bin[i],HM_bin[i+1]) for i in range(0,len(HM_bin)-1)]
    #return out
    centers=numpy.array(list(zip(*out))[0])
    counts=numpy.array(list(zip(*out))[1])
    dM=centers*numpy.diff(numpy.log(centers))[0]
    HMF=counts/dM/box_size_Mpc**3
    dHMF=numpy.sqrt(counts)/dM/box_size_Mpc**3
    return centers,HMF,dHMF


def BH_mass_function_AGN_fraction(bhmass,bolometric_luminosity,box_size,Nbins,log_HM_min,log_HM_max,log_lbol_cut):
    box_size_Mpc=box_size/1000.
    #print(HM)
    HM=bhmass
    def extract(HM_min,HM_max):
        mask=(HM>HM_min)&(HM<HM_max)
        mask2=mask&(bolometric_luminosity>10**log_lbol_cut)
        return (HM_min+HM_max)/2,len(HM[mask]),float(len(HM[mask2]))/(len(HM[mask])+0.00000001),numpy.sqrt(float(len(HM[mask2])))/(len(HM[mask])+0.00000001)

    HM_bin=numpy.logspace(log_HM_min,log_HM_max,Nbins,endpoint=True)
    out=[extract(HM_bin[i],HM_bin[i+1]) for i in range(0,len(HM_bin)-1)]
    #return out
    centers=numpy.array(list(zip(*out))[0])
    counts=numpy.array(list(zip(*out))[1])
    AGN_fraction=numpy.array(list(zip(*out))[2])
    d_AGN_fraction=numpy.array(list(zip(*out))[3])
    dM=centers*numpy.diff(numpy.log(centers))[0]
    HMF=counts/dM/box_size_Mpc**3
    dHMF=numpy.sqrt(counts)/dM/box_size_Mpc**3
    return centers,HMF,dHMF,AGN_fraction,d_AGN_fraction




def get_mass_function(category,object_type,desired_redshift,output_path,Nbins,log_mass_min,log_mass_max,list_all=True,dynamical_bh_mass=True):
    box_size=get_box_size(output_path)
      #print (box_)
    if (object_type=='group'):
        if (category=='total'):
            mass,output_redshift=get_group_property(output_path,'GroupMass',desired_redshift,list_all=list_all)
            centers,HMF,dHMF=mass_function(mass*1e10,box_size,Nbins,log_mass_min,log_mass_max)
            return centers,HMF,dHMF,output_redshift
        
    if (object_type=='subhalo'):
        if (category=='total'):
            mass,output_redshift=get_subhalo_property(output_path,'SubhaloMass',desired_redshift,list_all=list_all)
            centers,HMF,dHMF=mass_function(mass*1e10,box_size,Nbins,log_mass_min,log_mass_max)
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
                
            if ((p_type==5)&(dynamical_bh_mass==True)):
                mass,output_redshift=get_subhalo_property(output_path,'SubhaloBHMass',desired_redshift,list_all=list_all) 
            else:
                masstype,output_redshift=get_subhalo_property(output_path,'SubhaloMassType',desired_redshift,list_all=list_all)           
                mass=masstype[:,p_type]
            centers,HMF,dHMF=mass_function(mass*1e10,box_size,Nbins,log_mass_min,log_mass_max)
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
    p=4./3*scipy.pi*(rmax**3-rmin**3)/(BOXSIZE/1e3)**3
    return p 

    
def correlate_info(data, NBINS, RMIN, RMAX, BOXSIZE, WRAP):
    if data is not None:
        if RMAX is None:
            RMAX = BOXSIZE
        
        if WRAP:
            wrap_length = BOXSIZE/1e3
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


def get_particle_property_within_groups(output_path,particle_property,p_type,desired_redshift,subhalo_index,group_type='groups',list_all=True):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift,list_all)
    if (list_all):
        print('Below are the list of properties for ptype ',p_type)
        print(il.snapshot.loadSubset(output_path,output_snapshot,p_type).keys())
    
    requested_property=il.snapshot.loadSubset(output_path,output_snapshot,p_type)[particle_property]

    if (group_type=='groups'):              
        group_lengths,output_redshift=(get_group_property(output_path,'GroupLenType', desired_redshift,list_all=False))
        group_lengths=group_lengths[:,p_type] 
        group_offsets=numpy.array([sum(group_lengths[0:i]) for i in range(0,len(group_lengths))])
    elif (group_type=='subhalo'):              
        group_lengths,output_redshift=(get_subhalo_property(output_path,'SubhaloLenType', desired_redshift,list_all=False))
        group_lengths=group_lengths[:,p_type] 
        group_offsets=numpy.array([sum(group_lengths[0:i]) for i in range(0,len(group_lengths))])
    else:
        print("Error:Unidentified group type")
        return
    return requested_property[group_offsets[subhalo_index]:group_offsets[subhalo_index]+group_lengths[subhalo_index]],output_redshift

def get_particle_property_within_groups(output_path,particle_property,p_type,desired_redshift,subhalo_index,group_type='groups',list_all=True):

    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift,list_all=False)
    requested_property=il.snapshot.loadSubset(output_path,output_snapshot,p_type)[particle_property]
    #requested_property_parent_group=il.snapshot.loadSubset(output_path,output_snapshot,p_type)[particle_property]
    

    if (group_type=='groups'):              
        group_lengths,output_redshift=(get_group_property(output_path,'GroupLenType', desired_redshift,list_all=False))
        group_lengths=group_lengths[:,p_type] 
        group_offsets=numpy.array([sum(group_lengths[0:i]) for i in range(0,len(group_lengths))])                   
        group_particles=requested_property[group_offsets[subhalo_index]:group_offsets[subhalo_index]+group_lengths[subhalo_index]]
        return group_particles,output_redshift
    elif (group_type=='subhalo'):              

        group_lengths,output_redshift=(get_group_property(output_path,'GroupLenType', desired_redshift))
        group_lengths=group_lengths[:,p_type] 
        group_offsets=numpy.array([sum(group_lengths[0:i]) for i in range(0,len(group_lengths))])
        
        subhalo_lengths,output_redshift=(get_subhalo_property(output_path,'SubhaloLenType', desired_redshift))
        subhalo_lengths=subhalo_lengths[:,p_type] 
        subhalo_indices=numpy.arange(0,len(subhalo_lengths))
        
        
        subhalo_group_number,output_redshift=(get_subhalo_property(output_path,'SubhaloGrNr', desired_redshift,list_all=False));
        desired_group_number=subhalo_group_number[subhalo_index]  
        subhalo_lengths=subhalo_lengths[subhalo_group_number==desired_group_number]
        subhalo_offsets=numpy.array([sum(subhalo_lengths[0:i]) for i in range(0,len(subhalo_lengths))])
        
        mask=subhalo_group_number==desired_group_number
        #print(len(mask)),mask
        subhalo_indices=subhalo_indices[mask]  
        subhalo_final_indices=numpy.arange(0,len(subhalo_indices))
        group_particles=requested_property[group_offsets[desired_group_number]:group_offsets[desired_group_number]+group_lengths[desired_group_number]]   

        #subhalo_indices=subhalo_indices[subhalo_group_number==desired_group_number]
        final_index=(subhalo_final_indices[subhalo_indices==subhalo_index])[0]
        
        subhalo_particles=group_particles[subhalo_offsets[final_index]:subhalo_offsets[final_index]+subhalo_lengths[final_index]]
      
        return subhalo_particles,group_particles,output_redshift     
    else:
        print("Error:Unidentified group type")
        
        
def make_image(Coordinates,Coordinates_for_COM,plane,obj,boxsize,NBINS,colormap='Blues_r',opacity=1):
    x_pos=Coordinates[:,0]
    y_pos=Coordinates[:,1]
    z_pos=Coordinates[:,2]

    x_pos_COM=Coordinates_for_COM[:,0]
    y_pos_COM=Coordinates_for_COM[:,1]
    z_pos_COM=Coordinates_for_COM[:,2]
    
    #if (CENTER_OF_MASS):   
    COM_x=numpy.median(x_pos_COM)
    COM_y=numpy.median(y_pos_COM)
    COM_z=numpy.median(z_pos_COM)

    def min_dis(median_position, position,box_size):
        pos_1=position-median_position
        pos_2=position-median_position+boxsize
        pos_3=position-median_position-boxsize

        new_position_options=numpy.array([pos_1,pos_2,pos_3])
        get_minimum_distance=numpy.argmin(numpy.abs(new_position_options))
        #print(new_position_options)

        #print(get_minimum_distance)
        return new_position_options[get_minimum_distance]

    vectorized_min_dis = numpy.vectorize(min_dis)
    x_pos_wrapped=vectorized_min_dis(COM_x,x_pos,boxsize)
    y_pos_wrapped=vectorized_min_dis(COM_y,y_pos,boxsize)
    z_pos_wrapped=vectorized_min_dis(COM_z,z_pos,boxsize)

    plane='xz'


    if (plane=='xy'):
        obj.hist2d(x_pos_wrapped,y_pos_wrapped, bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap=colormap,alpha=opacity);
    if (plane=='yz'):
        obj.hist2d(y_pos_wrapped,z_pos_wrapped, bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap=colormap,alpha=opacity);
    if (plane=='xz'):
        obj.hist2d(x_pos_wrapped,z_pos_wrapped, bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap=colormap,alpha=opacity);

        
def get_merger_events(output_path):

    output_file_names=os.listdir(output_path+'blackhole_mergers/')
    snapshot_space=[]
    redshift_space=[]

    file_id_complete=numpy.array([],dtype=int)
    scale_fac_complete=numpy.array([])

    BH_id1_complete=numpy.array([],dtype=int)
    BH_mass1_complete=numpy.array([])
    BH_id2_complete=numpy.array([],dtype=int)
    BH_mass2_complete=numpy.array([])

    N_empty=0

    for name in output_file_names[:]:
        data=numpy.loadtxt(output_path+'blackhole_mergers/'+name)
        #data=numpy.transpose(data)
        try:
            file_id=data[:,0].astype(int)
            scale_fac=data[:,1]

            BH_id1=data[:,2].astype(int)
            BH_mass1=data[:,3]
            BH_id2=data[:,4].astype(int)
            BH_mass2=data[:,5]

            file_id_complete=numpy.append(file_id_complete,file_id)
            scale_fac_complete=numpy.append(scale_fac_complete,scale_fac)
            BH_id1_complete=numpy.append(BH_id1_complete,BH_id1)
            BH_mass1_complete=numpy.append(BH_mass1_complete,BH_mass1)    
            BH_id2_complete=numpy.append(BH_id2_complete,BH_id2)
            BH_mass2_complete=numpy.append(BH_mass2_complete,BH_mass2) 
        except IndexError:
            N_empty+=1
            aaa=1
    mass_tuple=list(zip(BH_mass1_complete,BH_mass2_complete))
    id_tuple=list(zip(BH_id1_complete,BH_id2_complete))

    #primary_mass=numpy.array([numpy.amax([dat[0],dat[1]]) for dat in mass_tuple])
    #secondary_mass=numpy.array([numpy.amin([dat[0],dat[1]]) for dat in mass_tuple])

    primary_index=numpy.array([numpy.argmax([dat[0],dat[1]]) for dat in mass_tuple])
    secondary_index=numpy.array([numpy.argmin([dat[0],dat[1]]) for dat in mass_tuple])

    primary_mass=numpy.array([mass_t[index] for (mass_t,index) in list(zip(mass_tuple,primary_index))])
    secondary_mass=numpy.array([mass_t[index] for (mass_t,index) in list(zip(mass_tuple,secondary_index))])

    primary_id=numpy.array([id_t[index] for (id_t,index) in list(zip(id_tuple,primary_index))])
    secondary_id=numpy.array([id_t[index] for (id_t,index) in list(zip(id_tuple,secondary_index))])
    return scale_fac_complete,primary_mass,secondary_mass,primary_id,secondary_id,file_id_complete,N_empty

def get_blackhole_history_high_res(output_path,desired_id):
    def parse_id_col(BH_ids_as_string):
        return numpy.int(BH_ids_as_string[3:])
    vec_parse_id_col=numpy.vectorize(parse_id_col)
    
    output_file_names=os.listdir(output_path+'blackhole_details/')

    BH_ids_for_id=numpy.array([],dtype=int)
    scale_factors_for_id=numpy.array([])
    BH_masses_for_id=numpy.array([])
    BH_mdots_for_id=numpy.array([])
    rhos_for_id=numpy.array([])
    sound_speeds_for_id=numpy.array([])
    
    
    merging_time,primary_mass,secondary_mass,primary_id,secondary_id,file_id_complete=get_merger_events(output_path)
    extract_events_as_secondary_BH=secondary_id==desired_id
    extract_events_as_primary_BH=primary_id==desired_id
    merging_partner_ids=numpy.append(primary_id[extract_events_as_secondary_BH],secondary_id[extract_events_as_primary_BH])    
    merging_times=numpy.append(merging_time[extract_events_as_secondary_BH],merging_time[extract_events_as_primary_BH])    

    total_desired_ids=numpy.append(numpy.array([desired_id]),merging_partner_ids)
    
    for output_file_name in output_file_names[:]:
        if ('blackhole_details' in output_file_name):
            full_data=numpy.loadtxt(output_path+'blackhole_details/'+output_file_name,dtype='str')
            BH_ids=vec_parse_id_col(full_data[:,0])
            scale_factors=(full_data[:,1]).astype('float')
            BH_masses=(full_data[:,2]).astype('float')
            BH_mdots=(full_data[:,3]).astype('float')
            rhos=(full_data[:,3]).astype('float')
            sound_speeds=(full_data[:,3]).astype('float')

            final_extract=numpy.array([False]*len(sound_speeds))
            for d_id in total_desired_ids:
                extract_id=(d_id==BH_ids)
                final_extract=final_extract+extract_id
            #print(len(BH_ids),len(BH_ids[extract_id]))
            BH_ids_for_id=numpy.append(BH_ids_for_id,BH_ids[final_extract])
            scale_factors_for_id=numpy.append(scale_factors_for_id,scale_factors[final_extract])
            BH_masses_for_id=numpy.append(BH_masses_for_id,BH_masses[final_extract])
            BH_mdots_for_id=numpy.append(BH_mdots_for_id,BH_mdots[final_extract])
            rhos_for_id=numpy.append(rhos_for_id,rhos[final_extract])
            sound_speeds_for_id=numpy.append(sound_speeds_for_id,sound_speeds[final_extract])
            
    return BH_ids_for_id,scale_factors_for_id,BH_masses_for_id,BH_mdots_for_id,rhos_for_id,sound_speeds_for_id,merging_times


def get_blackhole_history_high_res_all_progenitors(output_path,desired_id):
    def parse_id_col(BH_ids_as_string):
        return numpy.int(BH_ids_as_string[3:])
    vec_parse_id_col=numpy.vectorize(parse_id_col)
    
    output_file_names=os.listdir(output_path+'blackhole_details/')

    BH_ids_for_id=numpy.array([],dtype=int)
    scale_factors_for_id=numpy.array([])
    BH_masses_for_id=numpy.array([])
    BH_mdots_for_id=numpy.array([])
    rhos_for_id=numpy.array([])
    sound_speeds_for_id=numpy.array([])
    
  
    total_desired_ids,merging_times=get_progenitors_and_descendants(output_path,desired_id)
    ii=0
    for output_file_name in output_file_names[:]:
        if ('blackhole_details' in output_file_name):
            print(ii)
            ii+=1
            full_data=numpy.loadtxt(output_path+'blackhole_details/'+output_file_name,dtype='str')
            BH_ids=vec_parse_id_col(full_data[:,0])
            scale_factors=(full_data[:,1]).astype('float')
            BH_masses=(full_data[:,2]).astype('float')
            BH_mdots=(full_data[:,3]).astype('float')
            rhos=(full_data[:,3]).astype('float')
            sound_speeds=(full_data[:,3]).astype('float')

            final_extract=numpy.array([False]*len(sound_speeds))
            for d_id in total_desired_ids:
                extract_id=(d_id==BH_ids)
                final_extract=final_extract+extract_id
            #print(len(BH_ids),len(BH_ids[extract_id]))
            BH_ids_for_id=numpy.append(BH_ids_for_id,BH_ids[final_extract])
            scale_factors_for_id=numpy.append(scale_factors_for_id,scale_factors[final_extract])
            BH_masses_for_id=numpy.append(BH_masses_for_id,BH_masses[final_extract])
            BH_mdots_for_id=numpy.append(BH_mdots_for_id,BH_mdots[final_extract])
            rhos_for_id=numpy.append(rhos_for_id,rhos[final_extract])
            sound_speeds_for_id=numpy.append(sound_speeds_for_id,sound_speeds[final_extract])
            
    return BH_ids_for_id,scale_factors_for_id,BH_masses_for_id,BH_mdots_for_id,rhos_for_id,sound_speeds_for_id,merging_times

        
        
def get_progenitors_and_descendants(output_path,desired_id,MAX_ITERATION=100):

    BH_ids_for_id=numpy.array([],dtype=int)
    scale_factors_for_id=numpy.array([])
    BH_masses_for_id=numpy.array([])
    BH_mdots_for_id=numpy.array([])
    rhos_for_id=numpy.array([])
    sound_speeds_for_id=numpy.array([])
    merging_time,primary_mass,secondary_mass,primary_id,secondary_id,file_id_complete,N_empty=get_merger_events(output_path)

    progenitor_ids=numpy.array([desired_id],dtype=int)
    final_merging_times=numpy.array([])
    progenitor_ids_before_update=numpy.array([],dtype=int)
    i=0
    while (len(progenitor_ids_before_update)<len(progenitor_ids)):
            
        progenitor_ids_before_update=progenitor_ids+False
        
        extract_events_as_primary_BH=numpy.array([False]*len(secondary_id))
        extract_events_as_secondary_BH=numpy.array([False]*len(secondary_id))
        for ids in progenitor_ids:
            extract_events_as_secondary_BH=extract_events_as_secondary_BH+(secondary_id==ids)
            extract_events_as_primary_BH=extract_events_as_primary_BH+(primary_id==ids)
        merging_partner_ids=numpy.append(primary_id[extract_events_as_secondary_BH],secondary_id[extract_events_as_primary_BH])    
        merge_times=numpy.append(merging_time[extract_events_as_secondary_BH],merging_time[extract_events_as_primary_BH])    

        progenitor_ids=numpy.append(progenitor_ids,merging_partner_ids)
        final_merging_times=numpy.append(final_merging_times,merge_times)
         
        progenitor_ids_before_distinct=progenitor_ids+0
        progenitor_ids=numpy.unique(progenitor_ids)
        
        #temp_merging_times=[]
        #for prog_id in progenitor_ids:
        #    extract_id=progenitor_ids_before_distinct==prog_id
        #    temp=final_merging_times[extract_id]
        #    if (len(temp)==0):
        #        print("Warning! Merging time missing")
        #    else:
        #        temp_merging_times.append(temp[0])
        final_merging_times=numpy.unique(final_merging_times)
        i+=1
        if (i>MAX_ITERATION):
            print("Maximum number of iterations reached! Aborting")
            break
    return progenitor_ids,final_merging_times


def generate_group_ids(output_path,desired_redshift,p_type,save_output_path='./',group_type='groups',create=False):    
    global complete_particle_ids
    particle_property='ParticleIDs'
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)
  
    if not os.path.exists(save_output_path):
        print("Making directory for storing groupids")
        os.makedirs(save_output_path)
    
    if ((os.path.exists(save_output_path+'group_ids_%d.npy'%output_snapshot))&(create==False)):
        print("File exists!! group ids exist already")
        group_ids=numpy.load(save_output_path+'group_ids_%d.npy'%output_snapshot)
        return group_ids
        
    else:
        print("Constructing_group_ids")
        
    
    
    complete_particle_ids,output_redshift=get_particle_property(output_path, particle_property, p_type, desired_redshift)
    #print(len(complete_particle_ids))
    group_ids=numpy.array([-1]*len(complete_particle_ids))
    def find_index(id_to_be_searched):
        #print(numpy.where(complete_particle_ids==id_to_be_searched))
        
        return (numpy.where(complete_particle_ids==id_to_be_searched))[0]
    vec_find_index=numpy.vectorize(find_index)
    
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift,list_all=False)
    requested_property=il.snapshot.loadSubset(output_path,output_snapshot,p_type)[particle_property]

    print('Reading group lengths')
    if (group_type=='groups'):              
        group_lengths,output_redshift=(get_group_property(output_path,'GroupLenType', desired_redshift,list_all=False))
        group_lengths=group_lengths[:,p_type] 
        #group_offsets=numpy.array([sum(group_lengths[0:i]) for i in range(0,len(group_lengths))])   
        print('Constructing offsets')
        group_offsets=numpy.append(0,(numpy.cumsum(group_lengths))[:-1])
       #return group_particles,output_redshift
    elif (group_type=='subhalo'):     
        
        print("Warning!! This is going to be very slow! Still under development")

        group_lengths,output_redshift=(get_group_property(output_path,'GroupLenType', desired_redshift,list_all=False))
        group_lengths=group_lengths[:,p_type] 
        group_offsets=numpy.array([sum(group_lengths[0:i]) for i in range(0,len(group_lengths))])
        
        subhalo_lengths,output_redshift=(get_subhalo_property(output_path,'SubhaloLenType', desired_redshift,list_all=False))
        subhalo_lengths=subhalo_lengths[:,p_type] 
        subhalo_indices=numpy.arange(0,len(subhalo_lengths))
        
        
        subhalo_group_number,output_redshift=(get_subhalo_property(output_path,'SubhaloGrNr', desired_redshift,list_all=False));

 
        #subhalo_indices=subhalo_indices[subhalo_group_number==desired_group_number]
        #final_index=(subhalo_final_indices[subhalo_indices==subhalo_index])[0]
        
        #group_particles=group_particles[subhalo_offsets[final_index]:subhalo_offsets[final_index]+subhalo_lengths[final_index]]
      
        #return subhalo_particles,group_particles,output_redshift     



    #MAX_ITERATIONS=1000000
    for subhalo_index in range(0,len(group_offsets)): 
        if(subhalo_index%10==0):
            aaa=1
            #print(subhalo_index)
            #print(len(group_ids[group_ids==-1]))
            #print("-------------")
        if (group_type=='groups'):
            #print(len(group_offsets),subhalo_index)
            complete_particle_ids_group_wise=requested_property[group_offsets[subhalo_index]:group_offsets[subhalo_index]+group_lengths[subhalo_index]]
        if (group_type=='subhalo'):     
            desired_group_number=subhalo_group_number[subhalo_index]  
            subhalo_lengths_for_group=subhalo_lengths[subhalo_group_number==desired_group_number]
            subhalo_offsets_for_group=numpy.array([sum(subhalo_lengths[0:i]) for i in range(0,len(subhalo_lengths))])

            mask=subhalo_group_number==desired_group_number
            #print(len(mask)),mask
            subhalo_indices_for_group=subhalo_indices[mask]  
            subhalo_final_indices=numpy.arange(0,len(subhalo_indices_for_group))
            group_particles=requested_property[group_offsets[desired_group_number]:group_offsets[desired_group_number]+group_lengths[desired_group_number]]   
            #print(subhalo_indices_for_group,subhalo_index)
            final_index=(subhalo_final_indices[subhalo_indices_for_group==subhalo_index])[0]
            complete_particle_ids_group_wise=group_particles[subhalo_offsets_for_group[final_index]:subhalo_offsets_for_group[final_index]+subhalo_lengths_for_group[final_index]]

        if (len(complete_particle_ids_group_wise)==0):
            continue
        #print(complete_particle_ids_group_wise)
        indices_to_be_assigned=vec_find_index(complete_particle_ids_group_wise)
        
        group_ids[indices_to_be_assigned]=subhalo_index
        print(len(group_ids[group_ids==-1]))
        if(len(group_ids[group_ids==-1])==0):
            break
    
    #if (subhalo_index==MAX_ITERATIONS-1):
    #    print("Warning! Maximum number of iterations reached!")
    print(output_redshift,output_snapshot)
    numpy.save(save_output_path+'group_ids_%d.npy'%output_snapshot,group_ids)
    return group_ids
            


def generate_subhalo_ids(output_path,desired_redshift,p_type,save_output_path='./',create=False):
    output_redshift,output_snapshot=desired_redshift_to_output_redshift(output_path,desired_redshift)   
    if ((os.path.exists(save_output_path+'group_ids_%d.npy'%output_snapshot))&(create==False)):
        print("File exists!! subhalo ids exist already")
        subhalo_ids=numpy.load(save_output_path+'subhalo_ids_%d.npy'%output_snapshot)
        return subhalo_ids
    
    
    group_ids=generate_group_ids(output_path,desired_redshift,p_type,save_output_path=save_output_path)
    BH_Pos,output_redshift=get_particle_property(output_path,'Coordinates',p_type,desired_redshift,list_all=False)
    ParticleIDs,output_redshift=get_particle_property(output_path,'ParticleIDs',p_type,desired_redshift,list_all=False)
    boxsize=get_box_size(output_path)

    SubhaloGrNr,output_redshift=get_subhalo_property(output_path,'SubhaloGrNr',desired_redshift,list_all=False)
    SubhaloPos,output_redshift=get_subhalo_property(output_path,'SubhaloPos',desired_redshift,list_all=False)
    SubhaloMass,output_redshift=get_subhalo_property(output_path,'SubhaloMass',desired_redshift,list_all=False)
    SubhaloMass=SubhaloMass*1e10
    Subhalo_Indices=numpy.arange(len(SubhaloGrNr))
    mask=SubhaloMass==SubhaloMass
    #mask=SubhaloMass>=1e11
    Subhalo_Indices_cut=Subhalo_Indices[mask]
    SubhaloPos_cut=SubhaloPos[mask]
    SubhaloGrNr_cut=SubhaloGrNr[mask]

    def min_dis(median_position, position,box_size):
            pos_1=position-median_position
            pos_2=position-median_position+boxsize
            pos_3=position-median_position-boxsize

            new_position_options=numpy.array([pos_1,pos_2,pos_3])
            get_minimum_distance=numpy.argmin(numpy.abs(new_position_options))
            #print(new_position_options)

            #print(get_minimum_distance)
            return new_position_options[get_minimum_distance]


    def get_subhalo_id(blackhole_info):
        #print(blackhole_info)
        blackhole_group_id=blackhole_info[0]
        blackhole_position=blackhole_info[1]
        extract_ids_within_the_parent_FOF=blackhole_group_id==SubhaloGrNr_cut
        Subhalo_Indices_current_BH=Subhalo_Indices_cut[extract_ids_within_the_parent_FOF]
        #print(Subhalo_Indices_current_BH)
        SubhaloPos_current_BH=SubhaloPos_cut[extract_ids_within_the_parent_FOF]
        x_pos_Subhalo=SubhaloPos_current_BH[:,0]
        y_pos_Subhalo=SubhaloPos_current_BH[:,1]
        z_pos_Subhalo=SubhaloPos_current_BH[:,2]

        vectorized_min_dis = numpy.vectorize(min_dis)
        try:
            x_dis=vectorized_min_dis(blackhole_position[0],x_pos_Subhalo,boxsize)
            y_dis=vectorized_min_dis(blackhole_position[1],y_pos_Subhalo,boxsize)
            z_dis=vectorized_min_dis(blackhole_position[2],z_pos_Subhalo,boxsize)

            distance_sq=x_dis**2+y_dis**2+z_dis**2
            min_distance_sq=numpy.amin(distance_sq)
            subhalo_id=(Subhalo_Indices_current_BH[distance_sq==min_distance_sq])[0]
            return subhalo_id,min_distance_sq
        except ValueError:
            return -1,-1


    blackhole_info_space=list(zip(group_ids,BH_Pos))
    data=[get_subhalo_id(blackhole_info) for blackhole_info in blackhole_info_space]
    subhalo_ids=(numpy.array(data))[:,0]
    distance_from_subhalo_center=(numpy.array(data))[:,1]
    numpy.save(save_output_path+'subhalo_ids_%d.npy'%output_snapshot,subhalo_ids)
    numpy.save(save_output_path+'distance_from_subhalo_center_%d.npy'%output_snapshot,distance_from_subhalo_center)    
    return subhalo_ids,distance_from_subhalo_center
        
def mean_plot(x,y,xscl,yscl,nbins):
    #nbins = 5
    if(yscl==True):
        y=log10(y)
    if(xscl==True):
        x=log10(x)
    n, _ = numpy.histogram(x, bins=nbins)
    sy, _ = numpy.histogram(x, bins=nbins, weights=y)
    sy2, _ = numpy.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    std = numpy.sqrt(sy2/n - mean*mean)
    #std=1/numpy.sqrt(n)
    #plt.plot(x, y, 'bo')
    #plt.errorbar((_[1:] + _[:-1])/2, mean,std, color='blue', label = 'z = 8')
    #mean= savitzky_golay(mean, 11, 3)
    #print (_[1:] + _[:-1])/2
    #print mean
    mask=(std/mean)*100<100.
    
    x=((_[1:] + _[:-1])/2)[mask]
    y=mean[mask]
    yul=y+std[mask]
    yll=y-std[mask]
    return x,y,yul,yll#,plt.errorbar(x,y,color=colour,linewidth=3,label=labl),plt.fill_between(x,yll, yul,color=colour,alpha=0.5)
 

def median_plot(x,y,xscl,yscl,nbins):
    #nbins = 5
    if(yscl==True):
        y=log10(y)
    if(xscl==True):
        x=log10(x)
    x_space=numpy.linspace(numpy.amin(x),numpy.amax(x),nbins)
    y_space_med=[numpy.median(y[(x>x_space[i])&(x<x_space[i+1])]) for i in range(0,len(x_space)-1)]
    x_space_med=[(x_space[i]+x_space[i+1])/2 for i in range(0,len(x_space)-1)]
    return x_space_med,y_space_med
    
    
    
    
def luminosity_function(HM,box_size,log_HM_min,log_HM_max,Nbins):
        def extract(HM_min,HM_max):
            mask=(HM>HM_min)&(HM<HM_max)
            #print len(HM[mask])
            return (HM_min+HM_max)/2,len(HM[mask])
        
        HM_bin=numpy.logspace(log_HM_min,log_HM_max,Nbins,endpoint=True)
        out=[extract(HM_bin[i],HM_bin[i+1]) for i in range(0,len(HM_bin)-1)]
        centers=numpy.array(list(zip(*out))[0])
        counts=numpy.array(list(zip(*out))[1])
        #print counts
        dlogM=numpy.diff(numpy.log10(centers))[0]
        HMF=counts/dlogM/box_size**3
        dHMF=numpy.sqrt(counts)/dlogM/box_size**3
        return centers,HMF,dHMF        
    

    
    
        
        
        
        
        
        
        
        
        
        
        







