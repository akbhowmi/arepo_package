import arepo_package

global redshift_cut,mass_cut,mass_ratio_cut,total_mass_cut,apply_host_cuts,HDF5,ZOOM




def get_merger_distributions(basePath,nbins,min_value,max_value,mass_cut,linear=0,HDF5=0,path_to_hdf5='./',upto_redshift=0):
    if(HDF5):
        fobj=h5py.File(path_to_hdf5)
        scale_fac_complete_sorted=fobj.get('time')[:]
        primary_mass_sorted=fobj.get('mass_out')[:]
        secondary_mass_sorted=fobj.get('mass_in')[:]
    else:
        scale_fac_complete_sorted,primary_mass_sorted,secondary_mass_sorted,primary_id_sorted,secondary_id_sorted,file_id_complete_sorted,N_empty=arepo_package.get_merger_events_from_snapshot(basePath,upto_redshift,HOSTS=0)
    mass_ratio=secondary_mass_sorted/primary_mass_sorted
    merger_redshifts=1./scale_fac_complete_sorted-1.
    mask_primary=primary_mass_sorted>=mass_cut/1e10
    mask_secondary=secondary_mass_sorted>=mass_cut/1e10
    mask_total_mass=(primary_mass_sorted+secondary_mass_sorted)>=total_mass_cut/1e10
    mask_redshift=merger_redshifts>=redshift_cut
    mask_ratio_cut=mass_ratio>=mass_ratio_cut
    mask=((mask_primary&mask_secondary)&(mask_redshift&mask_ratio_cut))&mask_total_mass
    merger_redshifts=merger_redshifts[mask]
    bins, dist, dist_err,norm,counts_sum=arepo_package.get_probability_density(merger_redshifts,nbins,min_value,max_value,linear=linear) 
    dist*=counts_sum
    dist_err*=counts_sum
    dt=numpy.array([T(zz-norm/2,zz+norm/2) for zz in bins])
    dz=numpy.diff(bins)[0]
    dz_dt=numpy.array([dz/T(zz-dz/2,zz+dz/2) for zz in bins])   
    return bins,dist,dist_err,norm,dz_dt


def get_merger_distributions_hosts(basePath,nbins,min_value,max_value,mass_cut,linear=0,HDF5=0,path_to_hdf5='./',upto_redshift=0):
    if(HDF5):
        fobj=h5py.File(path_to_hdf5)
        scale_fac_complete_sorted=fobj.get('time')[:]
        primary_mass_sorted=fobj.get('mass_out')[:]
        secondary_mass_sorted=fobj.get('mass_in')[:]
    else:
        scale_fac_complete,BH_id_primary,BH_mass_primary,hosthalomass_primary,hosthalostellarmass_primary,hosthalogasmass_primary,hosthalodmmass_primary,BH_id_secondary,BH_mass_secondary,hosthalomass_secondary,hosthalostellarmass_secondary,hosthalogasmass_secondary,hosthalodmmass_secondary,file_id_complete=arepo_package.get_merger_events_hosts(basePath,HDF5=1)        #print(primary_mass_sorted,secondary_mass_sorted)
    mass_ratio=BH_id_secondary/BH_mass_primary
    mask_ratio_cut=mass_ratio>mass_ratio_cut
    merger_redshifts=1./scale_fac_complete-1.
    mask_primary=BH_mass_primary>mass_cut/1e10
    mask_secondary=BH_mass_secondary>mass_cut/1e10
    mask_total_mass=(BH_mass_primary+BH_mass_secondary)>total_mass_cut/1e10
    mask_redshift=merger_redshifts>redshift_cut
    mask=((mask_primary&mask_secondary)&(mask_redshift&mask_ratio_cut))&mask_total_mass 
    maskhost1=hosthalomass_primary+hosthalomass_secondary>14.2e10/1e10*h
    maskhost2=hosthalomass_primary>7.1e10/1e10*h
    maskhost3=hosthalodmmass_primary>2e9/1e10*h
    maskhost4=hosthalodmmass_secondary>2e9/1e10*h
    maskhost5=hosthalostellarmass_primary>1e8/1e10*h   
    maskhost6=hosthalostellarmass_secondary>1e8/1e10*h
    if(apply_host_cuts==0):
        maskhost=mask
    else:
        maskhost=((((maskhost1&maskhost2)&maskhost3)&maskhost4)&maskhost5)&maskhost6    
    merger_redshifts=merger_redshifts[mask&maskhost]   
    bins, dist, dist_err,norm,counts_sum=arepo_package.get_probability_density(merger_redshifts,nbins,min_value,max_value,linear=linear)  
    dt=numpy.array([T(zz-norm/2,zz+norm/2) for zz in bins])
    dz=numpy.diff(bins)[0]
    dz_dt=numpy.array([dz/T(zz-dz/2,zz+dz/2) for zz in bins])  
    return bins,dist,dist_err,norm,dz_dt


#f,axx=plt.subplots(1,4,figsize=(30,8),sharey=True,sharex=True)
def merger_rates(basePath,labl='.',width=8,HDF5=0,path_to_hdf5='.',box_length=1,mark='.',nbins=20,min_value=0, max_value=25):
    bins,dist,dist_err,norm,dz_dt=get_merger_distributions(basePath,nbins,min_value,max_value,mass_cut,linear=1,HDF5=HDF5,path_to_hdf5=path_to_hdf5)
    volume_of_redshift_shell=numpy.array([4*3.14*DC(0,z_bin)**2*dx_dz(z_bin) for z_bin in bins])
    if (HDF5):
        simulation_volume=box_length**3
    else:
        simulation_volume=(arepo_package.get_box_size(basePath)/1000)**3
    merger_rate=dist*1.
    merger_rate_err=dist_err*1.
    merger_rate*=dz_dt
    merger_rate_err*=dz_dt
    merger_rate/=simulation_volume
    merger_rate_err/=simulation_volume
    merger_rate*=volume_of_redshift_shell
    merger_rate_err*=volume_of_redshift_shell
    merger_rate/=(1+bins)
    merger_rate_err/=(1+bins)
    merger_rate2=dist*dz_dt/(1+bins)*volume_of_redshift_shell/simulation_volume
    merger_rate_err2=dist_err*dz_dt/(1+bins)*volume_of_redshift_shell/simulation_volume
    return bins, merger_rate2,merger_rate_err2

