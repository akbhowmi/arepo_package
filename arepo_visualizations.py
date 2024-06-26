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
import itertools
import functools

from scipy.ndimage import gaussian_filter,gaussian_filter1d

from multiprocessing import Pool

import arepo_package

def find_BH_position(basePath,desired_redshift,p_id,choose_most_massive=0): 
    bh_IDs,output_redshift=arepo_package.get_particle_property(basePath,'ParticleIDs',5,desired_redshift,list_all=False)
    bh_positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',5,desired_redshift,list_all=False)
 
    bh_masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',5,desired_redshift,list_all=False)
    #print(bh_positions,bh_masses)
    if (choose_most_massive):
        center=bh_positions[bh_masses==numpy.amax(bh_masses)]
    else:
        center=(bh_positions[bh_IDs==p_id])
    return center

def find_FOF_position(basePath,desired_redshift,choose_most_massive=0):
    fofpos,output_redshift=arepo_package.get_group_property(basePath,'GroupPos',desired_redshift,list_all=True)
    print(output_redshift)
    fofmass,output_redshift=arepo_package.get_group_property(basePath,'GroupMass',desired_redshift,list_all=False)
    if (choose_most_massive):
        maxpos=numpy.where(fofmass==numpy.amax(fofmass))
        center=fofpos[maxpos,:]
    return center

def orient_plane(positions,perpendicular_vector):
    def mag(vector):
        return numpy.sqrt(sum(vector**2))    
    A=numpy.array([1,1,1])
    unit_A=A/mag(A)
    unit_perpendicular_vector=perpendicular_vector/mag(perpendicular_vector)
    
    L=unit_perpendicular_vector      
    LxA=numpy.cross(L,unit_A)
    LxLxA=numpy.cross(L,LxA)
    
    LxA/=mag(LxA)
    LxLxA/=mag(LxLxA)
    
    zpos_new=numpy.dot(positions,L)
    
    xpos_new=numpy.dot(positions,LxA)
    
    ypos_new=numpy.dot(positions,LxLxA)

    new_positions=numpy.transpose(numpy.array([xpos_new,ypos_new,zpos_new]))
    
    return new_positions

def extract_slice(basePath,p_type,desired_center,desired_redshift,field,planeofsky_dimensions=(1000,1000),lineofsight_dimension=20,plane='xy',orient=0,get_full_box=0,only_zoom=0,HighResGasFractionThreshold=0.,angular_momentum_type=4,put_cut=0,cut=1):
    positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',p_type,desired_redshift,list_all=False)    
    if (only_zoom==1):
        HighResGasMass,output_redshift=arepo_package.get_particle_property(basePath,'HighResGasMass',p_type,desired_redshift,list_all=False)    
        masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
        HighResGasFraction=HighResGasMass/masses
    
    if (field=='Masses'):
        Masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
        Masses*=1e10 
        particle_property=Masses
        
    if (field=='Density'):
                    
        if(p_type==1):
            Potential,output_redshift=arepo_package.get_particle_property(basePath,'Potential',p_type,desired_redshift,list_all=False)

            Masses=numpy.array([arepo_package.load_snapshot_header(basePath,desired_redshift)['MassTable'][1]]*len(Potential))
        else:
            Masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
        Masses*=1e10 
        particle_property=Masses
        
    if ('LymanWernerIntensity' in field):
        Masses,output_redshift=arepo_package.get_particle_property(basePath,field,p_type,desired_redshift,list_all=False)
        particle_property=Masses
        
    if (field=='Metallicity'):
        metallicity,output_redshift=arepo_package.get_particle_property(basePath,'GFM_Metallicity',p_type,desired_redshift,list_all=False)  
        metallicity/=0.0127
        #print(metallicity)
        particle_property=metallicity     
        
    if (field=='SFR'):
        SFR,output_redshift=arepo_package.get_particle_property(basePath,'StarFormationRate',p_type,desired_redshift,list_all=False)  
        particle_property=SFR  
        
    if (field=='Temperature'):
        internal_energy,output_redshift=arepo_package.get_particle_property(basePath,'InternalEnergy',p_type,desired_redshift,list_all=False)
        electron_abundance,output_redshift=arepo_package.get_particle_property(basePath,'ElectronAbundance',p_type,desired_redshift,list_all=False)
        g_minus_1 = (5.0/3.0) - 1.0
        XH = 0.76
        mu=(4*ac.MP)/(1+3*XH+4*XH*electron_abundance)
        gas_temperature = g_minus_1*(internal_energy/ac.KB)*(10**10)*mu
        #print(gas_temperature)
        particle_property=gas_temperature
    if ((field=='Velocity Magnitude') | (field=='Velocities')):
        velocity,output_redshift=arepo_package.get_particle_property(basePath,'Velocities',p_type,desired_redshift,list_all=False) 
        if (field=='Velocity Magnitude'):
            xvel=velocity[:,0]
            yvel=velocity[:,1]
            zvel=velocity[:,2]
            mag_vel=numpy.sqrt((xvel**2)+(yvel**2)+(zvel**2))
            particle_property=mag_vel
        if (field=='Velocities'):
            particle_property=velocity
    print(output_redshift)
    positions_relative_to_center=positions-numpy.ravel(desired_center)
      
    if (orient):
        perpendicular_vector=stellar_ang_mom(basePath,desired_redshift,desired_center,p_type=angular_momentum_type)
        positions_relative_to_center=orient_plane(positions_relative_to_center,perpendicular_vector)
         
    if (plane=='xy'):
        planeofsky_pos1=positions_relative_to_center[:,0]
        planeofsky_pos2=positions_relative_to_center[:,1]
        lineofsight_pos=positions_relative_to_center[:,2]
        
    if (plane=='xz'):
        planeofsky_pos1=positions_relative_to_center[:,0]
        planeofsky_pos2=positions_relative_to_center[:,2]
        lineofsight_pos=positions_relative_to_center[:,1]
    
    if (plane=='yz'):
        planeofsky_pos1=positions_relative_to_center[:,1]
        planeofsky_pos2=positions_relative_to_center[:,2]
        lineofsight_pos=positions_relative_to_center[:,0]
        
    observable_positions=numpy.transpose(numpy.array([planeofsky_pos1,planeofsky_pos2,lineofsight_pos]))
    
    mask_planeofsky_1=numpy.abs(planeofsky_pos1)<planeofsky_dimensions[0]/2
    mask_planeofsky_2=numpy.abs(planeofsky_pos2)<planeofsky_dimensions[1]/2
    mask_lineofsight=numpy.abs(lineofsight_pos)<lineofsight_dimension/2
    mask=(mask_planeofsky_1&mask_planeofsky_2)&mask_lineofsight
    
    if (only_zoom):
        #print(HighResGasFraction)
        mask_zoom=HighResGasFraction>HighResGasFractionThreshold
        mask=mask&mask_zoom
        mask_lineofsight=mask_lineofsight&mask_zoom
    if (put_cut==1):
        mask=mask&cut
        
    
    if (get_full_box):
        return observable_positions[mask_lineofsight],particle_property[mask_lineofsight],output_redshift  
        
    return observable_positions[mask],particle_property[mask],output_redshift

def apply_mask(mask2,planeofsky_pos1,pixelsize_planeofsky1,final_property,pixel_volume,pixel_position1):
    mask1=(planeofsky_pos1>(pixel_position1-pixelsize_planeofsky1/2.0))&(planeofsky_pos1<(pixel_position1+pixelsize_planeofsky1/2.))
    mask=(mask1) & (mask2)
    return numpy.sum(final_property[mask])/pixel_volume

def apply_mask2(field,planeofsky_pos2,pixelsize_planeofsky2,planeofsky_pos1,pixelsize_planeofsky1,final_property,pixel_volume,planeofsky_tuple):
    pixel_position2=planeofsky_tuple[0]
    pixel_position1=planeofsky_tuple[1]
    mask1=(planeofsky_pos1>(pixel_position1-pixelsize_planeofsky1/2.0))&(planeofsky_pos1<(pixel_position1+pixelsize_planeofsky1/2.))
    mask2=(planeofsky_pos2>(pixel_position2-pixelsize_planeofsky2/2.0))&(planeofsky_pos2<(pixel_position2+pixelsize_planeofsky2/2.))
    mask=(mask1) & (mask2)
    if (field=='Density'):
        return numpy.sum(final_property[mask])/pixel_volume
    else:
        return numpy.sum(final_property[mask])/(len(final_property[mask])+1e-19)



def construct_grid(final_positions,final_property,number_of_pixels,field='.',PARALLEL=0,np=1):
    if (PARALLEL==0):
        print("property:",final_property[final_property>1e-12])    

    n_planeofsky1=number_of_pixels
    n_planeofsky2=number_of_pixels
    planeofsky_pos1=final_positions[:,0]
    planeofsky_pos2=final_positions[:,1]
    lineofsight_pos=final_positions[:,2]
    

    
    planeofsky1_grid=numpy.linspace(min(planeofsky_pos1),max(planeofsky_pos1),n_planeofsky1)
    planeofsky2_grid=numpy.linspace(min(planeofsky_pos2),max(planeofsky_pos2),n_planeofsky2)
    pixelsize_planeofsky1=numpy.diff(planeofsky1_grid)[0]
    pixelsize_planeofsky2=numpy.diff(planeofsky2_grid)[0]
    pixelsize_lineofsight=numpy.amax(lineofsight_pos)-numpy.amin(lineofsight_pos)

    proj_property=[]
    pixel_volume=1.
    if (field=='Density') | (field=='SFR'):
        pixel_volume=pixelsize_lineofsight*pixelsize_planeofsky1*pixelsize_planeofsky2
 
    First,Second=numpy.meshgrid(planeofsky1_grid,planeofsky2_grid)
    if (PARALLEL==1):
        proj_property=numpy.array(proj_property)
        p=Pool(np)
        for pixel_position2 in planeofsky2_grid:
            mask2=(planeofsky_pos2>(pixel_position2-pixelsize_planeofsky2/2.0))&(planeofsky_pos2<(pixel_position2+pixelsize_planeofsky2/2.))    
            proj_property=numpy.append(proj_property,p.map(functools.partial(apply_mask,mask2,planeofsky_pos1,pixelsize_planeofsky1,final_property,pixel_volume),planeofsky1_grid))
        p.close()

    elif (PARALLEL==2) :
        p=Pool(np)
        planeofsky_tuple_grid=list(itertools.product(planeofsky2_grid,planeofsky1_grid))
        proj_property=numpy.array(p.map(functools.partial(apply_mask2,field,planeofsky_pos2,pixelsize_planeofsky2,planeofsky_pos1,pixelsize_planeofsky1,final_property,pixel_volume),planeofsky_tuple_grid))
        p.close()
        
    else:
        for pixel_position2 in planeofsky2_grid:
            mask2=(planeofsky_pos2>(pixel_position2-pixelsize_planeofsky2/2.0))&(planeofsky_pos2<(pixel_position2+pixelsize_planeofsky2/2.))
            for pixel_position1 in planeofsky1_grid:
                mask1=(planeofsky_pos1>(pixel_position1-pixelsize_planeofsky1/2.0))&(planeofsky_pos1<(pixel_position1+pixelsize_planeofsky1/2.))

                mask=(mask1) & (mask2)
                if (field=='Density') | (field=='SFR'):
                    proj_property.append(numpy.sum(final_property[mask])/pixel_volume)
                elif(field=='Masses'): 
                    proj_property.append(numpy.sum(final_property[mask]))
                else:
                    proj_property.append(numpy.average(final_property[mask]))


 

    
    proj_property=numpy.asarray(proj_property)
    print(proj_property)
    #if (field=='Metallicit'):
    #    proj_property[(numpy.isnan(proj_property)) | (proj_property==0)]=min(proj_property[proj_property>0])
    #else:
    #    proj_property[(numpy.isnan(proj_property)) | (proj_property==0)]=1e-19
    #print(proj_property)
    Proj_property=proj_property.reshape(number_of_pixels,number_of_pixels)

   # proj_property[(numpy.isnan(proj_property)) | (proj_property==0)]=1e-19 
    return First,Second,Proj_property

def filter_property(First,Second,proj_property,field,apply_filter=1,sigma_filter=1,show_colorbar=1,colormap='Greys_r',valuemin=0.1,valuemax=1e6,alph=1,show_labels=1,show_ticks=1,cbar_orient='horizontal',variable_filter=0,neighborhood_size=100,characteristic_mass=1,pixel_size=1,variable_filter_property='',sigma_fac=5,no_periodic=0):
    def compute_local_sigma(arr, i, j,characteristic_mass=characteristic_mass):
        density = arr[i, j]
        size_of_cell=pow(characteristic_mass/density,1./3)
        size_of_cell_pixels=size_of_cell/pixel_size
        return size_of_cell_pixels
    if (apply_filter):
        if(variable_filter==0):
            Proj_property=gaussian_filter(proj_property,sigma=sigma_filter)
        else:



            # Create an array to store the filtered values
            filtered_property = numpy.zeros_like(proj_property)
            size_of_cell=pow(characteristic_mass/variable_filter_property,1./3) * 3./4./3.14 
            local_sigma_all_un=sigma_fac*size_of_cell/pixel_size
            local_sigma_all = gaussian_filter(local_sigma_all_un, sigma=0.3)
            # Iterate through each point in proj_property
            for i in range(proj_property.shape[0]):
                if(i%100==0):
                    print(i)
                for j in range(proj_property.shape[1]):
                    #if(j%1==0):
                        #print(j)
                    # Compute the local sigma based on the local density
                    local_sigma = local_sigma_all[i,j]
                    if(local_sigma>10):
                        local_sigma=10
                    if(local_sigma==0):
                        local_sigma=0.5
                    # Create a neighborhood around the current point
                    neighborhood = numpy.roll(proj_property,
                               (neighborhood_size//2 - i, neighborhood_size//2 - j), 
                               axis=(0, 1))[0:neighborhood_size, 0:neighborhood_size]
                    
                    if(no_periodic):
                        start_row = max(0, i - neighborhood_size // 2)
                        end_row = min(proj_property.shape[0], i + neighborhood_size // 2 + 1)
                        start_col = max(0, j - neighborhood_size // 2)
                        end_col = min(proj_property.shape[1], j + neighborhood_size // 2 + 1)

                        # Extract the non-periodic neighborhood
                        neighborhood = proj_property[start_row:end_row, start_col:end_col]
                    
                    



                    # Apply Gaussian filter with the computed sigma
                    filtered_value = gaussian_filter(neighborhood, sigma=local_sigma)

                    # Take the average of the filtered neighborhood
                    filtered_property[i, j] = filtered_value[neighborhood_size//2, neighborhood_size//2]
            Proj_property=filtered_property
    else:
         Proj_property=proj_property
    return First,Second,Proj_property
            
            
def visualize(First,Second,Proj_property,field,fig,ax,apply_filter=1,sigma_filter=1,show_colorbar=1,colormap='Greys_r',valuemin=0.1,valuemax=1e6,alph=1,show_labels=1,show_ticks=1,cbar_orient='horizontal',variable_filter=0,neighborhood_size=100,characteristic_mass=1,pixel_size=1):            
    if (field=='Density'):
        print("making density")
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(proj_property[proj_property>0]),vmax=Proj_property.max()),cmap=colormap)
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=mpl.colors.LogNorm(vmin=min(proj_property[proj_property>0]),vmax=max(proj_property[proj_property>0])),cmap=colormap)
        fig_object=ax.pcolor(First,Second,Proj_property,norm=mpl.colors.LogNorm(vmin=valuemin,vmax=valuemax),cmap=colormap,alpha=alph)
        if (show_colorbar):
            cbar=fig.colorbar(fig_object,ax=ax,orientation=cbar_orient)
            cbar.set_label(r'$\rho_{gas}$ ($M_{\odot}h^{3}kpc^{-3}$)',fontsize=30)
            cbar.ax.tick_params(labelsize=30) 
        #bh=plt.Circle((bhcoord[0],bhcoord[1]),0.5,color='black')
        #ax.add_artist(bh)

        
    if ('LymanWernerIntensity' in field):
        print('making LymanWernerIntensity')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='plasma')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=mpl.colors.LogNorm(vmin=valuemin,vmax=valuemax),cmap=colormap,alpha=alph)
        if (show_colorbar):
            cbar=fig.colorbar(fig_object,ax=ax,orientation=cbar_orient)
            if ('type2' in field):
                cbar.set_label('$J/J_{21}$ (Pop II)',fontsize=30)
            if ('type3' in field):
                cbar.set_label('$J/J_{21}$ (Pop III)',fontsize=30)
            cbar.ax.tick_params(labelsize=30)
        
        
    if (field=='Metallicity'):
        print('making metal')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='plasma')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=mpl.colors.LogNorm(vmin=valuemin,vmax=valuemax),cmap=colormap,alpha=alph)
        if (show_colorbar):
            cbar=fig.colorbar(fig_object,ax=ax,orientation=cbar_orient)
            cbar.set_label('$Z/Z_{\odot}$',fontsize=30)
            cbar.ax.tick_params(labelsize=30)

        #cbar.ax.set_ylim([cbar.norm(10e1),cbar.norm(10e-7)])

    if (field=='Temperature'):
        print('making temperature')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='hot')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=1e2,vmax=1e7),cmap='hot')
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('log gas temperature $K$',fontsize=15)
 
    if (field=='Velocity Magnitude'):
        print('making velocity mag')
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='viridis')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=valuemin,vmax=valuemax),cmap=colormap,alpha=alph)
        cbar=fig.colorbar(fig_object,ax=ax)
        cbar.set_label('SFR density \n ($M_{\odot}/\mathrm{yr}^{-1}/\mathrm{kpc}^{-3}$)',fontsize=15)
    if (field=='SFR'):
        print('making SFR')
        
        #fig_object=ax.pcolor(First,Second,Proj_property,norm=colors.LogNorm(vmin=min(final_property),vmax=Proj_property.max()),cmap='viridis')
        fig_object=ax.pcolor(First,Second,Proj_property,norm=mpl.colors.LogNorm(vmin=valuemin,vmax=valuemax),cmap=colormap,alpha=alph)
        if (show_colorbar):
            cbar=fig.colorbar(fig_object,ax=ax,orientation=cbar_orient)
            cbar.set_label('SFR density \n ($M_{\odot} \ \mathrm{yr}^{-1} \ \mathrm{kpc}^{-3}$)',fontsize=30)
            cbar.ax.tick_params(labelsize=30)    
   
    
    if (show_labels):
        ax.set_xlabel('Plane of sky 1 ($h^{-1}kpc$)',fontsize=30)
        ax.set_ylabel('Plane of sky 2 ($h^{-1}kpc$)',fontsize=30)
    if (show_ticks==0):
        ax.tick_params(labelleft=False,labelbottom=False)


    
def stellar_ang_mom(basePath,desired_redshift,desired_center,box_length=40,p_type=4):
    positions,output_redshift=arepo_package.get_particle_property(basePath,'Coordinates',p_type,desired_redshift,list_all=False)       
    masses,output_redshift=arepo_package.get_particle_property(basePath,'Masses',p_type,desired_redshift,list_all=False)
    velocities,output_redshift=arepo_package.get_particle_property(basePath,'Velocities',p_type,desired_redshift,list_all=False)
    masses*=1e10
    velocities/=3.0886e16
    
    xvel=velocities[:,0] #units switched from km/s to kpc/s
    yvel=velocities[:,1]
    zvel=velocities[:,2]

    xpos=positions[:,0]
    ypos=positions[:,1]
    zpos=positions[:,2]
    
    center=desired_center[0]
    center = numpy.ravel(center) 
    mask_x=numpy.abs(xpos-center[0])<box_length/2
    mask_y=numpy.abs(ypos-center[1])<box_length/2
    mask_z=numpy.abs(zpos-center[2])<box_length/2
    mask=(mask_x&mask_y)&mask_z
    xvel=xvel[mask]
    yvel=yvel[mask]
    zvel=zvel[mask]
    
    xpos=xpos[mask]
    ypos=ypos[mask]
    zpos=zpos[mask]

    masses=masses[mask]
    
    #Calculate angular momentum
    Lx=(ypos*zvel-zpos*yvel)*masses
    Lx=sum(Lx)
    Ly=(xpos*zvel-zpos*xvel)*masses
    Ly=sum(Ly)
    Lz=(xpos*yvel-ypos*xvel)*masses
    Lz=sum(Lz)    
    L_total=numpy.transpose(numpy.array([Lx,Ly,Lz]))
     
    p_particles=numpy.transpose(numpy.array([xvel*masses,yvel*masses,zvel*masses]))


    p_gal=numpy.sum(p_particles,axis=0)

    num=numpy.transpose(numpy.array([xpos*masses,ypos*masses,zpos*masses]))
    com=numpy.sum(num,axis=0)/numpy.sum(masses)
    

    L_orbital=numpy.cross(com,p_gal)

    L_spin=L_total-L_orbital


    return L_spin
