import sys
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import mdot_to_Lbol
import arepo_package
import scipy.interpolate
import numpy
import os

def extract_ids(ids_to_be_extracted,complete_ids):
    #print(complete_ids)
    #print(ids_to_be_extracted)
    return numpy.where(complete_ids==ids_to_be_extracted)[0]
    
    
vec_extract_ids=numpy.vectorize(extract_ids,excluded=['complete_ids'])

#run='L25_n256'
#basePath='/ufrc/lblecha/aklantbhowmick/arepo_runs_aklant/'+run+'/output/'
#save_output_path='/ufrc/lblecha/aklantbhowmick/richness_output_AGN_activity/'

run='L205n2500TNG'
#p_type=5
#basePath='/ufrc/lblecha/aklantbhowmick/arepo_runs_aklant/'+run+'/output/'

basePath='/n/hernquistfs3/IllustrisTNG/Runs/'+run+'/output/'
save_output_path='/n/boslfs/TRANSFER/aklantbhowmick/richness_output_AGN_activity/'
load_output_path='/n/home00/aklantbhowmick/illustris_processing/FOF_tag_outputs/'


if not os.path.exists(save_output_path):
        print("Making directory for storing richenss")
        os.makedirs(save_output_path)

radiative_efficiency=0.1
total_conv=mdot_to_Lbol.get_conversion_factor_arepo(radiative_efficiency)

redshift_space=[3.]
log_LINKING_LENGTH_space=numpy.linspace(-2,2,40)
log_luminosity_cuts_space=[0,42.0,43.0,44.0,45.0]
log_bhmass_cuts_space=[6]

#log_lum_cut=.
#log_luminosity_cuts_space=[42.5,42.0]


#log_bhmass_cut=[6.]

for redshift in redshift_space:
    p_type=5
    particle_property='BH_Mdot'
    desired_redshift=redshift
    BH_Mdot,output_redshift=arepo_package.get_particle_property(basePath,'BH_Mdot', p_type, desired_redshift, list_all=False)
    BH_Mass,output_redshift=arepo_package.get_particle_property(basePath,'BH_Mass', p_type, desired_redshift, list_all=False)

    ParticleIDs,output_redshift=arepo_package.get_particle_property(basePath,'ParticleIDs', p_type, desired_redshift, list_all=False)
    bolometric_luminosity=BH_Mdot*total_conv 
    BH_Mass=BH_Mass*1e10
    
    for log_bhmass_cut in log_bhmass_cuts_space:
#        for log_LINKING_LENGTH in log_LINKING_LENGTH_space:
        scale_of_interest_space=numpy.array([0.015,0.1,1.0])
        for scale_of_interest in scale_of_interest_space:
        #for log_LINKING_LENGTH in log_LINKING_LENGTH_space:
                log_scale_of_interest=log_LINKING_LENGTH_space[log_LINKING_LENGTH_space<numpy.log10(scale_of_interest)][-1]
                log_LINKING_LENGTH=log_scale_of_interest
                final_FOF_tag,subsampled_quantities=numpy.load(load_output_path+run+'_log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f.npy'%(log_bhmass_cut,redshift,log_LINKING_LENGTH),allow_pickle=True)
                #print(final_FOF_tag)
                position_vector_cut,velocity_vector_gadget_units_cut,blackhole_id_cut,subhalo_id_cut,hosthalo_id_cut=subsampled_quantities
                unique_FOF_tag=numpy.unique(final_FOF_tag)
                RICHNESS=[]
                RICHNESS_active_42=[]
                RICHNESS_active_43=[]
                RICHNESS_active_44=[]
                RICHNESS_active_45=[]
                RICHNESS_BH_70=[]
                RICHNESS_BH_80=[]
                RICHNESS_BH_90=[]
                RICHNESS_BH_100=[]

                average_log_bolometric_luminosity=[]
                average_log_BH_Mass=[]
                for tag in unique_FOF_tag:
                    extract_FOFs=(final_FOF_tag==tag)
                    ids_of_members=blackhole_id_cut[extract_FOFs]
                    bolometric_luminosity_members=bolometric_luminosity[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]
#                    extract_active_members=bolometric_luminosity_members>10**log_lum_cut                    


                    BH_Mass_members=BH_Mass[vec_extract_ids(ids_to_be_extracted=ids_of_members,complete_ids=ParticleIDs)]
                    
                    
                    RICHNESS.append(len(ids_of_members))
                    extract_active_members=bolometric_luminosity_members>10**42
                    RICHNESS_active_42.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**43
                    RICHNESS_active_43.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**44
                    RICHNESS_active_44.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=bolometric_luminosity_members>10**45
                    RICHNESS_active_45.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=BH_Mass_members>10**7
                    RICHNESS_BH_70.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=BH_Mass_members>10**8
                    RICHNESS_BH_80.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=BH_Mass_members>10**9
                    RICHNESS_BH_90.append(len(ids_of_members[extract_active_members]))

                    extract_active_members=BH_Mass_members>10**10
                    RICHNESS_BH_100.append(len(ids_of_members[extract_active_members]))

                                    
                RICHNESS=numpy.array(RICHNESS)
                RICHNESS_active_42=numpy.array(RICHNESS_active_42)   
                RICHNESS_active_43=numpy.array(RICHNESS_active_43)
                RICHNESS_active_44=numpy.array(RICHNESS_active_44)
                RICHNESS_active_45=numpy.array(RICHNESS_active_45)

                RICHNESS_BH_70=numpy.array(RICHNESS_BH_70)
                RICHNESS_BH_80=numpy.array(RICHNESS_BH_80)
                RICHNESS_BH_90=numpy.array(RICHNESS_BH_90)
                RICHNESS_BH_100=numpy.array(RICHNESS_BH_100)


                #average_log_bolometric_luminosity=numpy.array(average_log_bolometric_luminosity)
                #average_log_BH_Mass=numpy.array(average_log_BH_Mass)
                
                numpy.save(save_output_path+run+'log_bhmass_cut_%.1f_redshift_%.2f_log_LINKING_LENGTH_%.2f_all_info.npy'%(log_bhmass_cut,redshift,log_LINKING_LENGTH),[RICHNESS,RICHNESS_active_42,RICHNESS_active_43,RICHNESS_active_44,RICHNESS_active_45,RICHNESS_BH_70,RICHNESS_BH_80,RICHNESS_BH_90,RICHNESS_BH_100])
        #print(RICHNESS)
        #print(RICHNESS_active)
