#!/bin/bash
#SBATCH --export=NONE
##SBATCH --ntasks-per-node=1
#SBATCH --account=lblecha
#SBATCH --qos=lblecha-b
##SBATCH --reservation=bluetest
#SBATCH --time=40:00:00
#SBATCH --job-name=test
#SBATCH --output=cs_tagging.stdout
#SBATCH --error=cs_tagging.stderr
#SBATCH --mem-per-cpu=64000mb
##SBATCH --time=48:00:00
#SBATCH --mail-user=akbhowmi@andrew.cmu.edu
#SBATCH --mail-type=END,FAIL

#python get_subhalo_mass_luminosity_relations.py
#python connect_black_hole_systems_to_AGN_activity_and_massive_BH_new_max_lum_bound_median_all_eddington.py
#python finding_blackhole_systems_MBII.py

#python get_AGN_luminosity_functions.py
#python get_particle_history.py
#python generate_group_ids.py
#python generate_mass_functions_with_bolometric_cuts.py
#python convert_blackhole_data_to_hdf5.py
#python create_central_satellite_tags_max_lum.py
#python generate_maximum_gas_denisty_for_FOFs_with_BHs.py
#python connect_black_hole_systems_to_AGN_activity_and_massive_BH_new_max_lum_bound_median_all_eddington_mass_ratio_cut_random.py
#python get_metallicity_vs_redshift.py
#python get_SFmass_vs_redshift.py
#python get_halomass_SFRZmass_relation.py
#python create_merger_trees2.py
#python generate_DM_gas_spins.py
#python Analyse_merger_trees_4.py $1
#python Finalise_merger_trees.py $1
#python convert_subhalo_to_halo_merger_tree.py	
#python stack_halo_properties.py 
#python get_offsets.py $1
#python get_halo_environment.py $1
#python finding_gas_clumps_smuggle.py $1 
#python calculate_gas_clump_properties.py $1
echo "generating progenitors"
#python generate_progenitors.py
echo "applying delay"
python apply_delay.py
#python apply_delay_triple_induced_short.py
echo "generating the revised masses" 
#python generate_new_BH_masses.py
