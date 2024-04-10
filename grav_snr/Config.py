Mh,Msfmp=1000,5
path_to_output='/orange/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS/L12p5n512//'
run='/AREPO/' # name of the simulation runs
#output='/output_ratio%d_SFMFGM%d_seed3.19_bFOF/'%(Mh,Msfmp)
#basePath=path_to_output+run+output

output='/output_Msfmp%d_Mh%d_double_powerlaw_env_5_0.1_0.3_slope1.6/'%(Msfmp,Mh)
basePath=path_to_output+run+output
#path_to_output='/blue/lblecha/aklantbhowmick/CURRENT_RUNS2/L12p5n512/'
#run='/AREPO/'
#output= '/output_ratio1000_SFMFGM5_seed5.00_bFOF_LW10_spin_rich/'
#basePath=path_to_output+run+output
upto_redshift=0   #keep it at z=0 and it'll naturally choose the latest snapshot available

