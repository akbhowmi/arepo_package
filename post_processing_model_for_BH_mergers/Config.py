#path_to_output='/orange/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS///' # this is the folder containing the simulation run
#run='/L3p125n512/AREPO/' # name of the simulation runs
#output='/output_ratio1000_SFMFGM5_seed3.19_bFOF/' 
#basePath=path_to_output+run+output

Mh=1000
Msfmp=5
#path_to_output='/blue/lblecha/aklantbhowmick/CURRENT_RUNS2/L25n1024//' # this is the folder containing the simulation run
path_to_output='/orange/lblecha/aklantbhowmick/GAS_BASED_SEED_MODEL_UNIFORM_RUNS/L25n1024//'
run='/AREPO/' # name of the simulation runs
output='/output_Msfmp%d_Mh%d_double_powerlaw_env_5_0.1_0.3_slope1.6/'%(Msfmp,Mh)
basePath=path_to_output+run+output
    #make_plot(axx[0],axx[1],axx[2],style='solid',width=8,col='C1',plot1=True,plot2=True,plot3=True)


desired_redshift=7

new_directory = "postprocessing_mergers"

omk=0
oml=0.6911
omm=0.3089
c=3*10**5
h=0.6771

FIXED_DELAY=1
fixed_delay_Myr=100
triple_induced_delay=125


RANDOM_DELAY=0
min_delay=0
max_delay=0
	

