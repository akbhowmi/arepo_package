import sys
from Config import *
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
#sys.path.append('/home/aklantbhowmick/.local/lib/python3.12/site-packages/')
#sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.11//site-packages/')

import mdot_to_Lbol
import arepo_package
import scipy.interpolate
radiative_efficiency=0.2
total_conv=mdot_to_Lbol.get_conversion_factor_arepo(radiative_efficiency)
import h5py
import illustris_python as il
import numpy 
import pandas as pd
import scipy.integrate
from scipy.optimize import curve_fit
import os

title_fontsize=30
import warnings
warnings.filterwarnings("ignore")


from scipy.constants import golden_ratio
#import legwork

import astropy.constants as const
import time
import astropy.units as u

import gwent

import lal
import gwent.binary as binary
import gwent.detector as detector
import gwent.snr as snr
from gwent.snrplot import Plot_SNR

#Turn off warnings for tutorial
import warnings
warnings.filterwarnings('ignore')

load_directory = gwent.__path__[0] + '/LoadFiles/InstrumentFiles/'

def Initialize_Source(instrument,approximant='pyPhenomD',lalsuite_kwargs={}):
    """Initializes a source binary based on the instrument type and returns the source

    Parameters
    ----------
    instrument : object
        Instance of a gravitational wave detector class
    approximant : str, optional
        the approximant used to calculate the frequency domain waveform of the source.
        Can either be the python implementation of IMRPhenomD ('pyPhenomD', the default) given below,
        or a waveform modelled in LIGO's lalsuite's lalsimulation package.
    lalsuite_kwargs: dict, optional
        More specific user-defined kwargs for the different lalsuite waveforms
    """

    #q = m2/m1 reduced mass
    q = 1.0
    q_min = 1.0
    q_max = 18.0
    q_list = [q,q_min,q_max]

    #Chi = S_i*L/m_i**2, spins of each mass i
    chi1 = 0.0 #spin of m1
    chi2 = 0.0 #spin of m2
    chi_min = -0.85 #Limits of PhenomD for unaligned spins
    chi_max = 0.85
    chi1_list = [chi1,chi_min,chi_max]
    chi2_list = [chi2,chi_min,chi_max]

    #Redshift
    z_min = 1e-2
    z_max = 1e3

    if isinstance(instrument,detector.GroundBased):
        #Total source mass
        M_ground_source = [10.,1.,1e4]
        #Redshift
        z_ground_source = [0.1,z_min,z_max]

        source = binary.BBHFrequencyDomain(M_ground_source,
                                           q_list,
                                           z_ground_source,
                                           chi1_list,
                                           chi2_list,
                                           approximant=approximant,
                                           lalsuite_kwargs=lalsuite_kwargs)
    elif isinstance(instrument,detector.SpaceBased):
        M_space_source = [1e6,1.,1e10]
        z_space_source = [1.0,z_min,z_max]
        source = binary.BBHFrequencyDomain(M_space_source,
                                           q_list,
                                           z_space_source,
                                           chi1_list,
                                           chi2_list,
                                           approximant=approximant,
                                           lalsuite_kwargs=lalsuite_kwargs)
    elif isinstance(instrument,detector.PTA):
        M_pta_source = [1e9,1e8,1e11]
        z_pta_source = [0.1,z_min,z_max]
        source = binary.BBHFrequencyDomain(M_pta_source,
                                           q_list,
                                           z_pta_source,
                                           chi1_list,
                                           chi2_list,
                                           approximant=approximant,
                                           lalsuite_kwargs=lalsuite_kwargs)
    return source


def Initialize_aLIGO():
    #Observing time in years
    T_obs_ground_list = [4*u.yr,1*u.yr,10*u.yr]
    #aLIGO
    noise_dict_aLIGO = {'Infrastructure':
                  {'Length':[3995,2250,4160]},
                  'Laser':
                  {'Power':[125,10,1e3]},
                  'Seismic':
                  {'Gamma':[0.8,0.1,1.0]}}
    aLIGO = detector.GroundBased('aLIGO',T_obs_ground_list,noise_dict=noise_dict_aLIGO)

    return aLIGO


def Initialize_LISA():
    #Values taken from the ESA L3 proposal, Amaro-Seaone, et al., 2017 (https://arxiv.org/abs/1702.00786)
    T_obs_space_list = [4*u.yr,1*u.yr,10*u.yr]

    #armlength in meters
    L = 2.5e9*u.m
    L_min = 1.0e7*u.m
    L_max = 1.0e11*u.m
    L_list = [L,L_min,L_max]

    #Acceleration Noise Amplitude
    A_acc = 3e-15*u.m/u.s/u.s
    A_acc_min = 1e-16*u.m/u.s/u.s
    A_acc_max = 1e-14*u.m/u.s/u.s
    A_acc_list = [A_acc,A_acc_min,A_acc_max]

    #The Low Acceleration Noise Break Frequency
    f_acc_break_low = .4*u.mHz.to('Hz')*u.Hz
    f_acc_break_low_min = .1*u.mHz.to('Hz')*u.Hz
    f_acc_break_low_max = 1.0*u.mHz.to('Hz')*u.Hz
    f_acc_break_low_list = [f_acc_break_low,f_acc_break_low_min,f_acc_break_low_max]

    #The High Acceleration Noise Break Frequency
    f_acc_break_high = 8.*u.mHz.to('Hz')*u.Hz
    f_acc_break_high_min = 1.*u.mHz.to('Hz')*u.Hz
    f_acc_break_high_max = 10.*u.mHz.to('Hz')*u.Hz
    f_acc_break_high_list = [f_acc_break_high,f_acc_break_high_min,f_acc_break_high_max]

    #The Optical Metrology Noise Break Frequency
    f_IFO_break = 2.*u.mHz.to('Hz')*u.Hz
    f_IFO_break_min = 1.*u.mHz.to('Hz')*u.Hz
    f_IFO_break_max = 10.*u.mHz.to('Hz')*u.Hz
    f_IFO_break_list = [f_IFO_break,f_IFO_break_min,f_IFO_break_max]

    #Detector Optical Metrology Noise
    A_IFO = 10e-12*u.m
    A_IFO_min = 1.0e-13*u.m
    A_IFO_max = 1.0e-10*u.m
    A_IFO_list = [A_IFO,A_IFO_min,A_IFO_max]

    #Unresolved Galactic WD Background
    Background = False

    #Numerical Transfer Function
    T_type = 'N'

    LISA_prop1 = detector.SpaceBased('LISA_prop1',
                                     T_obs_space_list,L_list,A_acc_list,
                                     f_acc_break_low_list,f_acc_break_high_list,
                                     A_IFO_list,f_IFO_break_list,
                                     Background=Background,T_type=T_type)
    return LISA_prop1

def Initialize_NANOGrav():
    #NANOGrav calculation using 11.5yr parameters https://arxiv.org/abs/1801.01837
    #Observing time in years
    T_obs_ptas_list = [11.42*u.yr,5*u.yr,30*u.yr]
    #rms timing residuals in seconds
    sigma = 100*u.ns.to('s')*u.s
    sigma_min = 100*u.ns.to('s')*u.s
    sigma_max = 500*u.ns.to('s')*u.s
    sigma_list = [sigma,sigma_min,sigma_max]
    #Number of pulsars
    n_p = 34
    n_p_min = 18
    n_p_max = 200
    n_p_list = [n_p,n_p_min,n_p_max]
    #Avg observation cadence of 1 every 2 weeks in num/year
    cadence = 1/(2*u.wk.to('yr')*u.yr)
    cadence_min = 2/u.yr
    cadence_max = 1/(u.wk.to('yr')*u.yr)
    cadence_list = [cadence,cadence_min,cadence_max]

    #NANOGrav 11.4 yr WN only
    NANOGrav_WN = detector.PTA('NANOGrav_WN',n_p_list,T_obs=T_obs_ptas_list,sigma=sigma_list,cadence=cadence_list)
    return NANOGrav_WN

hubble=0.6771
scale_fac_complete_sorted,primary_mass_sorted,secondary_mass_sorted,primary_id_sorted,secondary_id_sorted,file_id_complete_sorted,N_empty=arepo_package.get_merger_events_from_snapshot(basePath,upto_redshift,HOSTS=0)

redshifts_complete = 1./scale_fac_complete_sorted-1
mass_tuple=numpy.transpose([primary_mass_sorted,secondary_mass_sorted])
ratios_complete=numpy.amax(mass_tuple,axis=1)/numpy.amin(mass_tuple,axis=1)
mass_total_complete=numpy.sum(mass_tuple,axis=1)*1e10/hubble


chi1=chi2=0
snr_space=[]
i=0
for Mtotal,q,z in list(zip(mass_total_complete,ratios_complete,redshifts_complete)):
        source=binary.BBHFrequencyDomain(Mtotal,q,z,chi1,chi2)
        source.h_f
        LISA_detector=Initialize_LISA()
        snr_space.append(snr.Calc_Chirp_SNR(source,LISA_detector))
        if(i%100==0):
            print(i)
        i+=1

snr_space=numpy.array(snr_space)

numpy.save(basePath + './SNRs_for_merger_events_zero_spin.npy',[mass_total_complete,ratios_complete,redshifts_complete,snr_space])




