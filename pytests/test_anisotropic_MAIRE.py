
import datetime as dt

import numpy as np
from anisotropicMAIRE import anisotropicMAIRE

def test_run_from_spectra():

    anisotropicMAIRE.run_from_spectra(proton_rigidity_spectrum=lambda x:1,
                     alpha_rigidity_spectrum=lambda x:1,
                     Kp_index=6,
                     date_and_time=dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0))
    
def test_run_from_spectra_two_locations():

    anisotropicMAIRE.run_from_spectra(proton_rigidity_spectrum=lambda x:1,
                     alpha_rigidity_spectrum=lambda x:1,
                     Kp_index=6,
                     date_and_time=dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0),
                     altitudes_in_km=np.append(np.array(range(0,13)) * 0.3048,11.28),
                     array_of_lats_and_longs=[[46.2,187.4],[-28.3,-92.7]])
    
def test_run_from_spectra_proton_only():

    anisotropicMAIRE.run_from_spectra(proton_rigidity_spectrum=lambda x:1,
                     Kp_index=6,
                     date_and_time=dt.datetime(year=2000,month=7,day=14,hour=10,minute=45,second=0))
    
def test_Mishev_spec():
    anisotropicMAIRE.run_from_power_law_gaussian_distribution(J0=256_000.0, gamma=3.41, deltaGamma=0.22, sigma=np.sqrt(0.19),
                                         reference_pitch_angle_latitude=-17.0,reference_pitch_angle_longitude=148.0,
                                         Kp_index=3, date_and_time=dt.datetime(2006, 12, 13, 3, 0),
                                         altitudes_in_km=np.append(np.array(range(0,13)) * 0.3048,11.28),
                                         )
    
def test_Mishev_spec_max_asymp_dir():

    array_of_lats_and_longs = np.array([[65.0,25.0],[-35.0,78.0]])
    array_of_zeniths_and_azimuths = np.array([(i,j) for i in np.linspace(0,20,5) for j in np.linspace(0,360,5)])

    outputted_dose = anisotropicMAIRE.run_from_power_law_gaussian_distribution(J0=256_000.0, gamma=3.41, deltaGamma=0.22, sigma=np.sqrt(0.19),
                                         reference_pitch_angle_latitude=-17.0,reference_pitch_angle_longitude=148.0,
                                         Kp_index=3, date_and_time=dt.datetime(2006, 12, 13, 3, 0),
                                         altitudes_in_km=np.append(np.array(range(0,13)) * 0.3048,11.28),
                                         array_of_lats_and_longs=array_of_lats_and_longs,
                                         array_of_zeniths_and_azimuths=array_of_zeniths_and_azimuths,
                                         )
    
    print(outputted_dose)


def test_DLR_spec():
    anisotropicMAIRE.run_from_DLR_cosmic_ray_model(OULU_count_rate_in_seconds=106.0,
                                         Kp_index=3, date_and_time=dt.datetime(2006, 12, 13, 3, 0),
                                         altitudes_in_km=np.append(np.array(range(0,13)) * 0.3048,11.28),
                                         )