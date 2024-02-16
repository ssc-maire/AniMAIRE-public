# AniMAIRE


**N.B. Currently this tool only runs on Linux-based machines (not on Windows)**

A Python toolkit for calculating dose rates in Earth's atmosphere based on any incoming proton or alpha particle spectra, with any pitch angle distribution. 

**If you use this software for scientific research, please reference AniMAIRE as C. S. W. Davis, F. Lei, F. Baird, K. Ryden and C. Dyer (2023). AniMAIRE version {version number}. https://github.com/ssc-maire/AniMAIRE-public , https://pypi.org/project/AniMAIRE/ . Surrey Space Center, University of Surrey.**

## Installation

To install this toolkit from this Github repository, first clone this repository to your local system, and then from the cloned respository, run

```
sudo pip3 install .
```

in the cloned directory. 

Note that there are quite a few sizeable data files within some of the dependencies for this package that get copied during installation (on the order of about several hundred megabytes in total) so installation may take a couple of minutes.

## Usage

After installation, to import the toolkit into a particular Python script, run 
```
from AniMAIRE import AniMAIRE
```
All of the main useful functions are contained within this `AniMAIRE` module, and all other modules contained in this toolkit are primarily intended to be accessed internally (although don't let that stop you from using or editing them for your own purposes if you wish).

## Calculating dose rates at any location in Earth's atmosphere

The basic function for performing a run to calculate dose rates in `AniMAIRE` is the `run_from_spectra` function, which has the format:

```
def run_from_spectra(
        proton_rigidity_spectrum=None,
        alpha_rigidity_spectrum=None,
        reference_pitch_angle_latitude=None,
        reference_pitch_angle_longitude=None,
        proton_pitch_angle_distribution=isotropicPitchAngleDistribution(),
        alpha_pitch_angle_distribution=isotropicPitchAngleDistribution(),
        altitudes_in_kft=[0,10,20] + [i for i in range(25, 61 + 1, 3)],
        altitudes_in_km=None,
        Kp_index=None,
        date_and_time=dt.datetime.now(),
        array_of_lats_and_longs=default_array_of_lats_and_longs,
        cache_magnetocosmics_run=True,
)
```

`run_from_spectra` performs a run at a single date and time and Kp index to calculate dose rates across Earth's atmosphere based on proton, alpha particle, or proton + alpha particle spectra. **Particle spectra here must be described in units of cm-2 s-1 sr-1 (GV/n)-1, and with respect to rigidity in units of GV**.

Particle spectra, as well as pitch angle distributions, can be set as any 'callable' object in Python, i.e. a function, as will be shown in examples below. At least one particle spectrum must be specified, as well as a Kp index, so this function to execute successfully. For runs designed to simulate dose rates during particular dates and times, the argument `date_and_time` must also be supplied with a Python `datetime` corresponding to the timestamp being investigated (by default, the function assumes that the current date and time should be used).

Note that while this function takes an alpha particle spectrum as an input, it actually interpolates the dose rates due to alpha particles to those of heavier ions too, so outputted dose rates due to an alpha particle spectrum are in fact the combined total of all ions heavier than protons.

`AniMAIRE` performs runs of the MAGNETOCOSMICS as part of dose rate calculations (using the [AsympDirsCalculator](https://github.com/ssc-maire/AsymptoticDirectionsCalculator) package), and these currently take up by far the majority of `AniMAIRE` runtime - on the order of over half an hour on the developer's computer versus less than 6 minutes for the rest of the program. Therefore if the `cache_magnetocosmics_run` argument is set to `True`, which it is by default, `AniMAIRE` will cache the results of MAGNETOCOSMICS simulations in the directory that `AniMAIRE` is run from in the generated `cachedMagnetocosmicsRunData` and `cacheAsymptoticDirectionOutputs` directories. This significantly speeds up any tasks where users wish to investigate a constant `Kp_index` and `date_and_time`, but wish to vary the spectrum and pitch angle distribution and investigate how dose rates are impacted, as magnetocosmics runs are only performed once.

### Simple isotropic runs and plotting

A basic run of the `run_from_spectra` function might look like this:

```
from AniMAIRE import AniMAIRE
import datetime as dt

test_isotropic_dose_rates = AniMAIRE.run_from_spectra(
        proton_rigidity_spectrum=lambda x:2.56*(x**-3.41),
        Kp_index=3,
        date_and_time=dt.datetime(2006, 12, 13, 3, 0),
)
```

in this example, the proton rigidity spectrum is set to be a power law with a normalisation factor of 2.56 cm-2 s-1 sr-1 (GV/n)-1, and a spectral index of 3.41, using the commonly used `lambda` approach to create a function within a single line. Kp index is set to be 3, and the date and time to simulate are set to be 13th of December 2006, 03:00. This function will likely take at least several minutes to run, depending on the speed of the machine and number of cores, and should output a Pandas DataFrame to `test_isotropic_dose_rates`, giving:

```
	latitude	longitude	altitude (km)	edose	adose	dosee	tn1	tn2	tn3	SEU	SEL
0	-90.0	0.0	0.0000	0.010442	0.012540	0.010010	0.004437	0.002729	0.001828	2.729229e-16	2.729229e-11
1	-90.0	0.0	3.0480	0.101786	0.117658	0.085755	0.051895	0.033617	0.022979	3.361684e-15	3.361684e-10
2	-90.0	0.0	6.0960	0.672702	0.742332	0.457695	0.326731	0.211853	0.145046	2.118530e-14	2.118530e-09
3	-90.0	0.0	7.6200	1.442377	1.541670	0.975436	0.665785	0.431261	0.295516	4.312608e-14	4.312608e-09
4	-90.0	0.0	8.5344	2.165860	2.249419	1.426324	0.964791	0.623291	0.426927	6.232913e-14	6.232913e-09
...	...	...	...	...	...	...	...	...	...	...	...
42619	90.0	355.0	14.9352	16.953199	14.033322	9.762795	5.138425	3.234975	2.164808	3.234975e-13	3.234975e-08
42620	90.0	355.0	15.8496	21.327714	16.796347	11.858681	5.949880	3.711156	2.463976	3.711156e-13	3.711156e-08
42621	90.0	355.0	16.7640	26.122787	19.802434	14.318360	6.759204	4.182604	2.757482	4.182604e-13	4.182604e-08
42622	90.0	355.0	17.6784	31.856083	23.278901	17.431863	7.562120	4.638132	3.029771	4.638132e-13	4.638132e-08
42623	90.0	355.0	18.5928	38.265701	26.680896	20.232603	8.340509	5.067928	3.275415	5.067928e-13	5.067928e-08
```

when `test_isotropic_dose_rates` is printed.

The outputted dose rate (or flux) labels represent the following dose rate/flux types:

|label | dose rate/flux type|
|------|--------------------|
|adose| ambient dose equivalent in µSv/hr |
|edose| effective dose in µSv/hr |
|dosee| dose equivalent in µSv/hr |
|tn1| >1 MeV neutron flux, in n/cm2/s |
|tn2| >10 MeV neutron flux, in n/cm2/s |
|tn3| >60 MeV neutron flux, in n/cm2/s |
|SEU| single event upset rate for an SRAM device in upsets/second/bit |
|SEL| single event latch-up rate for an SRAM device in latch-ups/second/device |

These dose rates are produced by default at every latitude and longitude corresponding to 5 by 5 degree intervals across Earth's surface, and for altitudes of 0 kilofeet, 10 kilofeet, 20 kilofeet and between 25 kilofeet and 61 kilofeet at intervals of 3 kilofeet. This can be altered using the `altitudes_in_kft`, `altitudes_in_km` and `array_of_lats_and_longs`. 

Any particular altitudes the user wants to use can be supplied to `altitudes_in_kft` or `altitudes_in_km` as a `list` or numpy array. 

If you want to perform calculations only at a specific set of latitudes and longitudes you should use the `array_of_lats_and_longs` argument, supplying it as a 2 dimensional `list` or numpy array, where the first column refers to latitudes and the second column refers to longitudes. All longitudes in this case should be specified in terms of longitude east (i.e. 0.00 degrees - 359.99 degrees). **Using the `array_of_lats_and_longs` argument significantly speeds up the running of `AniMAIRE` if you're only interested in a small number of coordinates, so its use is highly recommended in those situations.**

There are many ways you could plot this data. An example function, `create_single_dose_map_plot`, has been supplied in `AniMAIRE` that uses plotly to plot the dose rates across Earth (i.e. as a function of latitude and longitude) at a given altitude. Its specification is the following:

```
def create_single_dose_map_plot(DF_to_use,
                                selected_altitude_in_km)
```
where `DF_to_use` is the Pandas DataFrame outputted by a run of `AniMAIRE` and altitude is one of the altitudes in kilometers supplied to/outputted by the run.

To use this function to create a map of the isotropic situation as given as an example above, you could run 
```
isotropic_dose_rate_map = AniMAIRE.create_single_dose_map_plot(test_isotropic_dose_rates,
                            			      selected_altitude_in_km = 12.1920)
```

which should plot the following figure as a plotly plot:

![isotropic_test_plot](https://user-images.githubusercontent.com/16866485/223750107-402e56f6-bea1-48b4-83e9-2179be7f7449.png)

and assign the plot to the `isotropic_dose_rate_map` variable for the user to use as they wish.

### Anisotropic runs

`run_from_spectra` defaults to an isotropic spectrum if no pitch angle distribution is supplied for either protons or alpha particles by the user.  

To run an anisotropic spectrum, differential pitch angle distributions must be supplied to the `proton_pitch_angle_distribution` and/or `alpha_pitch_angle_distribution` arguments in `run_from_spectra` along with a reference location specified in terms of latitude and longitude in the `reference_pitch_angle_latitude` and `reference_pitch_angle_longitude` arguments respectively. The pitch angle distributions must be supplied as 2 dimensional functions, where the first argument is the pitch angle, and the second argument is particle rigidity (in many cases the pitch angle distribution might not depend on rigidity, but for programmatic reasons the function must at least take in rigidity as an argument although it does not need to have a dependence on it). The pitch angle here must be specified in units of **radians**.

`reference_pitch_angle_latitude` and `reference_pitch_angle_longitude` are the reference latitude and longitude in GEO coordinates representing a pitch angle of 0 in the supplied pitch angle distribution used. `AniMAIRE` currently makes the assumption that incoming particle distributions are cylindrically symmetric about this reference direction, and therefore that only the pitch angle with respect to this latitude and longitude are required to calculate dose rates anisotropically across Earth. Incoming solar particle events are frequently oriented near to the direction of the Interplanetary Magnetic Field (IMF), so you could specify pitch angles relative to the IMF here, and use the latitude and longitude of the IMF as the reference latitude and longitude.

The pitch angle distributions and rigidity spectrum must be specified in units normalised such that the product of the pitch angle distribution and rigidity spectrum multiplied together is in units of **cm-2 s-1 sr-1 (GV/n)-1**. 

The pitch angle distributions can be specified using the Python `lambda` as with the rigidity spectra, but with 2 dimensions rather than 1. For example, to specify a Gaussian pitch angle distribution you could use:

```
import numpy as np

sigma = np.sqrt(0.19)

test_pitch_angle_dist_function = lambda pitch_angle,rigidity:np.exp(-(pitch_angle**2)/(sigma**2))
```
where `sigma` here has arbitrarily been chosen to be the square root of 0.19 for example purposes. 

here `pitch_angle_dist_function` would be a viable input as a pitch angle distribution to `run_from_spectra`. While the function itself does not depend on `rigidity`, `rigidity` is specified as the second argument of the function nonetheless, as required for calculations to work.

An example of using this might be:

```
import numpy as np
from AniMAIRE import AniMAIRE
import datetime as dt

sigma = np.sqrt(0.19)
pitch_angle_reference_latitude = -17.0
pitch_angle_reference_longitude = 148.0

test_pitch_angle_dist_function = lambda pitch_angle,rigidity:np.exp(-(pitch_angle**2)/(sigma**2))

test_anisotropic_dose_rates = AniMAIRE.run_from_spectra(
        proton_rigidity_spectrum=lambda x:2.56*(x**-3.41),
        proton_pitch_angle_distribution=test_pitch_angle_dist_function,
        reference_pitch_angle_latitude=pitch_angle_reference_latitude,reference_pitch_angle_longitude=pitch_angle_reference_longitude,
        Kp_index=3,
        date_and_time=dt.datetime(2006, 12, 13, 3, 0),
)
```

when run this should output a Pandas DataFrame to `test_anisotropic_dose_rates` with the same general output format as given in the previously described isotropic dose rates case. 

In this case printing `test_anisotropic_dose_rates` should output:
```
	latitude	longitude	altitude (km)	edose	adose	dosee	tn1	tn2	tn3	SEU	SEL
0	-90.0	0.0	0.0000	1.976895e-07	2.027961e-07	3.222401e-07	1.286575e-08	7.877480e-09	5.401425e-09	7.877480e-22	7.877480e-17
1	-90.0	0.0	3.0480	4.838923e-07	4.869626e-07	7.032983e-07	7.429541e-08	4.866265e-08	3.410940e-08	4.866265e-21	4.866265e-16
2	-90.0	0.0	6.0960	1.453834e-06	1.411544e-06	2.045308e-06	2.472576e-07	1.611648e-07	1.136527e-07	1.611648e-20	1.611648e-15
3	-90.0	0.0	7.6200	2.352658e-06	2.382057e-06	3.304989e-06	3.720997e-07	2.420436e-07	1.708843e-07	2.420436e-20	2.420436e-15
4	-90.0	0.0	8.5344	2.970462e-06	2.791970e-06	4.201789e-06	4.473736e-07	2.919136e-07	2.063902e-07	2.919136e-20	2.919136e-15
...	...	...	...	...	...	...	...	...	...	...	...
42619	90.0	355.0	14.9352	1.062397e-06	7.975138e-07	5.902218e-07	2.766678e-07	1.764588e-07	1.216189e-07	1.764588e-20	1.764588e-15
42620	90.0	355.0	15.8496	1.289216e-06	9.183397e-07	6.964574e-07	3.063998e-07	1.945207e-07	1.338582e-07	1.945207e-20	1.945207e-15
42621	90.0	355.0	16.7640	1.525359e-06	1.040969e-06	8.062018e-07	3.339381e-07	2.112220e-07	1.452375e-07	2.112220e-20	2.112220e-15
42622	90.0	355.0	17.6784	1.777447e-06	1.163693e-06	9.297615e-07	3.568426e-07	2.248452e-07	1.543556e-07	2.248452e-20	2.248452e-15
42623	90.0	355.0	18.5928	2.036052e-06	1.275021e-06	1.040287e-06	3.755276e-07	2.357292e-07	1.614028e-07	2.357292e-20	2.357292e-15
```

which will produce the following plot when
```
anisotropic_dose_rate_map = AniMAIRE.create_single_dose_map_plot(test_anisotropic_dose_rates,
                            			      selected_altitude_in_km = 12.1920)
```
is run:

![anisotropic_test_plot](https://user-images.githubusercontent.com/16866485/223751057-5d0cff98-cf9e-4654-b71f-d1ae55c75602.png)

### Functions for running `AniMAIRE` for specific situations and for a past timestamp

In addition to the quite general `run_from_spectra` function, `AniMAIRE` currently contains several functions for running calculations for specific types of spectra and situations, to make it easier for users to perform runs without having to determine and feed in spectra to `run_from_spectra` themselves.

The `run_from_DLR_cosmic_ray_power_law` function allows users to run full atmospheric dose rate calculations (for cosmic ray only/'quiet' time periods only) from just a date and time, or alternatively from just a single OULU count rate, or just a value of 'W parameter', as well as Kp index. This function utilises the [CosRayModifiedISO package](https://github.com/ssc-maire/CosRayModifiedISO) to determine the spectra due to protons and alpha particles during cosmic ray only time periods, and then runs `run_from_spectra` using both of those spectra under isotropic conditions. The specifications of `run_from_DLR_cosmic_ray_power_law` are:

```
def run_from_DLR_cosmic_ray_model(OULU_count_rate_in_seconds=None,
                                      W_parameter=None,
                                      Kp_index=None,
				      date_and_time=None,
                                      **kwargs)
```

`**kwargs` here can be used to supply any arguments you wish to `run_from_spectra` as specified previously, such as the list of altitudes and list of coordinates to perform calculations for. Details on what `OULU_count_rate_in_seconds` and `W_parameter` mean can be found at https://github.com/ssc-maire/CosRayModifiedISO . If either `OULU_count_rate_in_seconds` or `W_parameter` are used, only one of them should be specified. Otherwise, `AniMAIRE` will determined their values using the `date_and_time` parameter supplied.

In addition to running from the DLR-ISO isotropic cosmic ray model, you can also run `AniMAIRE` from a combined power law rigidity spectrum and Gaussian pitch angle distribution. This can be done using the `run_from_power_law_gaussian_distribution` function:

```
def run_from_power_law_gaussian_distribution(J0, gamma, deltaGamma, sigma, 
                                             reference_pitch_angle_latitude, reference_pitch_angle_longitude, 
                                             Kp_index,date_and_time,
                                             **kwargs)
```

Here `J0`, `gamma`, `deltaGamma`, `sigma`, `reference_pitch_angle_latitude`, `reference_pitch_angle_longitude` are all defined as specified in the format of papers like [Mishev, A., Usoskin, I. Analysis of the Ground-Level Enhancements on 14 July 2000 and 13 December 2006 Using Neutron Monitor Data. Sol Phys 291, 1225–1239 (2016). https://doi.org/10.1007/s11207-016-0877-2](https://link.springer.com/article/10.1007/s11207-016-0877-2).






