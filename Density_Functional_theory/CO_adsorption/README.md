# This folder contains the DFT scripts to calculate CO adsorption energy on serveral different binding sites with different configurations and vibrational frequency of CO

### Objectives
The goal of this calculation is to get familiar with the DFT code for computing surface adsorption properties.
### Method
The adsorption energy of CO on Au can be calculated by the following equation
E<sub>ad</sub> = E<sub>Au*CO</sub> - E<sub>Au</sub> - E<sub>CO</sub>

A trilayer Au [111] slab was used as the surface with the bottom layer fixed. CO was allowed to approach the Au surface at three different sites, namely the bridge sight, the fcc site, and the hcp site. 

  <img src="https://user-images.githubusercontent.com/66216181/109054624-b5656500-76a3-11eb-9ae9-5388e5f35eb4.png" width="1000" height="300">

During the calculation only the positions of CO molecule and the top two layers of Au slab were allowed to relax.

### Resutls
After optimization,<br/>
* E<sub>CO</sub> = -14.80682643 eV<br/>
* E<sub>Au_3layer</sub> = -41.07678214 eV<br/>
* E<sub>Au_4layer</sub> = -56.74706518 ev

The calculated adsorption energies are summarized in the following table:


| Binding Site | Binding Atom | Binding Energy<sup>*</sup> |
| ------------ | ------------ | -------------- |
|    Bridge    |      C       |    -0.6903 eV  |
|     FCC      |      C       |    -0.6737 eV  |
|     HCP      |      C       |    -0.6199 eV  |
|    Bridge    |      O       |    -0.1117 eV  |
|     HCP      |      O       |    -0.1101 eV  |
|     Top      |      C       |    -0.5202 eV  |

\* Top configuration was calculated on a 4-layer Au slab, others were calcualted on a 3-layer Au slab.


Apparently, the bridge configuration with C approaching the Au surface is the most stable one. Using this configuration, the internal stretching frequency of CO was calculated to be 1890cm<sup>-1</sup>, which is quite different from experimental values. The internal stretching frequency for top configuration was calculated to be 2064 cm<sup>-1</sup>, also different from experimental values.

### Conclusion 
The adsorption energy and CO internal stretching frequency were calculated using PBE functional, the resutls are within acceptable ranges compared to DFT calculated values using PBE in the literature. However, PBE functional is known to be problematic for this type of calculation. A better functional to use is RPBE.

Reference: 
- The Journal of Physical Chemistry Letters, 2019, 10 (5), 1043-1047. DOI: 10.1021/acs.jpclett.9b00069
- Surface Science, 2014, 627, 57-69. DOI: 10.1016/j.susc.2014.04.012.
- Adsorption Science & Technology, 2013, 31 (9), 823-828. DOI: 10.1260/0263-6174.31.9.823
