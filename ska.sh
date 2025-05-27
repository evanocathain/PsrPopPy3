#!/bin/bash

#
# Once you know the spectral index distribution (!)
# This does the SKA surveys and estimates the 
# yields of _detected_ pulsars
#
# 2025/05/26 EFK
#

label=$1
si_mean=-1.45
si_sigma=0.15

#
# Make Populations
#
# MSPs
#python3.10 psrpoppy/populate.py -ldist pow -l 0.1 1000 -1.45 -n 30000 -pdist lorimer12 -z 0.5 -w 6 -o msp.model
echo "Simulating MSP Population"
python3.10 psrpoppy/populate.py -n 30000 -pdist lorimer12 -z 0.5 -w 6 -o msp.model -si $si_mean $si_sigma --nostdout 
# Slowies
pmps_slowies=1115 # psrcat v2.6.1
pmps_msps=43      # psrcat v2.6.1
echo "Simulating Slowie Population"
python3.10 psrpoppy/populate.py -n $pmps_slowies -w 6 -surveys PMSURV -o slow.model -si $si_mean $si_sigma --nostdout
 
#
# Do the surveys 
#
# Do AA4 and AA* separately

# AA4
# MSPs
echo "AA4: Searching for MSPs"
python3.10 psrpoppy/dosurvey.py -f msp.model -surveys SKA_low_allsky_AA4_2025 SKA_mid_allsky_band1_AA4_2025 SKA_mid_allsky_band2_AA4_2025 --asc --summary --nostdout
# mv the .summary, .det and .results
mkdir SKA_low_allsky_AA4_2025_msp_it$label/
mkdir SKA_mid_allsky_band1_AA4_2025_msp_it$label/
mkdir SKA_mid_allsky_band2_AA4_2025_msp_it$label/
mv SKA_low_allsky_AA4_2025.* SKA_low_allsky_AA4_2025_msp_it$label/
mv SKA_mid_allsky_band1_AA4_2025.* SKA_mid_allsky_band1_AA4_2025_msp_it$label/
mv SKA_mid_allsky_band2_AA4_2025.* SKA_mid_allsky_band2_AA4_2025_msp_it$label/

# Slowies
echo "AA4: Searching for Slowies"
python3.10 psrpoppy/dosurvey.py -f slow.model -surveys SKA_low_allsky_AA4_2025 SKA_mid_allsky_band1_AA4_2025 SKA_mid_allsky_band2_AA4_2025 --asc --summary --nostdout
# mv the .summary, .det and .results
mkdir SKA_low_allsky_AA4_2025_slow_it$label/
mkdir SKA_mid_allsky_band1_AA4_2025_slow_it$label/
mkdir SKA_mid_allsky_band2_AA4_2025_slow_it$label/
mv SKA_low_allsky_AA4_2025.* SKA_low_allsky_AA4_2025_slow_it$label/
mv SKA_mid_allsky_band1_AA4_2025.* SKA_mid_allsky_band1_AA4_2025_slow_it$label/
mv SKA_mid_allsky_band2_AA4_2025.* SKA_mid_allsky_band2_AA4_2025_slow_it$label/

# AA*
# MSPs
echo "AAstar: Searching for MSPs"
python3.10 psrpoppy/dosurvey.py -f msp.model -surveys SKA_low_allsky_AAstar_2025  SKA_mid_allsky_band1_AAstar_2025 SKA_mid_allsky_band2_AAstar_2025 --asc --summary --nostdout
# mv the .summary, .det and .results
mkdir SKA_low_allsky_AAstar_2025_msp_it$label/
mkdir SKA_mid_allsky_band1_AAstar_2025_msp_it$label/
mkdir SKA_mid_allsky_band2_AAstar_2025_msp_it$label/
mv SKA_low_allsky_AAstar_2025.* SKA_low_allsky_AAstar_2025_msp_it$label/
mv SKA_mid_allsky_band1_AAstar_2025.* SKA_mid_allsky_band1_AAstar_2025_msp_it$label/
mv SKA_mid_allsky_band2_AAstar_2025.* SKA_mid_allsky_band2_AAstar_2025_msp_it$label/

# Slowies
echo "AAstar: Searching for Slowies"
python3.10 psrpoppy/dosurvey.py -f slow.model -surveys SKA_low_allsky_AAstar_2025 SKA_mid_allsky_band1_AAstar_2025 SKA_mid_allsky_band2_AAstar_2025 --asc --summary --nostdout
# mv the .summary, .det and .results
mkdir SKA_low_allsky_AAstar_2025_slow_it$label/
mkdir SKA_mid_allsky_band1_AAstar_2025_slow_it$label/
mkdir SKA_mid_allsky_band2_AAstar_2025_slow_it$label/
mv SKA_low_allsky_AAstar_2025.* SKA_low_allsky_AAstar_2025_slow_it$label/
mv SKA_mid_allsky_band1_AAstar_2025.* SKA_mid_allsky_band1_AAstar_2025_slow_it$label/
mv SKA_mid_allsky_band2_AAstar_2025.* SKA_mid_allsky_band2_AAstar_2025_slow_it$label/

exit
