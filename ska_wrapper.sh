#!/bin/bash

#
# Run the SKA survey yields niter times
#
# 2025/05/26 EFK
#
if [ $# -ne 1 ]; then

  echo "This is a wrapper for running SKA survey scripts"
  echo ""
  echo "Usage: bash ska_wrapper.sh <num_iter>"
  exit

fi

niter=$1

for iteration in `seq 0 $niter`                    # Loop over niter iterations
do
  echo "Iteration "$iteration
  bash ska.sh $iteration
done

#AA4
echo "AA4 Low"
grep Detected SKA_low_allsky_AA4_2025_msp*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "MSPs: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'
grep Detected SKA_low_allsky_AA4_2025_slow*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "Slowies: "s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'

echo "AA4 Mid Band 1"
grep Detected SKA_mid_allsky_band1_AA4_2025_msp*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "MSPs: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'
grep Detected SKA_mid_allsky_band1_AA4_2025_slow*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "Slowies: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'

echo "AA4 Mid Band 2"
grep Detected SKA_mid_allsky_band2_AA4_2025_msp*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "MSPs: "s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'
grep Detected SKA_mid_allsky_band2_AA4_2025_slow*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "Slowies: "s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'

echo ""
#AAstar
echo "AAstar Low"
grep Detected SKA_low_allsky_AAstar_2025_msp*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "MSPs: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'
grep Detected SKA_low_allsky_AAstar_2025_slow*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "Slowies: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'

echo "AAstar Mid Band 1"
grep Detected SKA_mid_allsky_band1_AAstar_2025_msp*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "MSPs: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'
grep Detected SKA_mid_allsky_band1_AAstar_2025_slow*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print"Slowies: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'

echo "AAstar Mid Band 2"
grep Detected SKA_mid_allsky_band2_AAstar_2025_msp*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "MSPs: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'
grep Detected SKA_mid_allsky_band2_AAstar_2025_slow*/*summary | awk '{print $NF}' | awk '{s+=$1; ss+=$1*$1}END{print "Slowies: ",s/NR," +/- ",sqrt(ss/NR - (s/NR)*(s/NR))}'

exit
