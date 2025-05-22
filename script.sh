#!/bin/bash

if [ $# -ne 2 ]; then

  echo "This is a super script to run PsrPopPy3 loads of times"
  echo ""
  echo "Usage: bash script.sh <num_PSRs> <num_iter>"
  exit

fi

# Numbers from each of the 3 Parkes surveys are:
# psrcatdog -c survey | grep pks70 | wc -l
# 298, 1158, 18 # to be investigated: why is it 1160 for PMPS on web version of psrcat?!

pmsurv_number=$1
niter=$2

#for sp_mean in `seq -f %.2f -2.00 0.01 -1.00`
#for sp_mean in `seq -f %.2f -2.00 0.05 -1.00`
for sp_mean in `seq -f %.2f -1.95 0.05 -1.95`
do 

#  for sp_sigma in `seq -f %.2f 0.00 0.01 2.00`
  for sp_sigma in `seq -f %.2f 1.35 0.05 2.00`
  do

#    mkdir results
    for iteration in `seq 0 $niter`
    do

      echo "Generating population with Spectral index mean" $sp_mean "and spectral sigma" $sp_sigma
      python3.10 psrpoppy/populate.py -n $pmsurv_number -surveys PMSURV -si $sp_mean $sp_sigma --nostdout
      python3.10 psrpoppy/dosurvey.py -surveys PKS70 PMSURV MMB --asc --summary --accel --nostdout

      # Rename files to mark them per iteration
      # .det files
      mv PKS70.det results/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".det"
      mv PMSURV.det results/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".det"
      mv MMB.det results/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".det"
      # .det files
      mv PKS70.results results/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".results"
      mv PMSURV.results results/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".results"
      mv MMB.results results/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".results"
      # .results files
      mv PKS70.summary results/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".summary"
      mv PMSURV.summary results/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".summary"
      mv MMB.summary results/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".summary"

      pks70=`grep "Detected" results/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".summary" | awk '{print $2}'`
      pmsurv=`grep "Detected" results/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".summary" | awk '{print $2}'`
      mmb=`grep "Detected" results/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".summary" | awk '{print $2}'`

      echo $iteration $pks70 $pmsurv $mmb 
    done
    grep Detected results/PKS70"_"$sp_mean"_"$sp_sigma"_"*.summary | awk '{print $NF}' | awk -v m=$sp_mean -v r=$sp_sigma '{s+=$1; ss+=$1*$1}END{print m, r, s/NR, sqrt(ss/NR - (s/NR)*(s/NR))}' >> PKS70_vals
    grep Detected results/PMSURV"_"$sp_mean"_"$sp_sigma"_"*.summary | awk '{print $NF}' | awk -v m=$sp_mean -v r=$sp_sigma '{s+=$1; ss+=$1*$1}END{print m, r, s/NR, sqrt(ss/NR - (s/NR)*(s/NR))}' >> PMSURV_vals
    grep Detected results/MMB"_"$sp_mean"_"$sp_sigma"_"*.summary | awk '{print $NF}' | awk -v m=$sp_mean -v r=$sp_sigma '{s+=$1; ss+=$1*$1}END{print m, r, s/NR, sqrt(ss/NR - (s/NR)*(s/NR))}' >> MMB_vals

  done
done

exit
