#!/bin/bash

if [ $# -ne 2 ]; then

  echo "This is a super script to run PsrPopPy3 loads of times"
  echo "It is designed to figure out the yield for 3 Parkes"
  echo "surveys for the provided spectral index mean vallue."
  echo ""
  echo "Usage: bash script.sh <sp_mean> <num_iter>"
  exit

fi

# Numbers from each of the 3 Parkes surveys are:
# PSRCAT v2.6.1
#
# psrcatdog -c survey | grep pksmb | wc -l
# psrcatdog -c survey | grep pks70 | wc -l
# NOTE: MMB is not identified separately in PSRCAT
#
# 298, 1158, 18 # to be investigated: why is it 1160 for PMPS on web version of psrcat?!
# 
# For this we focus on the fraction that are slowies, so we take P>50ms
#
# psrcatdog -c "p0 survey" | grep pksmb | awk '{if ($1>0.050) print $0}' | wc -l
# psrcatdog -c "p0 survey" | grep pks70 | awk '{if ($1>0.050) print $0}' | wc -l
# NOTE: All MMB pulsars are slowies
#
# 279, 1115, 18 are the relevant target numbers for slowies for these surveys

pmsurv_number=1115
sp_mean=$1
niter=$2

for sp_sigma in `seq -f %.2f 0.00 0.05 2.00`           # Loop over sp sigma parameter
do
  results_dir="results_"$sp_mean"_"$sp_sigma
  mkdir $results_dir
  for iteration in `seq 0 $niter`                    # Loop over niter iterations
  do
    echo "Generating population with Spectral index mean" $sp_mean "and spectral sigma" $sp_sigma
    python3.10 psrpoppy/populate.py -n $pmsurv_number -surveys PMSURV -si $sp_mean $sp_sigma --nostdout
    python3.10 psrpoppy/dosurvey.py -surveys PKS70 PMSURV MMB --asc --summary --accel --nostdout

    # Rename and move files to mark them per iteration
    # .det files
    mv PKS70.det  $results_dir/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".det"
    mv PMSURV.det $results_dir/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".det"
    mv MMB.det    $results_dir/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".det"
    # .results files
    mv PKS70.results  $results_dir/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".results"
    mv PMSURV.results $results_dir/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".results"
    mv MMB.results    $results_dir/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".results"
    # .summary files
    mv PKS70.summary  $results_dir/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".summary"
    mv PMSURV.summary $results_dir/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".summary"
    mv MMB.summary    $results_dir/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".summary"

    pks70=`grep "Detected" $results_dir/PKS70"_"$sp_mean"_"$sp_sigma"_"$iteration".summary" | awk '{print $2}'`
    pmsurv=`grep "Detected" $results_dir/PMSURV"_"$sp_mean"_"$sp_sigma"_"$iteration".summary" | awk '{print $2}'`
    mmb=`grep "Detected" $results_dir/MMB"_"$sp_mean"_"$sp_sigma"_"$iteration".summary" | awk '{print $2}'`

    echo $iteration $pks70 $pmsurv $mmb 
  done
  # Now that all iterations are done for this sp_mean and sp_sigma values, 
  # work out mean and RMS of detected for each of the survey
  grep Detected $results_dir/PKS70"_"$sp_mean"_"$sp_sigma"_"*.summary | awk '{print $NF}' | awk -v m=$sp_mean -v r=$sp_sigma '{s+=$1; ss+=$1*$1}END{print m, r, s/NR, sqrt(ss/NR - (s/NR)*(s/NR))}' >> $results_dir/"PKS70_"$sp_mean"_"$sp_sigma
  grep Detected $results_dir/PMSURV"_"$sp_mean"_"$sp_sigma"_"*.summary | awk '{print $NF}' | awk -v m=$sp_mean -v r=$sp_sigma '{s+=$1; ss+=$1*$1}END{print m, r, s/NR, sqrt(ss/NR - (s/NR)*(s/NR))}' >> $results_dir/"PMSURV_"$sp_mean"_"$sp_sigma
  grep Detected $results_dir/MMB"_"$sp_mean"_"$sp_sigma"_"*.summary | awk '{print $NF}' | awk -v m=$sp_mean -v r=$sp_sigma '{s+=$1; ss+=$1*$1}END{print m, r, s/NR, sqrt(ss/NR - (s/NR)*(s/NR))}' >> $results_dir/"MMB_"$sp_mean"_"$sp_sigma

done

exit
