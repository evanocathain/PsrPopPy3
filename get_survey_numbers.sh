#!/bin/bash

psrcat='/Users/e.keane/PSRCAT/psrcat_tar/psrcat -db_file /Users/e.keane/PSRCAT/psrcat_tar/psrcat.db'
psrcatdog=$psrcat' -o short -nohead -nonumber'

# Numbers from each of the 3 Parkes surveys are:
# 
# First up get the numbers from the time of the Bates et al. (2013) work
# Number with P>50ms
#279, 1038, 18

# Next up get the most up to date numbers at time of writing
# psrcat v2.6.1 released 8/5/2025
$psrcatdog -v
pks70=`$psrcatdog -c survey | grep pks70 | awk '{if ($2>0.050) print $0}' | wc -l`
pksmb=`$psrcatdog -c survey | grep pksmb | awk '{if ($2>0.050) print $0}' | wc -l`
pksmmb=18 # hard-coded
# No period restriction: 298, 1158, 18
# P>50ms restriction:    279, 1115, 18

psrcatdog -c "name p0 survey" | grep pks70 | 
echo $pks70 $pksmb $pksmmb

exit

