Notes on Simulations
--------------------

If you want to simulate survey yields note that you CANNOT do this:

populate # to generate a population
loop over niter
  dosurvey

You must do the following:

loop over niter
  populate
  dosurvey

The latter takes A LOT longer but the variations you see in dosurvey
iterations of the same population model do not represent the variation
in the population as sampled by the various distributions.

So if you want to (say) work out the spectral index parameters a la
Bates et al. (2013) you need to do the following:

loop over sp mean
  loop over sp sigma
    loop over niter
      populate
      dosurvey

The script bates.sh, to be found in this directory, does this for a
given input sp mean, i.e. it is designed so that you parallelise sp
mean over your CPUs. To run it you do:

bash bates.sh $sp_mean

or maybe something like:

loop over sp_mean
  schedtool -a $sp_mean bates.sh $sp_mean




