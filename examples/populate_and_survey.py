#!/usr/bin/python

# import the populate and dosurvey modules
import populate
import dosurvey
import matplotlib.pyplot as plt

# run the populate.generate script, setting parameters as necessary.
# the resulting model is stored in the variable 'pop'
# Any unspecified variables use the default values (see populate
# for more details)

pop = populate.generate(0.1, 
               surveyList=['PMSURV'],
               radialDistType='lfl06',
               siDistPars=[-1.41, 0.96], # non-standard SI distribution
               duty_percent=6.,
               electronModel='lmt85',
               nostdout=True # switches off output to stdout
               )

# now run "dosurvey.run" on the model. Provide a list of surveys to use
surveyPopulations = dosurvey.run(pop, ['PMSURV'], nostdout=True)

# write out the survey results however you like
dosurvey.write(surveyPopulations, nores=False, summary=False, asc=False)

# write out the population model, if required
# pop.write(outf="populate.model")

# lists to store p/pdot
periods = [pulsar.period for pulsar in pop.population if not pulsar.dead]
pdots = [pulsar.pdot for pulsar in pop.population if not pulsar.dead]

# plot a scatter log-log plot of the p/pdot values
plt.loglog(periods, pdots, 'C0.')
plt.show()
# plt.savefig("test1.png")
