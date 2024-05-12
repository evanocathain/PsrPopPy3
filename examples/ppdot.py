#!/usr/bin/python

"""
**********************************************
python3 PsrPopPy/examples/ppdot.py <file_name>
**********************************************
1 argument --> file_name --> *.model file
"""

import sys
import pickle

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# open file, read in model
filename = sys.argv[1]
f = open(filename, 'rb')
pop = pickle.load(f)
f.close()

# lists to store p/pdot
periods = [pulsar.period for pulsar in pop.population if not pulsar.dead]
pdots = [pulsar.pdot for pulsar in pop.population if not pulsar.dead]

# plot a scatter log-log plot of the p/pdot values
plt.loglog(periods, pdots, 'k.')
plt.plot("test8.png")
plt.show()
