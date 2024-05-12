#!/usr/bin/env python

import sys
import argparse
import math
import random

import pickle

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, RadioButtons, CheckButtons

import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.signal import argrelmax

from psrpoppy.population import Population 
from psrpoppy.pulsar import Pulsar


class Visualize:
    """
    Class for plotting different parameters of a population file
    using the matplotlib module.

    """
    def __init__(self, popfile, frac):
        """Initialise the Visualize object."""

        # read in population file to population self.pop
        try:
            f = open(popfile, 'rb')
        except IOError:
            print(f"Could not open file {popfile}.")
            sys.exit()

        self.pop = pickle.load(f)
        f.close()

        # print out the population parameters 
        print(self.pop)

        self.textlabels = ['Period',
                      'DM',
                      'Gal X',
                      'Gal Y',
                      'Gal Z',
                      'W',
                      'Age',
                      'Pdot',
                      'SI',
                      'S1400',
                      'gl',
                      'gb',
                      'D',
                      'r0',
                      'B']
        self.axislabels = ['Period (ms)',
                      'DM (cm^-3 pc)',
                      'X (kpc)',
                      'Y (kpc)',
                      'Z (kpc)',
                      'Width (degrees)',
                      'Age [P/2Pdot] (years)',
                      'Pdot',
                      'Spectral Index',
                      'S1400 (mJy)',
                      'Galactic Longitude (degrees)',
                      'Galactic Latitude (degrees)',
                      'Distance (kpc)',
                      'GalacticRadius (kpc)',
                      'B [(PPdot)^(1/2)] (G)']

        if len(self.textlabels) != len(self.axislabels):
            print("Label list lengths not identical.")
            sys.exit()

        # for each of those labels, get the data from population into lists
        # I suppose most sense is to read it into arrays
        
        # create numpy array of right size
        dataArray = np.zeros((len(self.textlabels), self.pop.size()), float)
        # loop over pulsars and fill array
        npsr = 0

        # going to throw in a factor to reduce the number of plotted pulsars
        # to improve speed
        # make this an option
        if frac is None or frac > 1.0:
            frac = 1.0

        # print(self.pop.population[0].__dir__())

        for psr in self.pop.population:
            if random.random() < frac:
                dataArray[0][npsr] = psr.period
                dataArray[1][npsr] = psr.dm
                dataArray[2][npsr] = psr.galCoords[0]
                dataArray[3][npsr] = psr.galCoords[1]
                dataArray[4][npsr] = psr.galCoords[2]
                dataArray[5][npsr] = psr.width_degree
                dataArray[6][npsr] = psr.period / (2 * psr.pdot) / 31556952
                dataArray[7][npsr] = psr.pdot
                dataArray[8][npsr] = psr.spindex
                dataArray[9][npsr] = psr.s_1400()
                dataArray[10][npsr] = psr.gl
                dataArray[11][npsr] = psr.gb
                dataArray[12][npsr] = psr.dtrue 
                dataArray[13][npsr] = psr.r0
                dataArray[14][npsr] = np.sqrt(psr.period * psr.pdot) * 3.2e16
                npsr+=1

        self.dataArray = dataArray
        # delete population object, we don't need it now (I think!)
        # TBH not sure if this does anything, what with garbage collection
        del self.pop

    def display(self):
        """Method to create the plotting window and fill it with buttons."""

        # Create matplotlib window dressing
        self.fig = plt.figure(figsize=(10, 10), facecolor='lightgoldenrodyellow')
        # self.fig.canvas.set_window_title('PyPop: Visualization')
        self.fig.clear()

        # give titles
        self.fig.text(0.01, 0.97, 'PLOT SELECTION')
        self.fig.text(0.01, 0.80, 'X Axis')
        self.fig.text(0.01, 0.50, 'Y Axis')

        # arguments are [left, bottom, width, height] and label

        # radio buttons for x and y axis values
        radioXAxis = RadioButtons(plt.axes([0.070, 0.650, 0.150, 0.300]), self.textlabels,
                                    active=0, activecolor='blue')
        self.xindex = 0

        radioYAxis = RadioButtons(plt.axes([0.070, 0.350, 0.150, 0.300]), self.textlabels,
                                    active=1, activecolor='blue')
        self.yindex = 1

        # checkboxes to do log scale
        buttonXLog = CheckButtons(plt.axes([0.025, 0.225, 0.075, 0.050]), ['log x'],
                                    actives=[False])
        self.xlog = False

        buttonYLog = CheckButtons(plt.axes([0.125, 0.225, 0.075, 0.050]), ['log y'],
                                    actives=[False])
        self.ylog = False

        # button to do plot
        buttonPlot =      Button (plt.axes([0.075, 0.150, 0.075, 0.050]), 'Plot')

        # button to do histogram plot
        buttonHist =      Button (plt.axes([0.075, 0.050, 0.075, 0.050]), 'Hist')

        # if radio button is clicked, change active X/Y index
        radioXAxis.on_clicked(self._radioXClick)
        radioYAxis.on_clicked(self._radioYClick)

        # callbacks for the log plot switches
        buttonXLog.on_clicked(self._logClick)
        buttonYLog.on_clicked(self._logClick)

        # add callback to the plot button(s?)
        buttonPlot.on_clicked(self._scatterplot)
        buttonHist.on_clicked(self._histplot)

        plt.show()

    def _logClick(self, label):
        """Function to switch logarithm booleans"""
        if label == 'log x':
            self.xlog = not self.xlog
        elif label == 'log y':
            self.ylog = not self.ylog
        else:
            print("Weird log label!")
            sys.exit()

    def _radioXClick(self, label):
        self.xindex = self.textlabels.index(label)

    def _radioYClick(self, label):
        self.yindex = self.textlabels.index(label)

    def _scatterplot(self, event):
        # if there's already a scatter/histogram plot, delete it!
        try: self.fig.delaxes(self.histplot)
        except: pass

        try: self.fig.delaxes(self.scatterplt)
        except: pass

        # define axis position and dimensions
        self.scatterplt = self.fig.add_axes([0.32, 0.15, 0.65, 0.8])

        # axis labels
        self.scatterplt.set_xlabel(self.axislabels[self.xindex])
        self.scatterplt.set_ylabel(self.axislabels[self.yindex])

        # plot stuff
        self.scatterplt.plot(self.dataArray[self.xindex],
                             self.dataArray[self.yindex],
                             'r.')
        
        # if log switches are on, do a log plot
        if self.xlog:
            try:
                self.scatterplt.set_xscale('log')
            except ValueError:
                print("matplotlib refuses to do a log plot of the X data!")
                pass
        if self.ylog:
            try:
                self.scatterplt.set_yscale('log')
            except ValueError:
                print("matplotlib refuses to do a log plot of the Y data!")
                pass

        # finally, display the plot
        plt.show()

    def _histplot(self, event):
        # if there's already a scatter/histogram plot, delete it!
        try: self.fig.delaxes(self.scatterplt)
        except: pass

        try: self.fig.delaxes(self.histplot)
        except: pass

        # define axis position and dimensions
        self.histplot = self.fig.add_axes([0.32, 0.15, 0.65, 0.8])

        # axis labels
        self.histplot.set_xlabel(self.axislabels[self.xindex])

        # define data
        data = self.dataArray[self.xindex]

        if self.xindex in [2, 3, 10, 12]:
            # plotting two gaussians for better fit
            hist, bins_edges, _ = self.histplot.hist(data, bins=100, density=True, alpha=0.5, edgecolor="black")

            # Calculate bin centers
            bins_centers = (bins_edges[:-1] + bins_edges[1:]) / 2

            # Find peaks in histogram
            peak_indices = argrelmax(hist)[0]  # Find local maximums
            peak_values = hist[peak_indices]
            top_two_peaks = np.argsort(peak_values)[::-1][:2]  # Select the two highest peaks
            μ1, μ2 = bins_centers[peak_indices[top_two_peaks]]

            # Define the model function as a mixture of two Gaussians
            # with corresponding parameters (amplitude A, mean μ, standard deviation σ)
            def gaussian_mixture(x, A1, μ1, σ1, A2, μ2, σ2):
                return (A1 * norm.pdf(x, μ1, σ1) +
                        A2 * norm.pdf(x, μ2, σ2))

            # Initial guess for the parameters
            initial_guess = [1, μ1, 1, 1, μ2, 1]

            # Fit the model to the data [params (A1, μ1, σ1, A2, μ2, σ2), covariance matrix]
            params, cov = curve_fit(gaussian_mixture, bins_centers, hist, p0=initial_guess)

            # Plot the fitted curve [*params --> unpacking params]
            x = np.linspace(min(data), max(data), 1000)
            plt.plot(x, gaussian_mixture(x, *params), "r")
            plt.title("A1={0:.2f}, μ1={1:.2f}, σ1={2:.2f}, A2={3:.2f}, μ2={4:.2f}, σ2={5:.2f}".format(*params))

        else:
            # plot histogram along with gaussian
            self.histplot.hist(data, bins=100, density=True, alpha=0.5, edgecolor="black")

            μ, σ = norm.fit(data)
            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, 1000)
            p = norm.pdf(x, μ, σ)
            plt.plot(x, p, "k")
            plt.title(f"μ={μ:.2f}, σ={σ:.2f}")

        plt.show()

if __name__ == '__main__':
    """ 'Main' function for calling from command line"""
    parser = argparse.ArgumentParser(description='Visualize a population object')
    parser.add_argument('-f', metavar='filename', required=True,
                          help='file containing population model')

    parser.add_argument('-frac', metavar='fraction', type=float, default=None,
                          help='plot only this fraction of pulsars (default=None)')

    args = parser.parse_args()

    v = Visualize(popfile=args.f, frac=args.frac)

    v.display()
