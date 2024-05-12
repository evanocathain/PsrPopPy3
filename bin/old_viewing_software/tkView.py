#!/usr/bin/env python

import sys
import os
import argparse
import random
import warnings

import pickle
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.signal import argrelmax

import tkinter as tk
import matplotlib
matplotlib.use('TKAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import \
    FigureCanvasTkAgg, \
    NavigationToolbar2Tk as NavigationToolbar

class ViewException(Exception):
    pass

class DataObj:
    def __init__(self, name, textlabels, axislabels, dataarray, size):
        self.name = name
        self.textlabels = textlabels
        self.axislabels = axislabels
        self.dataArray = dataarray
        self.count = size

def makeDataObj(frac=None, extnList=['.model', '.results']):
    """
    Function to get different parameters of a population file

    """
    textlabels = ['Period',
                    'DM',
                    'Gal X',
                    'Gal Y',
                    'Gal Z',
                    'W',
                    'Pdot',
                    'Age',
                    'B',
                    'SI',
                    'S1400',
                    'gl',
                    'gb',
                    'D',
                    'r0',
                    'n']
    axislabels = ['Period (ms)',
                    'DM (cm^-3 pc)',
                    'X (kpc)',
                    'Y (kpc)',
                    'Z (kpc)',
                    'Width (°)',
                    'Pdot',
                    'Age [P/2Pdot] (yr)',
                    'B [(PPdot)^(1/2)] (G)',
                    'Spectral Index',
                    'S1400 (mJy)',
                    'Galactic Longitude (°)',
                    'Galactic Latitude (°)',
                    'Distance (kpc)',
                    'GalacticRadius (kpc)',
                    'Array Index']

    if len(textlabels) != len(axislabels):
        print("Label list lengths not identical.")
        sys.exit()

    dataObjList = []
    
    for filename in os.listdir("/PsrPopPy/"):
        for extn in extnList:
            if filename.endswith(extn):
                # read in population file to population self.pop
                # filename = "/PsrPopPy/populate.model"
                try:
                    f = open(filename, 'rb')
                except IOError:
                    print("Could not open file {0}.".format(filename))
                    sys.exit()

                pop = pickle.load(f)
                f.close()

                # for each of those labels, get the data from population into lists
                # I suppose most sense is to read it into arrays
                
                # create numpy array of right size
                dataArray = np.zeros((len(textlabels), pop.size()), float)
                # loop over pulsars and fill array
                npsr = 0

                # going to throw in a factor to reduce the number of plotted pulsars
                # to improve speed
                # make this an option
                if frac is None or frac > 1.0:
                    frac = 1.0

                for psr in pop.population:
                    if random.random() < frac:
                        dataArray[0][npsr] = psr.period
                        dataArray[1][npsr] = psr.dm
                        dataArray[2][npsr] = psr.galCoords[0]
                        dataArray[3][npsr] = psr.galCoords[1]
                        dataArray[4][npsr] = psr.galCoords[2]
                        dataArray[5][npsr] = psr.width_degree
                        dataArray[6][npsr] = psr.pdot
                        dataArray[7][npsr] = (psr.period / (2 * psr.pdot) / 31556952 ) if psr.pdot else 0
                        dataArray[8][npsr] = (np.sqrt(psr.period * psr.pdot) * 3.2e16) if psr.pdot else 0
                        dataArray[9][npsr] = psr.spindex
                        dataArray[10][npsr] = psr.s_1400()
                        dataArray[11][npsr] = psr.gl
                        dataArray[12][npsr] = psr.gb
                        dataArray[13][npsr] = psr.dtrue 
                        dataArray[14][npsr] = psr.r0
                        dataArray[15][npsr] = npsr
                        npsr+=1

                dataObjList.append(DataObj(filename, textlabels, axislabels, dataArray, pop.size()))

    if len(dataObjList) == 0:
        raise ViewException('No files found matching the extensions')

    # sort list from largest -> smallest so plotting is always done that way
    dataObjList.sort(key=lambda x: x.count, reverse=True)

    return dataObjList


class VisualizeFrame(tk.Frame):
    def __init__(self, dataObjList):

        self.dataObjList = dataObjList
        tk.Frame.__init__(self, master=root)

        self.colour_list = ['r.', 'g.', 'y.', 'm.', 'c.', 'b.', 'k.']
        self.xIndex = tk.IntVar()
        self.yIndex = tk.IntVar(value=1)
        self.logx_var = tk.BooleanVar()
        self.logy_var = tk.BooleanVar()
        self.grid_var = tk.BooleanVar()
        self.hist_var = tk.BooleanVar()
        self.gaus_var = tk.IntVar()
        self.plot_var = tk.BooleanVar(value=True)
        self.create_main_panel()


    def create_main_panel(self):
        self.panel = tk.Frame(self)
        self.panel.pack(fill=tk.BOTH, expand=True)  # Pack the panel to fill the entire space of the master window

        self.dpi = 100
        self.fig = Figure((5., 5.), dpi=self.dpi)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.panel)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)  # Pack the canvas to the top of the panel, filling both x and y directions

        self.toolbar = NavigationToolbar(self.canvas, self.panel)
        # self.toolbar.update()
        # self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)  # Pack the canvas again (it's the same canvas but now with toolbar), to the top of the panel, filling both x and y directions

        self.axes = self.fig.add_subplot(111)

        self.buttons1 = tk.Frame(self.panel)
        self.buttons1.pack(side=tk.TOP)

        self.buttons2 = tk.Frame(self.panel)
        self.buttons2.pack(side=tk.TOP)

        tk.Button(self.buttons1, text="Plot", command=self.on_draw_button).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Checkbutton(self.buttons1, text="logX", variable=self.logx_var).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Checkbutton(self.buttons1, text="logY", variable=self.logy_var).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Checkbutton(self.buttons1, text="grid", variable=self.grid_var).pack(side=tk.LEFT, padx=5, pady=5)

        tk.Checkbutton(self.buttons2, text="hist", variable=self.hist_var, command=self._hist).pack(side=tk.LEFT, padx=5, pady=5)
        self.gauss1 = tk.Radiobutton(self.buttons2, text="gauss1", variable=self.gaus_var, value=0, state=tk.DISABLED)
        self.gauss1.pack(side=tk.LEFT, padx=5, pady=5)
        self.gauss2 = tk.Radiobutton(self.buttons2, text="gauss2", variable=self.gaus_var, value=1, state=tk.DISABLED)
        self.gauss2.pack(side=tk.LEFT, padx=5, pady=5)

        # create list of population models. Set all "on" by default
        modelList = [d.name for d in self.dataObjList]
        self.modelCheckList = tk.Listbox(self.panel, selectmode=tk.MULTIPLE)
        for model in modelList:
            self.modelCheckList.insert(tk.END, model)
        self.modelCheckList.select_set(0, tk.END)  # Select all items by default
        self.modelCheckList.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.Y)

        # create LabelFrames for holding Radiobuttons
        self.radioBoxX = tk.LabelFrame(self.panel, text='X Axis')
        self.radioBoxX.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.X)

        self.radioBoxY = tk.LabelFrame(self.panel, text='Y Axis')
        self.radioBoxY.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.Y)

        # Create and pack Radiobuttons for pulsar props inside the LabelFrames
        labels = self.dataObjList[0].textlabels
        for i in range(len(labels)):
            tk.Radiobutton(self.radioBoxX, text=labels[i], variable=self.xIndex, value=i).pack(anchor=tk.W)
            tk.Radiobutton(self.radioBoxY, text=labels[i], variable=self.yIndex, value=i).pack(anchor=tk.W)

        # top horizontal panel - canvas and radio boxes
        self.radioBoxX.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.Y)
        self.radioBoxY.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.Y)
        self.canvas.get_tk_widget().pack(side=tk.LEFT, padx=5, pady=5, fill=tk.BOTH, expand=True)
        # self.modelCheckList.pack(side=tk.LEFT, padx=5, pady=5, fill=tk.Y)

        # add the matplotlib toolbar
        self.toolbar.pack(side=tk.TOP, fill=tk.X)

    def draw_figure(self):
        self.axes.clear()
        if not self.plot_var.get():
            self.axes.text(0.5, 0.5, "Not plottable")
            self.plot_var.set(True)
            return

        for dataObjIndex in self.modelCheckList.curselection():
            dataX = self.dataObjList[dataObjIndex].dataArray[self.xIndex.get()]
            dataY = self.dataObjList[dataObjIndex].dataArray[self.yIndex.get()]
            color = self.colour_list[dataObjIndex]

            if not self.hist_var.get():
                self.axes.plot(dataX, dataY, color,
                                label=self.dataObjList[dataObjIndex].name)
            else:
                # try:
                hist, bins_edges, _ = \
                self.axes.hist(dataX, color=color[0],
                                label=self.dataObjList[dataObjIndex].name,
                                bins=100, density=True, alpha=0.5, edgecolor='k')
                # except:
                #     self.plot_var.set(False)
                #     self.draw_figure()
                #     return
                    
                if not self.gaus_var.get():
                    # try:
                    μ, σ = norm.fit(dataX)
                    xmin, xmax = self.axes.get_xlim()
                    # warnings.simplefilter('error')
                    # warnings.warn("Not allowed")
                    x = np.linspace(xmin, xmax, 1000)
                    p = norm.pdf(x, μ, σ)
                    self.axes.plot(x, p, color[0])
                    # except:
                    #     self.plot_var.set(False)
                    #     self.draw_figure()
                    #     return
                else:
                    # try:
                    bins_centers = (bins_edges[:-1] + bins_edges[1:]) / 2
                    peak_indices = argrelmax(hist)[0]
                    peak_values = hist[peak_indices]
                    top_two_peaks = np.argsort(peak_values)[::-1][:2]
                    μ1, μ2 = bins_centers[peak_indices[top_two_peaks]]

                    def gaussian_mixture(x, A1, μ1, σ1, A2, μ2, σ2):
                        return (A1 * norm.pdf(x, μ1, σ1) +
                                A2 * norm.pdf(x, μ2, σ2))

                    initial_guess = [1, μ1, 1, 1, μ2, 1]
                    params, cov = curve_fit(gaussian_mixture, bins_centers, hist, p0=initial_guess)
                    # warnings.simplefilter('error')
                    # warnings.warn("Not allowed")
                    x = np.linspace(min(dataX), max(dataX), 1000)
                    self.axes.plot(x, gaussian_mixture(x, *params), color[0])
                    # except:
                    #     self.plot_var.set(False)
                    #     self.draw_figure()
                    #     return

                # self.axes.set_title(f"μ={μ:.2f}, σ={σ:.2f}")

            self.axes.legend(loc='upper center', bbox_to_anchor=(0.5, -0.075),
                              ncol=len(self.modelCheckList.curselection()),
                              )

        if len(self.modelCheckList.curselection()) > 0:
            self.axes.set_xlabel(self.dataObjList[0].axislabels[self.xIndex.get()],
                                    fontsize=10)
            self.axes.set_ylabel(self.dataObjList[0].axislabels[self.yIndex.get()],
                                    fontsize=10)

        if self.logx_var.get():
            self.axes.set_xscale('log')

        if self.logy_var.get():
            self.axes.set_yscale('log')

        # for label in self.axes.get_xticklabels():
        #     label.set_fontsize(7)
        # for label in self.axes.get_yticklabels():
        #     label.set_fontsize(7)

        self.axes.grid(self.grid_var.get())
        self.canvas.draw()

    def _hist(self):
        if self.hist_var.get():
            self.gauss1.config(state=tk.NORMAL)
            self.gauss2.config(state=tk.NORMAL)
        else:
            self.gauss1.config(state=tk.DISABLED)
            self.gauss2.config(state=tk.DISABLED)

    def on_draw_button(self):
        # redraw the canvas
        self.draw_figure()

if __name__ == '__main__':
    """ 'Main' function for calling from command line"""
    parser = argparse.ArgumentParser(description='Visualize a population object')
    parser.add_argument('-f', metavar='fname', default='populate.model',
                          help='file containing population model (def="populate.model")')

    parser.add_argument('-frac', nargs=1, type=float, default=None, 
                          help='plot only this fraction of pulsars')

    parser.add_argument('-extn', nargs='+', type=str,
                        default=['.results', '.model'],
                        help='extension(s) to look for when finding population models')
    args = parser.parse_args()

    dataObj = makeDataObj(frac=None, extnList=['.results', '.model'])

    root = tk.Tk()
    frame = VisualizeFrame(dataObj)
    frame.pack(fill=tk.BOTH, expand=True)
    root.mainloop()

