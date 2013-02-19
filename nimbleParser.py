#!/usr/bin/env python
"""
*NimbleGen tiling array file parser*

Andre F. Rendeiro <afrendeiro at gmail.com>
2013

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.
"""

import sys
import os
import linecache
import re
from optparse import OptionParser
import math as math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import statsmodels.api as sm #pandas, patsy


def main():
    # option parser
    usage = 'python nimbleparser.py [OPTIONS] <sample_key.txt>'
    parser = OptionParser(usage=usage)
    parser.add_option("-f", "--format",
    type="str", dest="format", default='nimblegen',
    help="Format of .ndf and .pair files, default='Nimblegen', read documentation for more.")
    parser.add_option("-c", "--onecolor",
    type="int", dest="color", default='2',
    help="Number of colours in the array experiment, default=2")
    parser.add_option("-v", "--verbose",
    dest="verbose", action="store_true",
    help="Set flag if verbose mode is to be activated.")

    # read arguments and options
    (options, args) = parser.parse_args()
    if len(args) > 1 or len(args) == 0:
        # return help message if argument number is incorrect
        print __doc__
        parser.print_help()
        sys.exit(0)

    global verbose
    verbose = True #options.verbose

    form = options.format

    # Start program
    if verbose: print __doc__

    sample_key = args[0]

    # Initialize Experiments and Samples
    if verbose: print "\nIdentifying and initializing experiments..."
    nsamples = 0
    experiments = initializeExperiments(sample_key)
    if verbose: print "Experiments are:"
    for exp_name, exp in experiments.iteritems():
        print "\nExperiment: %s" % (exp_name)
        print "  With samples:"
        for smp_name, sample in exp.samples.iteritems():
            print "\t" + str(smp_name)
            nsamples += 1

    # Initialize Array with probes
    if verbose: print "\nInitializing array object..."
    array = Array(nsamples, form)
    array.addProbes(array.filename)
    if verbose: print "  Array " + str(array) + " has " + str(array.probe_number) + " probes."
    if verbose: print "  Done."

    # Add samples to array
    if verbose: print "\nAdding samples to array..."
    for exp_name, exp in experiments.iteritems():
        for smp_name, sample in exp.samples.iteritems():
            sample.addArray(array)
    if verbose: print "  Done."

    # Read intensity for each sample
    if verbose: print "\nReading probe signal intensities..."
    for exp_name, exp in experiments.iteritems():
        if verbose: print "  Experiment %s:" % (exp_name)
        for sample_name, sample in exp.samples.iteritems():
            if verbose: print "\tReading sample %s '%s' probe signals" % (sample.number, sample.name)
            sample.addIntensities()
            if verbose: print "\tFinished reading sample '%s'.\n" % (sample.name)

    # Plot sample density function of intensities
    for exp_name, exp in experiments.iteritems():
        for sample_name, sample in exp.samples.iteritems():
            sample.plotDensity()
            sample.qqPlot()
    if verbose: print "Produced QC plots per sample"

    # Generate MAPlots for each experiment (replicates)
    for exp_name, exp in experiments.iteritems():
        exp.maPlot()
    if verbose: print "Produced MAPlots per group of replicates"


class Array(object):
    """Generic Class for any array in the experiment"""
    def __init__(self, nsamples, form):
        self.filename = self.findNDF()
        self.nsamples = nsamples
        self.format = form  # associate the provided format with the record array characteristics
        #self.x = 3 + self.nsamples # 3 for number of probe atributes to store
        self.y = self.arrayDimensions(self.filename)
        self.createArray(nsamples, form)
        self.probe_number = 0
        self.samples = {}

    def createArray(self, nsamples, form):
        # change to appropriate types latter
        dtype = [('probeID','a12'),('X','i4'),('Y','i4')]
        for sample in range(nsamples + 1):
            dtype.append((str(sample),'f8'))
        self.probes = np.zeros(shape=(self.y, ), dtype=dtype)

    def __str__(self):
        return self.filename

    def findNDF(self, filename='default'):
        """ looks in the working dir for the required .ndf file.
        If not found reads tile layout from the designname info on the pair files first line.
        Returns string.
        """
        found = False
        for filenm in os.listdir('.'):
            if '.ndf' in filenm:
                return filenm
        if not found:
            raise NotImplementedError("Feature not yet implemented. Put .ndf file in current directory")

    def arrayY(self, filename):
        with open(filename) as f:
            for i, l in enumerate(f):
                pass
        return i

    def arrayDimensions(self, filename, skip=1):
        """
        """
        y = self.arrayY(filename)
        return y

    def addProbes(self,  filename, skip=1):
        """Reads .ndf file and returns probes atributes"""
        i = 0
        fl = open(filename, 'r')
        for line in fl.readlines():
            if i < skip:
                i += 1
            else:
                line = line.rstrip().split('\t')
                probeID = line[12]
                X = line[15]
                Y = line[16]
                self.probes[i - skip]['probeID'] = probeID
                self.probes[i - skip]['X'] = X
                self.probes[i - skip]['Y'] = Y
                if i == 1: print self.probes[i - skip]
                self.probe_number += 1
                i += 1
        fl.close()


class Experiment(object):
    """Class to describe experiments (sets of sample replicates and their inputs"""
    def __init__(self, name):
        self.name = name
        self.samples = {}

    def subtract(self, ls):
        x = 0
        for i in range(len(ls)):
            if i == 0:
                x = ls.pop(0)
            else:
                x -= ls.pop(0)
            i += 1
        return x

    def computeMA(self, samples):
        """ Computes M [average(sum(log2(signals)))] and A [subtraction(log2(signals))] values per row.
            input: list of Sample objects
            Returns two lists of M and A values.
        """
        nsamples = len(samples)
        A = []
        M = []
        i = 0
        for row in range(samples[0].array.y):
            Xs = []
            for sample in samples:
                x = math.log(float(sample.array.probes[i][2 + sample.number]), 2)
                Xs. append(x)
            A.append(sum(Xs) / nsamples)
            M.append(self.subtract(Xs))
        return M, A


    def maPlot(self):
        """ Sorts samples per type and plots MA values.
            input: Experiment object
        """
        e = []
        i = []
        for smp_name, sample in self.samples.iteritems():
            if type(sample) == ChIP:
                e.append(sample)
            else:
                i.append(sample)

        for s in [e, i]:
            M, A = self.computeMA(s)
            plt.figure('en')
            plt.plot(A, M, 'ro')
            plt.xlabel('M')
            plt.ylabel('A')
            plt.title('')
            plt.savefig("%s_.png") % (self.name) # don't understand why name = obj


class Sample(object):
    """Generic Class for a Sample made in a hybridization experiment"""
    def __init__(self, filename, number):
        self.filename = filename
        self.name = re.split('.pair', self.filename)[0]
        self.number = number
        self.array = ''

    def addArray(self, array):
        self.array = array
        array.samples[self.name] = self

    def isFileStandard(self, filename, skip=2):
        standard_header = ['IMAGE_ID', 'GENE_EXPR_OPTION', 'SEQ_ID', 'PROBE_ID', 'POSITION', 'X', 'Y', 'MATCH_INDEX', 'SEQ_URL', 'PM', 'MM']
        if not ".pair" in filename:
            raise ValueError("File not in .pair format")
        elif not filename in os.listdir('.'):
            raise IOError("File not found")
        elif getHeader(filename, skip) != standard_header:
            raise ValueError("File not in standard .pair format")
        else:
            return True

    def addIntensities(self, skip=2):
        """ Reads Nimble file of Sample. Returns list of lists
        skip = starting line of data
        Adds relevant info to the probes.
        """
        if self.isFileStandard(self.filename, skip):
            i = 0
            fil = open(self.filename, 'r')
            for line in fil:
                if i < skip:
                    i += 1
                else:
                    line = line.rstrip().split('\t')
                    probeID = line[3]
                    PM = line[9]
                    print probeID
                    print np.where(self.array.probes[i - skip] == probeID) # Cant find it in array!!!!
                    x = np.where(self.array.probes[i - skip] == probeID)
                    self.array.probes[x][2 + self.number] = (PM)
            fil.close()

    def plotDensity(self):
        """ Plots sample signals density function"""
        data = self.array.probes[:, 2 + self.number]  # log2
        density = st.gaussian_kde(data)
        xs = np.linspace(0,8,200)
        density.covariance_factor = lambda : .25
        density._compute_covariance()
        plt.figure(self.number)
        plt.plot(xs,density(xs))
        plt.xlabel('Probe intensities')
        plt.ylabel('Density')
        plt.title('Probe intensities for %s. n= %s' % (self.name, str(self.array.probe_number)))
        plt.savefig("%s_density.png" % (self.name))

    #def qqPlot(self):
        #""" Plots sample signals against theorethical distribution"""
        #data = self.array.probes[:, 2 + self.number]
        #plt.figure(self.number)
        #osm, osr = st.probplot(data, fit=0, dist='norm')  # compute
        #plt.plot(osm, osr, '.')                              # plot
        ## Add fit line
        ##(osm, osr), (m, b, r) = st.probplot(data, dist='norm')  # compute
        ##osmf = osm.take([0, -1])  # endpoints
        ##osrf = m * osmf + b       # fit line
        ##ax.plot(osm, osr, '.', osmf, osrf, '-')
        #plt.xlabel('Theoretical quantiles')
        #plt.ylabel('Sample quantiles')
        #plt.title('Probe intensities for %s' % (self.name))
        #plt.savefig("%s_qqprob.png" % (self.name))

    def qqPlot(self):
        """ Plots sample signals against theorethical distribution"""
        import statsmodels.api as sm #pandas, patsy
        import matplotlib.pyplot as plt
        data = self.array.probes[:, 2 + self.number]  # add log2
        plt.figure(self.number)
        fig = sm.qqplot(data)
        plt.xlabel('Theoretical quantiles')
        plt.ylabel('Sample quantiles')
        plt.title('Probe intensities for %s' % (self.name))
        plt.savefig("%s_qqprob.png" % (self.name))


class ChIP(Sample):
    """ ChIP samples"""
    def __init__(self, filename, number):
        Sample.__init__(self, filename, number)

class Input(Sample):
    """ Input samples """
    def __init__(self, filename, number):
        Sample.__init__(self, filename, number)


def initializeExperiments(filename, skip=0):
    """ Creates Sample objects to be included in the analysis"""
    i = 0
    s = 1
    experiments = {}
    for line in open(filename, 'r'):
        if i < skip:
            i += 1
        else:
            line = line.rstrip().split('\t')
            filename = line[0]
            sname = re.split('.pair', filename)[0]
            number = s
            stype = line[1]
            name = line[2]
            if not name in experiments:
                experiments[name] = Experiment(name)  # make new experiment
            if stype == 'e':
                experiments[name].samples[sname] = ChIP(filename, number)
            elif stype == 'i':
                experiments[name].samples[sname] = Input(filename, number)
            else:
                raise ValueError("'Sample type' in sample/key file isn't 'e' or 'i'. See README.md")
            i += 1
            s += 1
    return experiments


def getHeader(filename, skip=2):
    """ Extracts the column headers from a file. Returns list
    filename = name of the file
    skip = number of lines before the header line
    """
    return linecache.getline(filename, skip).rstrip().split('\t')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nProgram canceled by user!\n")
        sys.exit(0)