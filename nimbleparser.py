#!/usr/bin/env python
"""
NimbleGen tiling array file parser

Andre F. Rendeiro <afrendeiro at gmail.com>
2013

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.
"""

import sys
import linecache
import os
from optparse import OptionParser
import numpy as np

def main():
    # option parser
    usage = 'python nimbleparser.py [OPTIONS] <sample_key.txt>'
    parser = OptionParser(usage=usage)
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

    # Start program
    if verbose: print __doc__

    sample_key = args[0]

    #initialize Samples (keep list of samples in experiment)
    samples = initializeSamples(sample_key)
    if verbose:
        print "Samples are:"
        for sample in samples: print sample.filename

    #initialize Array
    array = Array(len(samples))

    #add array to samples
    for sample in samples:
        sample.addArray(array)

    #populate array with probe values
    array.addProbes(array.filename)
    if verbose: print "Array " + str(array) + " has " + str(array.probe_number) + " probes."

    #read intensity for each sample
    for sample in samples:
        if verbose: print "Reading sample %s (%s) probe signals" % (sample.number, sample.filename)
        sample.addIntensities(sample.filename)
        if verbose: print "Finished reading %s.\n" % (sample.filename)


class Array(object):
    """Generic Class for any array in the experiment"""
    def __init__(self, lsamples):
        self.filename = self.findNDF()
        self.lsamples = lsamples
        self.x = 3 + self.lsamples # 3 for number of probe atributes to store
        self.y = self.arrayDimensions(self.filename)
        self.probes = np.zeros(shape=(self.y, self.x), dtype=('a10')) # change to appropriate type latter
        self.probe_number = 0
        self.samples = []

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
                self.probes[i - skip][:3] = (probeID, X, Y)
                self.probe_number += 1
                i += 1
        fl.close()


class Sample(object):
    """Generic Class for a Sample made in a hybridization experiment"""
    def __init__(self, filename, sample_number, sample_type, sample_experiment):
        self.filename = filename
        self.number = sample_number
        self.sample_type = sample_type  # consider making separate subclasses for IP and input Samples
        self.sample_experiment = sample_experiment
        self.array = ''

    def addArray(self, array):
        self.array = array
        array.samples.append(self)

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

    def addIntensities(self, filename, skip=2):
        """ Reads Nimble file of Sample. Returns list of lists
        skip = starting line of data
        Adds relevant info to the probes.
        """
        if self.isFileStandard(filename, skip):
            i = 0
            fil = open(filename, 'r')
            for line in fil:
                if i < skip:
                    i += 1
                else:
                    line = line.rstrip().split('\t')
                    probeID = line[3]
                    PM = line[9]
                    x = np.where(self.array.probes==probeID)[0][0]
                    self.array.probes[x][2 + self.number] = (PM)
            fil.close()


def initializeSamples(filename, skip=0):
    """ Creates Sample objects to be included in the analysis"""
    i = 0
    s = 1
    samples = []
    for line in open(filename, 'r'):
        if i < skip:
            i += 1
        else:
            line = line.rstrip().split('\t')
            sample_file = line[0]
            sample_number = s
            sample_type = line[1]
            sample_experiment = line[2]
            samples += [Sample(sample_file, sample_number, sample_type, sample_experiment)]
            i += 1
            s += 1
    return samples


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