#!/usr/bin/env python
"""
NimbleGen tiling array file parser

Andre F. Rendeiro <afrendeiro at gmail.com>
2012

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.
"""

import sys
import linecache
import os
from optparse import OptionParser

def main():
    # option parser
    usage = 'python nimbleparser.py [OPTIONS] <sample_key.txt>'
    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--onecolor",
    type="int", dest="c_option", default='2',
    help="Number of colours in the array experiment, default=2")
    parser.add_option("-v", "--verbose",
    type="str", dest="v_option", default='false',
    help="'F' if verbose mode is to be deactivated. Default='T'")

    # read arguments and options
    (options, args) = parser.parse_args()
    if len(args) > 1 or len(args) == 0:
        # return help message if argument number is incorrect
        raise ValueError("Incorrect number of arguments provided. type 'python nimbleparser.py -h' to see help")
        #print __doc__
        #parser.print_help()
        sys.exit(0)
    
    #if options[0] == 'F':
    #    verbose = False
    sample_key = args[0]

    #initialize Samples (keep list of samples in experiment)
    samples = initializeSamples(sample_key)

    #initialize Array (& probes)
    array = Array()

    #add array to samples
    for sample in samples:
        sample.addArray(array)

    #read intensity for each sample
    #for sample in samples:
        #if verbose:
        #    print "Reading ", sample.filename, " probe intensities..."
        #sample.readIntensities(sample.filename)
    samples[0].readIntensities(samples[0].filename)

    #Reporting...
    print "Samples are: "
    for sample in samples:
        print sample.filename
    print "Array ", array, " has ", array.getArrayLenght(), " probes."


    #Quering probe intensity values
    print "Probe 1 intensity values:"
    print array.probes[0].intensity


class Array(object):
    """Generic Class for any array in the experiment"""
    def __init__(self):
        self.name = self.getArrayLayout()
        self.probes = []
        self.addProbes(self.getArrayLayout())

    def __str__(self):
        """Returns a string representation of self"""
        return self.name

    def getArrayLenght(self):
        return len(self.probes)

    def getArrayLayout(self, filename='default'):
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
            #first_line = linecache.getline(filename, 0).rstrip().split('\t')
            #if '.ndf' in first_line:
                #find name + '.ndf'

    def addProbes(self, filename, skip=1):
        """Reads .ndf file and returns probes atributes"""
        i = 0
        for line in open(filename, 'r'):
            if i < skip:
                i += 1
            else:
                line = line.rstrip().split('\t')
                vars()['p' + str(line[12])] = Probe(self, line[1], line[5], line[12], line[15], line[16])
                self.probes.append(vars()['p' + str(line[12])])
                i += 1


class Probe(object):
    """Class to describe probes, their fixed atributes and intensity per sample"""
    def __init__(self, array, container, probe_seq, probeID, X, Y):
        self.array = array
        self.container = container
        self.probe_seq = probe_seq
        self.probeID = probeID
        self.X = X
        self.Y = Y
        self.intensity = {}  # a dic of {sample:probe_intensity}

    def __str__(self):
        """Returns a string representation of self"""
        return self.probeID

    def getprobeID(self):
        return self.probeID

    def addExperimentIntensity(self, sample, PM):
        self.intensity[sample] = PM


class Sample(object):
    """Generic Class for a Sample made in a hybridization experiment"""
    def __init__(self, filename, sample_type, sample_name):
        self.filename = filename
        self.sample_type = sample_type  # consider making separate subclasses for IP and input Samples
        self.sample_name = sample_name
        self.array = ''

    def getSampleName(self):
        """Returns a string representation of self"""
        return self.filename

    def addArray(self, array):
        self.array = array

    def readIntensities(self, filename, skip=2):
        """ Reads Nimble file of Sample. Returns list of lists
        skip = starting line of data
        Adds relevant info to the probes.
        """
        standard_header = ['IMAGE_ID', 'GENE_EXPR_OPTION', 'SEQ_ID', 'PROBE_ID', 'POSITION', 'X', 'Y', 'MATCH_INDEX', 'SEQ_URL', 'PM', 'MM']
        if not ".pair" in filename:
            raise ValueError("File not in .pair format")
        if not filename in os.listdir('.'):
            raise IOError("File not found")
        if getHeader(filename, skip) != standard_header:
            raise ValueError("File not in standard .pair format")
        i = 0
        for line in open(filename, 'r'):
            if i < skip:
                i += 1
            else:
                line = line.rstrip().split('\t')
                probeID = line[3]
                PM = line[9]
                probeInd = self.array.probes.index(probeID)
                self.array.probes[probeInd].addExperimentIntensity(self, PM)


def initializeSamples(filename, skip=1):
    """ Creates Sample objects to be included in the analysis"""
    i = 0
    samples = []
    for line in open(filename, 'r'):
        if i < skip:
            i += 1
        else:
            line = line.rstrip().split('\t')
            vars()[str(line[2])] = Sample(line[0], line[1], line[2])
            samples.append(vars()[line[2]])
            i += 1
    return samples


def getHeader(filename, skip=2):
    """ Extracts the column headers from a file. Returns list
    filename = name of the file
    skip = number of lines before the header line
    """
    return linecache.getline(filename, skip).rstrip().split('\t')





def readPairs(filenames):
    """read a set of .pairs files and returns perfect matches per probe"""
    table = []
    for fl in filenames:
        sample = readNimble(fl)
        table.append(sample)
    return table        
#readPairs('asd.pair')



if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("Program canceled by user!\n")
        sys.exit(0)