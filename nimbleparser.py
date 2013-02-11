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
        print __doc__
        parser.print_help()
        sys.exit(0)

    #if options[0] == 'F':
    #    verbose = False
    sample_key = args[0]

    #initialize Samples (keep list of samples in experiment)
    samples = initializeSamples(sample_key)
    #debug
    print "Samples are: "
    for sample in samples:
        print sample.filename

    #initialize Array (& probes)
    array = Array()
    #debug
    print "Array " + str(array) + " has " + str(array.probe_number) + " probes."

    #add array to samples
    for sample in samples:
        sample.addArray(array)

    #read intensity for each sample
    for sample in samples:
        print "Reading ", sample.filename, " probe intensities..."
        sample.readIntensities(sample.filename)

    ##debug
    #for probe in array.probes[:20]:
        #print probe.intensity

class Array(object):
    """Generic Class for any array in the experiment"""
    def __init__(self):
        self.name = self.getArrayLayout()
        self.probes = {}
        self.probe_number = 0
        self.readNDF(self.getArrayLayout())

    def __str__(self):
        return self.name

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
            #getNameFromFile(filename)

        def getNameFromFile(self, filename):
            """        """
            pass
            #first_line = linecache.getline(filename, 0).rstrip().split('\t')
            #if '.ndf' in first_line:
                #find name + '.ndf'

    def addProbe(self, probeID, probe_seq, X, Y):
        """Adds probe to array"""
        self.probes[probe_seq] = [probeID, X, Y]
        self.probe_number += 1

    def readNDF(self, filename, skip=1):
        """Reads .ndf file and returns probes atributes"""
        i = 0
        fl = open(filename, 'r')
        for line in fl.readlines():
            if i < skip:
                i += 1
            else:
                line = line.rstrip().split('\t')
                probeID = line[5]
                probe_seq = line[12]
                X = line[15]
                Y = line[16]
                self.addProbe(probeID, probe_seq, X, Y)
                i += 1
        fl.close()


class Sample(object):
    """Generic Class for a Sample made in a hybridization experiment"""
    def __init__(self, filename, sample_type, sample_experiment):
        self.filename = filename
        self.sample_type = sample_type  # consider making separate subclasses for IP and input Samples
        self.sample_experiment = sample_experiment
        self.array = ''

    def addArray(self, array):
        self.array = array

    def readIntensities(self, filename, skip=2):
        """ Reads Nimble file of Sample. Returns list of lists
        skip = starting line of data
        Adds relevant info to the probes.
        """
        def isFileStandard(filename, skip=2):
            standard_header = ['IMAGE_ID', 'GENE_EXPR_OPTION', 'SEQ_ID', 'PROBE_ID', 'POSITION', 'X', 'Y', 'MATCH_INDEX', 'SEQ_URL', 'PM', 'MM']
            if not ".pair" in filename:
                raise ValueError("File not in .pair format")
            elif not filename in os.listdir('.'):
                raise IOError("File not found")
            elif getHeader(filename, skip) != standard_header:
                raise ValueError("File not in standard .pair format")
            else:
                return True

        if isFileStandard(filename, skip):
            i = 0
            fil = open(filename, 'r')
            for line in fil:
                if i < skip:
                    i += 1
                else:
                    line = line.rstrip().split('\t')
                    probeID = line[3]
                    PM = line[9]
                    self.array.probes[probeID].append({self.filename:PM})
            fil.close()


def initializeSamples(filename, skip=0):
    """ Creates Sample objects to be included in the analysis"""
    i = 0
    samples = []
    for line in open(filename, 'r'):
        if i < skip:
            i += 1
        else:
            line = line.rstrip().split('\t')
            sample_file = line[0]
            sample_type = line[1]
            sample_experiment = line[2]
            samples += [Sample(sample_file, sample_type, sample_experiment)]
            i += 1
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