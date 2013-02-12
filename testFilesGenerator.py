#!/usr/bin/env python
"""
Generator of test files for "NimbleGen tiling array file parser"

Andre F. Rendeiro <afrendeiro at gmail.com>
2012

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License version 2 as
published by the Free Software Foundation.
"""

import sys
from optparse import OptionParser
import string
import random
import csv

def main():
    # option parser
    usage = 'python nimble_test_generator.py [OPTIONS]'
    parser = OptionParser(usage=usage)
    parser.add_option("-s", "--samples",
    type="int", dest="samples", default='2',
    help="Number of samples to be generated. Default=2")
    parser.add_option("-l", "--length",
    type="int", dest="length", default='10000',
    help="Number of probes in the array to be generated. Default=1000")

    (options, args) = parser.parse_args()

    if len(args) > 0:
        # return help message if argument number is incorrect
        print __doc__
        parser.print_help()
        sys.exit(0)

    # PARAMETERS TO BE TUNED:
    number_pair_files = options.samples
    length = options.length

    print "Generating 1 .nsf file, and %s .pair files with %s probes\n" % (number_pair_files, length)

    probeID = generatePairs(number_pair_files, length)
    generateNDF(probeID, length)
    generateSampleKey(number_pair_files)

    print "\nFinished!\n"


def generatePairs(number_pair_files, length):
    """ Generates pairs files all with same probeIDs"""
    probeID = []
    candidate_chars = string.ascii_letters+string.digits
    ncol = 11
    pair_header = ['#', 'software=NimbleScan', 'version=2.5.28', 'imagefile=...']
    pair_header2 = ['IMAGE_ID', 'GENE_EXPR_OPTION', 'SEQ_ID', 'PROBE_ID', 'POSITION', 'X', 'Y', 'MATCH_INDEX', 'SEQ_URL', 'PM', 'MM']

    # generate first pair file
    fl = open('sample1.pair', 'w+')
    wr = csv.writer(fl, dialect='excel-tab', quoting=csv.QUOTE_NONE)
    print "Started sample1.pair"
    for ln in range(length):
        line = []
        if ln == 0:
            wr.writerow(pair_header)
            wr.writerow(pair_header2)
        for col in range(ncol):
            if col == 3:
                probe = ''.join(random.choice(candidate_chars) for _ in range(10))
                line.append(probe)
                probeID.append(probe)
            elif col == 9:
                line.append(random.choice(range(0,65535)))
            else:
                line.append(''.join(random.choice(candidate_chars) for _ in range(10)))
        wr.writerow(line)
    print "Finished sample1.pair\n"
    fl.close()

    # generate extra .pair files
    for fil in range((number_pair_files - 1)):
        fl = open('sample' + str(fil + 2) + '.pair', 'w+')
        wr = csv.writer(fl, dialect='excel-tab', quoting=csv.QUOTE_NONE)
        print "Started file%s.pair" % (fil + 2)
        for ln in range(length):
            line = []
            if ln == 0:
                wr.writerow(pair_header)
                wr.writerow(pair_header2)
            for col in range(ncol):
                if col == 3:
                    line.append(probeID[ln])
                elif col == 9:
                    line.append(random.choice(range(0,65535)))
                else:
                    line.append(''.join(random.choice(candidate_chars) for _ in range(10)))
            wr.writerow(line)
        print "Finished sample%s.pair\n" % (fil + 2)
        fl.close()

    return probeID


def generateNDF(probeID, length):
    """ Generates a NDF file with the same probes as the pair files """
    ncol = 18
    candidate_chars = string.ascii_letters+string.digits
    ndf_header = ['PROBE_DESIGN_ID','CONTAINER','DESIGN_NOTE','SELECTION_CRITERIA','SEQ_ID','PROBE_SEQUENCE','MISMATCH','MATCH_INDEX','FEATURE_ID', 'ROW_NUM', 'COL_NUM', 'PROBE_CLASS', 'PROBE_ID','POSITION','DESIGN_ID', 'X','Y','DMD']

    fl = open('array.ndf', 'w+')
    wr = csv.writer(fl, dialect='excel-tab', quoting=csv.QUOTE_NONE)
    print "Started array.nsf"
    for ln in range(length):
        line = []
        if ln == 0:
            wr.writerow(ndf_header)
        for col in range(ncol):
            if col == 12:
                line.append(probeID[ln])
            else:
                line.append(''.join(random.choice(candidate_chars) for _ in range(10)))
        wr.writerow(line)
    print "Finished array.nsf\n"
    fl.close()


def generateSampleKey(number_pair_files):
    """ Generates a sample/key file describing the samples generated"""

    fl = open('sample_key.txt', 'w+')
    wr = csv.writer(fl, dialect='excel-tab', quoting=csv.QUOTE_NONE)
    for ln in range(number_pair_files):
        line = []
        line.append('sample' + str(ln + 1) + '.pair')
        if ln % 2 == 0:
            line.append('e')
        else:
            line.append('i')
        line.append('H3K27me3')
        wr.writerow(line)
    print "Finished sample_key.txt\n"
    fl.close()

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nProgram canceled by user!\n")
        sys.exit(0)