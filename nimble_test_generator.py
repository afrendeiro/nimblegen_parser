#!/usr/bin/env python

import string
import random
import csv



# PARAMETERS TO BE TUNED:
# length of files to be generated
number_pair_files = 2
length = 1000
ndf_files = 1


# DONT TOUCH!
probeID = []
CANDIDATE_CHARS = string.ascii_letters+string.digits

# for pair files
ncol = 11

# generate first pair file
fl = open('file1.pair', 'w+')
for ln in range(length):
    line = []
    for col in range(ncol):
        if col == 3:
            probe = ''.join(random.choice(CANDIDATE_CHARS) for _ in range(10))
            line.append(probe)
            probeID.append(probe)
        else:
            line.append(''.join(random.choice(CANDIDATE_CHARS) for _ in range(10)))
    wr = csv.writer(fl, dialect='excel-tab')
    wr.writerow(line)

fl.close()

# generate extra .pair files
for fil in range((number_pair_files - 1)):
    fl = open('file' + str(fil + 1) + '.pair', 'w+')
    for ln in range(length):
        line = []
        for col in range(ncol):
            if col == 3:
                line.append(probeID[ln])
            else:
                line.append(''.join(random.choice(CANDIDATE_CHARS) for _ in range(10)))
        wr = csv.writer(fl, dialect='excel-tab')
        wr.writerow(line)

    fl.close()



# for ndf files
ncol = 18

fl = open('file.ndf', 'w+')
for ln in range(length):
    line = []
    for col in range(ncol):
        if col == 12:
            line.append(probeID[ln])
        else:
            line.append(''.join(random.choice(CANDIDATE_CHARS) for _ in range(10)))
    wr = csv.writer(fl, dialect='excel-tab')
    wr.writerow(line)

fl.close()
