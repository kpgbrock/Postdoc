# We have a data file containing expression quantities and 
# localization parameters for different yeast genes.  We want to
# only take lines that have cytoplasm in the location and that have
# an expression quantity measured.

import sys

import shutil
import os
import time
import datetime
import math
import urllib
from array import array
import re

if len(sys.argv) != 4:
	print("Must include 3 arguments:  input file name with list of protein ids, output file name, end of command")
	print(len(sys.argv))
	

# read in input sequences and process one at a time, then write the output to the output file
inNameSL = sys.argv[1]
outName = sys.argv[2]
lastPartCmd = sys.argv[3]

infileSL = open(inNameSL,'r')
outfile = open(outName,'w')

firstPartCmd = "bash /groups/marks/pipelines/buildali4_realign/buildali.sh -p "
#lastPartCmd = " -c plmc -r full -m 30 -g -W 72"

for line in infileSL:
	
	line = line.strip('\n')
	outfile.write(firstPartCmd+line+"  "+lastPartCmd+"\n")
	
# close files at the end
infileSL.close()
outfile.close()
