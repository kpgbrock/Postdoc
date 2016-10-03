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

inNameSL = 'yeastPrionDomains.txt'
OutName = 'yeastPrionDomainsProcessed.txt'
    
# read in input sequences and process one at a time, then write the output to the output file
infileSL = open(inNameSL,'r')
outfile = open(OutName,'w')

for line in infileSL:
	
	line = line.strip('\n')
	line = line.split()
	outfile.write(line[0]+" -r "+line[1]+"-"+line[2]+"\n")
	
# close files at the end
infileSL.close()
outfile.close()
