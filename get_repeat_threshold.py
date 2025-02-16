#!/usr/bin/python

# Copyright 2015 Lina Sieverling

# This file is part of TelomereHunter.

# TelomereHunter is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# TelomereHunter is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with TelomereHunter.  If not, see <http://www.gnu.org/licenses/>.


############################################################################################
### Returns the repeat threshold from the input BAM file (treshold = read length*6/100). ###
### If the read lengths of the first 100 reads differ, it returns "n".                   ###
############################################################################################

import os
import sys
import pysam


def get_repeat_threshold(bam_file):

      # open input bam_file for reading
      bamfile = pysam.Samfile( bam_file, "rb" )


      ### check if read lengths of the first 100 non-supplementary or secondary alignments are the same
      cntr=0
      different_read_lengths=False

      for read in bamfile.fetch(until_eof=True):   

          if read.is_secondary:        #skips all secondary alignments
              continue

          if read.flag >= 2048:        # skip supplementary alignments
              continue
            
          cntr+=1
          
          if cntr == 1:
              read_length = len(read.seq)
          else:
              if read_length != len(read.seq):
                  different_read_lengths=True
                  repeat_threshold = "n"
                  print ("Read lengths in sample differ: repeat threshold will be set individually for each read.")
                  break
              
          if cntr == 100:    
              break


      #if they are the same, set threshold
      if different_read_lengths==False:
          repeat_threshold = int(round(float(read_length)*6/100))
          print ("Read length is " + str(read_length) + ". Repeat threshold is set to " + str(repeat_threshold) + ".")


      return repeat_threshold
    
    
