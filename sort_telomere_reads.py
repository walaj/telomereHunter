#!/usr/bin/python

# Copyright 2015 Lina Sieverling, Philip Ginsbach, Lars Feuerbach

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


import os
import sys
import re
import pysam
import subprocess
from telomerehunter.filter_telomere_reads import getReverseComplement

###############################################################################################################################
### sorts the reads from an input BAM file (needs to be name sorted!) containing telomere reads into 4 different fractions: ###
### intratelomeric, junction-spanning, subtelomeric, intrachromosomal. Additionally, a table containing the number of       ###
### telomere reads and the pattern occurrences per chromosome band is written.                                              ###
###############################################################################################################################

#def sort_telomere_reads(input_dir, band_file, out_dir, pid, mapq_threshold, t_type, c_type, g_type, j_type, n_type):
def sort_telomere_reads(input_dir, band_file, out_dir, pid, mapq_threshold, repeats):

        #####################
        ### get patterns  ###
        #####################

        patterns = []

        for repeat in repeats:
                patterns.append(repeat)
                patterns.append(getReverseComplement(repeat))


        # open input bam_file for reading
#        bamfile = pysam.Samfile( input_dir + "/" + pid + "_filtered_name_sorted.bam" , "rb" )
        print(f"Attempting to open BAM file: {input_dir}/{pid}_filtered_name_sorted.bam")
        try:
                bamfile = pysam.Samfile(input_dir + "/" + pid + "_filtered_name_sorted.bam", "rb")
    
                # Check if the BAM file was opened successfully
                print(f"BAM file {input_dir}/{pid}_filtered_name_sorted.bam opened successfully.")

        except Exception as e:
                print(f"Error opening BAM file: {e}")

        # open filtered files for writing
        intratelomeric_file = pysam.Samfile(out_dir + "/" + pid + "_filtered_intratelomeric.bam", "wb", template=bamfile)
        junctionspanning_file = pysam.Samfile(out_dir + "/" + pid + "_filtered_junctionspanning.bam", "wb", template=bamfile)
        subtelomeric_file = pysam.Samfile(out_dir + "/" + pid + "_filtered_subtelomeric.bam", "wb", template=bamfile)
        intrachromosomal_file = pysam.Samfile(out_dir + "/" + pid + "_filtered_intrachromosomal.bam", "wb", template=bamfile)


        ############################
        ### make chromosome list ###
        ############################

        # get references
        references = bamfile.references 

        # autodetect chromosome prefix in BAM file
        if references[0][0:3] == "chr":
                bam_chr_prefix = "chr"
        else:
                bam_chr_prefix = ''

        # Generate list of chromosome names.
        chromosome_list = [ str(i) for i in range(1,22+1) ] + ["X","Y"]
        chromosome_list_with_prefix = [bam_chr_prefix + chr for chr in chromosome_list]

        ####################################
        ### make band and spectrum lists ###
        ####################################

        # Read list of chromosome bands for each chromosome, strip chromosome prefixes.
        bands_list = { chr:{ "band_name":[], "end":[] } for chr in chromosome_list + ["unmapped"] }

        # make spectrum list which will be written to spectrum file
        spectrum_list = { chr:{} for chr in chromosome_list }

        # add junction_spanning fields for each chromosome to spectrum list
        for chr in chromosome_list:
                spectrum_list[chr]["junction1"] = {pattern:0 for pattern in patterns}
                spectrum_list[chr]["junction1"]["other"]=0.0
                spectrum_list[chr]["junction1"]["reads_with_pattern"]=0
                spectrum_list[chr]["junction2"] = {pattern:0 for pattern in patterns}
                spectrum_list[chr]["junction2"]["reads_with_pattern"]=0
                spectrum_list[chr]["junction2"]["other"]=0.0


        for line in open( band_file, "r" ):

                try:
                        line = line.rstrip().split()
                        
        #               start = line[1]
                        end = line[2]
                        band_name = line[3]

                        chrom_name = ""

                        if line[0][:3] == "chr":
                                chrom_name = line[0][3:]
                        else:
                                chrom_name = line[0]

                        bands_list[chrom_name]["band_name"] += [band_name]
                        bands_list[chrom_name]["end"] += [int(end)]
                        
                        spectrum_list[chrom_name][band_name] = {pattern:0 for pattern in patterns}
                        spectrum_list[chrom_name][band_name]["other"]=0.0
                        spectrum_list[chrom_name][band_name]["reads_with_pattern"]=0
                        

                except:

                        print( "Invalid line in banding file: \'"+" ".join( line )+"\'" )


        # add unmapped to spectrum list
        spectrum_list["unmapped"] = {"unmapped": {pattern:0 for pattern in patterns} }
        spectrum_list["unmapped"]["unmapped"]["other"]=0.0
        spectrum_list["unmapped"]["unmapped"]["reads_with_pattern"]=0

        bands_list["unmapped"]["band_name"] += ["unmapped"]
        bands_list["unmapped"]["end"] += [0]

        #flag to break out of double while loop
        break_flag = False


        #################################################################
        ### go through name sorted BAM file containing telomere reads ###
        ### => always look at mate pairs (if possible)                ###
        #################################################################

        for read1 in bamfile:
                try:
                        # Try to get the next read
                        read2 = next(bamfile)
                except StopIteration:
                        # This handles the case where there is no second read available
                        sort_reads_without_mate(read1, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
                                 intratelomeric_file, subtelomeric_file, intrachromosomal_file)
                        break


                # Reads without mates
                read1_first_flag = True  # read 1 is first read

                while read1.qname != read2.qname:
                        if read1_first_flag:
                                sort_reads_without_mate(read1, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
                                                        intratelomeric_file, subtelomeric_file, intrachromosomal_file)
            
                                try:
                                        read1 = next(bamfile)  # Get next read
                                except StopIteration:
                                        # Handle case where no more reads are available
                                        sort_reads_without_mate(read2, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
                                                                intratelomeric_file, subtelomeric_file, intrachromosomal_file)
                                        break  # Exit loop after processing read2
                                read1_first_flag = False  # read2 is now first read
                        else:
                                sort_reads_without_mate(read2, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
                                     intratelomeric_file, subtelomeric_file, intrachromosomal_file)

                                try:
                                        read2 = next(bamfile)  # Get next read
                                except StopIteration:
                                        # Handle case where no more reads are available
                                        sort_reads_without_mate(read1, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
                                                                intratelomeric_file, subtelomeric_file, intrachromosomal_file)
                                        break  # Exit loop after processing read1
                                read1_first_flag = True


                if break_flag == True:                  # break out of outer while loop if the end of the input BAM file is reached
                        break

                ### reads with mates

                read1_chromosome, read1_band, read1_is_unmapped, read1_junctionspanning, read1_junction = read_check (read1, references, bands_list, chromosome_list, mapq_threshold)
                read2_chromosome, read2_band, read2_is_unmapped, read2_junctionspanning, read2_junction = read_check (read2, references, bands_list, chromosome_list, mapq_threshold)

                # INTRATELOMERIC: both mates are unmapped
                if read1_is_unmapped==True and read2_is_unmapped==True:

                        sort_read_with_mate(read=read1, fraction_file=intratelomeric_file, chromosome="unmapped", band="unmapped", spectrum_list=spectrum_list, patterns=patterns)
                        sort_read_with_mate(read=read2, fraction_file=intratelomeric_file, chromosome="unmapped", band="unmapped", spectrum_list=spectrum_list, patterns=patterns)


                # JUNCTION SPANNING: one mate is unmapped and the other mate is mapped to first or last band of chromosome
                elif read1_is_unmapped==True and read2_junctionspanning==True or read2_is_unmapped==True and read1_junctionspanning==True:      

                        if read1_is_unmapped==True:
                                chromosome = read2_chromosome
                                band = read2_junction
                        else:
                                chromosome = read1_chromosome
                                band = read1_junction

                        sort_read_with_mate(read=read1, fraction_file=junctionspanning_file, chromosome=chromosome, band=band, spectrum_list=spectrum_list, patterns=patterns)
                        sort_read_with_mate(read=read2, fraction_file=junctionspanning_file, chromosome=chromosome, band=band, spectrum_list=spectrum_list, patterns=patterns)  
                        
        

                # SUBTELOMERIC: both mates are mapped to first or last band or chromosome (and have similar mapping positions)
                elif read1_junctionspanning==True and read2_junctionspanning==True:     

                        sort_read_with_mate(read=read1, fraction_file=subtelomeric_file, chromosome=read1_chromosome, band=read1_band, spectrum_list=spectrum_list, patterns=patterns)
                        sort_read_with_mate(read=read2, fraction_file=subtelomeric_file, chromosome=read2_chromosome, band=read2_band, spectrum_list=spectrum_list, patterns=patterns)  


                
                # SUBTELOMERIC/INTRACHROMOSOMAL: one read is subtelomeric, the other is intra-chromosomal => count and sort individually
                elif read1_junctionspanning==True or read2_junctionspanning==True:
                  
                        if read1_junctionspanning==True:
                                read1_file = subtelomeric_file
                                read2_file = intrachromosomal_file
                        else:
                                read2_file = subtelomeric_file
                                read1_file = intrachromosomal_file

                        sort_read_with_mate(read=read1, fraction_file=read1_file, chromosome=read1_chromosome, band=read1_band, spectrum_list=spectrum_list, patterns=patterns)
                        sort_read_with_mate(read=read2, fraction_file=read2_file, chromosome=read2_chromosome, band=read2_band, spectrum_list=spectrum_list, patterns=patterns)         

                
                
                # INTRACHROMOSOMAL: one read is intra-chromosomal, the other is unmapped => count both to position of mapped read
                elif read1_is_unmapped==True and read2_junctionspanning==False or read2_is_unmapped==True and read1_junctionspanning==False:
                        
                        if read1_is_unmapped==True:
                                chromosome = read2_chromosome
                                band = read2_band
                        else:
                                chromosome = read1_chromosome
                                band = read1_band

                        sort_read_with_mate(read=read1, fraction_file=intrachromosomal_file, chromosome=chromosome, band=band, spectrum_list=spectrum_list, patterns=patterns)
                        sort_read_with_mate(read=read2, fraction_file=intrachromosomal_file, chromosome=chromosome, band=band, spectrum_list=spectrum_list, patterns=patterns)  

                

                # INTRACHROMOSOMAL: both mapped
                else:

                        sort_read_with_mate(read=read1, fraction_file=intrachromosomal_file, chromosome=read1_chromosome, band=read1_band, spectrum_list=spectrum_list, patterns=patterns)
                        sort_read_with_mate(read=read2, fraction_file=intrachromosomal_file, chromosome=read2_chromosome, band=read2_band, spectrum_list=spectrum_list, patterns=patterns)      






                
        
         # while True:                     

         #        try:
         #                 read1 = bamfile.next()
         #        except Exception as e:
         #                print(f"Error while reading next BAM entry: {e}")
         #                break;
         #        try:
         #                read2 = bamfile.next()
         #        except: 
         #                sort_reads_without_mate (read1, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
         #                                        intratelomeric_file,
         #                                        subtelomeric_file,
         #                                        intrachromosomal_file )
         #                break

                

         #        ### reads without mates
         #        read1_first_flag = True # read 1 is first read

         #        while read1.qname != read2.qname:
                        
         #                if read1_first_flag == True:
                                
         #                        sort_reads_without_mate (read1, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
         #                                                intratelomeric_file,
         #                                                subtelomeric_file,
         #                                                intrachromosomal_file )

         #                        try:
         #                                read1 = bamfile.next()
         #                        except: 
         #                                sort_reads_without_mate (read2, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
         #                                                intratelomeric_file,
         #                                                subtelomeric_file,
         #                                                intrachromosomal_file )
         #                                break_flag = True
         #                                break                # break out of inner while loop
                                
         #                        read1_first_flag = False # read2 is now first read

                                
         #                else:
         #                        sort_reads_without_mate (read2, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
         #                                                intratelomeric_file,
         #                                                subtelomeric_file,
         #                                                intrachromosomal_file )
                                
         #                        try:
         #                                read2 = bamfile.next()
         #                        except: 
         #                                sort_reads_without_mate (read1, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
         #                                                intratelomeric_file,
         #                                                subtelomeric_file,
         #                                                intrachromosomal_file )
         #                                break_flag = True
         #                                break                 # break out of inner while loop
                                        
         #                        read1_first_flag = True

                


                        

        ###########################
        ### write spectrum file ###
        ###########################

        # open spectrum file for writing
        spectrum_file = open(out_dir + "/" + pid + ".spectrum", "w")

        # header
        spectrum_file.write( "chr\tband\treads_with_pattern")

        for pattern in patterns:
                spectrum_file.write( "\t" + pattern )
                
        spectrum_file.write( "\tother\n" )

        #write line for each chromosome band
        for chromosome in chromosome_list:

                for band in ["junction1"] + bands_list[chromosome]["band_name"] + ["junction2"]:
  
                        spectrum_file.write( "%s\t%s\t%i" % (chromosome, band, spectrum_list[chromosome][band]["reads_with_pattern"]))

                        for pattern in patterns:
                                spectrum_file.write("\t" + str(spectrum_list[chromosome][band][pattern]))

                        spectrum_file.write("\t" + str(int(round(spectrum_list[chromosome][band]["other"]))) + "\n" )


        spectrum_file.write( "unmapped\tunmapped\t%i" % (spectrum_list["unmapped"]["unmapped"]["reads_with_pattern"]))

        for pattern in patterns:
                spectrum_file.write("\t" + str(spectrum_list["unmapped"]["unmapped"][pattern]))

        spectrum_file.write("\t" + str(int(round(spectrum_list["unmapped"]["unmapped"]["other"]))) + "\n" )


        ##########################
        ### close file handles ###
        ##########################

        bamfile.close()
        intratelomeric_file.close()
        junctionspanning_file.close()
        subtelomeric_file.close()
        intrachromosomal_file.close()
        spectrum_file.close()



# function for checking reads without mate and sorting
def read_check (read, references, bands_list, chromosome_list, mapq_threshold):
  
        # get reference
        read_tid = read.tid
        read_ref = ''
        if read_tid != -1:
                read_ref = references[read_tid]
                read_ref = read_ref.replace("chr", "")
        
        # get position
        read_pos = read.pos
                
        # get mapping quality
        read_mapq = read.mapq

        
        read_junctionspanning = False
        read_junction = ""
        
        # check if read is considered unmapped
        if read.is_unmapped or read_ref not in chromosome_list or read.mapq < mapq_threshold:
                read_is_unmapped=True
                chromosome = "unmapped"
                band = "unmapped"
        else:
                read_is_unmapped=False
                
                chromosome = read_ref
                
                i=0
                while read_pos > bands_list[chromosome]["end"][i] and i<(len(bands_list[chromosome]["end"])-1):
                        i += 1
                else:
                        band = bands_list[chromosome]["band_name"][i]
                        
                        
                        # check if read is mapped to first or last chromosome band
                        if i==0:
                                read_junctionspanning = True
                                read_junction = "junction1"
                        elif i==(len(bands_list[chromosome]["end"])-1):
                                read_junctionspanning = True
                                read_junction = "junction2"

        return ( chromosome, band, read_is_unmapped, read_junctionspanning, read_junction )


# function for sorting reads without a mate into correct fraction and add counts to spectrum list
def sort_reads_without_mate (read, references, bands_list, chromosome_list, mapq_threshold, spectrum_list, patterns,
                             intratelomeric_file,
                             subtelomeric_file,
                             intrachromosomal_file ):

        # do read check
        chromosome, band, read_is_unmapped, read_junctionspanning, read_junction = read_check (read, references, bands_list, chromosome_list, mapq_threshold)

        # write read to correct fraction
        if read_is_unmapped:
                intratelomeric_file.write(read)
        elif read_junctionspanning:
                subtelomeric_file.write(read)
        else:
                intrachromosomal_file.write(read)

        # add read counts to spectrum list
        spectrum_list[chromosome][band]["reads_with_pattern"]+=1

        read_total_pattern_count = 0
        
        for pattern in patterns:
                spectrum_list[chromosome][band][pattern] += read.seq.count( pattern )
                read_total_pattern_count += read.seq.count( pattern )


        spectrum_list[chromosome][band]["other"] += float(len(read.seq))/6 - read_total_pattern_count


# function for writing reads with mates to fraction and adding counts to spectrum list
def sort_read_with_mate (read, fraction_file, chromosome, band, spectrum_list, patterns):

        fraction_file.write(read)
        
        spectrum_list[chromosome][band]["reads_with_pattern"]+=1

        read_total_pattern_count = 0

        #count occurrence of pattern types      
        for pattern in patterns:
                spectrum_list[chromosome][band][pattern] += read.seq.count( pattern )
                read_total_pattern_count += read.seq.count( pattern )

        spectrum_list[chromosome][band]["other"] += float(len(read.seq))/6 - read_total_pattern_count

