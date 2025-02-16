#!/usr/bin/env python

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



import os
#from PyPDF2 import PdfFileReader, PdfFileMerger

from PyPDF2 import PdfMerger  # Import PdfMerger instead of PdfFileMerger
from PyPDF2 import PdfReader  # PdfReader is also used in place of PdfFileReader

###############################################################################
### merge all PDF plots produced as TelomereHunter output into one PDF file ###
###############################################################################

def mergeTelomereHunterPDFs(pid, outdir):

    # all possible pdf file names in correct order for merged pdf
    chromosomes = [str(i) for i in range(1, 22 + 1)] + ["X", "Y"]
    possible_file_names = [pid + '_telomere_content', pid + '_summary', pid + '_hist_telomere_repeats_per_intratelomeric_read', pid + '_gc_content'] + [pid + "_" + i for i in chromosomes] 
    possible_pdf_names = [i + '.pdf' for i in possible_file_names]


    # check which of the possible pdf files exist and sort
    files_dir = os.path.join(outdir, pid, "plots")
    pdf_files = [f for f in os.listdir(files_dir) if f.endswith("pdf")]
    pdf_files_ordered = [pdf for pdf in possible_pdf_names if pdf in pdf_files]


    # merge files using PdfMerger
    merger = PdfMerger()  # Use PdfMerger instead of PdfFileMerger

    for filename in pdf_files_ordered:
        with open(os.path.join(files_dir, filename), "rb") as file:
            reader = PdfReader(file)  # Use PdfReader instead of PdfFileReader
            merger.append(reader)  # Append the PdfReader object

    # Write the merged PDF
    output_pdf = os.path.join(outdir, pid, pid + "_merged_plots.pdf")
    with open(output_pdf, "wb") as output_file:
        merger.write(output_file)

# def mergeTelomereHunterPDFs(pid, outdir):

#     # all possible pdf file names in correct order for merged pdf
#     chromosomes = [ str(i) for i in range(1,22+1) ] + ["X","Y"]
#     possible_file_names = [pid + '_telomere_content',pid + '_summary', pid + '_hist_telomere_repeats_per_intratelomeric_read', pid + '_gc_content'] + [ pid + "_" + i for i in chromosomes] 
#     possible_pdf_names = [ i + '.pdf' for i in possible_file_names ]


#     # check which of the possible pdf files exist and sort
#     files_dir = os.path.join(outdir, pid, "plots")
#     pdf_files = [f for f in os.listdir(files_dir) if f.endswith("pdf")]
#     pdf_files_ordered =  [pdf for pdf in possible_pdf_names if pdf in pdf_files]


#     # merge files
#     merger = PdfFileMerger()

#     for filename in pdf_files_ordered:
#         merger.append(PdfFileReader(os.path.join(files_dir, filename), "rb"))

#     merger.write(os.path.join(outdir, pid, pid + "_merged_plots.pdf"))


