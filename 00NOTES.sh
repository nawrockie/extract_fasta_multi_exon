# EPN, Wed Mar  4 10:03:48 2015
# 
# Example of running idfetch, id_fasta.pl and extract_fa_multi_exon to
# get specific regions (potentially mult-exonic) of specific
# sequences.
#
# Based on:
# /panfs/pan1/dnaorg/programs/15_0213_extract_fasta_richa/
# which contains a 00NOTES.sh that derived from notes
# created by Alejandro Schaffer.
#
# The major difference between this directory and
# /panfs/pan1/dnaorg/programs/15_0213_extract_fasta_richa/ is that
# 'extract_fasta' has been renamed and modified here to
# 'extract_fasta_multi_exon'. The modifications are to allow fetching
# of multiple exon files. Also a test script that checks the
# error-checking and some normal run cases of this program are
# included in this dir too.
#
#######################
# More information
#######################
# See /home/nawrocke/notebook/15_0303_dnaorg_extract_fasta_multiple_exons/00LOG.txt
# for notes on development and testing of this program
# 
#######################
# Example command
#######################
# a single line command that runs all three programs:
idfetch -t 5 -c 1 -G sample.1.in | id_fasta.pl | extract_fasta_multi_exon sample.2.in > sample.out
#
######################
# Programs
# (idfetch and id_fasta.pl sections copied from ../programs/15_0213_extract_fasta_richa/00NOTES.sh)
######################
#
# idfetch
# type:     C program
# author:   ?
# location: system-wide; /netopt/ncbi_tools64/bin/idfetch
#
# Important flags for our purposes:
# idfetch: 
# '-t 5' produces FASTA
# '-c 1' is the complexity of the sequences ('get the bioseq of interest')
# '-G'   input file of accessions, GenBank IDs, one per line
# 
# For more info: 'idfetch -'
#
######################
#
# id_fasta.pl 
# type:     Perl script
# author:   ?
# location: local, in this directory
# 
# Renames fasta sequences so the IDs are Genbank IDs,
# modifying only the name (first token attached to the '>')
# and regurgitating everything else. Original name is 
# included as the first token in the sequence description,
# followed by the original description.
# 
# Usage: id_fasta.pl <fasta file, from 'idfetch'>
#
######################
#
# extract_fasta_multi_exon
# type:     C program
# location: local, in this directory
# author:   Richa Agarwala [extract_fasta]
#           Eric Nawrocki [modified original by Richa]
#
# Extracts specified intervals/strands from DNA sequences in a
# provided fasta file. Those intervals can be multi-exonic (this is
# the main difference between extract_fasta_multi_exon and
# extract_fasta).
#
# Usage:
# > ./extract_fasta_multi_exon
# extract_fasta_multi_exon <interval_list> [<fa_file>]
# Format of <interval_list>, each line must look like this:
#
# <accession/id> <num-pieces (n)> <start_1> <end_1> <start_2> <end_2> ... <start_n> <end_n> <strand>
#
# For all n <start> must be <= <end> and 
# <start_j> must be > <end_i> for all i < j.
#
# Finally, if a specific <accession/id> occurs multiple times, none of
# the intervals for that <accession/id> (defined as
# <start_1>..<end_n>) can overlap (by >= 1 nt) with any of the other
# intervals on either strand.
#
# Compilation:
#  > gcc -o extract_fasta_multi_exon extract_fasta_multi_exon.c
#
# Source of 'extract_fasta_multi_exon':
# /home/nawrocke/notebook/15_0303_dnaorg_extract_fasta_multiple_exons/
# (notes in the 00LOG.txt file in the above dir).
#
# Source of original 'extract_fasta':
# ~agarwala/programs/extract_fasta.c
#
######################
# Input files
#
# Note: sample.2.in would actually work as input for both 'idfetch'
#       and 'extract_fasta_multi_fetch' because 'idfetch' ignores all
#       but the first white-space delimited token on each line (at
#       least with the options/args we're using). But that won't
#       always be the case, because it's important that each sequence
#       name be listed exactly once in the idfetch input, even though
#       it may occur more than once in the extract_fasta_multi_exon
#       input.
# 
######################
## sample.1.in
# JOQW01000002.1
# KN275973.1
#
######################
## sample.2.in
# JOQW01000002.1 1 101803 101810 +
# KN275973.1 3 1 50 70 100 130 150 -
#
######################
# Output files
######################
#
## sample.out
#>JOQW01000002.1:101803_101810:+
#ATGACAGA
#>KN275973.1:1_50:70_100:130_150:-
#TTATAAATAAATTAAATATAAAAAGTTTATAAAATTATTAAATATAAAATATTAATTAACTATAATACTTTAAAATTGAA
#TTAATAATATATAGAAAAATCG
#
######################
# Test script for extract_fasta_multi_exon
######################
#
# test_extract_fasta_multi_exon.pl: 
# 
# Tests extract_fasta_multi_exon for about 12 error cases (in which
# the program is expected to fail) and a few normal cases where the
# program is not expected to fail. For the non-failures, output is
# checked against expected output.
#
# test-files/: this directory includes the input and output files
#              that test_extract_fasta_multi_exon.pl uses to test 
#              extract_fasta_multi_exon.
# 
# Usage: perl test_extract_fasta_multi_exon.pl <executable-dir> <testfile-dir> <'1' for verbose mode, '0' for quiet mode>
# 
# Running it (change '1' to '0' for quiet mode)
perl test_extract_fasta_multi_exon.pl . test-files 1
# 
# If all tests pass, final line of output is:
# PASS [11 tests failed as expected and 12 tests succeeded as expected]
#
# last updated [EPN, Mon Mar  9 09:37:55 2015]
