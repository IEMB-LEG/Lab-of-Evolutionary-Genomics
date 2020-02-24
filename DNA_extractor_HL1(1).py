#! /bin/python
										###############################################
										#          Scripted by Hongan Long 6/10/2013  #
										#            Indiana University               #
										#           longhongan@gmail.com              #
										#   To extract target DNA sequences using     #
										#          start and end coordinates          #
										###############################################
# module load biopythn first
# argv1 is the genome fasta file, scaffold name is the fasta header of the scaffold/chromosome, pos_start/end are the start and end cooridinates of he sequences to be extracted
import sys
from Bio import SeqIO
from sys import argv
in_file=argv[1] #genome file
scaffold_name=argv[2] # scaffold/chromosome name in the  fasta header, without ">" 
pos_start=int(argv[3]) # start of the target sequence
pos_end=int(argv[4]) # end of the target sequence, start and end coordinates can be identical to get a single nucleotide position
handle_in=open(in_file,"rU") 
Extracted_seq=[]                                # create blank list to contain the 1.5kb sequences to be extracted
Seq_ref=SeqIO.parse(handle_in,"fasta")          # parse the genome fasta file
seq_sp=[]
for record in Seq_ref:				# iterate over each scaffold of the reference genome
	if record.id==scaffold_name:		# if the reference record has the same name as the target scaffold/chromosome
		print ">"+scaffold_name+":"+str(pos_start)+":"+str(pos_end)	# append the scaffold name into the final sequence file
		seq_sp=list(str(record.seq))		# turn the sequence of the scaffold into a list containing single-letter elements
		Extracted_seq.append('%s\n' % ''.join(seq_sp[pos_start-1:pos_end]))
for i in Extracted_seq:
		print i
handle_in.close()

