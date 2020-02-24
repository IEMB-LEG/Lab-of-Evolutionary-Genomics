#! /bin/env python
# Only one file to convert, can be updated using loop
# module load biopython in shell first
from Bio import SeqIO
handle=open("/N/dc/projects/longh/pseduomonas/130320_I1135_FCD24LLACXX_L8_CCAPEI13010029_1.fq","rU")
outputhandle=open("/N/dc/projects/longh/pseduomonas/Megablast/rawfq1.fasta","w")
standardfna=[]
for rec in SeqIO.parse(handle,"fastq"):
	standardfna.append('>%s\n' % rec.id)
	standardfna.append('%s\n' % rec.seq.tostring())
for i in standardfna:
	outputhandle.write(i)
handle.close()
outputhandle.close()
