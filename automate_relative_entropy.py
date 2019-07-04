import os
import numpy as np
genomes=os.listdir("/home/fauzul/Videos/project_zika_virus/entropy_zkv_any_nucl")
len(genomes)
number_genomes=len(genomes)
os.chdir("/home/fauzul/Videos/project_zika_virus/entropy_zkv_any_nucl")
from Bio import SeqIO
      # This is the same as count = count + 1
for seq_record in SeqIO.parse(genomes[1], "fasta"):
	print(seq_record.seq.count('A'))
	print(seq_record.seq.count('T'))
	print(seq_record.seq.count('G'))
        print(seq_record.seq.count('C'))
	nucleotide=['A','T','G','C']
	genome_length=len(seq_record)
	codons=['ATT', 'ATC', 'ATA','CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG','GTT', 'GTC', 'GTA', 'GTG','TTT', 'TTC','ATG','TGT', 'TGC','GCT', 'GCC', 'GCA', 'GCG','GGT', 'GGC', 'GGA', 'GGG','CCT', 'CCC', 'CCA', 'CCG','ACT', 'ACC', 'ACA', 'ACG','TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC','TAT', 'TAC','TGG','CAA', 'CAG','AAT', 'AAC','CAT', 'CAC','GAA', 'GAG','GAT', 'GAC','AAA', 'AAG','CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG','TAA', 'TAG', 'TGA']
	print(nucleotide)
	print(codons)
	#calculate frequency of nucleotide in fasta sequence#
	
	

nucleotide_frequency = []
for i in range(len(nucleotide)):
	freq_n=seq_record.seq.count(nucleotide[i])/float(genome_length)
	nucleotide_frequency.append(freq_n)
	#calculate frequency of codons in fasta sequence#
	


codons_frequency = []
for j in range(len(codons)):
	freq_c=seq_record.seq.count(codons[j])
	codons_frequency.append(freq_c)


	
total_codons=sum(codons_frequency)
Dkli = []
for k in range(len(codons)):
	a=seq_record.seq.count(list(codons[k])[0])/float(genome_length)
	b=seq_record.seq.count(list(codons[k])[1])/float(genome_length)
	c=seq_record.seq.count(list(codons[k])[2])/float(genome_length)
	d=codons_frequency[k]/float(total_codons)
	#print(a)
	#print(b)
	#print(c)
	#print(d)
	ent= (d*np.log2(d/(a*b*c)))
	Dkli.append(ent)
	#print (Dkli)

print (sum(Dkli))

