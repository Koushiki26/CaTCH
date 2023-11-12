import re
from Bio import SeqIO
in_file = open("intron_seq.fa", "r")
out_file = open("branchpoint.txt","w")

sequences = [seq for seq in SeqIO.parse(in_file, "fasta")]
for each_line in sequences:
    gene_id = each_line.id
    length = len(each_line.seq)
    string = each_line.seq
    Id = gene_id.split('_')
    out_file.write(str(Id[6]) + "_" + str(Id[7]) + "\t" + str(string) + "\t" + str(length) + "\n")


