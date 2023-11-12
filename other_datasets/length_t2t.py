from Bio import SeqIO
in_file = open(input("Enter the name of the fasta file: "), "r")
out_file = open(input("Enter the name of output file: "), "w")

sequences = [seq for seq in SeqIO.parse(in_file, "fasta")]
for each_line in sequences:
    gene_id = each_line.id
    length = len(each_line.seq)
    Id = gene_id.split('_')
    out_file.write(Id[6] + "_" + Id[7] + "\t" + str(length) + "\n")

