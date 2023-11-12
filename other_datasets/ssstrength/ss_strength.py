file1 = open("chr1.fa", "r")
Genome_sequence1 = ""
i1=0
for each_line in file1:
    n1=each_line.strip()
    if i1>0:
        Genome_sequence1 = Genome_sequence1 + n1
    i1=i1+1
file2 = open("chr2.fa", "r")
Genome_sequence2 = ""
i2=0
for each_line in file2:
    n2=each_line.strip()
    if i2>0:
        Genome_sequence2 = Genome_sequence2 + n2
    i2=i2+1
file3 = open("chr3.fa", "r")
Genome_sequence3 = ""
i3=0
for each_line in file3:
    n3=each_line.strip()
    if i3>0:
        Genome_sequence3 = Genome_sequence3 + n3
    i3=i3+1
file4 = open("chr4.fa", "r")
Genome_sequence4 = ""
i4=0
for each_line in file4:
    n4=each_line.strip()
    if i4>0:
        Genome_sequence4 = Genome_sequence4 + n4
    i4=i4+1
file5 = open("chr5.fa", "r")
Genome_sequence5 = ""
i5=0
for each_line in file5:
    n5=each_line.strip()
    if i5>0:
        Genome_sequence5 = Genome_sequence5 + n5
    i5=i5+1
file6 = open("chr6.fa", "r")
Genome_sequence6 = ""
i6=0
for each_line in file6:
    n6=each_line.strip()
    if i6>0:
        Genome_sequence6 = Genome_sequence6 + n6
    i6=i6+1
file7 = open("chr7.fa", "r")
Genome_sequence7 = ""
i7=0
for each_line in file7:
    n7=each_line.strip()
    if i7>0:
        Genome_sequence7 = Genome_sequence7 + n7
    i7=i7+1
file8 = open("chr8.fa", "r")
Genome_sequence8 = ""
i8=0
for each_line in file8:
    n8=each_line.strip()
    if i8>0:
        Genome_sequence8 = Genome_sequence8 + n8
    i8=i8+1
file9 = open("chr9.fa", "r")
Genome_sequence9 = ""
i9=0
for each_line in file9:
    n9=each_line.strip()
    if i9>0:
        Genome_sequence9 = Genome_sequence9 + n9
    i9=i9+1
file10 = open("chr10.fa", "r")
Genome_sequence10 = ""
i10=0
for each_line in file10:
    n10=each_line.strip()
    if i10>0:
        Genome_sequence10 = Genome_sequence10 + n10
    i10=i10+1
file11 = open("chr11.fa", "r")
Genome_sequence11 = ""
i11=0
for each_line in file11:
    n11=each_line.strip()
    if i11>0:
        Genome_sequence11 = Genome_sequence11 + n11
    i11=i11+1
file12 = open("chr12.fa", "r")
Genome_sequence12 = ""
i12=0
for each_line in file12:
    n12=each_line.strip()
    if i12>0:
        Genome_sequence12 = Genome_sequence12 + n12
    i12=i12+1
file13 = open("chr13.fa", "r")
Genome_sequence13 = ""
i13=0
for each_line in file13:
    n13=each_line.strip()
    if i13>0:
        Genome_sequence13 = Genome_sequence13 + n13
    i13=i13+1
file14 = open("chr14.fa", "r")
Genome_sequence14 = ""
i14=0
for each_line in file14:
    n14=each_line.strip()
    if i14>0:
        Genome_sequence14 = Genome_sequence14 + n14
    i14=i14+1
file15 = open("chr15.fa", "r")
Genome_sequence15 = ""
i15=0
for each_line in file15:
    n15=each_line.strip()
    if i15>0:
        Genome_sequence15 = Genome_sequence15 + n15
    i15=i15+1
file16 = open("chr16.fa", "r")
Genome_sequence16 = ""
i16=0
for each_line in file16:
    n16=each_line.strip()
    if i16>0:
        Genome_sequence16 = Genome_sequence16 + n16
    i16=i16+1
file17 = open("chr17.fa", "r")
Genome_sequence17 = ""
i17=0
for each_line in file17:
    n17=each_line.strip()
    if i17>0:
        Genome_sequence17 = Genome_sequence17 + n17
    i17=i17+1
file18 = open("chr18.fa", "r")
Genome_sequence18 = ""
i18=0
for each_line in file18:
    n18=each_line.strip()
    if i18>0:
        Genome_sequence18 = Genome_sequence18 + n18
    i18=i18+1
file19 = open("chr19.fa", "r")
Genome_sequence19 = ""
i19=0
for each_line in file19:
    n19=each_line.strip()
    if i19>0:
        Genome_sequence19 = Genome_sequence19 + n19
    i19=i19+1
file20 = open("chr20.fa", "r")
Genome_sequence20 = ""
i20=0
for each_line in file20:
    n20=each_line.strip()
    if i20>0:
        Genome_sequence20 = Genome_sequence20 + n20
    i20=i20+1
file21 = open("chr21.fa", "r")
Genome_sequence21 = ""
i21=0
for each_line in file21:
    n21=each_line.strip()
    if i21>0:
        Genome_sequence21 = Genome_sequence21 + n21
    i21=i21+1
file22 = open("chr22.fa", "r")
Genome_sequence22 = ""
i22=0
for each_line in file22:
    n22=each_line.strip()
    if i22>0:
        Genome_sequence22 = Genome_sequence22 + n22
    i22=i22+1
fileX = open("chrX.fa", "r")
Genome_sequenceX = ""
iX=0
for each_line in fileX:
    nX=each_line.strip()
    if iX>0:
        Genome_sequenceX = Genome_sequenceX + nX
    iX=iX+1
fileY = open("chrY.fa", "r")
Genome_sequenceY = ""
iY=0
for each_line in fileY:
    nY=each_line.strip()
    if iY>0:
        Genome_sequenceY = Genome_sequenceY + nY
    iY=iY+1

output1 = open("human_5ss.fasta", "w")
output2 = open("human_3ss.fasta", "w")
coor_file = open("coordinate.txt", "r")
trns_id = []
chr_num = []
fivess = []
threess = []
for each_line in coor_file:
    m = each_line.split()
    trns_id.append(m[0])
    chr_num.append(m[1])
    fivess.append(m[2])
    threess.append(m[3])

for j in range(0, len(trns_id)):
    if chr_num[j] == "chr1":
        sequence5 = Genome_sequence1[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence1[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr2":
        sequence5 = Genome_sequence2[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence2[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr3":
        sequence5 = Genome_sequence3[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence3[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr4":
        sequence5 = Genome_sequence4[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence4[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr5":
        sequence5 = Genome_sequence5[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence5[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr6":
        sequence5 = Genome_sequence6[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence6[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr7":
        sequence5 = Genome_sequence7[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence7[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr8":
        sequence5 = Genome_sequence8[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence8[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr9":
        sequence5 = Genome_sequence9[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence9[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr10":
        sequence5 = Genome_sequence10[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence10[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr11":
        sequence5 = Genome_sequence11[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence11[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr12":
        sequence5 = Genome_sequence12[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence12[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr13":
        sequence5 = Genome_sequence13[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence13[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr14":
        sequence5 = Genome_sequence14[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence14[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr15":
        sequence5 = Genome_sequence15[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence15[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr16":
        sequence5 = Genome_sequence16[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence16[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr17":
        sequence5 = Genome_sequence17[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence17[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr18":
        sequence5 = Genome_sequence18[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence18[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr19":
        sequence5 = Genome_sequence19[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence19[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr20":
        sequence5 = Genome_sequence20[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence20[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr21":
        sequence5 = Genome_sequence21[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence21[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chr22":
        sequence5 = Genome_sequence22[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequence22[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chrX":
        sequence5 = Genome_sequenceX[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequenceX[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
    elif chr_num[j] == "chrY":
        sequence5 = Genome_sequenceY[int(fivess[j])-4:int(fivess[j])+5]
        output1.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output1.write(sequence5 + "\n")
        sequence3 = Genome_sequenceY[int(fivess[j])-21:int(fivess[j])+2]
        output2.write(">" + str(trns_id[j]) + "_" + str(chr_num[j]) + "\n")
        output2.write(sequence3 + "\n")
        
        
        
