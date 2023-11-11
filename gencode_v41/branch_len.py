import re
from Bio import SeqIO

in_file1 = open("branchpoint.txt", "r")
out_file1 = open("branch_count_length.txt","w")
list0 = []
list1 = []

for each_line in in_file1:
    n = each_line.split()
    list0.append(n[0])
    list1.append(n[1])

for i in range(0, len(list1)):
    ssr = len(re.findall("((T|A)T(T|A|G|C)A(T|A))", list1[i]))
    if ssr!=0:
        count = 0
        len_list = []
        for match in re.finditer(r'((T|A)T(T|A|G|C)A(T|A))', list1[i]):
            count += 1
            len_list.append(int(len(list1[i]) - match.end()))
        out_file1.write(str(list0[i]) + "\t" + str(ssr) + "\t" + str(min(len_list)) + "\n")

