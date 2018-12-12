import random
from Bio import SeqIO


def readsmapling(read1, read2, percentage):
    dic = {}
    cnt = 0
    cnt2 = 0
    tmp = "sampling_1_tmp{}.fasta".format(percentage)
    with open(tmp, 'w') as output_file, \
    open ("sampling_{}.fasta".format(percentage), 'w') as output_file2:
        for seq in SeqIO.parse(read1,'fastq'):
            if random.randrange(1,101) <= percentage:
                # print(seq.id)
                cnt+=1
                dic[seq.id] = 0
                SeqIO.write(seq,output_file,'fasta')
        ind = SeqIO.index(tmp, "fasta")
        for seq in SeqIO.parse(read2,'fastq'):
            if seq.id in dic:
                seq1 = ind[seq.id]
                seq1.id += "/1"
                seq1.description = ""

                seq.id += "/2"
                seq.description = ""

                SeqIO.write(seq1,output_file2,'fasta')
                SeqIO.write(seq, output_file2, 'fasta')
                cnt2+=1

    print('Forward reads: ' + str(cnt))
    print('Reverse reads: ' + str(cnt2))

import sys

if __name__ == '__main__':
    readsmapling(sys.argv[1], sys.argv[2], float(sys.argv[3]))