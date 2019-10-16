"""
This script takes two fastq files and prepare final fasta file for RE run. The following commands can be carry out
1. fastq2fasta
2. random pairs
3. add prefix
4. rename paired reads (_f and _r are added)


# final file with INTERLACED reads looks like:
>Species1_read1_f
....
>Species1_read1_r

"""
import sys
import random
from Bio import SeqIO
class RERP():
    def __init__(self, fq1, fq2, species_name, random=0):
        self.fq1 = fq1
        self.fq2 = fq2
        self.species = species_name
        self.random = random
        self.OUT = self.species + "prepared_forRE.fasta"
        self.main()

    def writeFile(self, inFileFq, outFileFq, ids_dic):
        with open(outFileFq, 'w') as outFile:
            for seq in SeqIO.parse(inFileFq, "fastq"):
                if seq.id in ids_dic:
                    SeqIO.write(seq, outFile, 'fastq')
        print('  ', outFileFq,'has been created.')

    def getRandomReadNames(self):
        random_list_ids = {i:0 for i in random.sample([i.id for i in SeqIO.parse(self.fq1, "fastq")], self.random)}
        print("Random list of {} ids have been generated".format(len(random_list_ids)))

        new_name_fq1 = self.fq1 + "_random{}.fq".format(self.random)
        new_name_fq2 = self.fq2 + "_random{}.fq".format(self.random)

        self.writeFile(self.fq1, new_name_fq1, random_list_ids)
        self.writeFile(self.fq2, new_name_fq2, random_list_ids)

        self.fq1 = new_name_fq1
        self.fq2 = new_name_fq2

    def main(self):
        if self.random:
            self.getRandomReadNames()
        ind = SeqIO.index(self.fq2, "fastq")
        with open(self.OUT, "w") as outFile:
            for i, seq in enumerate(SeqIO.parse(self.fq1, "fastq")):
                seq1 = seq
                seq2 = ind[seq.id]
                seq1.id = self.species + "read{}".format(i)+ "_f"
                seq1.description = ""
                seq2.id = self.species + "read{}".format(i) + "_r"
                seq2.description = ""
                SeqIO.write(seq1, outFile, 'fasta')
                SeqIO.write(seq2, outFile, 'fasta')

            print("Success! File {0} has been created. It contains {1} interlaced reads with prefix {2}".format(self.OUT, i+1, self.species))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='A tool to create a single file for RE usage from two fq files')

    parser.add_argument("--randomN", default=0, help="The number of reads to be randomly selected", type=int)
    parser.add_argument("fq1",  help="path to fastq file 1")
    parser.add_argument("fq2",  help="path to fastq file 2")
    parser.add_argument("prefix", help="prefix top be added to the read names")

    args = parser.parse_args()

    RERP(args.fq1, args.fq2, args.prefix, random=args.randomN)