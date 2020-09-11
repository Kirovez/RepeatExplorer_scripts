
from Bio import SeqIO
import random
"""
Script to select unmapped reads after mapping the parired reads to the AcCent1K and AfCen1K sequences
This script will return two fasta files and table with position of mapped reads on corresponding sequences
"""
def _getMfromCigar(cigar):
    sp = cigar.split("M")
    total_M = 0
    for i,words in enumerate(sp):
        if i != len(sp) - 1:
            M_num = ''
            for i in words[::-1]:
                try:
                    int(i)
                    M_num = i + M_num
                except ValueError:
                    break
            total_M += int(M_num)
    return total_M

def _getPosition(sam):
    d = {} #read id:[pos1,pos1+50]
    for lines in sam:
        sp = lines.split("\t")
        #if sp[2] == target_chromosome:
        if not lines.startswith("@"):
            if _getMfromCigar(sp[5]) > 100: # or (int(sp[3]) < 300 or int(sp[3]) > 900):
                d[sp[0]] = 0
    return d

def _intersectGetTable(sam1_mapped, sam2_mapped):

    """
    :param sam1_mapped: dictionary with read ids from fq1 that were mapped to the target
    :param sam2_mapped: dictionary with read ids from fq2 that were mapped to the target
    :return: dictionary with two keys "sam1" and "sam2" containing the read id for which second mate was not mapped
    """
    cnt = 0
    cnt_neg = 0
    dic = {'sam1':{}, 'sam2':{}} # ['sam1': 'read id': ['pos1 /t pos1+50', 'pos2 /t pos2+50']
    for hits in sam1_mapped:
        if hits in sam2_mapped:
            cnt += 1
        else:
            cnt_neg += 1
            dic["sam1"][hits] = 0

    for hits in sam2_mapped:
        if hits not in sam1_mapped:
            cnt_neg += 1
            dic["sam2"][hits] = 0



    print("Read in pairs", cnt)
    print("Not in pairs from file 1", cnt_neg)

    return dic

def getReads(read_id_dic, fq):
    "dictionary with readids and fq file whereread have to be selected from"
    dic_in_file = {} # read ids that were wrote in fq file
    with open(fq + "_unmapped_to_target.fq", "w") as outFile:
        cnt = 0
        for seq in SeqIO.parse(fq, "fastq"):
            if seq.id in read_id_dic:
                SeqIO.write(seq, outFile, "fastq")
                dic_in_file[seq.id] = 0
                cnt += 1

        print("Number of sequences selected from {} is".format(fq), cnt)
        return dic_in_file

def getUnmappedReads(sam1, sam2, fq1, fq2):
    sam1 = open(sam1)
    sam2 = open(sam2)
    sam1_mapped = _getPosition(sam1)
    sam2_mapped = _getPosition(sam2)
    sam1.close()
    sam2.close()
    unmapped_reads = _intersectGetTable(sam1_mapped, sam2_mapped)
    written_reads_sam1 = getReads(unmapped_reads["sam1"] , fq2)
    written_reads_sam2 = getReads(unmapped_reads["sam2"] , fq1)
    #
    # ##just write the position of mapped read on target
    # with open("read_positions_mapped_to_target.tab", "w") as outTab:
    #     for dics in [sam1_mapped, sam2_mapped]:
    #         for reads in dics:
    #
    #             outTab.write("\t".join([reads, str(dics[reads][0]),str(dics[reads][1])]) + "\n")

# getUnmappedReads("Afistulosum_vs_AF_centromere_RE_contigs_1.sam",
#                  "Afistulosum_vs_AF_centromere_RE_contigs_2.sam",
#                  r"/home/lik/raw_reads_Allium_BGI/Afistulosum/FCH3YY7DSXX_L1_wHAXPI081155-83_1.fq",
#                   r"/home/lik/raw_reads_Allium_BGI/Afistulosum/FCH3YY7DSXX_L1_wHAXPI081155-83_2.fq")

getUnmappedReads("Acepa_vs_AC_centromere_RE_contigs_1.sam",
                 "Acepa_vs_AC_centromere_RE_contigs_2.sam",
                 r"/home/lik/raw_reads_Allium_BGI/Acepa/FCH3YY7DSXX_L1_wHAXPI081197-82_1.fq",
                  r"/home/lik/raw_reads_Allium_BGI/Acepa/FCH3YY7DSXX_L1_wHAXPI081197-82_2.fq")