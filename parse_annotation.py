with open(r"D:\Allium\2018_Allium\READS_BGI\forRE_filtered_sampled\AC_AR_AF_RE2\CLUSTER_TABLE.csv") as infile, \
    open("annotation_AC_AR_AF.tab", "w") as outfile:

    for i,lines in enumerate(infile):
        if i!= 0:
            sp = lines.rstrip().split("\t")
            if len(sp) > 3:
                cluster = sp[0]
                annotation = sp[-1].replace('"','')
                if annotation:
                    print(annotation)
                    percent = annotation.split(' ')[0][:-1]
                    if float(percent)> 5.0:
                        annot = annotation.split(' ')[1]
                    else:
                        annot = "Unknown"
                else:
                    percent = "0"
                    annot = "Unknown"

                outfile.write("\t".join([cluster,annot]) + "\n")

