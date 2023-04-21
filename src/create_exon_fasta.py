import csv
import math
from os.path import exists
import random
import argparse

def parse_fasta(query_file):
    f = open(query_file,"r")
    temp = ''
    seq = {}
    name = "placeholder"
    for line in f:
        if line.startswith(">"):
            seq[name] = temp
            name = line.split()[0][1:]
            name = chr_conversion(name)
            temp = ''
        else:
            temp += line.replace('\n','') #remove whitespaces in the sequence

    seq[name] = temp
    seq.pop("placeholder") #remove the first key put into the dict as a placeholder                                                                    
    f.close()
    return seq

def chr_conversion(nc_name):
    if nc_name[:3] == "chr":
        return nc_name
    num=int(nc_name.split('.')[0][3:])
    if int(num) < 23:
        return "chr"+str(num)
    elif int(num)==23:
        return "chrX"
    elif int(num)==24:
        return "chrY"
    else:
        return "UnknownChr"

def extract_exon_seq(seq_dict,gff_file,output_file,neg_file):
    exon_seq_list=[]
    headers_list=[]
    a_string="ACTGactg"
    conversion=a_string.maketrans("ACTGactg","TGACtgac")
    neg_dir_exons = []
    first_gene_or_transcript = False
    with open(gff_file,'r') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for row in reader:
            if len(row) != 9:
                continue
            if ("Parent" not in row[8]) and ("gene" in row[8]) and (row[2] in "transcript mRNA gene"):
                first_gene_or_transcript = True
                direction = row[6]
                seq_name =  chr_conversion(row[0])
                if "gene=" in row[8]:
                    gene_name = row[8].split(sep="gene=")[1].split(sep=";")[0]
                elif "gene_name" in row[8]:
                    gene_name = row[8].split(sep="gene_name=")[1].split(sep=";")[0]
                elif "geneID" in row[8]: #The stringtie assembly uses geneID instead of gene_name, but in hg38c_protein_and_lncRNA.gff geneID is meaningless
                    gene_name = row[8].split(sep="geneID=")[1].split(sep=";")[0]
                else:
                    print("INVALID GFF FORMAT?")
                    return
            else:
                #exon rows
                if first_gene_or_transcript == False or seq_name not in seq_dict.keys(): #eg. chromosome M
                    continue
                start = int(row[3])
                end=int(row[4])
                full_seq=seq_dict[seq_name]
                exon_seq = full_seq[start-1:end]
                if (direction=='-'):
                    neg_dir_exons.append("{}_{}".format(seq_name,gene_name))
                    exon_seq=exon_seq.translate(conversion)[::-1]
                exon_seq_list.append(exon_seq)
                headers_list.append("{}_{}_{}_{}".format(seq_name,gene_name,start,end))
            #1-indexed
            
    header_exon_seq_dict = {headers_list[i]: exon_seq_list[i] for i in range(len(headers_list))}
    #random.shuffle(headers_list)
    with open(output_file,'w') as of:
        for seqname in list(set(headers_list)):
            of.write(">{}\n".format(seqname))
            new_seq = split_output(header_exon_seq_dict[seqname],60)
            for l in new_seq:
                of.write(l+"\n")
    
    
    neg_dir_exons=list(set(neg_dir_exons))
    with open(neg_file,'w') as of:
        of.write('\n'.join(neg_dir_exons))

def split_output(seq, num_per_line=60): #make a new line after num_per_line bases
    lines = math.ceil(len(seq)/num_per_line)
    output = []
    for i in range(lines):
        output.append(seq[num_per_line*i:num_per_line*(i+1)])
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--reads")
    parser.add_argument("-g","--gff")
    parser.add_argument("-o","--out")
    parser.add_argument("-n","--neg")
    args = parser.parse_args()
    seq_dict=parse_fasta(args.reads)
    extract_exon_seq(seq_dict,args.gff,args.out,args.neg)

if __name__ == '__main__':
    main()
