import csv
import sys
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
            temp = ''
        else:
            temp += line.replace('\n','') #remove whitespaces in the sequence

    seq[name] = temp
    seq.pop("placeholder") #remove the first key put into the dict as a placeholder                                                                    
    f.close()
    return seq


def extract_exon_seq(seq_dict,gff_file,gtf_file,output_file):
    exon_seq_list=[]
    headers_list=[]
    a_string="ACTGactg"
    conversion=a_string.maketrans("ACTGactg","TGACtgac")
    neg_dir_exons = []
    gene_name = ""
    if (gff_file == None and gtf_file == None) or (gff_file != None and gtf_file != None):
        sys.stderr.write("Please supply a reference annotation in GTF format (using —-gtf) or in GFF format (using —-gff)")
        return
    elif gff_file != None:
        with open(gff_file,'r') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                if len(row) != 9:
                    continue
                if row[2].endswith("transcript") or row[2].endswith("mRNA"):
                    if "gene=" in row[8]:
                        gene_name = row[8].split(sep="gene=")[1].split(sep=";")[0]
                    elif "gene_name" in row[8]:
                        gene_name = row[8].split(sep="gene_name=")[1].split(sep=";")[0]
                    elif "geneID" in row[8]:
                        gene_name = row[8].split(sep="geneID=")[1].split(sep=";")[0]
                    elif "gene_id" in row[8]:
                        gene_name = row[8].split(sep="gene_id=")[1].split(sep=";")[0]
                    else:
                        sys.stderr.write("INVALID GFF FORMAT? It is recommended to clean up your GFF file so that the third field is transcript,mRNA,or exon, and the ninth field of transcript or mRNA entries contain the gene info (in format of gene, gene_name or geneID).")
                        sys.stderr.write(row)
                        return

                if "exon" in row[2]: 
                    seq_name =  row[0]
                    if seq_name not in seq_dict.keys():
                        continue
                    direction = row[6]
                    start = int(row[3])
                    end=int(row[4])
                    full_seq=seq_dict[seq_name]
                    if "gene=" in row[8]:
                        gene_name = row[8].split(sep="gene=")[1].split(sep=";")[0]
                    elif "gene_name" in row[8]:
                        gene_name = row[8].split(sep="gene_name=")[1].split(sep=";")[0]
                    elif "geneID" in row[8]:
                        gene_name = row[8].split(sep="geneID=")[1].split(sep=";")[0]
                    elif "gene_id" in row[8]:
                        gene_name = row[8].split(sep="gene_id=")[1].split(sep=";")[0]
                    
                    exon_seq = full_seq[start-1:end]
                    orientation='F'
                    if (direction=='-'):
                        orientation = 'R'
                        exon_seq=exon_seq.translate(conversion)[::-1]
                    exon_seq_list.append(exon_seq)
                    headers_list.append("{}-{}-{}-{}-{}".format(seq_name,gene_name,orientation,start,end))
    else: #GTF
        with open(gtf_file,'r') as csvfile:
            reader = csv.reader(csvfile, delimiter='\t')
            for row in reader:
                if len(row) != 9:
                    continue
                if row[2].endswith("transcript") or row[2].endswith("mRNA"):
                    if "gene \"" in row[8]:
                        gene_name = row[8].split(sep="gene ")[1].split(sep=";")[0]
                    elif "gene_name" in row[8]:
                        gene_name = row[8].split(sep="gene_name ")[1].split(sep=";")[0]
                    elif "geneID" in row[8]:
                        gene_name = row[8].split(sep="geneID ")[1].split(sep=";")[0]
                    elif "gene_id" in row[8]:
                        gene_name = row[8].split(sep="gene_id ")[1].split(sep=";")[0]
                    else:
                        sys.stderr.write("INVALID GFF FORMAT? It is recommended to clean up your GFF file so that the third field is transcript,mRNA,or exon, and the ninth field of transcript or mRNA entries contain the gene info (in the format of gene, gene_name or geneID).")
                        sys.stderr.write(row)
                        return
                    gene_name =	gene_name.replace("\"","")
                    
                if "exon" in row[2]: 
                    seq_name =  row[0]
                    if seq_name not in seq_dict.keys():
                        continue
                    direction = row[6]
                    start = int(row[3])
                    end=int(row[4])
                    full_seq=seq_dict[seq_name]
                    if "gene \"" in row[8]:
                        gene_name = row[8].split(sep="gene ")[1].split(sep=";")[0]
                    elif "gene_name" in row[8]:
                        gene_name = row[8].split(sep="gene_name ")[1].split(sep=";")[0]
                    elif "geneID" in row[8]:
                        gene_name = row[8].split(sep="geneID ")[1].split(sep=";")[0]
                    elif "gene_id" in row[8]:
                        gene_name = row[8].split(sep="gene_id ")[1].split(sep=";")[0]
                    gene_name =	gene_name.replace("\"","")
                    
                    exon_seq = full_seq[start-1:end]
                    orientation='F'
                    if (direction=='-'):
                        orientation = 'R'
                        exon_seq=exon_seq.translate(conversion)[::-1]
                    exon_seq_list.append(exon_seq)
                    headers_list.append("{}-{}-{}-{}-{}".format(seq_name,gene_name,orientation,start,end))
       
        
                      
    header_exon_seq_dict = {headers_list[i]: exon_seq_list[i] for i in range(len(headers_list))}
    with open(output_file,'w') as of:
        for seqname in list(set(headers_list)):
            of.write(">{}\n".format(seqname))
            new_seq = split_output(header_exon_seq_dict[seqname],60)
            for l in new_seq:
                of.write(l+"\n")
    
    

def split_output(seq, num_per_line=60): #make a new line after num_per_line bases
    lines = math.ceil(len(seq)/num_per_line)
    output = []
    for i in range(lines):
        output.append(seq[num_per_line*i:num_per_line*(i+1)])
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--reads")
    parser.add_argument("--gff",default=None)
    parser.add_argument("--gtf",default=None)
    parser.add_argument("-o","--out",default="output.exons.fna")
    args = parser.parse_args()
    seq_dict=parse_fasta(args.reads)
    extract_exon_seq(seq_dict,args.gff,args.gtf,args.out)

if __name__ == '__main__':
    main()
