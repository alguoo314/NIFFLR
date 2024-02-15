import argparse
import sys
from collections import OrderedDict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--assembled",default="sorted_good_output.gtf",help="Path to the sorted gtf file containing the assembled exon paths that have a l\
ow gap penalty")
    parser.add_argument("-r","--ref",default="reference.gff",help="The GFF file containing the reference")
    parser.add_argument("-o","--outp",default="read_count_added_reference.gff",help="Path to the output GFF file with the number of reads added to the referenc\
e transcripts")
    parser.add_argument("-c","--chr_names",default="chr_names.txt",help="A file containing the chromosome names, aka the first field of the reference gff file")
    args = parser.parse_args()
    global ref_file
    global assembled_file
    global pre_output_text
    global chr_names_orders
    global cov_info
    cov_info = True
    chr_names_orders = {}
    with open(args.chr_names, "r") as chrfile:
        i=0
        for line in chrfile:
            chrom = line.strip()
            chr_names_orders[chrom]=i
            i+=1

    ref_file = open(args.ref,'r')
    assembled_file = open(args.assembled,'r')
    output_file = open(args.outp,'w')
    output_text = ""
    pre_output_text = OrderedDict()
    ref_line = ref_file.readline()
    assembled_line= assembled_file.readline()
    
    while ref_line[0] == "#" or ref_line[0] == "":
        output_text +=ref_line
        ref_line = ref_file.readline()
    
    while assembled_line[0] == "#":
        assembled_line = assembled_file.readline()
    
    first_transcript = True
    assembled_line_fields = assembled_line.split('\t')
    ref_line_fields = ref_line.split('\t')
    assembled_exon_chain = ""
    global ref_exon_chain_record
    ref_exon_chain_record = {}
    global coverage_record
    coverage_record = {}
    global transcript_reads_record
    transcript_reads_record = {}
    global transcript_support_record
    transcript_support_record = {}
    global intron_match_record
    intron_match_record = {}
    num_reads = 0
    transcript_len = 0
    max_read_len = 0
    
    output_file.write(output_text)
    while ref_line and assembled_line:
        if ref_line[0] == "#":
            continue
        if assembled_line_fields[2] == "transcript":
            chr_name,coord_starts,coord_ends,assembled_exon_chain,first_exon_end,num_reads,assembled_line_fields,assembled_line,single_exon,max_read_len,num_junc = process_assembled_transcript(assembled_line_fields)
            #now the pointer is at the next transcript of the assembled gtf file
            if first_transcript == False:
                calc_read_proportions(assembled_exon_chain,num_reads,single_exon,max_read_len,transcript_len,max(coord_starts,coord_ends),num_junc)
            
            prev_ref_line = ref_line
            
            ref_line_fields,ref_line,transcript_len = process_ref_transcript(ref_line,ref_line_fields,coord_starts,first_exon_end,chr_name,transcript_len)
            
            first_transcript = False
            if prev_ref_line != ref_line: #if the pointer for ref file moved
                calc_read_proportions(assembled_exon_chain,num_reads,single_exon,max_read_len,transcript_len,max(coord_starts,coord_ends),num_junc)
                
            #now the pointer is at the next transcript of the ref gff file
        
        if len(pre_output_text) > 2000: #keep the size of the dict small
            it = 0
            while it < 1900:
                k,v =pre_output_text.popitem(last=False)
                intron_matched_frac=intron_match_record.pop(k,[0,0])
                output_file.write(v[0].strip('\n')+";read_num={};transcript_support={};covered_junctions={}/{}\n".format(round(transcript_reads_record.pop(k,0),3),round(transcript_support_record.pop(k,[0,0,0])[2],3),intron_matched_frac[0],intron_matched_frac[1]))
                output_file.write(''.join(v[1:]))
                it +=1
                
          
                
                
      #lastly
    if first_transcript == False:
        for k,v in pre_output_text.items():
            intron_matched_frac=intron_match_record.pop(k,[0,0])
            output_file.write(v[0].strip('\n')+";read_num={};transcript_support={};covered_junctions={}/{}\n".format(round(transcript_reads_record.pop(k,0),3),round(transcript_support_record.pop(k,[0,0,0])[2],3),intron_matched_frac[0],intron_matched_frac[1]))
            output_file.write(''.join(v[1:]))

        #if there are any ref transcripts remaining with no reads
        last_buffer = []
        exon_num = 0
        for line in ref_file.read().splitlines():
            if line[0] == "#":
                continue
            if line.split('\t')[2] == "transcript":
                if exon_num > 0 and len(last_buffer)>0:
                    output_file.write(last_buffer[0]+str(exon_num-1)+'\n'+last_buffer[1])
                last_buffer = [line+";read_num=0;transcript_support=0;covered_junctions=0/",'']
                exon_num = 0
            else:
                if len(last_buffer)>0:
                    last_buffer[1]+=str(line+'\n')
                if line.split('\t')[2] == "exon":
                    exon_num +=1

        if exon_num > 0 and len(last_buffer) > 0:
            output_file.write(last_buffer[0]+str(exon_num-1)+'\n'+last_buffer[1])

        
    ref_file.close()
    assembled_file.close()
    output_file.close()
        
        
 
    
 
def process_assembled_transcript(assembled_line_fields):
    global ref_file
    global assembled_file
    single_exon = False
    chr_name = assembled_line_fields[0]
    coord_starts = int(assembled_line_fields[3])
    coord_ends = int(assembled_line_fields[4])
    assembled_exon_chain = ""
    num_reads = len((assembled_line_fields[8].split("source_reads=")[1]).split(","))
    max_read_len = int(assembled_line_fields[8].split("longest_mapped_read_len=")[1])
    assembled_line = assembled_file.readline()
    if assembled_line !="":
        assembled_line_fields = assembled_line.split('\t')

    exon_nums = 0    
    while assembled_line_fields[2] == "exon":
        exon_nums+=1
        assembled_exon_chain+=str(assembled_line_fields[3])+"-"+str(assembled_line_fields[4])+"-"
        assembled_line = assembled_file.readline()
        if assembled_line != "":
            assembled_line_fields = assembled_line.split('\t')
        else:
            break
    first_exon_end = int(assembled_exon_chain.split("-")[1])
    intron_chain = assembled_exon_chain.split("-")[1:-2]
    if len(intron_chain) > 0:
        assembled_exon_chain = "-".join(intron_chain)
    else:
        assembled_exon_chain=assembled_exon_chain.split("-")
        single_exon = True

    num_junc=exon_nums-1
    return chr_name,coord_starts,coord_ends,assembled_exon_chain,first_exon_end,num_reads,assembled_line_fields,assembled_line,single_exon,max_read_len,num_junc



def process_ref_transcript(ref_line,ref_line_fields,coord_starts,first_exon_end,chr_name,transcript_len):
    global output_text
    global coverage_record
    global pre_output_text
    global ref_file
    global ref_exon_chain_record
    global intron_match_record
    global cov_info
    transcript_id = None
    if  chr_num_conversion(ref_line_fields[0]) > chr_num_conversion(chr_name):
         return ref_line_fields,ref_line,transcript_len

    exon_num=0
    while (ref_line_fields[2] != "transcript" and ref_line_fields[2] != "mRNA") or int(ref_line_fields[4]) < coord_starts or chr_num_conversion(ref_line_fields[0]) < chr_num_conversion(chr_name):
        if ref_line_fields[2] != "transcript" and ref_line_fields[2] != "mRNA":
            if transcript_id != None:
                pre_output_text[transcript_id].append(ref_line)
            if ref_line_fields[2]=="exon":
                exon_num +=1
        else:
            if exon_num > 0:
                intron_match_record[transcript_id] = [0,exon_num-1]
                exon_num = 0
            transcript_id = ref_line_fields[8].split(";")[0]
            pre_output_text[transcript_id]=[ref_line]
            
        ref_line = ref_file.readline()
        if ref_line != "":
            ref_line_fields = ref_line.split('\t')
        else:
            return ref_line_fields,ref_line,transcript_len
        #linear scan
    if exon_num > 0:
        intron_match_record[transcript_id] = [0,exon_num-1]
    while int(ref_line_fields[3]) <  first_exon_end and ref_line_fields[0] == chr_name:
        if ref_line_fields[2] == "transcript" or ref_line_fields[2] == "mRNA":
            field_8 = ref_line_fields[8]
            transcript_id = field_8.split(";")[0]
            if cov_info == False:
                coverage = 1.0
            elif cov_info== True and "cov=" not in field_8:
                coverage = 1.0
                cov_info = False
                sys.stderr.write("Transcript Coverage is not listed as 'cov=' in the attribute column of the GFF file. Assuming equal coverage for all transcripts.\n")
            else:
                coverage = float((field_8.split("cov=")[1]).split(";")[0])
            coverage_record[transcript_id] = coverage
            transcript_len = 0
            pre_output_text[transcript_id]=[ref_line.strip('\n')]
            ref_line = ref_file.readline()
            ref_line_fields = ref_line.split('\t')
            ref_exon_chain = ""
        
        exon_num = 0
        while ref_line_fields[2] == "exon":
            exon_num +=1
            pre_output_text[transcript_id].append(ref_line)
            transcript_len += int(ref_line_fields[4])-int(ref_line_fields[3])+1
            ref_exon_chain+=str(ref_line_fields[3])+"-"+str(ref_line_fields[4])+"-"
            ref_line = ref_file.readline()
            if ref_line !="":
                ref_line_fields = ref_line.split('\t')
            else:
                 break
        
        #now we have reached a new transcript
            
        if len(ref_exon_chain_record) >0 and int(ref_exon_chain_record[next(reversed(ref_exon_chain_record))].split("-")[-2]) < int(ref_exon_chain.split("-")[0]):
            #reduce the total size of ref_exon_chain_record
            #is this too time consuming?
            ref_exon_chain_record = {}
            
        ref_exon_chain_record[transcript_id] = ref_exon_chain
        
        intron_match_record[transcript_id] = [0,exon_num-1]
   
    return ref_line_fields,ref_line,transcript_len

 

def calc_read_proportions(assembled_exon_chain,num_reads,single_exon,max_read_len,transcript_len,end_pos,num_junc): #end_pos is the position of the 3' end of the transcript
    global transcript_reads_record
    global transcript_support_record
    global intron_match_record
    coverages = []
    transcript_ids = []
    
    if single_exon:
        left = assembled_exon_chain[0]
        right = assembled_exon_chain[1]
        for k in ref_exon_chain_record.keys():
            if '-'+left in ref_exon_chain_record[k] or right+'-' in ref_exon_chain_record[k]:
                coverages.append(coverage_record[k])
                transcript_ids.append(k)
                transcript_support_info = transcript_support_record.get(k,[0,0,0])
                if num_junc > transcript_support_info[0] or (num_junc == transcript_support_info[0] and end_pos > transcript_support_info[1]) or (num_junc == transcript_support_info[0] and end_pos == transcript_support_info[1] and max_read_len/transcript_len > transcript_support_info[2]):
                    transcript_support_record[k]=[num_junc,end_pos,max_read_len/transcript_len]
                
                intron_match_record[k][0]=max(0,intron_match_record[k][0])
                
    else:
        for k in ref_exon_chain_record.keys():
            if assembled_exon_chain in ref_exon_chain_record[k]:
                coverages.append(coverage_record[k])
                transcript_ids.append(k)
                transcript_support_info = transcript_support_record.get(k,[0,0,0])
                if num_junc > transcript_support_info[0] or (num_junc == transcript_support_info[0] and end_pos > transcript_support_info[1]) or (num_junc == transcript_support_info[0] and end_pos == transcript_support_info[1] and max_read_len/transcript_len > transcript_support_info[2]):
                    transcript_support_record[k]=[num_junc,end_pos,max_read_len/transcript_len]
                intron_match_record[k][0] = max(intron_match_record[k][0],num_junc) #this is the number of junctions in the assembled transcript
                
                
    coverage_proportions = [num_reads*x / sum(coverages) for x in coverages]
    
    
    for i in range(len(transcript_ids)):
        transcript_id = transcript_ids[i]
        if transcript_id in transcript_reads_record.keys():
            transcript_reads_record[transcript_id] += coverage_proportions[i]
        else:
            transcript_reads_record[transcript_id] =coverage_proportions[i]
    return


def chr_num_conversion(chro):
    return chr_names_orders[chro]
    
if __name__ == '__main__':
    main()
