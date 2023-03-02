import argparse
from collections import OrderedDict

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--assembled",default="sorted_good_output.gtf",help="Path to the sorted gtf file containing the assembled exon paths that have a l\
ow gap penalty")
    parser.add_argument("-r","--ref",default="reference.gff",help="The GFF file containing the reference")
    parser.add_argument("-o","--outp",default="read_count_added_reference.gff",help="Path to the output GFF file with the number of reads added to the referenc\
e transcripts")
    args = parser.parse_args()
    global ref_file
    global assembled_file
    global pre_output_text
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
    num_reads = 0
    
    output_file.write(output_text)
    while ref_line and assembled_line:
        if assembled_line_fields[2] == "transcript":
            chr_name,coord_starts,coord_ends,assembled_exon_chain,first_exon_end,num_reads,assembled_line_fields,assembled_line,single_exon = process_assembled_transcript(assembled_line_fields)
            #now the pointer is at the next transcript of the assembled gtf file
            if first_transcript == False:
                calc_read_proportions(assembled_exon_chain,num_reads,single_exon)
            
            prev_ref_line = ref_line
            
            ref_line_fields,ref_line = process_ref_transcript(ref_line,ref_line_fields,coord_starts,first_exon_end,chr_name)
            
            first_transcript = False
            if prev_ref_line != ref_line: #if the pointer for ref file moved
                calc_read_proportions(assembled_exon_chain,num_reads,single_exon)
                
            #now the pointer is at the next transcript of the ref gff file
        
        if len(pre_output_text) > 2000: #keep the size of the dict small
            it = 0
            while it < 1900:
                k,v =pre_output_text.popitem(last=False)
                output_file.write(v[0].strip('\n')+";read_num="+str(round(transcript_reads_record.pop(k,0),3))+'\n')
                output_file.write(''.join(v[1:]))
                it +=1
                
          
                
                
      #lastly
    if first_transcript == False:
        for k,v in pre_output_text.items():
            output_file.write(v[0].strip('\n')+";read_num="+str(round(transcript_reads_record.pop(k,0),3))+'\n')
            output_file.write(''.join(v[1:]))
        #if there are any ref transcripts remaining with no reads
        for line in ref_file.read().splitlines():
            if line.split('\t')[2] == "transcript":
                output_file.write(line+";read_num=0\n")
            else:
                output_file.write(line+'\n')
        
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
    assembled_line = assembled_file.readline()
    if assembled_line !="":
        assembled_line_fields = assembled_line.split('\t')
        
    while assembled_line_fields[2] == "exon":
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
    return chr_name,coord_starts,coord_ends,assembled_exon_chain,first_exon_end,num_reads,assembled_line_fields,assembled_line,single_exon



def process_ref_transcript(ref_line,ref_line_fields,coord_starts,first_exon_end,chr_name):
    global output_text
    global coverage_record
    global pre_output_text
    global ref_file
    global ref_exon_chain_record
    transcript_id = None
    if  chr_num_conversion(ref_line_fields[0]) > chr_num_conversion(chr_name):
         return ref_line_fields,ref_line
    
    while ref_line_fields[2] != "transcript" or int(ref_line_fields[4]) < coord_starts or chr_num_conversion(ref_line_fields[0]) < chr_num_conversion(chr_name):
        if ref_line_fields[2] != "transcript":
            pre_output_text[transcript_id].append(ref_line)
        else:
            transcript_id = ref_line_fields[8].split(";")[0]
            pre_output_text[transcript_id]=[ref_line]
        
        ref_line = ref_file.readline()
        if ref_line != "":
            ref_line_fields = ref_line.split('\t')
        else:
            return ref_line_fields,ref_line
        #linear scan
    
    while int(ref_line_fields[3]) <  first_exon_end and ref_line_fields[0] == chr_name:
        if ref_line_fields[2] == "transcript":
            field_8 = ref_line_fields[8]
            transcript_id = field_8.split(";")[0]
            coverage = float((field_8.split("cov=")[1]).split(";")[0])
            coverage_record[transcript_id] = coverage
            pre_output_text[transcript_id]=[ref_line.strip('\n')]
            ref_line = ref_file.readline()
            ref_line_fields = ref_line.split('\t')
            ref_exon_chain = ""
            
        while ref_line_fields[2] == "exon":
            pre_output_text[transcript_id].append(ref_line)
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
   
    return ref_line_fields,ref_line

 

def calc_read_proportions(assembled_exon_chain,num_reads,single_exon):
    global transcript_reads_record
    coverages = []
    transcript_ids = []
    
    
    if single_exon:
        left = assembled_exon_chain[0]
        right = assembled_exon_chain[1]
        for k in ref_exon_chain_record.keys():
            if '-'+left in ref_exon_chain_record[k] or right+'-' in ref_exon_chain_record[k]:
                coverages.append(coverage_record[k])
                transcript_ids.append(k)
    else:
        for k in ref_exon_chain_record.keys():
            if assembled_exon_chain in ref_exon_chain_record[k]:
                coverages.append(coverage_record[k])
                transcript_ids.append(k)
    coverage_proportions = [num_reads*x / sum(coverages) for x in coverages]

    for i in range(len(transcript_ids)):
        transcript_id = transcript_ids[i]
        if transcript_id in transcript_reads_record.keys():
            transcript_reads_record[transcript_id] += coverage_proportions[i]
        else:
            transcript_reads_record[transcript_id] =coverage_proportions[i]
    return


def chr_num_conversion(chro):
    num = chro.split("chr")[1]
    if num.isnumeric():
        return int(num)
    elif num == "X":
        return 23
    elif num=="Y":
        return 24
    else: #for example M
        return 25
    
if __name__ == '__main__':
    main()

