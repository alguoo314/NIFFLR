#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 19:26:56 2022

@author: alibenbar
"""
import argparse
import collections

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inp",default='best_paths.txt',help="Path to the file containing the best paths")
    parser.add_argument("-g","--good",default='good_output.gtf',help="Path to the output gtf file containing transcripts formed from the best paths")
    args = parser.parse_args()
    inputfile = open(args.inp,'r')
    Lines = inputfile.readlines()
    inputfile.close()
    outputgoodfile = open(args.good,'w')
    outputgoodfile.close()
    
    
    first_line = False
    score = None
    transcript_entries_count_dict_good = {}
    transcript_entries_count_dict_bad = {}
    list_of_good_entries = []
    list_of_bad_entries=[]
    num_lines = 0
    mapped_read_coord = []
    big_dictionary_good = {} #A big dictionary. key is the sequence of exons, value is the lines to output in gtf file (both transcript and exon)
    for l in Lines:
        if l[0]=='>':
            first_line = True   
            if num_lines > 0:
                #process the last read
                mapped_read_info_final=mapped_read_info[:-1].copy()
                mapped_read_info_final.append(read_mapped_len)
                total_mapped = sum(mapped_read_info_final) #mapped_len 
                transcript_entries_count_dict_good,big_dictionary_good=add_gtf_lines(list_of_good_entries,transcript_entries_count_dict_good,big_dictionary_good,total_mapped,list_of_bad_entries)
                
            read_info = l.strip().split('\t')
            read_name = read_info[0][1:]
            score = float(read_info[1])
            maxscore = float(read_info[2])
            if score==-1:
                list_of_bad_entries.append(read_name)
            list_of_good_entries=[read_name,score,maxscore]

            num_lines=0
            mapped_read_info = []
            
  
        else:
            num_lines+=1
            splitted = l.split('\t')
            read_mapped_len = int(splitted[3])-int(splitted[2])+1
            if num_lines == 1:
                mapped_read_info.append(read_mapped_len)
            else:
                mapped_read_info.append(int(splitted[6]))
            if first_line == True:
                first_line = False
                gene_seg = splitted[1]
                gene_name = '-'.join(gene_seg.split('-')[:-2])
                if gene_name[-1] == 'R':
                    direction = '-'
                else:
                    direction = '+'
                gene_name=gene_name[:-2]
                list_of_good_entries.extend([gene_name,direction,gene_seg])
            else:
               gene_seg = splitted[1]
               list_of_good_entries.append(gene_seg)
               
    #finally
    mapped_read_info_final=mapped_read_info[:-1].copy()
    mapped_read_info_final.append(read_mapped_len)
    total_mapped = sum(mapped_read_info_final) #mapped_len
    transcript_entries_count_dict_good,big_dictionary_good=add_gtf_lines(list_of_good_entries,transcript_entries_count_dict_good,big_dictionary_good,total_mapped,list_of_bad_entries)
    
    
    with open(args.good,'a') as of:
        for k, v in big_dictionary_good.items():
           transcript_first_8 = v[0] 
           transcript_attributes_col = v[1]
           
           transcript_attributes_col = "gene_id \"{}\"; transcript_id \"{}\"; source_reads \"{}\"; longest_mapped_read_len \"{}\"; best_matched_reads_avg_penality_score \"{}\"; best_matched_reads_max_penality_score \"{}\";".format(transcript_attributes_col[0],transcript_attributes_col[1],transcript_attributes_col[2],transcript_attributes_col[3],transcript_attributes_col[4],transcript_attributes_col[5])
           of.write(transcript_first_8+transcript_attributes_col+'\n')
           for w in range(2,len(v)):
               if w %2 == 0:
                   exon_first_8 = v[w]
               else:
                   exon_attributes_col = v[w]
                   exon_attributes_col = "gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\; exon_id \"{}\;".format(exon_attributes_col[0],exon_attributes_col[1],exon_attributes_col[2],exon_attributes_col[3])
                   of.write(exon_first_8+transcript_attributes_col.split(" source_reads")[0]+'\n')




def add_gtf_lines(list_of_entries,transcript_entries_count_dict,big_dictionary,total_mapped,list_of_bad_entries):
    seqname = list_of_entries[3].split('-')[0]
    gene_name = "-".join(list_of_entries[3].split('-')[1:])
    read_name = list_of_entries[0]
    source = 'jguo54'
    starts = []
    ends= []
    for i in range(5,len(list_of_entries)):
        s = int(list_of_entries[i].split('-')[-2])
        e = int(list_of_entries[i].split('-')[-1])
        starts.append(s)
        ends.append(e)
    
    if starts[0] > ends[0]: #reversed mapping, index 3 are ends, 4 are starts
        transcript_start = min(ends)
        transcript_end = max(starts)
    else:
        transcript_start = min(starts)
        transcript_end = max(ends)
    
    score = list_of_entries[1]
    maxscore = list_of_entries[2]
    strand = list_of_entries[4]
    frame = '.'
    transcript_id = '-'.join([list_of_entries[3],str(transcript_start),str(transcript_end)])
    if transcript_id in transcript_entries_count_dict.keys():
        transcript_entries_count_dict[transcript_id] +=1
    else:
        transcript_entries_count_dict[transcript_id]=1
        
    transcript_attributes = []
    unduplicated_transcript_id = transcript_id+'-'+str(transcript_entries_count_dict[transcript_id])
    big_dictionary_key = ",".join(list_of_entries[5:])
    reversed_big_dictionary_key = ",".join(list_of_entries[5:][::-1])
    
    if big_dictionary_key in big_dictionary.keys():
        if read_name not in list_of_bad_entries:
            big_dictionary[big_dictionary_key][1][2] = big_dictionary[big_dictionary_key][1][2]+","+read_name
            big_dictionary[big_dictionary_key][1][3]=max(big_dictionary[big_dictionary_key][1][3],total_mapped)
            if big_dictionary[big_dictionary_key][1][4] > score:
                big_dictionary[big_dictionary_key][1][5] = maxscore
                big_dictionary[big_dictionary_key][1][4]=score
        
        return transcript_entries_count_dict,big_dictionary #duplication found, no need to rewrite the exon lines below

    elif reversed_big_dictionary_key in big_dictionary.keys():
        if read_name not in list_of_bad_entries:
            big_dictionary[reversed_big_dictionary_key][1][2] = big_dictionary[reversed_big_dictionary_key][1][2]+","+read_name
            big_dictionary[reversed_big_dictionary_key][1][3]=max(big_dictionary[reversed_big_dictionary_key][1][3],total_mapped)
            if big_dictionary[reversed_big_dictionary_key][1][4] > score:
                big_dictionary[reversed_big_dictionary_key][1][5] = maxscore
                big_dictionary[reversed_big_dictionary_key][1][4]=score
        return transcript_entries_count_dict,big_dictionary

    else:
        if read_name in list_of_bad_entries:
            return transcript_entries_count_dict,big_dictionary
        transcript_attributes= [gene_name,unduplicated_transcript_id,read_name,total_mapped,score,maxscore] #gene_id, transcript_id, source_read,read_len,avg_gap_ovrelap_penalty,max_overlap_gap_penalty
        transcript_entry_first_8 = '\t'.join([seqname,source,'transcript',str(transcript_start),str(transcript_end),'.',strand,frame])+'\t'
        big_dictionary[big_dictionary_key] = [transcript_entry_first_8,transcript_attributes]
        #add exon lines below                    
        for j in range(5,len(list_of_entries)):
            exon_first_8 = '\t'.join(str(entry) for entry in [seqname,source,'exon',min(starts[j-5],ends[j-5]),max(starts[j-5],ends[j-5]),'.',strand,frame])+'\t'
            exon_attributes = [gene_name,unduplicated_transcript_id, str(j-3), unduplicated_transcript_id+'_'+str(j-4)] #gene_id, transcript_id, exon_number, exon_id
            big_dictionary[big_dictionary_key].extend([exon_first_8,exon_attributes])
        return transcript_entries_count_dict,big_dictionary
    


if __name__ == '__main__':
    main()
