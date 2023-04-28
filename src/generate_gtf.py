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
    parser.add_argument("-g","--good",default='good_output.gtf',help="Path to the output gtf file containing transcripts with penalty score <= 5")
    parser.add_argument("-b","--bad",default='bad_output.gtf',help="Path to the output gtf file containing transcripts with penalty score > 5")
    parser.add_argument("-n","--neg",default='negative_direction_exons.csv',help="Path to the csv file containing the negative direction exons")
    args = parser.parse_args()
    inputfile = open(args.inp,'r')
    Lines = inputfile.readlines()
    inputfile.close()
    outputgoodfile = open(args.good,'w')
    outputgoodfile.close()
    outputbadfile = open(args.bad,'w')
    outputbadfile.close()
    neg_exons_list_file = open(args.neg,'r')
    neg_exons_list=neg_exons_list_file.read().splitlines()
    neg_exons_list_file.close()
    first_line = False
    score = None
    transcript_entries_count_dict_good = {}
    transcript_entries_count_dict_bad = {}
    list_of_good_entries = []
    list_of_bad_entries = []
    num_lines = 0
    mapped_read_coord = []
    bad_entries = True
    big_dictionary_good = {} #A big dictionary. key is the sequence of exons, value is the lines to output in gtf file (both transcript and exon)
    big_dictionary_bad = {}
    for l in Lines:
        if l[0]=='>':
            first_line = True   
            if num_lines > 0:
                #process the last read
                mapped_read_info_final=mapped_read_info[:-1].copy()
                mapped_read_info_final.append(read_mapped_len)
                total_mapped = [sum(mapped_read_info_final),num_lines-1] #mapped_len,junctions covered 
                if bad_entries:
                    transcript_entries_count_dict_bad,big_dictionary_bad=add_gtf_lines(list_of_bad_entries,transcript_entries_count_dict_bad,big_dictionary_bad,total_mapped)
                else:
                    transcript_entries_count_dict_good,big_dictionary_good=add_gtf_lines(list_of_good_entries,transcript_entries_count_dict_good,big_dictionary_good,total_mapped)
            bad_entries = False
            read_info = l.strip().split('\t')
            read_name = read_info[0][1:]
            score = float(read_info[1])
            
            if score > 5 or score==-1:
                bad_entries = True
                list_of_bad_entries = [read_name,score]
            else:
                list_of_good_entries = [read_name,score]
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
                if gene_name in neg_exons_list:
                    direction = '-'
                else:
                    direction = '+'
                if bad_entries:
                    list_of_bad_entries.extend([gene_name,direction,gene_seg])
                else:
                    list_of_good_entries.extend([gene_name,direction,gene_seg])
            else:
               gene_seg = splitted[1] 
               if bad_entries:
                   list_of_bad_entries.append(gene_seg)
               else:
                   list_of_good_entries.append(gene_seg)
               
    #finally
    mapped_read_info_final=mapped_read_info[:-1].copy()
    mapped_read_info_final.append(read_mapped_len)
    total_mapped = [sum(mapped_read_info_final),num_lines-1] #mapped_len,junctions covered
    if bad_entries == False:
        transcript_entries_count_dict_good,big_dictionary_good=add_gtf_lines(list_of_good_entries,transcript_entries_count_dict_good,big_dictionary_good,total_mapped)
    elif bad_entries == True:
         transcript_entries_count_dict_bad,big_dictionary_bad=add_gtf_lines(list_of_bad_entries,transcript_entries_count_dict_bad,big_dictionary_bad,total_mapped)
    
    with open(args.good,'a') as of:
        for k, v in big_dictionary_good.items():
           transcript_first_8 = v[0] 
           transcript_attributes_col = v[1]
           transcript_attributes_col = "gene_id \"{}\"; transcript_id \"{}\"; source_reads \"{}\"; longest_mapped_read_len \"{}\";".format(transcript_attributes_col[0],transcript_attributes_col[1],transcript_attributes_col[2],transcript_attributes_col[3])
           of.write(transcript_first_8+transcript_attributes_col+'\n')
           for w in range(2,len(v)):
               if w %2 == 0:
                   exon_first_8 = v[w]
               else:
                   exon_attributes_col = v[w]
                   exon_attributes_col = "gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\; exon_id \"{}\;".format(exon_attributes_col[0],exon_attributes_col[1],exon_attributes_col[2],exon_attributes_col[3])
                   of.write(exon_first_8+transcript_attributes_col.split(" source_reads")[0]+'\n')


    with open(args.bad,'a') as of:
        for k, v in big_dictionary_bad.items():
           transcript_first_8 = v[0]
           transcript_attributes_col = v[1]
           transcript_attributes_col = "gene_id \"{}\"; transcript_id \"{}\"; source_reads \"{}\"; longest_mapped_read_len \"{}\";".format(transcript_attributes_col[0],transcript_attributes_col[1],transcript_attributes_col[2],transcript_attributes_col[3])
           of.write(transcript_first_8+transcript_attributes_col+'\n')
           for w in range(2,len(v)):
               if w %2 == 0:
                   exon_first_8 = v[w]
               else:
                   exon_attributes_col = v[w]
                   exon_attributes_col = "gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\; exon_id \"{}\;".format(exon_attributes_col[0],exon_attributes_col[1],exon_attributes_col[2],exon_attributes_col[3])
                   of.write(exon_first_8+transcript_attributes_col.split(" source_reads")[0]+'\n')


def add_gtf_lines(list_of_entries,transcript_entries_count_dict,big_dictionary,total_mapped):
    seqname = list_of_entries[2].split('-')[0]
    gene_name = "-".join(list_of_entries[2].split('-')[1:])
    read_name = list_of_entries[0]
    source = 'jguo54'
    starts = []
    ends= []
    for i in range(4,len(list_of_entries)):
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
    strand = list_of_entries[3]
    frame = '.'
    transcript_id = '-'.join([list_of_entries[2],str(transcript_start),str(transcript_end)])
    if transcript_id in transcript_entries_count_dict.keys():  
        transcript_entries_count_dict[transcript_id] +=1
    else:
        transcript_entries_count_dict[transcript_id]=1
        
    transcript_attributes = []
    unduplicated_transcript_id = transcript_id+'-'+str(transcript_entries_count_dict[transcript_id])
    big_dictionary_key = ",".join(list_of_entries[4:])
    reversed_big_dictionary_key = ",".join(list_of_entries[4:][::-1])

    if big_dictionary_key in big_dictionary.keys():
        big_dictionary[big_dictionary_key][1][2] = big_dictionary[big_dictionary_key][1][2]+","+read_name
        if total_mapped[1] >= big_dictionary[big_dictionary_key][1][4]:
            big_dictionary[big_dictionary_key][1][3]=max(big_dictionary[big_dictionary_key][1][3],total_mapped[0])
            big_dictionary[big_dictionary_key][1][4]=total_mapped[1]
        return transcript_entries_count_dict,big_dictionary #duplication found, no need to rewrite the exon lines below
    elif reversed_big_dictionary_key in big_dictionary.keys():
        big_dictionary[reversed_big_dictionary_key][1][2] = big_dictionary[reversed_big_dictionary_key][1][2]+","+read_name
        if total_mapped[1] >= big_dictionary[reversed_big_dictionary_key][1][4]:
            big_dictionary[reversed_big_dictionary_key][1][3]=max(big_dictionary[reversed_big_dictionary_key][1][3],total_mapped[0])
            big_dictionary[reversed_big_dictionary_key][1][4]=total_mapped[1]
        return transcript_entries_count_dict,big_dictionary

    else:
        updated_reads_list = read_name
        transcript_attributes= [gene_name,unduplicated_transcript_id,updated_reads_list,total_mapped[0],total_mapped[1]] #gene_id, transcript_id, source_read,read_len, covered_junctions
        transcript_entry_first_8 = '\t'.join([seqname,source,'transcript',str(transcript_start),str(transcript_end),str(score),strand,frame])+'\t'
        big_dictionary[big_dictionary_key] = [transcript_entry_first_8,transcript_attributes]
        #add exon lines below                    
        for j in range(4,len(list_of_entries)):
            exon_first_8 = '\t'.join(str(entry) for entry in [seqname,source,'exon',min(starts[j-4],ends[j-4]),max(starts[j-4],ends[j-4]),'.',strand,frame])+'\t'
            exon_attributes = [gene_name,unduplicated_transcript_id, str(j-3), unduplicated_transcript_id+'_'+str(j-3)] #gene_id, transcript_id, exon_number, exon_id
            big_dictionary[big_dictionary_key].extend([exon_first_8,exon_attributes])
        return transcript_entries_count_dict,big_dictionary
    


if __name__ == '__main__':
    main()
