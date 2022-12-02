import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inp",default="sorted_one_line.txt",help="Path to the input file")
    parser.add_argument("-o","--out",default=None,help="Path to the output file")
    parser.add_argument("-n","--neg",default=None,help="Path to the file containing the names of exons on the negative strand")
    args = parser.parse_args()
    
    neg_exons_list_file = open(args.neg,'r')
    neg_exons_list=neg_exons_list_file.read().splitlines()
    neg_exons_list_file.close()
    
    #majority voting and sorting by position of exon in the genome
    output_file=open(args.out, 'w')
    output_file.close()
    prev_read = "NaN"
    weight_dict={}
    exon_dict = {}
    with open(args.inp,'r') as f:
        for line in f:
            if line[0] != '>':
                continue
            info = line.split()
            read_name = info[1]
            if read_name != prev_read and prev_read != "NaN":
                write_exons(weight_dict,args.out,neg_exons_list,exon_dict,prev_read)
                weight_dict={}
                exon_dict = {}
            prev_read = read_name
            exon_name_and_pos=info[0][1:]
            exon_name='_'.join(exon_name_and_pos.split('_')[0:2])
            exon_len = int(info[2])
            exon_seg_start_index = int(exon_name_and_pos.split('_')[2])
            if exon_name not in weight_dict.keys():
                weight_dict[exon_name] = 0
                exon_dict[exon_name]=[]
            
           
            order ='+'
            matched_len = int(info[5])-int(info[4])+1
            exon_start = int(info[4])
            exon_end = int(info[5])
            match_start = int(info[6])
            match_end = int(info[7])
            if match_start > match_end:
                temp=match_start
                match_start = match_end
                match_end=temp
                order = '-'
                overhang_added_match_start = match_start-(exon_len-exon_end) #1 based
                overhang_added_match_end = match_end+exon_start-1
            else:
                overhang_added_match_start =	match_start-(exon_start-1)                                         
                overhang_added_match_end = match_end+(exon_len-exon_end)  
           
            exon_dict[exon_name].append([order,exon_seg_start_index,exon_name_and_pos,match_start,match_end,overhang_added_match_start,overhang_added_match_end,exon_len])
            matched_len_minus_errors = matched_len-int(info[8])
            weight_dict[exon_name]+= matched_len_minus_errors
          
       
    write_exons(weight_dict,args.out,neg_exons_list,exon_dict,prev_read)         
    return

def write_exons(weight_dict,out,neg_exons_list,exon_dict,read):
    if bool(weight_dict):
        best_exon=max(weight_dict, key=weight_dict.get)
        exon_data=exon_dict[best_exon]
        if len(exon_data) > 0:
           if (exon_data[0][0]=='+' and best_exon not in neg_exons_list) or (exon_data[0][0]=='-' and best_exon in neg_exons_list):
              exon_data.sort(key = lambda x: int(x[1]))
           else:
              exon_data.sort(key = lambda x: int(x[1]),reverse=True)
           with open(out, 'a') as of:    
               of.writelines(">"+read+'\n')
               n=0
               for line in exon_data:
                  n+=1
                  of.writelines('exon'+str(n)+' '+' '.join(str(item) for item in line[2:])+'\n')
            
    return

         
   
if __name__ == '__main__':
    main()
