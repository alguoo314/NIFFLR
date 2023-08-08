import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g","--aftergffread",default='check_discarded.combined.gtf',help="Path to gtf file as the product of gffread -T --cluster-only")
    parser.add_argument("-a","--anno",default='gffcmp.annotated.gtf',help="Path to the annotated.gtf file as the product of gffcompare -r reference.gff combined.gtf (produced by gffcompare -D)")
    parser.add_argument("-u","--unfiltered",default='final_output.gtf',help="Path to the unfiltered gtf file")
    parser.add_argument("-o","--output",default='reads_num_added_annotated.gtf',help="Output file: annotated.gtf with read number added for each transcript")
    
    
    args = parser.parse_args()
    id_reads_num_dict = match_transcript_id_and_read_count(args.unfiltered)
    add_reads_num_to_anno(args.anno,args.output,id_reads_num_dict)


def match_transcript_id_and_read_count(unfiltered):
    unfiltered_file = open(unfiltered,'r')
    unfiltered_lines = unfiltered_file.readlines()
    unfiltered_file.close()
    
    id_reads_num_dict = {}
    for l in unfiltered_lines:
        if l.split('\t')[2] == "transcript":
            attributes = l.split('\t')[-1]
            ID = attributes.split('\"')[3]
            reads_num =  len(attributes.split('\"')[5].split(','))
            id_reads_num_dict[ID]=reads_num
        else:
            continue
    
    
        
    return id_reads_num_dict
    
    
def add_reads_num_to_anno(anno,output,id_reads_num_dict):
    annofile = open(anno,'r')
    anno_lines = annofile.readlines()
    annofile.close()
    
    with open(output,'w') as of:
        for l in  anno_lines:
          if l.split('\t')[2] == "transcript":
              attributes =  l.split('\t')[-1]
              transcript_id = attributes.split('\"')[1]
              reads_num = id_reads_num_dict[transcript_id]
              of.write(l.strip()+" read_count \"{}\";".format(reads_num)+'\n')
          else:
              of.write(l)
    return    
    
if __name__ == '__main__':
    main()
