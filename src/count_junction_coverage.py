import argparse

def parse_gff_for_transcripts(gff_file):
    transcripts = {}
    with open(gff_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if fields[2] == 'transcript':
                    attributes = fields[8].split(';')
                    
                    key, value = attributes[2].split('=')
                    read_count = len(value.split(','))
                    transcript_id = attributes[0].split('=')[1]
                    transcripts[transcript_id]={'exons':[],'reads':read_count}

                if fields[2] == 'exon':
                    read_count = None
                    attributes = fields[8].split(';')
                    if attributes[0].split('=')[1] == transcript_id:
                        chromosome = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        strand = fields[6]
                    
                    if len(transcripts[transcript_id]['exons'])>0:
                        transcripts[transcript_id]['exons'].append([start, end])
                    else:
                        transcripts[transcript_id]['exons']=[[chromosome, start, end, strand]]
    return transcripts

def generate_junction_chains(transcripts):
    junctions = {}
    for transcript_id, data in transcripts.items():
        exons = data['exons']
        
        #exons.sort(key=lambda x: x[1])  # Sort exons by start position
        if len(exons) < 2:
            continue

        def_info = exons[0][0]+exons[0][-1]
        if def_info+str(exons[0][2])+","+str(exons[1][0]) in junctions:
            junctions[def_info+str(exons[0][2])+","+str(exons[1][0])] += data['reads']
        else:
            junctions[def_info+str(exons[0][2])+","+str(exons[1][0])] = data['reads']
        for i in range(1,len(exons) - 1):
            if def_info+str(exons[i][1])+","+str(exons[i+1][0]) in junctions:
                junctions[def_info+str(exons[i][1])+","+str(exons[i+1][0])] += data['reads']
            else:
                junctions[def_info+str(exons[i][1])+","+str(exons[i+1][0])]= data['reads']
    return junctions

def write_junction_counts(junctions, output_file):
    with open(output_file, 'w') as f:
        f.write("Junction\tRead Count\n")
        for junction, count in junctions.items():
            f.write(f"{junction}\t{count}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input",default="output.sorted.good_output.gff",help="Path to the sorted gff file containing the assembled exon paths")
    parser.add_argument("-o","--output",default='exon_junction_counts.csv',help="Path to the csv file containing the read counts per junction")
    args = parser.parse_args()
    gff_file = args.input
    output_file = args.output
    transcripts = parse_gff_for_transcripts(gff_file)
    junctions = generate_junction_chains(transcripts)
    write_junction_counts(junctions, output_file)
