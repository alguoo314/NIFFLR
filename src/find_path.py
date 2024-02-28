#!/usr/bin/env python
from collections import defaultdict
import argparse
import sys
import csv
import numpy as np

def main():
    sys.setrecursionlimit(2500)
    parser = argparse.ArgumentParser()
    parser.add_argument("-o","--outp",default='best_paths.txt',help="Path to the output file")
    args = parser.parse_args()
    outputfile = open(args.outp,'w')
    outputfile.close()
    outp=args.outp
    exons = []
    read_name = None
    exon_set = None
    same_exons_record=None
    exon_index_record = None
    i=0
    score_recorder = []
    to_be_written = []
    read_counter = 0
    for l in sys.stdin:
        if l[0]=='>':
            if read_counter == 1000: #write 1000 reads to file
                read_counter = 0
                with open(outp,'a') as of:
                    of.write("".join(to_be_written))
                to_be_written = []
            if exons != []:
                #if only one exon
                if exons[0][1].split("_rePlicate")[0] ==  exons[-1][1].split("_rePlicate")[0] or (exons[0][4] >=  exons[-1][4]) or exons[0][5] >=  exons[-1][5] or exons[-1][4] < 0 or exons[0][2] == exons[-1][2] or exons[0][3] >=  exons[-1][3]:
                    exons.sort(key = lambda x: (int(x[5])-int(x[4])-int(x[7])),reverse=True)
                    #just output the longest mapping with min overhang penalty as possible
                    score_recorder.append(0)
                    exons[0][1]=exons[0][1].split("_rePlicate")[0]
                    to_be_written.append(str(read_name+'\t'+str(0)+'\n'))
                    to_be_written.append('\t'.join(map(str,exons[0][:-1]))+'\n')
                    read_counter +=1 
                    
                else:
                    chr_and_gene_name = '-'.join(exons[0][1].split('-')[:-2])
                    if (chr_and_gene_name[-1] == 'R')+(mapped_ori=='-') == 1: #only 1 is true
                        neg = True
                    else:
                        neg = False
                    score_recorder,to_be_written,read_counter= construct_shortest_path(read_name,exons,outputfile,same_exons_record,exon_index_record,outp,score_recorder,to_be_written,read_counter,neg)
            read_name_and_mapped_ori = l.strip().split('\t')
            read_name = read_name_and_mapped_ori[0]
            mapped_ori=read_name_and_mapped_ori[1]
            exons = []
            same_exons_record = {}
            exon_index_record={}
            exon_set = set()
            i=0
        else:
            exon_info = l.split()
            if exon_info[1] in exon_set:
                exon_name = exon_info[1]+"_rePlicate_"+str(i) #in case the exactly same exon is mapped multiple times at slightly different places on the reads 
                same_exons_record[exon_info[1]].append(exon_name)
                i+=1
            else:
                exon_name = exon_info[1]
                exon_set.add(exon_info[1])
                same_exons_record[exon_info[1]]=[exon_name]
            start = int(exon_info[4])
            end = int(exon_info[5])
            #if start > 0:
                #overhang = int(exon_info[2])-start+end-int(exon_info[3])
            #else: # do not penalize the overhang that goes before the start of the read
                #overhang = end-int(exon_info[3])
            overhang = int(exon_info[2])-start+end-int(exon_info[3]) #feb 2024 edit: now we penalize start/end overhangs
            overhangs_penalty = max(0,(overhang-2)*0.1)
            exons.append([exon_info[0],exon_name,int(exon_info[2]),int(exon_info[3]),start,end,int(exon_info[6]),float(overhangs_penalty)])
            exon_index_record[exon_name] = exon_info
            
     #final
    
    if exons[0][1].split("_rePlicate")[0] ==  exons[-1][1].split("_rePlicate")[0] or (int(exons[0][4]) >=  int(exons[-1][4])) or int(exons[0][5]) >=  int(exons[-1][5]) or int(exons[-1][4]) < 0 or int(exons[0][2]) == int(exons[-1][2]) or int(exons[0][3]) >=  int(exons[-1][3]):
        exons.sort(key = lambda x: (int(x[5])-int(x[4])-int(x[7])),reverse=True)
        #just output the longest mapping with min overhang penalty as possible
        with open(outp,'a') as of:
            of.write("".join(to_be_written))
            score_recorder.append(0)
            exons[0][1]=exons[0][1].split("_rePlicate")[0]
            of.write(read_name+'\t'+str(0)+'\n')
            of.write('\t'.join(map(str,exons[0][:-1]))+'\n')             
    else:
        chr_and_gene_name = '-'.join(exons[0][1].split('-')[:-2])
        if (chr_and_gene_name[-1] == 'R')+(mapped_ori=='-') == 1:
            neg = True                 
        else:
            neg = False
        score_recorder,to_be_written,read_counter=construct_shortest_path(read_name,exons,outputfile,same_exons_record,exon_index_record,outp,score_recorder,to_be_written,read_counter,neg)
        with open(outp,'a') as of:
            of.write("".join(to_be_written))
    
    with open('scores.csv','w') as file: #this is mostly for statsitical calculation of the percentage of each score range before determining the cutoff was 5. Can remove before publishing
        csvwriter = csv.writer(file,delimiter='\n')
        csvwriter.writerow(score_recorder)
         

def extract_w(lst):
    return [item[2] for item in lst]
 

          
def construct_shortest_path(read_name,exons,outputfile,same_exons_record,exon_index_record,outp,score_recorder,to_be_written,read_counter,neg):
    exon_names_list = [x[1] for x in exons]
    overhang_penalties = [x[-1] for x in exons]
    overhang_pen_dict = dict(map(lambda i,j : (i,j) , exon_names_list,overhang_penalties))
    g = Graph(exon_names_list)
    edges_candidates = []
    unconnected_nodes = set()
    for i in range(len(exons)):
        node_1_info = exons[i]
        node_1_name =  node_1_info[1]
        node_1_coords = [int(node_1_name.split('-')[-2]),int(node_1_name.split('-')[-1])]
        unconnected_nodes.add(node_1_name)
        for j in range(i+1,len(exons)):
            node_2_info = exons[j]
            node_2_name =  node_2_info[1]
            diff = node_2_info[4]-node_1_info[5]
            node_2_coords = [int(node_2_name.split('-')[-2]),int(node_2_name.split('-')[-1])]
            
            #nonneg weight (gap or overlap)
            if diff > 0: 
                length = abs(diff)-1
            else:
                length = abs(diff)+1

            if (neg == False and node_1_coords[0] > node_2_coords[0]) or (neg and node_1_coords[0] < node_2_coords[0]):
                length = float("inf")
                #the order of exons in the gff reference file cannot be violated here 
            
            if check_overlap([node_1_coords,node_2_coords]):
                overlap_forbidden = float("inf")
            else:
                overlap_forbidden = length
             
            edges_candidates.append([node_1_name, node_2_name, length,overlap_forbidden])

    weights = extract_w(edges_candidates)
    thres = (np.min(weights)+10)*2
    for e in edges_candidates:
        g.addEdge(e[0], e[1], e[2],e[3]) #still add the edge even if the gap/overlap is above the threshold. Prevent segmentation in the middle
        if e[2]<=thres:
            unconnected_nodes.discard(e[0])
            unconnected_nodes.discard(e[1])
    
    while exons[0][1] in unconnected_nodes:
        del exons[0]
    while exons[-1][1] in unconnected_nodes:
        del exons[-1]
    
    origin = exons[0]
    destination=exons[-1]
    
    origin_candidate_nodes = set()
    origin_candidate_nodes.add(origin[1])
    destination_candidate_nodes = set()
    destination_candidate_nodes.add(destination[1])
    for ex in range(1,len(exons)):
        if int(exons[ex][4]) <= 0 or int(exons[ex][4]) <= int(exons[0][4]) or int(exons[ex][2]) <= int(exons[0][2]): 
            name = exons[ex][1]
            if name not in unconnected_nodes:
                origin_candidate_nodes.add(exons[ex][1])
            #add all exons that start with the very beginning of the read
    
    for ex2 in reversed(range(len(exons)-1)):
        if int(exons[ex2][5]) >= int(exons[-1][5]) or int(exons[ex2][3]) >= int(exons[-1][3]):
            name = exons[ex2][1]
            if name not in unconnected_nodes:
                destination_candidate_nodes.add(name)
    
    if origin[1].split("_rePlicate")[0] != destination[1].split("_rePlicate")[0]:              
        if ("rePlicate" in origin[1] or origin[1] in same_exons_record.keys()):
            for name in same_exons_record[origin[1].split("_rePlicate")[0]]:
                if name not in unconnected_nodes:
                    origin_candidate_nodes.add(name)
        if ("rePlicate" in destination[1] or destination[1] in same_exons_record.keys()):
            for name in same_exons_record[destination[1].split("_rePlicate")[0]]:
                if name not in unconnected_nodes:
                    destination_candidate_nodes.add(name)

        possible_paths = []
        possible_paths_no_overlaps = []
        intersec=origin_candidate_nodes.intersection(destination_candidate_nodes)
        if len(intersec) > 0:
            for node in intersec:
                possible_paths.append([0,overhang_pen_dict[node],[node]])
                possible_paths_no_overlaps.append([0,overhang_pen_dict[node],[node]])
        else:            
            for node1 in origin_candidate_nodes:
                for node2 in destination_candidate_nodes:
                    dist, overhang_penalized_dist,path = g.shortestPath(node1,node2,overhang_pen_dict,True) #allow overlap
                    dist_over, overhang_penalized_dist_no_over,path_no_over = g.shortestPath(node1,node2,overhang_pen_dict,False) #not allow overlap
                    if dist != None:
                        possible_paths.append([dist, overhang_penalized_dist, path])
                    if dist_over != None:
                        possible_paths_no_overlaps.append([dist_over, overhang_penalized_dist_no_over,path_no_over])
            
            possible_paths.sort(key = lambda x: (int(x[1])/(len(x[2])-1),-1*(int(exon_index_record[x[2][-1]][3]) - int(exon_index_record[x[2][0]][2]))))
            possible_paths_no_overlaps.sort(key = lambda x: (int(x[1])/(len(x[2])-1),-1*(int(exon_index_record[x[2][-1]][3]) - int(exon_index_record[x[2][0]][2]))))
        # if there is a tie of score between two paths, choose the path that spans the most of the reads (actual match, not overhangs)
        
        best = possible_paths[0]
        best_dist = best[0]
        best_path = best[2]
        best_path_no_overlap = None
        if len(possible_paths_no_overlaps) > 0:
            best_no_overlap = possible_paths_no_overlaps[0]
            best_dist_no_overlap =best_no_overlap[0]
            best_path_no_overlap = best_no_overlap[2]

        
        if best_path[0].split("_rePlicate")[0] ==  best_path[-1].split("_rePlicate")[0]:
            score = 0
            score_recorder.append(score)
            to_be_written.append(str(read_name+'\t'+str(0)+'\n'))
            exons.sort(key = lambda x: (int(x[5])-int(x[4])-int(x[7])),reverse=True)
            to_be_written.append('\t'.join(map(str,exons[0]))+'\n')
        else:
            
            #check for overlap: It refers to the overlap between exon positions in reference file, not where it aligns to the reads            
            if len(possible_paths_no_overlaps)==0 or best_path == best_path_no_overlap or (best_dist_no_overlap/(len(best_path_no_overlap)-1) > 5):
                if best_path != best_path_no_overlap:
                    score =-1
                else:
                    score = round(best_dist/(len(best_path)-1),3)
                to_be_written.append(str(read_name+'\t'+str(score)+'\n'))
                for node in best_path:
                    exon_info = exon_index_record[node]
                    if "rePlicate" in exon_info[1]:
                        exon_info[1]=exon_info[1].split("_rePlicate")[0]
                    exon_info_splitted_again=exon_info[1].split('-')
                    exon_info[1] = str(exon_info[1])
                    exon_info[2]=str(exon_info[2])
                    to_be_written.append(str('\t'.join(exon_info)+'\n'))
                
                score_recorder.append(score)
            
            elif best_dist_no_overlap/(len(best_path_no_overlap)-1) <= 5:
                score = round(best_dist_no_overlap/(len(best_path_no_overlap)-1),3)
                to_be_written.append(str(read_name+'\t'+str(score)+'\n'))
                for node in best_path_no_overlap:
                    exon_info = exon_index_record[node]
                    if "rePlicate" in exon_info[1]:
                        exon_info[1]=exon_info[1].split("_rePlicate")[0]
                    exon_info_splitted_again=exon_info[1].split('-')
                    exon_info[1] = str(exon_info[1])
                    exon_info[2]=str(exon_info[2])
                    to_be_written.append(str('\t'.join(exon_info)+'\n'))
                
                score_recorder.append(score)
                
    else:
        read_counter +=1
        score_recorder.append(0)
        exons[0][1]=exons[0][1].split("_rePlicate")[0]
        to_be_written.append(str(read_name+'\t'+str(0)+'\n'))
        to_be_written.append('\t'.join(map(str,exons[0]))+'\n')
        
            
    return score_recorder,to_be_written,read_counter
            
        
    
 
def check_overlap(intervals):
    intervals.sort() 
    if intervals[0][1] >= intervals[1][0]:
        return True
    return False

class Graph:
    def __init__(self,vertices):
 
        self.V = len(vertices) # No. of vertices
        self.verts = vertices
            
 
        # dictionary containing adjacency List
        self.graph = defaultdict(list)
 
    # function to add an edge to graph
    def addEdge(self,u,v,w,w_no_overlap):
        self.graph[u].append((v,[w,w_no_overlap]))
 
 
    # A recursive function used by shortestPath
    def topologicalSortUtil(self,v,visited,stack):
 
        # Mark the current node as visited.
        visited[v] = True
 
        # Recur for all the vertices adjacent to this vertex
        if v in self.graph.keys():
            for node,weight in self.graph[v]:
                if visited[node] == False:
                    self.topologicalSortUtil(node,visited,stack)
 
        # Push current vertex to stack which stores topological sort
        stack.append(v)
 
 
    ''' The function to find shortest paths from given vertex.
        It uses recursive topologicalSortUtil() to get topological
        sorting of given graph.'''
    def shortestPath(self,s,s2,overhang_pen_dict,overlap_allowed):
        visited = {}
        # Mark all the vertices as not visited
        for v in self.verts:
            visited[v] = False
        stack =[]
 
        # Call the recursive helper function to store Topological
        # Sort starting from source vertices
        for i in (self.verts):
            if visited[i] == False:
                self.topologicalSortUtil(s,visited,stack)
 
        # Initialize distances to all vertices as infinite and
        # distance to source as 0
        dist ={}
        dist_overhang_penalized = {}
        for v in self.verts:
            dist[v] = float("Inf")
            dist_overhang_penalized[v] = float("Inf")
        dist[s] = 0
        dist_overhang_penalized[s] = overhang_pen_dict[s]
        prevs = {}
        prevs[s] = None
        
 
        # Process vertices in topological order
        while stack:
 
            # Get the next vertex from topological order
            i = stack.pop()
 
            # Update distances of all adjacent vertices
            if overlap_allowed:
                for node,weight in self.graph[i]:
                    if dist_overhang_penalized[node] > dist_overhang_penalized[i] + weight[0] + overhang_pen_dict[node]:
                        dist[node] = dist[i] + weight[0]
                        dist_overhang_penalized[node] = dist_overhang_penalized[i] + weight[0] + overhang_pen_dict[node]
                        prevs[node]=i
            else:
                for node,weight in self.graph[i]:
                    if dist_overhang_penalized[node] > dist_overhang_penalized[i] + weight[1] + overhang_pen_dict[node]:
                        dist[node] = dist[i] + weight[1]
                        dist_overhang_penalized[node] = dist_overhang_penalized[i] + weight[1] + overhang_pen_dict[node]
                        prevs[node]=i
        
        final_dist = dist[s2]
        final_dist_overhang_penalized = dist_overhang_penalized[s2]
        path = []
        temp= s2
        path.append(temp)
        
        if temp not in prevs.keys():
            return None,None,None
        while prevs[temp] != None:
            prev=prevs[temp]
            path.append(prev)
            temp=prev
            if temp not in prevs.keys():
                return None,None,None
            
        path.reverse()
        return final_dist,final_dist_overhang_penalized,path
    
    
    
if __name__ == '__main__': 
    main() 
