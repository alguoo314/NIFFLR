from collections import defaultdict
import argparse
import sys
import csv
import numpy as np

def main():
    sys.setrecursionlimit(2500)
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--inp",default='majority_voted.fa',help="Path to the file after majority voting and sorting")
    parser.add_argument("-o","--outp",default='best_paths.txt',help="Path to the output file")
    args = parser.parse_args()
    inputfile = open(args.inp,'r')
    Lines = inputfile.readlines()
    inputfile.close()
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
    for l in Lines:
        if l[0]=='>':
            if exons != []:
                #if only one exon
                if exons[0][1].split("_rePlicate")[0] ==  exons[-1][1].split("_rePlicate")[0] or (exons[0][4] >=  exons[-1][4]) or exons[0][5] >=  exons[-1][5] or exons[-1][4] < 0:
                    exons.sort(key = lambda x: (int(x[5])-int(x[4])-int(x[7])),reverse=True)
                    #just output the longest mapping with min overhang penalty as possible
                    with open(outp,'a') as of:
                        score_recorder.append(0)
                        exons[0][1]=exons[0][1].split("_rePlicate")[0]
                        of.write(read_name+'\t'+str(0)+'\n')
                        of.write(' '.join(map(str,exons[0][:-1]))+'\n')
                        #of.write('The total distance = 0\n')
                        
                        
                                    
                else:        
                    score_recorder= construct_shortest_path(read_name,exons,outputfile,same_exons_record,exon_index_record,outp,score_recorder)
            read_name = l.strip()
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
            overhang = int(exon_info[2])-start-1+end-int(exon_info[3])
            overhangs_penalty = max(0,(overhang-2)*0.1)
            exons.append([exon_info[0],exon_name,exon_info[2],exon_info[3],start,end,exon_info[6],overhangs_penalty])
            exon_index_record[exon_name] = exon_info
            
     #final
    if exons[0][1].split("_rePlicate")[0] ==  exons[-1][1].split("_rePlicate")[0] or exons[0][4] >=  exons[-1][4] or exons[0][5] >=  exons[-1][5] or  exons[-1][4] < 0:
        exons.sort(key = lambda x: (int(x[5])-int(x[4])-int(x[7])),reverse=True)
        #just output the longest mapping with min overhang penalty as possible
        with open(outp,'a') as of:
            score_recorder.append(0)
            exons[0][1]=exons[0][1].split("_rePlicate")[0]
            of.write(read_name+'\t'+str(0)+'\n')
            of.write(' '.join(map(str,exons[0][:-1]))+'\n')
            #of.write('The total distance = 0\n')             
    else:        
        construct_shortest_path(read_name,exons,outputfile,same_exons_record,exon_index_record,outp,score_recorder)
    
    
    with open('scores.csv','w') as file:   
        csvwriter = csv.writer(file,delimiter='\n')
        csvwriter.writerow(score_recorder)
         

def extract_w(lst):
    return [item[2] for item in lst]
 

          
def construct_shortest_path(read_name,exons,outputfile,same_exons_record,exon_index_record,outp,score_recorder):
    exon_names_list = [x[1] for x in exons]
    overhang_penalties = [x[-1] for x in exons]
    overhang_pen_dict = dict(map(lambda i,j : (i,j) , exon_names_list,overhang_penalties))
    g = Graph(exon_names_list)
    edges_candidates = []
    unconnected_nodes = set()
    #path_len_dict ={} #negative for overlaps
    for i in range(len(exons)):
        node_1_info = exons[i]
        node_1_name =  node_1_info[1]
        unconnected_nodes.add(node_1_name)
        for j in range(i+1,len(exons)):
            node_2_info = exons[j]
            node_2_name =  node_2_info[1]
            diff = node_2_info[4]-node_1_info[5]
            #nonneg weight (gap or overlap)
            if diff > 0: 
                length = abs(diff)-1
            else:
                length = abs(diff)+1
            edges_candidates.append([node_1_name, node_2_name, length])

    weights = extract_w(edges_candidates)
    thres = (np.min(weights)+10)*2
    for e in edges_candidates:
        g.addEdge(e[0], e[1], e[2]) #still add the edge even if the gap/overlap is above the threshold. Prevent segmentation in the middle
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
        if exons[ex][4] <= 0 or (exons[ex][5] <= exons[0][5] and exons[ex][4] < exons[ex-1][5]): #the right part to the "and" is to prevent starting in the middle of an exon chain
            name = exons[ex][1]
            if name not in unconnected_nodes:
                origin_candidate_nodes.add(exons[ex][1])
            #add all exons that start with the very beginning of the read
    
    for ex2 in reversed(range(len(exons)-1)):
        if exons[ex2][5] >= exons[-1][5] or exons[ex2][4]== exons[-1][4]:
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
        for node1 in origin_candidate_nodes:
            for node2 in destination_candidate_nodes:
                if node1 != node2:
                    dist, overhang_penalized_dist,path = g.shortestPath(node1,node2,overhang_pen_dict)
                    if dist != None:
                        possible_paths.append([dist, overhang_penalized_dist, path])

        possible_paths.sort(key = lambda x: int(x[1]))
        best = possible_paths[0]
        best_dist = best[0]
        best_path = best[2]
        with open(outp,'a') as of:
            if best_path[0].split("_rePlicate")[0] ==  best_path[-1].split("_rePlicate")[0]:
                score = 0
                score_recorder.append(score)
                of.write(read_name+'\t'+str(0)+'\n')
                exons.sort(key = lambda x: (int(x[5])-int(x[4])-int(x[7])),reverse=True)
                of.write(' '.join(map(str,exons[0]))+'\n')
            else:
                score = round(best_dist/(len(best_path)-1),3)
                score_recorder.append(score)
                of.write(read_name+'\t'+str(score)+'\n')
                for node in best_path:
                    to_be_written = exon_index_record[node]
                    if "rePlicate" in to_be_written[1]:
                        to_be_written[1]=to_be_written[1].split("rePlicate")[0]
                        
                    to_be_written[1] = str(to_be_written[1])
                    to_be_written[2]=str(to_be_written[2])
                    of.write(' '.join(to_be_written)+'\n')
       
    else:
        with open(outp,'a') as of:
            score_recorder.append(0)
            exons[0][1]=exons[0][1].split("_rePlicate")[0]
            of.write(read_name+'\t'+str(0)+'\n')
            of.write(' '.join(map(str,exons[0]))+'\n')
    
    return score_recorder
            
        
    
 


class Graph:
    def __init__(self,vertices):
 
        self.V = len(vertices) # No. of vertices
        self.verts = vertices
            
 
        # dictionary containing adjacency List
        self.graph = defaultdict(list)
 
    # function to add an edge to graph
    def addEdge(self,u,v,w):
        self.graph[u].append((v,w))
 
 
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
    def shortestPath(self,s,s2,overhang_pen_dict):
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
            for node,weight in self.graph[i]:
                if dist_overhang_penalized[node] > dist_overhang_penalized[i] + weight + overhang_pen_dict[node]:
                    dist[node] = dist[i] + weight
                    dist_overhang_penalized[node] = dist_overhang_penalized[i] + weight + overhang_pen_dict[node]
                    prevs[node]=i
 
        # Print the calculated shortest distances
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
