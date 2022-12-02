import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inp",default=None,help="Path to the delta file")
args = parser.parse_args()



directory = '/dev/shm/temp_alina/'

#inputfile = open('/ccb/salz4-3/Alina/22fall/out.delta', 'r')
inputfile = open(args.inp,'r')
output_file= None
Lines = inputfile.readlines()
two_line_header = []
i=1
read_name = None
read_names =  set()
for l in Lines:
   if i < 3:
     two_line_header.append(l)  
     i+=1
   else:
      if l[0]==">":
           if (i!=3 and read_name!=l.split()[1]):
               output_file.close()
               read_name = l.split()[1]
               output_file = open(directory+read_name+'.delta', 'a')
           if (i==3): #first read
              read_name = l.split()[1]
              output_file = open(directory+read_name+'.delta', 'a')
           i+=1
           if (read_name not in read_names):
              #print(read_name)
              output_file.writelines(two_line_header)
              read_names.add(read_name)
           output_file.write(l)
      else:
           output_file.write(l)
    
output_file.close()            
        
