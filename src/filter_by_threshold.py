import sys

threshold = float(sys.argv[1])

count = {}
linecount = {}
linesupport = {}
linecov = {} 
min_junc_count=2
lines = []

for line in sys.stdin:
    line = line.strip()
    lines.append(line)

for j in range(len(lines)):
    fields = lines[j].split("\t")
    if fields[2] == "transcript":
        f = fields[8].split(";")
        for item in f:
            if item.startswith("read_num="):
                count[int(item.split("=")[1])] = count.get(int(item.split("=")[1]), 0) + 1
                linecount[j] = int(item.split("=")[1])
            elif item.startswith("transcript_support="):
                linesupport[j] = float(item.split("=")[1])
            elif item.startswith("least_junction_reads_coverage="):
                linecov[j]=float(item.split("=")[1])

thresh = sum(i * count[i] for i in count) * threshold
min_count = 0
n = 0

for i in count:
    n += i * count[i]
    if n > thresh:
        min_count = i
        break

print("#gff\n#produced by NIFFLR\n#min read count =", min_count)

for j in range(len(lines)):
    if linecount.get(j, 0) > min_count or linesupport.get(j, 0) > 0.85:
        if linecov.get(j, 0) >= min_junc_count: 
            print(lines[j])
