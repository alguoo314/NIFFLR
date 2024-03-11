import sys

def convert_fastq_to_fasta():
    for line in sys.stdin:
        if line.startswith("@"):  # type is fastq
            print(">" + line[1:].strip()) 
            seq = ""
            nlines = 0
            for line in sys.stdin:
                line = line.strip()
                if line.startswith("+"):
                    break
                seq += line
                nlines += 1
            print(seq)
            slines = 0
            for line in sys.stdin:
                slines += 1
                if slines == nlines:
                    break
        elif line.startswith(">"):  # type is fasta, output unchanged
            print(line.strip())
            for line in sys.stdin:
                print(line.strip())

if __name__ == "__main__":
    convert_fastq_to_fasta()
