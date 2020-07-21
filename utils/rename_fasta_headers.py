import sys
import re
import os

i=1
new_to_old = {}
with open(sys.argv[1]) as fasta:
    for line in fasta:
        if line.startswith(">"):
            oldname = re.sub("^>","", line).strip()
            newname = f"contig{i}"
            new_to_old[newname] = oldname
            print (f">contig_{i}")
            i+=1

        else:
            print (line.strip())

with open(f"{sys.argv[1]}.rename_map", 'w') as out:
    for n,o in new_to_old.items():
        out.write(f"{n}\t{o}\n")
