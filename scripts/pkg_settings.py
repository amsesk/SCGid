SPDB_OS_REGEXP_PATTERN = "OS=(.+?)(?:OX=.+|GN=.+|PE=.+|SV=.+|$)"
BLAST_OUTFMT = [
    "6",
    "qseqid", 
    "sseqid", 
    "pident",  
    "length", 
    "mismatch", 
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send", 
    "evalue", 
    "bitscore", 
    "staxids"
    ]
BLAST_OUTFMT_STR = ' '.join(BLAST_OUTFMT)
BLAST_HEADERS = {v: k for k,v  in enumerate(BLAST_OUTFMT[1::])}
