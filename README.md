## ***scgid***, a python-based tool for scaffold binning and genome prediction from single-cell sequencing libraries
####version 0.1b
-----

### What is ***scgid*** for?
*scgid* is a python-based tool aimed at extracting a draft genome sequence for a single organism from a mildly metagenomic sequencing library resulting from single-cell genomic amplifications. The only thing that you need to start is your assebmly in FASTA format and *scgid* will do the rest.

*scgid* takes your assebmly and subjects it to three binning methods each based on a different sequence signature. It takes the output of each, draws consensus based on majority-rule, and yields one final consensus-genome draft at the intersection of all three methods. 

As of this version, *scgid* is only compatible with genome assemblies generated with SPAdes (http://bioinf.spbau.ru/en/spades_for_remove) that have their sequence fasta headers **unchanged** (this is very important - just rename them after you're done with *scgid*). 

### Dependencies
* python 2.7 (***scgid* is not compatbile with python 3.x.x**)
	* pandas 1.15.0
	* numpy 0.23.4
	* ETE3 toolkit 3.1.1
	` pip install --user pandas numpy ete3`
	or
	`conda install pandas numpy ete3`
	or you can download and build these packages from source.
* R 3.4.1
	* ape 
	`install.packages("ape")`
	* Rscript (included with most R distributions)
* NCBI Blast+
* Augustus (http://bioinf.uni-greifswald.de/augustus/)
* ClaMS-CLI (http://clams.jgi-psf.org/)
    Scroll to the bottom of the page a and click on "Download"
    Accept the license agreement
    Download the "Command lin verion" of ClaMS
* Java
* Databionic ESOM (http://www.databionin-esom.sourceforge.net)
* BASH

### Downloading and Installing ***scgid***
* Download the tar.gz file from (LINK) or clone the repository at (GITHUB LINK)
* If you downloaded the .tar.gz file, decompress and untar the archive with...
`tar -xvf scgid-0.1b.tar.gz`
* Navigate to `bin` within the newly created directory
`cd scgid-0.1b/bin`

#### Automatic Setup (recommended)
* Type `./scgid init`
	* Follow the prompts to define some package-wide settings and download the necessary databases.
	* This script **will try** to install the required python packages, if they are not already available, using `conda` or `pip`. Make sure you know which package manager is being used on your system as installation via `pip` has been known to break some conda installations.
	* You are responsible for having downloaded and installed all of the third-party dependencies listed above. You will be asked where some of them are.
	* This script requires an internet connection in order to download the swissprot databases.
	
* Modify your `.bashrc` or `.bash_profile` file and add `scgid-0.1b/bin` to your enviornmental $PATH variable. For instance, add a line like... `export PATH=$PATH:/path/to/scgid-0.1b/bin`
* Ensure that other stand-alone dependencies (i.e. BLAST and Augustus) have also been added to $PATH.

#### Manual Setup (if auto isn't working)
* Download and decompress a copy of the most recent swissprot database.
* Edit `settings.py` to reflect the locations of ESOM, ClaMS, and Rscript, and your swissprot database, as follows:
```
esom_path ="<path_to_folder_containing_bin>"`
clams_path="<path_to_folder_containing_clams-cli_jar>"
path_to_Rscript="<path_to_folder_containing_Rscript>"
path_to_spdb="<path_to_swissprot_database>"
spdb_version="dd-Mon-yy" #eg 21-Jul-18
```
**Note** Don't forget the quotes! This file is interpereted by python.

* Modify your `.bashrc` or `.bash_profile` file and add `scgid-0.1b/bin` to your enviornmental $PATH variable. For instance, add a line like... 
	`export PATH=$PATH:/path/to/scgid-0.1b/bin`
* Ensure that other stand-alone dependencies (i.e. BLAST and Augustus) have also been added to $PATH.

### Running ***scgid***
Each module of *scgid* is designed to be run separately in a bash command line. To take full advantage of *scgid*'s consensus-based approach, run all three binning algorithms (gc-cov, kmers, codons) prior to running `scgid consesnsus` to determine your final genome draft. The basic workflow for a *scgid* run is as follows...
```
scgid gc-cov [args...]  
scgid codons [args...]
scgid kmers train [args...]
scgid kmers annotate [args...]
scgid kmers extract [args...]
```
**Note** Because of the need to manually select a region from the self-organizing kmer map, the `scgid kmers...` portion of the workflow is divided into three separate commands. You will select your region of interest using the `esomana` GUI in between calls to `scgid kmers annotate` and `scgid kmers extract`.

Finally, run the consensus portion of `scgid` to draw consensus between all three binning algorithms.

`scgid consensus [args...]`

### What is a "swissprot-style database" and how do I know I have one?
**Preface** - It is unlikely that you will ever run into an issue with this as long as you're working with the database downloaded by `./scgid init` and only update/edit it and its associated taxonomy database with the included scripts. So don't worry about this too much.

**If you like to worry** - *scigd* is currently only compatible with a swissprot-style protein database. What that means is that the fasta headers (ie everything after '>' in the fasta file) have to share some elements with the standard swissprot format. Namely, each fasta header must contain a unique sequence identifier or USID (first) and a description (second). The USID and the description must be separated by a space. The USID cannot have any spaces in it, but the description can have as many as you like. Furthermore, the description must have a species designation in the format `OS=<species>. The species name can have spaces but must occur as the end of the description (and header). *scgid* uses this species designation to link sequences in the protein database with lineage information in the taxonomy database, so this is very important. You will be stopped if you try to use a misformatted database, but the source of these errors can be a pain to locate and correct. 

So, it is recommended that you let the included scripts do the work when making updates or edits to the swissprot-style database. The only time that you should be aware of these formatting requirements is when you want to manually add protein sequences to the database. So, here are some examples of GOOD and BAD ways to format your fasta headers for inclusion in the *scgid*-linked swissprot-style database.

```
GOOD HEADERS
----------------------------------------
>USID[:space:]description OS=species
>unique|species|identifier|without|spaces description can have as many spaces as you like OS=Homo sapiens
>myprotdb|X763901|mb243yy8 this is a protein that exists in a bacterium and does something OS=Escherichia coli
>cell76|MG7584x|8756|p8 Manganese peroxidase m2 splice variant a OS=Alicyclobacillus sp. strain G47
>x7645361AB2 hypothetical protein OS=Zooin

BAD HEADERS
----------------------------------------
>A765|6475882|p18 Cyctochrome C subunit from Medicago 
>M1 G74 X19 p55 Vaculoar protein of unknown function
>Russula_brevipes_contig1
>cell18|X6547a|p4 OS=Bacillus sp. some kinda of lactase
```

###Expanding your swissprot-style database
scgid was developed while working with cryptic and uncultured early-diverging fungi that are under-represented in sequence databases, including swissprot. Because of this, we saw it beneficial to expand the swissprot database to include published draft genomes of other early-diverging fungi that weren't included in the curated database prior to running scgid. While there are a couple of potential issues with this, it allowed us to be more confident in the ability to taxonomically confirm the identities of our target contigs. *scgid* allows you to do this too!

**You will need two files:**
1) FASTA files of the proteomes or sets of proteins you want to add the database. 

*Note* You can either append the applicable "OS=<species>" portion of the sequence description to each header manually (see above section "What is a "swissprot-style database" and how do I know I have one?") or use the provided `bin/reformat_my_headers.py` script to append them for you (see `reformat_my_headers.py -h` for details)

2) A two-column, tab-separated list of the OS's (column 1) and semicolon-separated lineage information. Here's an example:

```
Stylopage hadra	Eukaryota;Fungi;Zoopagomycota;Zoopages;Stylopage
Escherichia coli	Bacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia
Sacce	Fungi;Ascomycota;Saccharomycetales;Saccharomycetaceae;Sacharomyces
etc..
```

### Output Directories and Content
*scgid* makes a separate folder for each of its modules, all contained within a shared working directory named with the `-f|--prefix` command line option. This means that if you want to keep the outputs of each module in the same working directory (which is what you should do), make sure to specify the same `-f|--prefix` for each and every call to scgid when working on a your workflow. 

For your reference, I'm going to go through the content of the output folders for each module (excluding the prefices).The names are intended to be intuitive.

```
<prefix>_scgid_output/
	./gc-cov/
		1) aug.out.fasta** -> protein sequences of gene models called by augustus in FASTA format
		
		2) aug.out.gff3 -> gene models called by augustus 
		
		3) final_genome.fasta -> the final draft genome based on the gc-coverage-based genome selection process in FASTA format
		
		4) final_nontarget_bin.fasta -> all of the contigs **thrown out** based on the gc-coverage-based genome selection process in FASTA format
		
		5) gc_coverage_plots.pdf -> PDF file containing static graphical representations of the gc-coverage plots and the best window used to select unclassified contigs
		
		6) info_table.tsv -> Tab-deliminted table of information on contigs and their proteins that got an annotation based on blastp searches of the swissprot-style database. The columns of this file are as follows: 
    (i) contig_name, (ii) protein_length, (iii) contig_coverage, (iv) contig_gc, (v) protein_id, (vi) top_hit_species_id, (vii) top_hit_lineage, (viii) evalue, (ix) parsed_taxonomy [ie target/nontarget]
    
		7) unclassified_info_table.tsv -> Tab delimited table of information on contigs that did not get an annotation based on blastp searches of the swissprot-style database. The columns of this file are as follows:
    (i) contig_name, (ii) contig_coverage, (iii) contig_gc, (iv) taxonomy [ie Unclassified]
    
		8) spdb.blast.out** -> the output file resulting from a blastp search of the swussprot-style database, tab-delimited 
		9) spdb.blast.out.best -> the filtered "best" blast hits from the blast output file, tab delimited
		
		10) spdb.blast.out.parsed -> the "parsed, best" blast hits from the best blast output file, annotated with their swissprot-style fasta headers
		
		11) windows.all.out -> Tab-delimited file containing information of all twelve gc-coverage selection windows generated by scgid gc-cov. The columns of this file are as follows:
    (i) method, (ii) p_target, (iii) p_nontarget, (iv) gc_window, (v) gc_window_width, (vi) cov_window (vii) cov_window_width
    
		12) exluded_by_sp_taxonomy -> a list of the contigs that were excluded by scgid gc-cov based solely on their taxonomy hit

	./kmers/
		<put file descriptions here>
	./codons/
		<put file descriptions here>
```

### Command Line Arguments, explained

### Citations
