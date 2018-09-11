## ***scgid***, a python-based tool for scaffold binning and genome prediction from single-cell sequencing libraries
### version 0.1b
-----

### What is ***scgid*** for?
*scgid* is a python-based tool aimed at extracting a draft genome sequence for a single organism from a mildly metagenomic sequencing library resulting from single-cell genomic amplifications. The only thing that you need to start is your assebmly in FASTA format and *scgid* will do the rest.

*scgid* takes your SPAdes assebmly and subjects it to three binning methods each based on a different sequence signature. It takes the output of each, draws consensus based on majority-rule, and yields one final consensus-genome draft at the intersection of all three methods. 

As of this version, *scgid* is only compatible with genome assemblies generated with SPAdes (http://bioinf.spbau.ru/en/spades_for_remove) that have their sequence fasta headers **unchanged** (this is very important - just rename them after you're done with *scgid*). 

Please note that this version of *scgid* constitutes an early-release beta version. Some aspects of it and its documentation may be incomplete and/or under development. We actively support *scgid* and are working to expand and test it. If you encounter an error or obstacle when running *scgid*, please open an issue on this github repositiory so that we can address it and get you up and running again. See the next section below for some more notes about this version of *scgid*.

### A couple of notes about the beta version, 0.1b
This pre-publication release version of *scgid* is still under active development and thus updates to the GitHub repository are being made quite consistently. To ensure that you are able to take advantage of new features and bug fixes, make sure that you are updating your local repo consistently. We've tried to make this easy for you.

Type `scgid update` and confirm to automatically pull the current source from this repo.

One of the most obvious inconsistencies that we still need to correct involve some naming differences in between the module calls and the file architecture of the output folder. This will be fixed in the next version, but in the mean time, here are the synonyms:  
```
scgid kmers [args...]		Outputs to <prefix>_scgid_output/esom
scgid gc-cov [args...]		Outputs to <prefix>_scgid_output/blob
scgid codons [args...]		Outputs to <prefix>_scgid_output/rscu
scgid consensus [args...]	Outputs to <prefix>_scgid_output/consensus
```

In this version of *scgid*,  
	(i) `scgid kmers` can only be annotated with blastn searches of the NCBI nt database. This will be expanded in later versions. `scgid kmers` will perform this BLAST search for you.  
	(ii) `scgid gc-cov` results can be annotated with blastp searches against a local swissprot-style protein database (see below section for what exactly this means. `scgid gc-cov` will perform this BLAST search for you.  
	(iii) `scgid codons` results can be annotated with either blastn searches against the NCBI nt database (default, done for you) OR blastp searches against a local swissprot-style protein database. This latter feature is still in the process of being fully integrated into *scgid*. We found that we got a higher proportion of tree leaf annotations using blastp searches. For the time being, `scgid codons` will not perform the blastp search for you. Instead, run `scgid gc-cov` on your assembly first and then provide the `blob/<prefix>_info_table.tsv` file to `scgid codons` using the `-i|--infotable` option and specify `--mode blastp`. Doing so will allow `scgid codons` to use the results of the blastp search against the local swissprot-style database.  

This version of *scgid* has been tested on MacOS X and several Linux distributions (CentOS, Ubuntu, and Arch).  

Please report any and all errors and issues that you run into while using *scgid*. Its development is a top priority for me at this time and I want it to work for you! Please open issues here on the GitHub page or [email me directly](amsesk@umich.edu).  

### Dependencies
* python 2.7 (*scgid* is not compatbile with python 3.x.x)
	* pandas 1.15.0
	* numpy 0.23.4
	* ETE3 toolkit 3.1.1  
**Note** *scgid* will try to install these python dependencies if using automatic setup (recommended) below.
* R 3.4.1
	* ape 
	`install.packages("ape")`
	* Rscript (included with most R distributions)
* NCBI Blast+
* Augustus (http://bioinf.uni-greifswald.de/augustus/)
* ClaMS-CLI (http://clams.jgi-psf.org/)  
    (i) Scroll to the bottom of the page a and click on "Download",  
    (ii) Accept the license agreement,  
    (iii) Download the "Command line verion" of ClaMS  
* Java
* Databionic ESOM (http://www.databionin-esom.sourceforge.net)
    (i) Download the ESOM Installer .jar file and follow the instructions in the GUI installer.
* BASH

### Downloading and Installing ***scgid***
* Clone the repository at [https://github.com/amsesk/scgid.git](https://github.com/amsesk/scgid.git)
* Navigate to `bin` within the newly created directory
`cd scgid/bin`

#### Automatic Setup (recommended)
* Type `./scgid init`
	* Follow the prompts to define some package-wide settings and download the necessary databases.
	* This script **will try** to install the required python packages, if they are not already available, using `conda` or `pip`. Make sure you know which package manager is being used on your system as installation via `pip` has been known to break some conda installations. If the flavor of python package manager on your system is up to you, I strongly recommend using `conda` since it will also install important module dependencies for you. You can download miniconda (here)[https://conda.io/miniconda.html].
	* You are responsible for having downloaded and installed all of the other third-party dependencies listed above. You will be asked where some of them are.
	* This script requires an internet connection in order to download the swissprot protein and taxonomy databases.
	
* Modify your `.bashrc` or `.bash_profile` file and add `scgid/bin` to your enviornmental $PATH variable. For instance, add a line like... `export PATH=$PATH:/path/to/scgid/bin`
* Ensure that other stand-alone dependencies (i.e. BLAST and Augustus) have also been added to $PATH.

#### Manual Setup (if auto isn't working)
* Download and decompress a copy of the most recent swissprot database.
* Download and save the file located [here](https://www.uniprot.org/taxonomy/?query=*&format=tab) and then run `./scgid buildtaxdb [args...]`. See `./scgid buildtaxdb -h` for more information.
* Edit `settings.py` to reflect the locations of ESOM, ClaMS, and Rscript, and your swissprot databases, as follows:
```
esom_path ="<path_to_folder_containing_bin>"`
clams_path="<path_to_folder_containing_clams-cli_jar>"
path_to_Rscript="<path_to_folder_containing_Rscript>"
path_to_spdb="<path_to_swissprot_database>"
path_to_taxdb="<path_to_taxdb>"
taxonomy_all_tab="<path_to_taxonomy_all_tab>"
spdb_version="dd-Mon-yy" #eg 21-Jul-18
```
**Note** Don't forget the quotes! This file is interpereted by python.

* Modify your `.bashrc` or `.bash_profile` file and add `scgid-0.1b/bin` to your enviornmental $PATH variable. For instance, add a line like... 
	`export PATH=$PATH:/path/to/scgid-0.1b/bin`
* Ensure that other stand-alone dependencies (i.e. BLAST and Augustus) have also been added to $PATH.

### Running ***scgid***
To run *scgid*, all you need is a SPAdes assembly (or at least an assembly with SPAdes-style fasta headers). In its current version, SPAdes-style fasta headers are a requirement for *scgid*. This means that each fasta header contains identification, length, and coverage information for each contig in the format `NODE_XXX_length_XXX_cov_XXX.XXX`. If this is an issue for you please open a new issue and we'll try to expand compatibility in future versions.

Each module of *scgid* is designed to be run separately in a bash command line. To enumerate the command-line arguments and their descriptions, merely type `scgid <module> -h` or `scgid <module> --help`. Try running `scgid -h` to get descriptions of the available module commands to see where to get started. In most cases, *scgid* will try to pull options that you don't specify from your `settings.py` file, so keep in mind what is specified there (you set-up this file when you ran `./scgid init`). You only need to explicitly specify these options when you want to use databases NOT pointed to in your `settings.py` file.

**IMPORTANT** Output directories for runs are determined by the `-f|--prefix` option supplied to each call to the module. So, if you would like the outputs of all modules to be located in the same output head directory, make sure you are in the parent directory of the `<prefix>_scgid_output` folder and that you specify the same prefix for each call to *scgid*. This will ensure that the outputs of each module for each run (on a particular assembly) are located in the same output directory. Changing the prefix of the call to *scgid* will create a new output directory. Further, calling *scgid* from a different directory will create a new output folder in that directory. Also note that because of this, results can be overwritten. For instance, if you call `scgid gc-cov [args...]` from the same directory twice with the same prefix, the outputs contained within the `<prefix>_scgid_output/gc-cov` directory will be overwritten.

To take full advantage of *scgid*'s consensus-based approach, run all three binning algorithms (gc-cov, kmers, codons) prior to running `scgid consesnsus` to determine your final genome draft. The basic workflow for a *scgid* run is as follows...
```
scgid gc-cov [args...] 
scgid codons [args...]
scgid kmers train [args...]
scgid kmers annotate [args...]
scgid kmers extract [args...]
```
**Note** Because of the need to manually select a region from the self-organizing kmer map, the `scgid kmers...` portion of the workflow is divided into three separate commands. You will select your region of interest using the `esomana` GUI in between calls to `scgid kmers annotate` and `scgid kmers extract`. See below section for more information.

Finally, run the consensus portion of `scgid` to draw consensus between all three binning algorithms.

`scgid consensus [args...]`

The final consensus genome draft is now located at `<prefix>_scgid_output/consensus/<prefix>_final_genome.fasta` 

### I just ran `scgid kmers annotate [args...]`, now what?
Completion of this command means that you have successfully trained and taxonomically-annotated the ESOM topology for your metagenome. Now it is time to look at the annotated "map" and decide which sector you want to carve-out as your target organism. Follow these steps to do so (**pictures coming**):   
1. Locate and open `esomana` in `path/to/ESOM/bin/`  
2. From the "File" drop-down menu, select `Load *.wts`  
3. Navigate to the kmers output folder in the current scgid run, i.e. `path/to/<prefix>_scgid_output/kmers/` and open the `.wts` file that you find there.
4. You should now see an ESOM topology with colored dots. This is your annotated map.  
5. In the bottom of the window, click on the Classes tab. This should show you which of your taxonomic groups are represented by each color in the map.  
6. When you are ready to select a region of the map, use the mouse cursor to click a shape around your region of interest - this can take a while depending on the size of your map.  
7. Once you are satisfied with your region selection, close off the region by clicking both mouse buttons at the same time. **This will create a NEW class that incorporates all of the colored dots within your region selection.**  
8. Notice that your new class has a number associated with it. Remember this number.   
9. From the "File" drop-down menu, select `Save *.cls` and save your new class file.  
10. Return to command line and run `scgid kmers extract [args...]` making sure to specify the .cls file you just created in step 9 for `-c|--cls` and the number of your new class for `-cid|--classnum`.  
11. Now you should have a final draft genome from the kmers module in your `path/to/<prefix>_scgid_output/kmers` directory.

### What is a "swissprot-style database" and how do I know I have one?
**Preface** - It is unlikely that you will ever run into an issue with this as long as you're working with the database downloaded by `./scgid init` and only update or edit it with `scgid updatedb` or `scgid spexpand` respectively. So don't worry about this too much.

**If you like to worry** - *scigd* is currently only compatible with a swissprot-style protein database. What that means is that the fasta headers (ie everything after '>' in the fasta file) have to share some elements with the standard swissprot format. Namely, each fasta header must contain a unique sequence identifier or USID (first) and a description (second). The USID and the description must be separated by a space. The USID cannot have any spaces in it, but the description can have as many as you like. Furthermore, the description must have a species designation in the format `OS=<species>`. The species name can have spaces but should occur at the end of the description (and fasta header). *scgid* uses this species designation to link sequences in the protein database with lineage information in the taxonomy database, so this is very important. You will be stopped if you try to use a misformatted database, but the source of these errors can be a pain to locate and correct. 

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

### Expanding your swissprot-style database
scgid was developed while working with cryptic and uncultured early-diverging fungi that are under-represented in sequence databases, including swissprot. Because of this, we saw it beneficial to expand the swissprot database to include published draft genomes of other early-diverging fungi that weren't included in the curated database prior to running scgid. While there are a couple of potential issues with this, it allowed us to be more confident in the ability to taxonomically confirm the identities of our target contigs. *scgid* allows you to do this too!

Use `scgid spexpand [args...]` to expand your swissprot protein and taxonomy databases as you see fit.

**You will need two files:**
1) FASTA files of the proteomes or sets of proteins you want to add the database. 

*Note* You can either append the applicable "OS=<species>" portion of the sequence description to each header manually (see above section "What is a "swissprot-style database" and how do I know I have one?") or use the provided `bin/reformat_my_headers.py` script to append them for you (see `reformat_my_headers.py -h` for details)

2) A two-column, tab-separated list of the OS's (column 1) and semicolon-separated lineage information. Here's an example:

```
Stylopage hadra	Eukaryota;Fungi;Zoopagomycota;Zoopagales;Stylopage
Escherichia coli	Bacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia
Sacce	Fungi;Ascomycota;Saccharomycetales;Saccharomycetaceae;Sacharomyces
etc..
```

Now all you have to do is run `scgid spexpand [args...]` remembering to provide File #1 as `-p|--proteins` and File #2 as `-l|--lineages` and you're all set to go.

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
**\<In progress\>**  
For now, type `scgid <module> -h` or `scgid <module> --help` into your command line for descriptions of supported command-line arguments.


### Citations
Altschup, S. F., Gish, W., Pennsylvania, T., & Park, U. (1990). Basic Local Alignment Search Tool 2Department of Computer Science, 403–410.  

Amses, K. R., Davis, W. J., & James, T. Y. scgid, a python-based tool for scaffold binning and genome prediction from single-cell sequencing libraries, (in prep).
  
Dick, G. J., Andersson, A. F., Baker, B. J., Simmons, S. L., Thomas, B. C., Yelton, A. P., & Banfield, J. F. (2009). Community-wide analysis of microbial genome sequence signatures, 10(8). https://doi.org/10.1186/gb-2009-10-8-r85  

Huerta-Cepas, J., Serra, F., and Bork, P. ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Mol Biol Evol 2016; doi: 10.1093/molbev/msw046
  
Kumar, S., Jones, M., Koutsovoulos, G., Clarke, M., & Blaxter, M. (2013). Blobology : exploring raw genome data for contaminants , symbionts , and parasites using taxon-annotated GC-coverage plots, 4(November), 1–12. https://doi.org/10.3389/fgene.2013.00237

McInerney, J. O. (1998). GCUA: general codon usage analysis. Bioinformatics, 14(4), 372–373. https://doi.org/10.1093/bioinformatics/14.4.372  
  
Mikhailov, K. V, Simdyanov, T. G., Aleoshin, V. V, & Belozersky, A. N. (2016). Genomic survey of a hyperparasitic microsporidian Amphiamblys sp. (Metchnikovellidae) Genome Biology and Evolution Advance Access. Genome Biology and Evolution, 9(3), 454–467. https://doi.org/10.1093/gbe/evw235  
  
Pati, A., Heath, L. S., Krypides, N. C., & Ivanova, N. (2011). ClaMS: A Classifier for Metagenomic Sequences. Standards in Genomic Sciences, 5, 248–253. https://doi.org/10.4056/sigs.2075298  

Stanke, M., & Morgenstern, B. (2005). AUGUSTUS: A web server for gene prediction in eukaryotes that allows user-defined constraints. Nucleic Acids Research, 33(SUPPL. 2), 465–467. https://doi.org/10.1093/nar/gki458
