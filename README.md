#### ***SCGid*** was updated to version 0.9b on April 22, 2020

-----

## SCGid: a consensus approach to contig filtering and genome prediction from single-cell sequencing libraries of uncultured eukaryotes
-----
version 0.9b

[Link to SCGid Publication](https://academic.oup.com/bioinformatics/article-abstract/36/7/1994/5640497)

### What is ***SCGid*** for?
*SCGid* is a python-based tool aimed at extracting a draft genome sequence for a single organism from a mildly metagenomic sequencing library resulting from single-cell genomic amplifications. The only thing that you need to start is your assebmly in FASTA format and *SCGid* will do the rest.

*SCGid* takes your nucleotide assebmly and subjects it to three binning methods each based on a different sequence signature. It takes the output of each, draws consensus based on majority-rule, and yields one final consensus-genome draft at the intersection of all three methods. 

*SCGid* is now compatible with any nucleotide assembly (SPAdes or other). If your input is a non-SPAdes assembly or a SPAdes asembly that has had the sequence headers modified, you will have to supply a contig coverage matrix to `scgid gct`. See `scgid gct -h` for more information.

Please note that this version of *scgid* constitutes an early-release beta version. Some aspects of it and its documentation may be incomplete and/or under development. We actively support *SCGid* and are working to expand and test it. If you encounter an error or obstacle when running *SCGid*, please open an issue on this github repositiory so that we can address it and get you up and running again. See the next section below for some more notes about this version of *SCGid*.

In this version of *SCGid*,  
	(i) `scgid kmers` can only be annotated with blastn searches of the NCBI nt database. This will be expanded in later versions. `scgid kmers` will perform this BLAST search for you.  
	(ii) `SCGid gc-cov` results can be annotated with blastp searches against a local swissprot-style protein database (see below section for what exactly this means. `scgid gc-cov` will perform this BLAST search for you.  
	(iii) `scgid codons` results can be annotated with either blastn searches against the NCBI nt database (default, done for you) OR blastp searches against a local swissprot-style protein database. We found that we got a higher proportion of tree leaf annotations using blastp searches. For the time being, `scgid codons` will not perform the blastp search for you, but can use the outputs of this search generated by `scgid gct` when using `--mode blastp` in `scgid codons`.

*SCGid* has been tested on MacOS X and several Linux distributions (CentOS, Ubuntu, and Arch).  

Please report any and all errors and issues that you run into while using *SCGid*. Its development is a top priority for me at this time and I want it to work for you! Please open issues here on the GitHub page.

### Dependencies
* python3 (*SCGid* is not compatbile with python 2.7)
	* pandas >= 0.23.4
	* numpy >= 1.15.0
	* ETE3 toolkit >= 3.1.1  
	* plotly
	* pyyaml
	* urllib3
	* importlib_metadata
	
**Note** *SCGid* will install the above python dependencies during installation.

* R 3.4.1
	* ape 
	* dplyr
    	* phytools
	* RColorBrewer
    	* Rscript (included with most R distributions)
* NCBI Blast+ 
* Augustus (http://bioinf.uni-greifswald.de/augustus/)
* ClaMS-CLI  
    (i) ClaMS-CLI is apparently no longer available vis the JGI website. I've forked the source in its original form with the original copyright notice and liscence to https://github.com/amsesk/ClaMS-CLI-fork.
    (ii) Clone it with `git clone https://github.com/amsesk/ClaMS-CLI-fork.git` 
* Java
* Databionic ESOM (http://www.databionin-esom.sourceforge.net)  
    (i) Download the ESOM Installer .jar file and follow the instructions in the GUI installer.

### Downloading and Installing ***SCGid***
```
git clone https://github.com/amsesk/SCGid.git
cd SCGid

# Running in a virtual environment is recommended
python -m venv scgidenv
source /path/to/scgidenv/bin/activate

# Use develop instead of install to make updates easier ahead of release version
python setup.py develop 

# Set some package-wide settings and download the necessary databases
scgid init

# Another pre-release workaround to make updates easier
cp config.yaml config.yaml.local

# Make sure that augustus and BLAST are installed, working, and available in $PATH by adding them to .bashrc
export PATH=$PATH:/path/to/augustus
export PATH=$PATH:/path/to/blast

# Ready to run SCGid
scgid --help
```

### Running ***SCGid***
To run *SCGid*, all you need is a nucleotide assembly.

Each module of *scgid* is designed to be run separately in a bash command line. To enumerate the command-line arguments and their descriptions, merely type `scgid <module> -h` or `scgid <module> --help`. Try running `scgid -h` to get descriptions of the available module commands to see where to get started. In most cases, *SCGid* will try to pull options that you don't specify from `config.yaml`, so keep in mind what is specified there. You only need to explicitly specify some options when you want to use databases NOT pointed to in `config.yaml`.

**IMPORTANT** Output directories for runs are determined by the `-f|--prefix` option supplied to each call to the module. So, if you would like the outputs of all modules to be located in the same output head directory, make sure you are in the parent directory of the `<prefix>_scgid_output` folder and that you specify the same prefix for each call to *SCGid*. This will ensure that the outputs of each module for each run (on a particular assembly) are located in the same output directory. Changing the prefix of the call to *scgid* will create a new output directory. Further, calling *SCGid* from a different directory will create a new output folder in that directory. Also note that because of this, results can be easily overwritten. For instance, if you call `scgid gct [args...]` from the same directory twice with the same prefix, the outputs contained within the `<prefix>_scgid_output/gct` directory will be overwritten.

To take full advantage of *SCGid*'s consensus-based approach, run all three binning algorithms (gc-cov, kmers, codons) prior to running `scgid consesnsus` to generate a consensus-filtered draft. The basic workflow for a *scgid* run is as follows...
```
scgid gct [args...] 
scgid codons [args...]
scgid kmers train [args...] (a good reference for which options to specify for ESOM training: https://github.com/tetramerFreqs/Binning)
scgid kmers annotate [args...]
scgid kmers extract [args...]
scgid consensus [args...]
```
**Note** Because of the need to manually select a region from the ESOM map, the `scgid kmers...` portion of the workflow is divided into three separate CLI calls. You will select your region of interest using the `esomana` GUI from Databioincs ESOM in between calls to `scgid kmers annotate` and `scgid kmers extract`. See below section for more information.

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
**Preface** - It is unlikely that you will ever run into an issue with this as long as you're working with the database downloaded by `./scgid init` and only update or edit it with `scgid spdbup` or `scgid spexpand` respectively. So don't worry about this too much.

**If you like to worry** - *SCGid* is currently only compatible with a swissprot-style protein database. What that means is that the fasta headers (ie everything after '>' in the fasta file) have to share some elements with the standard swissprot format. Namely, each fasta header must contain a unique sequence identifier or USID (first) and a description (second). The USID and the description must be separated by a space. The USID cannot have any spaces in it, but the description can have as many as you like. Furthermore, the description must have a species designation in the format `OS=<species>`. The species name can have spaces but should occur at the end of the description (and fasta header). *SCGid* uses this species designation to link sequences in the protein database with lineage information in the taxonomy database, so this is very important.

So, it is recommended that you let the included modules do the work when making updates or edits to the swissprot-style database. The only time that you should be aware of these formatting requirements is when you want to manually add protein sequences to the database. So, here are some examples of GOOD and BAD ways to format your fasta headers for inclusion in the *scgid*-linked swissprot-style database.

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
1) FASTA file of the proteomes or sets of proteins you want to add the database. 

2) A two-column, tab-separated list of the OS's (column 1) and semicolon-separated lineage information. Here's an example:

```
Stylopage hadra	Eukaryota;Fungi;Zoopagomycota;Zoopagales;Stylopage
Escherichia coli	Bacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia
Sacce	Fungi;Ascomycota;Saccharomycetales;Saccharomycetaceae;Sacharomyces
etc..
```

Now all you have to do is run `scgid spexpand [args...]` remembering to provide File #1 as `-p|--proteins` and File #2 as `-l|--lineages`.

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
