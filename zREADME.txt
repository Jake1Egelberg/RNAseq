***INSTALLATION INSTRUCTIONS***

1) Install R (ideally version 4.0.0) and R Studio

2) Install the BiocManager, stringr, gplots, and RColorBrewer packages. Use BiocManager to install the Rsubread, edgeR, limma, Glimma,  org.Mm.eg.db, org.Hs.eg.db, fgsea, and reactome.db; see http://bioconductor.org/packages. For package installation, can run the 0PACKAGE_INSTALLATION.bat file in the "scripts" subfolder.

3) Download this RNAseq folder onto your hardrive (filepath = "D:") and rename folder to "RNAseq" (remove version identifiers)

4) Configure filepaths in R. Open all three .R scripts in the "scripts" subfolder in R Studio and ensure that the filepath at the top of each script matches where this "RNAseq" folder is saved.

5) Configure the .bat files in this "RNAseq" folder. Ensure that .bat opens the version of R installed on your computer and correctly references where your R scripts are saved.

------------------------------------------------------------

GENERAL STEPS OF RNA-SEQ:
1) Build index of reference genome
2) Align RNA reads to index
3) Measure differential expression between high quality reads

------------------------------------------------------------

1) Identify sequences to analyze from NCBI SRA (Sequence Read Archive)
	https://www.ncbi.nlm.nih.gov/sra/?term=hg19
*Filter by  Access=Public, Source=RNA, Library Layout=<Your choice>, File Type= .fastq OR .bam
**Order results by taxon to according to which organism you are investigating (human/hg19 or mouse/mm10)

2) Download FASTQ/BAM from European Nucleotide Archive (downloading BAM speeds up analysis)
	https://www.ebi.ac.uk/ena/browser/view
*Use "View" function to search by SRR from SRA

3) Import FASTQ/BAM to 1fastqfiles folder in workflow path
	See history folder for old FASTQ files
	For paired data, ensure one set is added to pair folder in same order

4) Identify reference genome and import .fa.gz file to 2genome folder
	https://hgdownload.soe.ucsc.edu/downloads.html
*For taxon="Homo sapiens," use genome hg19. For taxon="Mus musculus," use genome mm10

5) Record sample metadata and configure design matrix

6) Configure parms per metadata and design matrix

7) Run the following .bat files in order:
	1) 1BUILD_INDEX.bat *SKIP THIS IF YOU ALREADY HAVE A BUILT INDEX*
	2) 2ALIGN_READS.bat *SKIP THIS IF YOU DOWNLOADED BAM FILES*
	3) 3DE_ANALYSIS.bat
*THIS CAN TAKE A LONG TIME

8) Monitor program progress as follows:
	a) See "progress" folder for current processes
	b) Monitor the buildindex and 1fastqfiles folders to see addition of new files
	c) Check task manager to ensure that resources are being directed towards RNA-seq program (shows up as Windows Command Processor)
	d) Ensure that analysis is complete in "progress" folder before the command window closes automatically
	e) If analysis stops on 3DE_ANALYSIS, re-configure threshold and sample values to include more genes


------------------------------------------------------------

EXAMPLE FASTQ SEQUENCE ACCESSION #S TO START:
SRR1552444
SRR1552445
SRR1552446
SRR1552447
SRR1552448
SRR1552449
SRR1552450
SRR1552451
SRR1552452
SRR1552453
SRR1552454
SRR1552455

Fu, N.Y., Rios, A., Pal, B., Soetanto, R., Lun, A.T.L., Liu, K., Beck, T., Best, S., Vaillant, F., Bouillet, P., Strasser, A., Preiss, T., Smyth, G.K., Lindeman, G., , Visvader, J.: EGF-mediated induction of Mcl-1 at the switch to lactation is essential for alveolar cell survival. Nature Cell Biology 17(4), 365–375 (2015)

------------------------------------------------------------

ANNOTATED PARMS FILE

#Specify reference genome file name (remove ending, .fa.gz implied), file to be indexed
index.file<-"chr1_mm10"

#Set true if data is paired-end
paired.end.status<-FALSE

#Specify reference genome of samples and index file
  #Mouse (mm10) and human (hg19) genome assemblies
ref.genome<-"mm10"

#Specifies whether to use existing raw feature counts or not
use.existing.counts: FALSE

#Configure design matrix
  #Groups as column names (factor)
  #Rows are samples (row number should correspond to order of FASTQ)
  #Fill columns with 1 or 0 to indicate status (level)

#Specify group to compare (should be column name in design matrix)
  #Can only compare 2 levels within 1 factor
interest.group<-"typeluminal"

#Specify threshold and sample values for removing lowly expressed genes
thresh.value<-6 #CPM<X selected
sample.value<-2 #Across at least X samples

------------------------------------------------------------

CITE THIS WORKFLOW:
https://zenodo.org/badge/latestdoi/396408760



