# GRAVE (Genomic Resistance Annotation and Visualization Environment)

Indian Copyright (Diary number) - 34038/2024-CO/SW

## Genomic Resistance Annotation and Visualization Environment
Current Version: 1.0 Tested on: Ubuntu 20.04.6 LTS (Linux 5.15.0-124-generic)

Antimicrobial resistance (AMR) represents a critical global challenge, necessitating advanced computational tools for the analysis and visualization of antimicrobial resistance genes (ARGs). While existing methods can process large-scale datasets, there remains a shortage of Linux-based solutions capable of consistently analyzing multiple samples while integrating statistical insights into visual outputs. To bridge this gap, we have developed the GRAVE-Genomic Resistance Annotation and Visualization Environment, a Linux-based platform designed to streamline the assembly, annotation, and visualization of ARGs from both shotgun metagenomic and whole-genome sequencing (WGS) data. GRAVE integrates a suite of scripts tailored for various stages of analysis. It employs metaSPAdes for assembling shotgun metagenomics reads and Unicycler for WGS reads, ensuring a standardized workflow. The tool translates nucleotide sequences into amino acid sequences using Prodigal, facilitating subsequent protein BLAST searches against a customized CARD database to identify genes, antibiotic classes, and mechanisms. GRAVE quantifies the abundances of these ARG elements and consolidates the data into comprehensive visualizations, including heatmaps, Principal Coordinates Analysis (PCoA) plots, correlation matrices, alpha diversity indices, and core gene representations. Furthermore, GRAVE generates stacked bar plots to illustrate the distribution of ARG classes and mechanisms based on metadata, and uses spider web plots to compare mechanisms across samples. By incorporating statistical analyses such as Welch’s t-test, GRAVE highlights significant variations between sample groups, enhancing interpretability of complex datasets. Its robust capacity for handling extensive datasets, combined with its ability to produce detailed statistical and graphical outputs, positions GRAVE as an invaluable tool for AMR research, enabling deeper insights into ARG dynamics across diverse environmental and clinical settings.

GRAVE developed by three developers equally
1) Parthkumar Prajapati
2) Neelam Nathani
3) Chandrashekar Mootapally


1_GRAVE_nASS.py

This script is designed for the assembly of multiple paired-end raw reads from sequencing data. For shotgun metagenomics data, it utilizes the metaSPAdes tool (Nurk et al., 2017), while Unicycler (Wick et al., 2017) is employed for whole-genome sequencing (WGS) assemblies. The input is provided through a manifest.csv (Manifest-1) file with three columns: paths for R1 reads, paths for R2 reads, and path of intended output directories. All raw reads are processed and assembled with a standardized command, ensuring consistent output across samples. Also, providing log file with list of assembled samples and its respective folder size as given.
   
   Dependencies required

 1) unicycler v0.5.0 (conda install -c bioconda unicycler=0.5.0) 
 2) spades v3.15.5 (conda install -c bioconda spades=3.15.5)

   Inputs required

manifest.csv [first column of the file containing path of R1 (Forward reads), second column containing path of R2 (Reverse reads), third column containing intended output path]
   
        python3 1_GRAVE_nAss.py --manifest_file MANIFEST_FILE --threads THREADS --type {WGS or metagenome} --log_file LOG_FILE  
            
2_GRAVE_keeper.py
          
After assembly, user will have to rename all output fasta files according to their sample IDs and consolidate them into a single directory manually. This script then translates nucleotide sequences (.fasta) into amino acid sequences (.faa) using Prodigal (Hyatt et al., 2010). The script accepts a manifest.csv (Manifest-2) file, where the first column specifies input file paths and the second column specifies output paths. It supports both “meta” mode for metagenomic assemblies and “single” mode for WGS assemblies, allowing translation with a single command. Also providing log file contain list of samples processed along with status, size of input fasta file and size of output faa file.
      
Dependencies required

1) prodigal V2.6.3 (conda install -c bioconda prodigal=2.6.3) 

Inputs required 
        
manifest.csv [first column of the file containing name of (Contig) FASTA file, second column containing name of required FAA file (Predicted protein assembly), no header required]

         python3 2_GRAVE_keeper.py --manifest MANIFEST_FILE --type {WGS or metagenome} --input_dir INPUT_DIR --output_dir OUTPUT_DIR 

3_GRAVE_yard.py
          
GRAVE_yard.py processes translated protein sequences for further analysis. It requires a manifest.csv (Manifest-3) file, input and output directories, metadata file (in .csv), and database paths for ARG, class, and mechanism. Using DIAMOND blastP (Buchfink et al., 2015), it searches against modified CARD 3.3.0 database(s) with a 60% identity threshold and an e-value cutoff of 0.00001. The resulting BLAST outputs of ARG, Class, and Mechanism are quantified as abundance hits for individual ARGs, class, and mechanism, producing SampleID_abund.csv files for each of the samples, which are further compiled into a combined_arg/class/mech.tsv file. This consolidated file is used for constructing heatmaps of top 50 ARGs with high abundance, class and dominant mechanisms and their respective tables in .csv file. Correlation coefficient matrix plot of ARGs and class using Karl pearson’s method will be generated along with the raw table of r value (correlation coefficient) for ARGs and Class; Alpha diversity will be computed as Shannon index for ARG and Class, PCoA (PERMANOVA) for the predicted ARGs and Class and a heatmap of core ARGs and Class will also be generated for visualization.
      
Dependencies required
     
1) diamond v2.1.8.162 (conda install -c bioconda diamond=2.1.8.162) 

Inputs required 

1) manifest.csv [first column of the file containing name of (Predicted protein assembly)   FAA file, second column contain name of output file (SampleID_Blast.csv), no header required]
2) metadata.csv [first column containing SampleID, and other columns with respective metadata of each sample, header required]
3) gene_db (CARD_arg_102024.dmnd), 
   class_db (CARD_class_102024.dmnd), 
   mech_db (CARD_mech_102024.dmnd), all three databases are trained database to be used as input in this script
     
        python3 3_GRAVE_yard.py --manifest_file MANIFEST_FILE --input_dir INPUT_DIR --input_dir --output_dir OUTPUT_ARG_DIR --class_output_dir CLASS_OUTPUT_DIR --mech_output_dir MECH_OUTPUT_DIR --metadata_file METADATA_FILE --gene_db GENE_DB --class_db CLASS_DB --mech_db MECH_DB --threads THREADS

4_GRAVE_filler.py
          
This script is written to process BLAST output files containing either class or mechanism data. In some instances, certain ARG classes and mechanisms are categorized in the database as exhibiting more than two antibiotic class name in one class (for classes) or employing multiple mechanisms strategies (for mechanisms). This categorization results in longer names that may be challenging to accommodate within the generated visualizations. To address this issue, the GRAVE_filler tool can be employed. This tool takes the blast.csv file as input and standardizes class or mechanism names that contain semicolons (";") by converting them into a specified type, such as "Multiclass_resistance" or "multi_mode". Subsequently, users can utilize the GRAVE_yard tool, specifying the output directory from GRAVE_filler as the input for mechanism processing. By using the --skip_abundance command argument, users can re-generate all visualizations with the updated files.

Inputs file required

manifest.csv [first column containing names of blast output file (SampleID_Blast.csv) file, second column containing same file names as mentioned in column one, no header required]
     
        python3 4_GRAVE_filler.py --manifest_file MANIFEST_FILE --replace_option {multiclass_resistant or multimode} --input_dir INPUT_DIR --output_dir OUTPUT_DIR

5_GRAVE_digger.py
          
This script constructs stacked bar plots to represent the composition of ARG classes and resistance mechanisms within sample groups. It uses metadata to group the samples, normalizes mean abundance values, converts them into percentages, and generates visual representations of the ARG distributions.

Inputs file required 

1) metadata.csv [first column containing SampleID, and other columns with respective metadata of each sample, header required]
2) combined_abundance.tsv (output of class and mechanisms of GRAVE_yard)
     
        python3 5_GRAVE_digger.py --metadata_file METADATA_FILE --type {Class or Mechanism} --compiled_abundance_file COMPILED_ABUNDANCE_FILE --output_dir OUTPUT_DIR --colour COLOUR[optional]

6_GRAVE_web.py
          
GRAVE_web.py uses the given combined_abundance.tsv file of either class or mechanism to generate spider web plots, visualizing the abundance of each class or mechanism of ARG across samples. It converts abundance values into percentages, creating a graphical display that highlights the trend in resistance mechanisms between samples.

Inputs file required 

1) combined_abundance.tsv (output of mechanisms of GRAVE_yard)

        python3 6_GRAVE_web.py --compiled_abundance_file COMPILED_ABUNDANCE_FILE --output_dir OUTPUT_DIR


References

For Unicycler (WGS) - 

      Wick RR, Judd LM, Gorrie CL, Holt KE (2017) Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 13(6): e1005595. https://doi.org/10.1371/journal.pcbi.1005595

For Prodigal (Protein prediction) - 

      Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119

For metaSPAdes (Shotgun) - 

      Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research, 27(5), 824-834.

For CARD database -

      Alcock, B. P., Raphenya, A. R., Lau, T. T., Tsang, K. K., Bouchard, M., Edalatmand, A., ... & McArthur, A. G. (2020). CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. Nucleic acids research, 48(D1), D517-D525.

For DIAMOND (BlastP) - 

      Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
