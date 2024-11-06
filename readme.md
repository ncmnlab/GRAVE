GRAVE (Genomic Resistance Annotation and Visualization Environment)

Copyright Diary Number - 34038/2024-CO/SW

Description
Antimicrobial resistance (AMR) represents a critical global challenge, necessitating advanced computational tools for the analysis and visualization of antimicrobial resistance genes (ARGs). While existing methods can process large-scale datasets, there remains a shortage of Linux-based solutions capable of consistently analyzing multiple samples while integrating statistical insights into visual outputs. To bridge this gap, we have developed the GRAVE-Genomic Resistance Annotation and Visualization Environment, a Linux-based platform designed to streamline the assembly, annotation, and visualization of ARGs from both shotgun metagenomic and whole-genome sequencing (WGS) data. GRAVE integrates a suite of scripts tailored for various stages of analysis. It employs metaSPAdes for assembling shotgun sequencing reads and Unicycler for WGS data, ensuring a standardized workflow. The tool translates nucleotide sequences into amino acid sequences using Prodigal, facilitating subsequent protein BLAST searches against a customized CARD database to identify genes, resistance classes, and mechanisms. GRAVE quantifies the abundances of these ARG elements and consolidates the data into comprehensive visualizations, including heatmaps, Principal Coordinates Analysis (PCoA) plots, correlation matrices, alpha diversity indices, and core gene representations. Furthermore, GRAVE generates stacked bar plots to illustrate the distribution of ARG classes and mechanisms based on metadata, and uses spider web plots to compare mechanisms across samples. By incorporating statistical analyses such as Welch’s t-test, GRAVE highlights significant variations between sample groups, enhancing interpretability of complex datasets. Its robust capacity for handling extensive datasets, combined with its ability to produce detailed statistical and graphical outputs, positions GRAVE as an invaluable tool for AMR research, enabling deeper insights into ARG dynamics across diverse environmental and clinical settings.

GRAVE developed by three developers equally- 
1) Parthkumar Prajapati 
2) Neelam Nathani
3) Chandrashekar Mootapally

1_GRAVE_nASS.py
          This script is designed for the assembly of multiple paired-end raw reads from sequencing data. For shotgun metagenomic data, it utilizes the metaSPAdes tool (Nurk et al., 2017), while Unicycler (Wick et al., 2017) is employed for whole-genome sequencing (WGS) assemblies. The input is provided through a manifest.csv (Manifest-1) file with three columns: paths for R1 reads, paths for R2 reads, and path of intended output directories. All raw reads are processed and assembled with a standardized command, ensuring consistent output across samples. Also, providing log file with list of assembled samples and its respective folder size as given (Annexure 2.1).
     Dependencies required
     1) import os
     2) import pandas as pd
     3) unicycler v0.5.0 (conda install -c bioconda unicycler=0.5.0) 
     4) spades v3.15.5 (conda install -c bioconda spades=3.15.5)

     Inputs required 
            manifest.csv [first column of the file containing path of R1 (Forward reads), second column contain path of R2 (Reverse reads), third column contain intended output path]
     Command
        python3 1_GRAVE_nAss.py --manifest_file MANIFEST_FILE --threads THREADS --type {WGS or metagenome} --log_file LOG_FILE  
            
2_GRAVE_keeper.py
          After assembly, user will have to rename all output fasta files according to their sample IDs and consolidates them into a single directory manually. This script then translates nucleotide sequences (.fasta) into amino acid sequences (.faa) using Prodigal (Hyatt et al., 2010). The script accepts a manifest.csv (Manifest-2) where the first column specifies input file paths and the second column specifies output paths. It supports both “meta” mode for metagenomic assemblies and “single” mode for WGS assemblies, allowing translation with a single command. Also providing log file contain list of samples processed along with status, size of input fasta file and size of output faa file (Annexure 2.2).
     Dependencies required
     1) import os
     2) import subprocess
     3) prodigal V2.6.3 (conda install -c bioconda prodigal=2.6.3) 

     Inputs required 
            manifest.csv [first column of the file containing name of (Contig) FASTA file, second column contain name of required FAA file (Predicted protein assembly), no header require]
     Command
        python3 2_GRAVE_keeper.py --manifest MANIFEST_FILE --type {WGS or metagenome} --input_dir INPUT_DIR --output_dir OUTPUT_DIR 

3_GRAVE_yard.py
          GRAVE_yard.py processes translated protein sequences for further analysis. It requires a manifest.csv (Manifest-3), input and output directories, metadata, and database paths for ARG, classes, and mechanisms. Using DIAMOND blastP (Buchfink et al., 2015), it searches against modified CARD 3.3.0 databases with a 60% identity threshold and an e-value cutoff of 0.00001. The resulting BLAST outputs of ARG, Class, and Mechanisms (Annexure 2.3A, 2.3B, 2.3C, respectively) are quantified ARGs, class, and mechanism abundances, producing SampleID_abund.csv (Annexure 2.3D, 2.3E, 2.3F, respectively) files, which are compiled into a combined_arg/class/mech.tsv (Annexure 2.4A, 2.4B, 2.4C) file. This consolidated file is used for constructing heatmaps of top 50 ARGs with high abundance (Annexure 2.5A), class (Annexure 2.5B) and dominant mechanisms (Annexure 2.5C) and their respective table in .csv file (Annexure 2.5a,2.5b, 2.5c, respectively) and correlation coefficient matrix plot of ARGs (Annexure 2.6A) and class (Annexure 2.6B) by Karl pearson’s method filter with p < 0.05 and also giving table for ARGs and Class (Annexure 2.6C, 2.6D respectively) of with r value (correlation coefficient) without any filter, Alpha diversity by Shannon index of ARG (Annexure 2.7A) and Class (Annexure 2.7B), PCoA (PERMANOVA) of ARG (Annexure 2.7A) and Class (Annexure 2.7B) and core ARGs (Annexure 2.8A)  and Class (Annexure 2.8B) for visualization.
     Dependencies required
     01) import os
     02) import subprocess
     03) import csv
     04) from collections import Counter
     05) import glob
     06) import pandas as pd
     07) import seaborn as sns
     08) import matplotlib.pyplot as plt
     09) import numpy as np
     10) from matplotlib.patches import Patch
     11) import plotly.express as px
     12) import plotly.figure_factory as ff
     13) import plotly.graph_objects as go
     14) from skbio.stats.ordination import pcoa
     15) from skbio.stats.distance import permanova, DistanceMatrix
     16) from scipy.spatial.distance import pdist, squareform
     17) import scipy.stats as stats
     18) from scipy.stats import pearsonr
     19) from skbio.diversity import alpha_diversity
     20) from scipy.stats import ttest_ind
     21) import time
     22) import sys
     23) diamond v2.1.8.162 (conda install -c bioconda diamond=2.1.8.162) 

     Inputs required 
        1) manifest.csv [first column of the file containing name of (Predicted protein assembly)   FAA file, second column contain name of output file (SampleID_Blast.csv), no header require]
        2) metadata.csv [first column contain SampleID, and other columns with respective metadata of each sample, header require]
        3) gene_db (CARD_arg_102024.dmnd), 
           class_db (CARD_class_102024.dmnd), 
           mech_db (CARD_mech_102024.dmnd), all three databases are trained database to be used as input in this script
     Command
        python3 3_GRAVE_yard.py --manifest_file MANIFEST_FILE --input_dir INPUT_DIR --input_dir --output_dir OUTPUT_ARG_DIR --class_output_dir CLASS_OUTPUT_DIR --mech_output_dir MECH_OUTPUT_DIR --metadata_file METADATA_FILE --gene_db GENE_DB --class_db CLASS_DB --mech_db MECH_DB --threads THREADS

4_GRAVE_filler.py
          This script is written to process BLAST output files containing either class or mechanism data. In some instances, certain ARG classes and mechanisms are categorized in the database as exhibiting more than two antibiotic class name in one class (for classes) or employing multiple mechanisms strategies (for mechanisms). This categorization results in longer names that may be challenging to accommodate within the generated visualizations. To address this issue, the GRAVE_filler tool can be employed. This tool takes the blast.csv file (Annexure 2.2A) as input and standardizes class or mechanism names that contain semicolons (";") by converting them into a specified type, such as "Multiclass_resistance" or "multi_mode" (Annexure 2.11A). Subsequently, users can utilize the GRAVE_yard tool, specifying the output directory from GRAVE_filler as the input for mechanism processing. By using the --skip_abundance command argument, users can generate all visualizations without encountering issues related to the display of lengthy class or mechanism names.
     Dependencies required
     1) import os
     2) import pandas as pd
     3) import matplotlib.pyplot as plt
     4) import seaborn as sns
     5) from matplotlib.patches import Patch
     6) import time 

     Inputs file required 
            manifest.csv [first column column contain name of blast file (SampleID_Blast.csv) file, second column contain same file names as mentioned in column one, no header require]
     Command
        python3 4_GRAVE_filler.py --manifest_file MANIFEST_FILE --replace_option {multiclass_resistant or multimode} --input_dir INPUT_DIR --output_dir OUTPUT_DIR

5_GRAVE_digger.py
          This script constructs stacked bar plots to represent the composition of ARG classes (Annexure 2.9A) and resistance mechanisms (Annexure 2.9B) within sample groups. It uses metadata to group the samples, normalizes mean abundance values, converts them into percentages, and generates visual representations of the ARG distributions.
     Dependencies required
     1) import os
     2) import csv 

     Inputs file required 
     1) metadata.csv [first column contain SampleID, and other columns with respective metadata of each sample, header require]
     2) combined_abundance.tsv (output of class and mechanisms of GRAVE_yard)
     Command
        python3 5_GRAVE_digger.py --metadata_file METADATA_FILE --type {Class or Mechanism} --compiled_abundance_file COMPILED_ABUNDANCE_FILE --output_dir OUTPUT_DIR --colour COLOUR[optional]

6_GRAVE_web.py
          GRAVE_web.py uses the given combined_abundance.tsv file of either class or mechanism to generate spider web plots, visualizing the abundance of each class or mechanism of ARG across samples. It converts abundance values into percentages, creating a graphical display (Annexure 2.10) that highlights the variations in resistance mechanisms between samples.
     Dependencies required
     1) import pandas as pd
     2) import numpy as np
     3) import matplotlib.pyplot as plt
     4) from math import pi

     Inputs file required 
     1) combined_abundance.tsv (output of mechanisms of GRAVE_yard)
     Command
        python3 6_GRAVE_web.py --compiled_abundance_file COMPILED_ABUNDANCE_FILE --output_dir OUTPUT_DIR


References
For Unicycler (WGS) - 

	Wick RR, Judd LM, Gorrie CL, Holt KE (2017) Unicycler: Resolving bacterial genome assemblies from short and long sequencing reads. PLoS Comput Biol 13(6): e1005595. https://doi.org/10.1371/journal.pcbi.1005595

For Prodigal (Protein prediction) - 

	Hyatt, D., Chen, GL., LoCascio, P.F. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics 11, 119 (2010). https://doi.org/10.1186/1471-2105-11-119

For metaSPAdes (Shotgun)-

	Nurk, S., Meleshko, D., Korobeynikov, A., & Pevzner, P. A. (2017). metaSPAdes: a new versatile metagenomic assembler. Genome research, 27(5), 824-834.
