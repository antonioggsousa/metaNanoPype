### Table of contents

   + [Introduction](#introduction)

   + [Installation](#installation)

   + [Tutorial](#tutorial)
   
   + [Support or Contact](#support-or-contact)

   + [Acknowledgement](#acknowledgement)


<br>

---

<br>

### Introduction

**metaNanoPype** is a modular python based pipeline for reproducible analysis of nanopore metabarcoding data with the following modules: 

   + (I) demultiplexing (*not implemented yet*); 

   + (II) quality-assessment (*implemented*); 
   
   + (III) quality-filtering and trimming (*implemented*);  
   
   + (IV) polishing/read correction (*not implemented yet*); 
   
   + (V) taxonomic classification (*under development*); 
   
   + (VI) diversity analyses (alpha- and beta-diversity) (*not implemented yet*). 
   
Depending on the step, the module can include more than one option, i.e., tool/algorithm, to allow flexibility to the user. During each step is generated a log file with relevant information to build a report (an additional option) in html format describing the results, versions of software used as well scripts and references of the software.

<br>

---

<br>

### Installation

**metaNanoPype** is a modular pipeline relying on python scripts that are wrappers to other algorithm/tools that need to be installed in your OS. 

In general, it requires: 

   + **python** >= v.3.6.9 (*[installation](https://www.python.org/downloads/)*)
   
   + **metaNanoPype** python scripts: `git clone https://github.com/antonioggsousa/metaNanoPype.git`  

<br>

Then, for each module step it requires several third-party, stand-alone tools installed in your system: 

   + (I) demultiplexing (*not implemented yet*); 

   + (II) quality-assessment (*implemented*); 
      
      + **fastqc** v.0.11.5 (*[installation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)*)
      
      + **multiqc** v.1.9 (*[installation](https://multiqc.info/docs/#installation)*)
      
      + **NanoPlot** v.1.29.1 (*[installation](https://github.com/wdecoster/NanoPlot#installation)*)
   
   + (III) quality-filtering and trimming (*implemented*);  

      + **porechop** v.0.2.4 (*[installation](https://github.com/rrwick/Porechop#installation)*)
      
      + **NanoFilt** v.2.7.1 (*[installation](https://github.com/wdecoster/nanofilt#installation-and-upgrading)*)
   
   + (IV) polishing/read correction (*not implemented yet*); 
   
   + (V) taxonomic classification (*under development*); 

      + **kraken2** v.2.1.1 (*[installation](https://github.com/DerrickWood/kraken2/wiki/Manual#installation)*) 
   
   + (VI) diversity analyses (alpha- and beta-diversity) (*not implemented yet*). 


<br>

---

<br>

### Tutorial

For this quick start tutorial it will be used the publicly available data of nanopore full-length 16S rRNA amplicon sequences published in Low *et al*. 2021: [Evaluation of full-length nanopore 16S sequencing for detection of pathogens in microbial keratitis](https://peerj.com/articles/10778/). 

<br>

>This tutorial was performed in a server with the OS Ubuntu 18.04.5 LTS and 20 cores and ~120 Gb RAM. Generally many of the commands used can be reproduced in any Linux distribution or even UNIX OS. It assumes that the several dependencies are already installed and available in your `PATH`.    

<br>

0. Create the directory structure to reproduce the tutorial and scripts:

    Change into the directory where you want to reproduce the tutorial and create the directory structure:

        mkdir bin scripts data results report

    Download **metaNanoPype** scripts:

        cd bin 

        git clone https://github.com/antonioggsousa/metaNanoPype.git 

    Next, add the **metaNanoPype** bin folder to your **PATH**:

        export PATH="$PWD/metaNanoPype/bin:$PATH"

<br>

<br>

1. Download the full-length 16S rRNA amplicon fastq files: 

    The fastq files were deposited in ENA (European Nucleotide Archive) under the project accession number: [PRJEB37709](https://www.ebi.ac.uk/ena/browser/view/PRJEB37709). 

    The files can be downloaded from several ways. A convenient way is by using the ENA toolkit **enaBrowserTools**, such as: `enaGroupGet` command. Download them from [github](https://github.com/enasequence/enaBrowserTools) and download the data (follow the steps below).

        git clone https://github.com/enasequence/enaBrowserTools.git

    Next, add the **enaBrowserTools** python3 folder to your **PATH**:

        export PATH="$PWD/enaBrowserTools/python3:$PATH"

    Download the PRJEB37709 project fastq sequencing data (*comment - this will take a while*): 

        cd ../data #change into data folder first

        enaGroupGet -g read -f fastq PRJEB37709 

    Create a new folder with all individual fastq files inside and delete the previous folder (with each fastq file inside of a sample specific folder) to work more conveniently with the files: 

        mkdir fastq

        mv */*.fastq.gz ./fastq

        rm -rf ./PRJEB37709

<br>

<br>

2. Assess the quality of the nanopore full-length 16S rRNA amplicon sequences with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [NanoPlot](https://github.com/wdecoster/NanoPlot) and [multiqc](https://multiqc.info/) in one command-line with the `fastqc-py` **metaNanoPype** script: 

    Change directory to scripts to save the `log` under the script folder:

        cd ../script

        fastqc-py --help # display options

        fastqc-py -f ../data/fastq -n True -t 10 -o ../results/QC

    The command above will give as input all the fastq files in the folder (`-f` option) `../data/fastq`, it will run `NanoPlot` (`-n True`), with 10 threads (`-t 10`) and the output result (`-o` option) will be saved in the folder `../results/QC`.

    You can inspect the quality of the nanopore 16S sequences by looking into the individual html files produced by `fastqc` or the aggregated html report produced by `multiqc` as well the report produced by `Nanoplot` at: `../results/QC`. Through this way, you can have a good picture about the quality of your data. 

<br>

<br>

3. The next step consists in filtering and trimming bad quality nanopore sequences using [porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt) in one command-line with the `filter_fastq-py` **metaNanoPype** script: 

        filter_fastq-py --help # display options
        
        filter_fastq-py -f ../data/fastq -t 10 -min_len 1000 -max_len 1700 -qs 10 -o ../data/trim    

    The command above is filtering/trimming the fastq files given in the folder (`-f` option) `../data/fastq`, using 10 threads (`-t 10` - passed to `Porechop`) and discarding reads shorter than 1000 bp (`-min_len 1000`) or longer than 1700 bp (`-max_len 1700`) as well as reads with a quality-score lower than 10 (`-qs 10`). The output good-quality full-length nanopore reads are saved at `../data/trim` (with the option `-o ../data/trim`). Find more options with the `filter_fastq-py --help` command. 


<br>

<br>

4. Re-assess the quality of the reads after filtering/trimming them with the `fastqc-py` **metaNanoPype** script: 

        fastqc-py -f ../data/trim -n True -t 10 -o ../results/QC_trim

<br>

<br>

5. Perform taxonomic assignment with [kraken2](https://ccb.jhu.edu/software/kraken2/), inclusive the download of 16S rRNA reference databases indexed, in one command-line with the `tax_assign-py` **metaNanoPype** script:

        tax_assign-py --help # display options

        # create variables to simplify
        $FASTQ=$(ls ../data/trim/*.fastq.gz | xargs echo | sed 's/ /,/g')
        $SAMPLES=$(echo $FASTQ | sed 's/..\/data\/trim\/\|.fastq.gz//g')
        $DB=../data/SILVA_DB

        tax_assign-py -i $FASTQ -t 10 -db $DB -db_down silva -s $SAMPLES -r -o ../results/tax

    The command above will download the 16S rRNA database indexed `silva` (`-db_down silva` option) into the directory `../data/SILVA_DB` (`-db $DB` option) which will be the reference to map nanopore 16S reads provided at `../data/trim/*.fastq.gz` (`-i $FASTQ` with the option - comma-separated list of fastq file directories) with 10 threads (`-t 10`). The output folder with the results will be `../results/tax` (`-o ../results/tax` option) with files named based on sample names provided at `-s $SAMPLES` (a comma-separated list of samples/files name in the same order as provided in the input fastq files). The ouput file name will be: `<sample_name>.out`. Since the option `-r` was provided, a report kraken2 file will be created also (with the extension `<sample_name>.report`). 

    Inspect the **kraken2** taxonomic assignments at: `../results/tax`

<br>

---

<br>

### Support or Contact


Please open an [issue](https://github.com/antonioggsousa/metaNanoPype/issues) for **support or contact**.


<br>

---

<br>

### Acknowledgement

This project is being developed under the scope of the [Open Life Science 3](https://openlifesci.org/ols-3/projects-participants/) program.

<br>

**Project lead**: [António Sousa](https://openlifesci.org/ols-3/projects-participants/#antonioggsousa)

**Mentor**: [Hans-Rudolf Hotz](https://openlifesci.org/ols-3#hrhotz)

**Advice & support**: [Ricardo Ramiro](https://scholar.google.com/citations?user=bZkDsqEAAAAJ&hl=en)
