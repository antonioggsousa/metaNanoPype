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
   
   + (VI) diversity analyses (alpha- and beta-diversity) (*not implemented yet*). 


<br>

---

<br>

### Tutorial





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
