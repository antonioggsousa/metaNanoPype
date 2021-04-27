#!/usr/bin/env python3

#-------------------------------------------------------------------------------------------------------
# author: António Sousa
# written: 26/12/2020
# updated: 28/03/2021
#
# programs: NCBI BLAST+ command-line suites (v.2.6.0) and python3 packages (re, os, sys, tqdm, pandas, 
# biopython, subprocess)
#
# TODO: 
# 
# - 'blast_remote_fasta()' substitute the option '-max_target_seqs' which can produce erroneous 
# results (check scientific literature).
# - 'fastq2fasta()' include the possibility to provide gzip fastq files. 
# 
#-------------------------------------------------------------------------------------------------------

# import modules
import re
import os
import sys
import subprocess
import datetime
#from tqdm import tqdm

def fastq2fasta(fastq, fasta):

    # TODO: include the possibility to provide gzip fastq files

    '''
    'fastq2fasta()' converts fastq to fasta format.

     ---

    Parameters: 
   
    'fastq' (mandatory): file name (str) of the input fastq file. 

    'fasta' (mandatory): file name (str) of the output fasta file.
    
    'fastq' and 'fasta' are mandatory input arguments.
    Each fastq sequence should have the following structure: 
    (1) fastq header which starts with "@" character (1 line); 
    (2) fasta sequence (>= 1 line(s)); 
    (3) plus "+" sign/character (1 line); 
    (4) quality sequence with Phred Q scores (>= 1 line(s)).
    Apart from the header and the line with the plus sign, 
    the fasta and quality sequences can have more than one line. 
    '''

    # count the no. of lines in fastq file to show progess bar w/ 'tqdm'
    total_lines = len(open(fastq, "r").readlines())

    # open files to read/write I/O
    in_file = open(fastq, "r")
    out_file = open(fasta, "w")

    for read in tqdm(in_file, total = total_lines): # loop over each line
        if read.startswith("@"): # header
            read_parsed = re.sub("@", ">", read)
            out_file.write(read_parsed)
            isFasta = True
        elif read.startswith("+"): # '+' line
            isFasta = False
        elif not read.startswith("@") and isFasta: # fasta seq   
            out_file.write(read)
        else: # quality seq
            continue

    # close files 
    in_file.close()
    out_file.close()

#-------------------------------------------------------------------------------------------------------

def build_kraken2_db(db_dir, db_name, threads = 1):

    # written: 25/04/2021
    # updated: 25/04/2021

    '''
    'build_kraken2_db()': download and build a special
    kraken2 database of 16S reference sequences. One 
    of 'silva', 'greengenes', 'rdp'.

    ---

    Parameters:

    'db_dir' (mandatory): database directory name to 
    save/write the database. 

    'db_name' (mandatory): database name to download 
    and build/index. One of 'silva', 'greengenes', 'rdp'.

    'threads' (mandatory): number of threads to use. 
    By default is 1 core.  
    '''

    ## Check input
    db_available: ["silva", "greengenes", "rdp"]
    if db_name not in db_available: 
        sys.exit("The " + db_name + " 'db_name' provided is not among the databases available!\n\
                 Please provide as db_name 'silva', 'greengenes' or 'rdp'.\n\
                 Aborting...")

    ## Check if database files already exist; if they exist, abort
    if os.path.exists(db_dir):
        db_dir_files = os.listdir(db_dir) # list the files of the db dir
        db_files = ["hash.k2d", "opts.k2d", "taxo.k2d"] # database files necessary to run kraken2 
        check_db_files = all(f in db_dir_files for f in db_files) # return 'True' if all databases already exist
        if check_db_files: 
            sys.exit("All the three database files required were found in 'db_dir':\n\
                    'hash.k2d', 'opts.k2d', 'taxo.k2d'!\n\
                    Please provide a different 'db_dir' if you want to download\n\
                    and build again the database!")  

    ## Download database 
    check_tool("kraken2-build") # check if tool is in your PATH
    print("Downloading " + db_name + " ...")
    subprocess.run(["kraken2-build", "--standard", "--db", db_dir, "--special", db_name, 
                    "--threads", str(threads)])
    print(db_name + " downloaded at: " + db_dir)

#-------------------------------------------------------------------------------------------------------

## Import packages
#

def map_reads_with_kraken2(fastq, db_dir, samples, report = False, threads = 1, out_dir = "./"):

    # written: 25/04/2021
    # updated: 25/04/2021

    '''
    'map_reads_with_kraken2()': map fastq/fasta nanopore reads
    against a reference database provided with kraken2.

    ---

    Parameters:

    'fastq' (mandatory): fastq/fasta file directories (str).
    A string separated by comma for the different file 
    directories. For example: 
    'sample_A.fastq.gz,sample_B.fastq.gz,sample_C.fastq.gz'.  

    'db_dir' (mandatory): database directory name where the  
    database indexed kraken2 files are saved. 

    'samples' (mandatory): name of the samples to give as 
    kraken2 output. It will be added the suffix '.out'.  

    'report' (mandatory): if 'True' report is created with the 
    name provided in 'samples' and the suffix '.report'. 
    By default is 'False' - it is not created.

    'threads' (mandatory): number of threads to use. 
    By default is 1 core.  

    'out_dir' (optional): output directory name. If not given 
    it will use the default that is current working directory.
    '''

    ## Check fastq
    fastq = str.split(fastq, ",")
    check_fastq_exists = all(os.path.isfile(fast) for fast in fastq)
    if not check_fastq_exists: 
        sys.exit("The 'fastq' files provided were not found!\n\
                 Please check if all the paths to the input fastq/fasta\n\
                 files provided exist or if you provide them in the right format:\n\
                 e.g., 'sample_A.fastq.gz,sample_B.fastq.gz,sample_C.fastq.gz'.\n\
                 Aborting...")
    samples = str.split(samples, ",")
    if len(fastq) != len(samples): 
        sys.exit("The number of fastq/fasta files provided in 'fastq'\n\
                 is different from 'samples'.\n\
                 Please provide a comma-separated list of 'fastq' and\n\
                 respective 'samples' of the same length")

    ## Check 'out_dir'
    if not os.path.exists(out_dir): # if folder not exists, create it
        os.mkdir(out_dir) 
    if not out_dir.endswith("/"): # check if path ends with '/', otherwise add it
        out_dir = out_dir + "/"
    samples_dir = [ out_dir + samp + ".out" for samp in samples ] # add extension '.out' to 'samples'

    ## Add report 
    if report: 
        report_name = [ out_dir + samp + ".report" for samp in samples ]

    ## Map reads with kraken2: 
    no_fastq = len(fastq)
    for f in range(no_fastq): 
        check_tool("kraken2") # check if tool is in your PATH
        if report: 
            print("Mapping sample " + samples[f] + " corresponding to the fastq/fasta " + fastq[f] + "...")
            subprocess.run(["kraken2", "--db", db_dir, fastq[f], 
                            "--threads", str(threads), "--output", samples_dir[f], 
                            "--report", report_name[f]])
            print("Sample " + samples[f] + " mapped and output created at: " + samples_dir[f] + "...")
            print("Report for sample " + samples[f] + " was saved at: " + report_name[f] + "...")
            print("")
            print("..........................................................................................")
            print("")
        else: 
            print("Mapping sample " + samples[f] + " corresponding to the fastq/fasta " + fastq[f] + "...")

            subprocess.run(["kraken2", "--db", db_dir, fastq[f], 
                            "--threads", str(threads), "--output", samples_dir[f]])     
            print("Sample " + samples[f] + " mapped and output created at: " + samples_dir[f] + "...")       
            print("")
            print("..........................................................................................")
            print("")

#-------------------------------------------------------------------------------------------------------

#import os
#import re
#import sys
#import wget # download files
#import subprocess
#requires: 'md5sum' GNU/Linux utility
#requires: 'update_blastdb' from NCBI BLAST+ command-line suites (v.2.6.0)

def get_blast_db(db_folder = None, ncbi_db = None, url = None, 
                 show_ncbi_blast_db = False, check_ncbi_md5 = True, 
                 passive = False, force = False, quiet = False, 
                 decompress = False, check_md5 = None): 

    # written: 26/12/2020
    # updated: 28/03/2021

    '''
    'get_blast_db()': it allows to download data from the NCBI BLAST
    databases repository - https://ftp.ncbi.nlm.nih.gov/blast/db/ -
    with 'update_blastdb.pl' perl script from NCBI BLAST+ utilities 
    or any other online database by providing a ftp/https url. It 
    allows to check for the md5 hash sum to assess if the database 
    was properly downloaded and decompression of the same database 
    (if in the format of: '.tar.gz', '.gz', '.zip'). 

    ---

    Parameters:

    'db_folder' (optional): database folder name/path (str) to where the 
    database will be downloaded. If 'None' (default) it will be downloaded 
    to the current working directory. 

    'ncbi_db' (optional): by default 'None'. Given a name of a NCBI database 
    (str) among the options provided (run 'get_blast_db(show_ncbi_blast_db = True)'
    to see the possible options), it will download the database in the 'db_folder'
    provided with the 'update_blastdb' command-line tool from NCBI BLAST+ 
    command-line suites (v.2.6.0).

    'url' (optional): download an online database, other than NCBI, by 
    providing the url (ftp/https hyperlink) (str) to the 'db_folder' folder
    provided. By default is 'None'. 

    'show_ncbi_blast_db' (optional): logical, by default is 'False'. Print 
    the NCBI databases names that are available to download and exit. 

    'check_ncbi_md5' (optional): logical, by deafult is 'True'. Only applied 
    if a NCBI database is downloaded. It uses the md5 hash sum file downloaded
    with the NCBI database to check if the hashes match or not. Based on the 
    'md5sum' GNU/Linux utility. 

    'passive' (optional): by default is 'False'. If 'passive' should be passed to the  
    'update_blastdb' command-line tool from NCBI BLAST+ command-line suites (v.2.6.0).

    'force' (optional): by default is 'False'. If 'force' should be passed to the
    'update_blastdb' command-line tool from NCBI BLAST+ command-line suites (v.2.6.0).
    
    'quiet' (optional): by default is 'False'. If 'quiet' should be passed to the
    'update_blastdb' command-line tool from NCBI BLAST+ command-line suites (v.2.6.0).

    'decompress' (optional): by default is 'False'. If 'True' it will decompress the 
    database downloaded by checking its extension. Formats supported are: '.tar.gz', 
    '.gz', '.zip'.  
    
    'check_md5' (optional): by deafult is 'None'. Only applied to a custom database, 
    or database other than NCBI. If 'True' (logical) it will attempt to download a 
    file with the md5 hash sum that is available in many repositories with the 
    same link as the database but with the extension ".md5". If a string is provided 
    instead, it will compare this string with the md5 hash sum of the database 
    downloaded. In both cases the database md5 hash sum is checked with the 
    'md5sum' GNU/Linux utility. 

    '''
    
    # check conditions
    # only one can be different from 'None': 'ncbi_db = None'; 'url = None'
    if (ncbi_db is not None) and (url is not None): 
        sys.exit("You can not download a NCBI database simultaneously with a custom database.\n \
           Only one database can be downloaded at a time.")

    # current wd
    wd = os.getcwd()

    # NCBI databases options to download/update
    ncbi_blast_dbs = ['16S_ribosomal_RNA','18S_fungal_sequences','28S_fungal_sequences',
                      'Betacoronavirus','ITS_RefSeq_Fungi','ITS_eukaryote_sequences',
                      'LSU_eukaryote_rRNA','LSU_prokaryote_rRNA','SSU_eukaryote_rRNA',
                      'cdd_delta','env_nr','env_nt','human_genome','landmark','mito',
                      'mouse_genome','nr','nt','pataa','patnt','pdbaa','pdbnt',
                      'ref_euk_rep_genomes','ref_prok_rep_genomes','ref_viroids_rep_genomes',
                      'ref_viruses_rep_genomes','refseq_protein','refseq_rna','refseq_select_prot',
                      'refseq_select_rna','swissprot','taxdb','tsa_nr','tsa_nt']

    # parse options to 'update_blastdb' perl script
    parse_update_blastdb = {'passive': '--passive', 
                            'force': '--force', 
                            'quiet': '--quiet'} # match/compatible with BLAST+
    update_blastdb_names = ['passive', 'force', 'quiet'] # names to be parsed
    update_blastdb_opts = [passive, force, quiet] # options to be parsed

    # show blast databases options to download/update
    if show_ncbi_blast_db: 
        print("NCBI blast databases:\n")
        [print(i) for i in ncbi_blast_dbs]
        sys.exit(0) # exit with success 

    try: 
        # mkdir and change directory 
        if db_folder is not None: 
            if not os.path.exists(db_folder): 
                os.makedirs(db_folder)
            os.chdir(db_folder)

        # download from NCBI repo: 
        if ncbi_db in ncbi_blast_dbs: 
            check_tool("update_blastdb") # check if tool is in your PATH
            opts_2_add = ["update_blastdb", ncbi_db]
            for opt in range(len(update_blastdb_opts)): # update and parse options to pass to 'update_blastdb'
                if update_blastdb_opts[opt]: 
                    #print(update_blastdb_names[i])
                    opts_2_add.append(parse_update_blastdb[update_blastdb_names[opt]])
            print("Downloading " + ncbi_db + " ...")
            subprocess.run(opts_2_add) # run 'update_blastdb'
            print(ncbi_db + " database downloaded!")

            # check md5sum: integrity
            if check_ncbi_md5:
                md5file = [i for i in os.listdir() if re.search("^"+ncbi_db+"(.*?)\.md5$", i)] # get hash md5 sum
                md5file = md5file[0]
                read_hash_from_file = open(md5file, 'r')
                md5good = [line.split(' ')[0] for line in read_hash_from_file]  
                md5good = md5good[0] # get hash md5sum 
                ncbi_db_file = re.sub(".md5", "", md5file) # get db file name compressed
                md5db = subprocess.check_output(["md5sum", ncbi_db_file]) # get db file name compressed
                md5db = md5db.decode('UTF-8').split(" ")[0] # convert bytes to unicode strings
                read_hash_from_file.close()
                database_name = ncbi_db_file
                if md5good == md5db: 
                    print("The database " + database_name + " was properly downloaded!")     
                else: # if the download fails, it will fail above (this message will never be printed - remove it)
                    sys.exit("The database " + database_name + " is corrupted!\n\
                        Please check your internet connection and download it again!") 
        else: 
            sys.exit("'ncbi_db' is not among the valid options of NCBI databases.\n \
            Please provide a valid database to download. Check valid\n \
            databases (run): get_blast_db(show_ncbi_blast_db = True)")
        
        # download other database that is not from NCBI
        if url is not None: 
            database_name = os.path.basename(url)
            print("Downloading " + database_name + " ...")
            wget.download(url) # download db
            print(database_name + " database downloaded!")
            if check_md5 is not None: # logical or string
                if check_md5: # if logical attempt to download the md5 hash sum with the database 
                    # by adding to the db file name the extension '.md5' - if the text file exists in the db repo
                    # it'll be downloaded and the first column & row imported as a hash string 
                    wget.download(url + ".md5") # attempt to download db md5 hash file (it it exists in the repo)
                    md5file = [i for i in os.listdir() if re.search("^"+database_name+"(.*?)\.md5$", i)] # get hash md5 sum
                    md5file = md5file[0]
                    read_hash_from_file = open(md5file, 'r')
                    md5good = [line.split(' ')[0] for line in read_hash_from_file]  
                    md5good = md5good[0] # get hash md5sum 
                    read_hash_from_file.close()
                else: # if the md5 hash string itself is provided
                    md5good = check_md5
                md5db = subprocess.check_output(["md5sum", database_name]) # md5sum hash str: get db file name compressed
                md5db = md5db.decode('UTF-8').split(" ")[0] # convert bytes to unicode strings
                if md5good == md5db: # check if hashes match
                    print("The database " + database_name + " was properly downloaded!")     
                else: # if the download fails, it will fail above (this message will never be printed - remove it)
                    sys.exit("The database " + database_name + " is corrupted!\n\
                        Please check your internet connection and download it again!") 

        # decompress database 
        if decompress: 
            print("Decompressing " + database_name + " database!")
            decompress_file(file_name = database_name)
            print(database_name + " database decompressed!")

    finally:
        # get back to the starting point dir
        os.chdir(wd)

#-------------------------------------------------------------------------------------------------------

def index_db_2_blast(db_path = None, db_index = None, 
                     dbtype = "nucl", title = None, 
                     parse_seqids = True, force = False):

    # written: 28/03/2021
    # updated: 28/03/2021

    # TODO: check if indexed database already exists for 'dbtype = "prot"'.
    # Currently, it only checks for nucleic databases.

    '''
    'index_db_2_blast()': index a local database to blast by using 
    the 'makeblastdb' command-line utility from NCBI BLAST+ 
    command-line suites (v.2.6.0).

    ---

    Parameters: 

    'db_path' (mandatory): the path to the database fasta file name (str) 
    to index by 'makeblastdb' command-line utility from NCBI BLAST+ 
    command-line suites (v.2.6.0).

    'db_index' (optional): the path (str) given will host the output of 
    indexed files and the prefix the name of the index files. For instance 
    for the path - "/path/out_db/prefix" - the output will be created under 
    the directory '/path/out_db' and the database name/prefix will be 'prefix'.
    If not given the prefix will be the name of the fasta file from the 
    database and the directory fiven in 'db_path'.

    'dbtype' (mandatory): one of 'prot' (for protein databases) or 'nucl' 
    (for nucleic databases) (str). By default 'nucl'.

    'title' (optional): by default 'None' (str).

    'parse_seqids' (optional): parse sequence ids (logical). By default is 
    'True'.

    'force' (mandatory): by default 'False'. The directory that will host 
    the indexed files is checked to verify if the files already exist. If
    any of the indexed files expected to be created already exists, it exists.
    To force it to overwrite the files, choose 'force = True'.
    '''

    # add index if 'db_index = None'
    if db_index is None:
        db_index = re.sub('.fasta$|.fa$|.fna$', "", db_path)
    db_name = os.path.basename(db_index)
    db_dir = os.path.dirname(db_index)

    # check if files already exist in the indexed directory:
    exts = [ db_name + '.nog', db_name + '.nsd', db_name + '.nsi', 
             db_name + '.nhr', db_name + '.nin', db_name + '.nsq']
    if force is False: 
        for file in os.listdir(db_dir): 
            if file in exts: 
                sys.exit("File " + file + " already exists under the " + db_dir + " \n \
           database directory provided!\n \
           Use the option 'force' if you want to overwrite it!\n \
           Aborting!")

    # parse options
    parse_update_opts = {'db_path': '-in', 
                            'db_index': '-out', 
                            'dbtype': '-dbtype', 
                            'title': '-title', 
                            'parse_seqids': '-parse_seqids'} # match/compatible with makeblastdb BLAST+
    update_opts = [db_path, db_index, dbtype, title, parse_seqids] # options to be parsed   
    update_opts_names = ['db_path', 'db_index', 'dbtype', 'title', 'parse_seqids']
    opts_2_add = ["makeblastdb"]
    for opt in range(len(parse_update_opts)): # update and parse options to pass to 'update_blastdb'
        if (type(update_opts[opt]) == str): 
            param_name = parse_update_opts[update_opts_names[opt]]
            param_opt = update_opts[opt]
            opts_2_add.append(param_name)
            opts_2_add.append(param_opt)
        if update_opts[opt] is True: 
            param_name = parse_update_opts[update_opts_names[opt]]
            opts_2_add.append(param_name)
        
    # check tool and run it with subprocess
    check_tool("makeblastdb") # check if tool is in your PATH
    print("Indexing database " + db_name + " at " + db_dir + " ...")
    subprocess.run(opts_2_add)
    print("Database " + db_name + " indexed at " + db_dir)

#-------------------------------------------------------------------------------------------------------

#import subprocess

def blast_remote_fasta(fasta, blast_out, blast_fmt = 6, blast_tool = "blastn", 
                       db_name = "nt", db_type = "remote", hits_per_seq = 1, 
                       evalue = 10):

    '''
    'blast_remote_fasta()': blast remotely (multi-)fasta file with BLAST+  
    command-line suites. If the option 'blast_fmt' is set to '6', the 
    blast output result is automaticaly parsed with the function 
    'parsing_blast_outfmt_six()' to keep one hit per sequence based on
    'bitscore' (1st) and 'pident' (2nd) (see more with: help(parsing_blast_outfmt_six))

    --- 

    Parameters: 

    'fasta' (mandatory): (multi-)fasta file (str) given as input to BLAST.

    'blast_out' (mandatory): file name to save blast results (str).

    'blast_fmt': blast output format (int). By default '6'. One from 0-18.
    
    'blast_tool': one from the BLAST+ command-line suites. By default 
    "blastn" (str).
    
    'db_name': database name. By default nucleotide (str): 'nt'. 
    
    'db_type': blast remotely (str): "remote". 
    
    'hits_per_seq': maximum hits per sequence (int). By default '1'.
                       
    'evalue' = e-value threshold (int). By default is '10'.
    '''

    # check if 'blast_tool' is in your PATH
    check_tool(tool = blast_tool)

    # blast remotely
    if db_type == "remote": 
        subprocess.run([blast_tool, "-db", db_name, "-" + db_type, 
                        "-query", fasta, "-out", blast_out, 
                        "-outfmt", str(blast_fmt), "-max_target_seqs", str(hits_per_seq), 
                        "-evalue", str(evalue)])

    # parse blast output if 'blast_fmt'/'-outfmt == 6' 
    if blast_fmt == str(6): 
        # define file name of the parsed output
        parse_ext = "_parsed.tsv"
        if "." in blast_out: 
            parsed_blast_out = re.sub("\..+", parse_ext, blast_out)
        else: 
            parsed_blast_out = blast_out + parse_ext
        parsing_blast_outfmt_six(blast_out, parsed_blast_out) # run parsing function
    else: 
        print('WARNING: "blast_fmt" parameter (corresponding to the "-outfmt" original option in blast) \
              it is different from "6". Exiting without performing parsing to the blast output.')

#-------------------------------------------------------------------------------------------------------

#import pandas as pd

def parsing_blast_outfmt_six(blast_out, parsed_blast_out):

    # TODO: check if there is a need to group by 'qseqid' and sort by 'bitscore'.
    # Blast results seem grouped by 'qseqid' and sorted by 'bitscore'.

    '''
    'parsing_blast_outfmt_six()': parses the blast result with output format 6 by
    sorting the columns 'bitscore' (1st) and 'pident' (2nd) by descending order 
    and grouping the hits by query ('qseqid'). The top hit by query based on the 
    columns 'bitscore' (1st) and 'pident' (2nd) is retrieved. 

    ---

    Parameters: 

    'blast_out' (mandatory): blast result (str file name) with output format 6.

    'parsed_blast_out' (mandatory): blast result parsed (str file name) to save.
    '''

    # import blast results data with pandas
    data = pd.read_table(blast_out, header = 0, 
                         names = ["qseqid", "sseqid", "pident", 
                                  "length", "mismatch", "gapopen", 
                                  "qstart", "qend", "sstart", "send", 
                                  "evalue", "bitscore"]) 

    # for each query seq ("qseqid"), pick up the hit/target with the highest 'bitscore' and 'pident'  
    data2save = data.sort_values(["bitscore", "pident"], 
                                  ascending = [False, False]).groupby("qseqid").head(1).reset_index(drop=True) 

    # save the parsed blast result
    data2save.to_csv(parsed_blast_out, sep = "\t", index=False)

#-------------------------------------------------------------------------------------------------------

#import os
#import sys
#import pandas as pd
#from tqdm import tqdm
#from Bio import Entrez

def get_taxonomy_from_NCBI(acc_ids_name, file_name = None,
                           col_ids = None, header = None, sep = '\t'):

    '''
    'get_taxonomy_from_NCBI()': get taxonomy information from NCBI based on a list 
    of nucleotide NCBI accession numbers using the biopython Entrez efetch function.
    The list of accession numbers can be provided directly or through a file. The 
    result is a table with 2 columns: 'Accession_no' (accession numbers used) and 
    'Taxonomy' (full taxonomy).
    
    ---

    Parameters: 

    'acc_ids_name' (mandatory): list or a file name (str) with a list of accession 
    numbers. If a file name is provided, the table will be imported with Pandas into 
    a data frame and later converted into a list.  

    'file_name': by default ('None') the results are returned to shell. Give a 
    file name (str) to save the table into a file.

    'col_ids': if a list of accession numbers is imported from a file ('acc_ids_name')
    you can use 'col_ids' to specify the column of the file that contains the accession
    numbers. It can be a number (starts from 0) or the name of the column (in case you
    have one - change the next parameter to 'header = 0' - if the header is the 1st row). 

    'header': if a list of accession numbers is imported from a file ('acc_ids_name'), 
    it assumes by default that the file does not have an header. Change it to 'header = 0'
    if the 1st row is the header. Otherwise set the header to the respective row number. 

    'sep': if a list of accession numbers is imported from a file ('acc_ids_name'), 
    by default assumes that the text file is a tab-separated file  with separator 
    field set to '\t'. If you have a 'csv' (comma-separated values/file) change it to 
    'sep = ',''.
    '''

    # import accession numbers from file; 
    # otherwise needs to be a list of accession numbers
    if type(acc_ids_name)==str and os.path.isfile(acc_ids_name): # check if the file exists
        acc_ids = pd.read_table(acc_ids_name, # import df
                                header = header, 
                                sep = sep)
        if col_ids != None: # select col from df
            if type(col_ids) == int:
                if not col_ids < acc_ids.shape[1]: 
                    sys.exit("Please provide a 'col_ids' number within the no. of columns \n \
           that the file '" + acc_ids_name + "' contains. The file \n \
           contains " + str(acc_ids.shape[1]) + " columns, column to select is higher: " + str(col_ids) + ".\n \
           Do not forget that python starts to count from 0!" )  
                acc_ids = acc_ids.iloc[:,[col_ids]]   
            if type(col_ids) == str:
                if not col_ids in acc_ids.columns: 
                    sys.exit("Please provide a 'col_ids' name within the column names\n \
           of the the file '" + acc_ids_name + "'. The column name provided '" + str(col_ids) + "'\n \
           is not among the column names comprised by the file: " + str(acc_ids.columns.values) + ".\n \
           If column names are numeric than you probably forgot to set the parameter 'header' to:\n \
           'header = 0' (assuming the 1st row as header).")              
                acc_ids = acc_ids.loc[:,[col_ids]]  
        if acc_ids.shape[1] != 1: 
            sys.exit("The text file provided '" + acc_ids_name + "' does not have one column.\n \
           Please make sure you provide a file with one column of accession numbers or \n \
           a file with multiple columns but in that case specify the column with \n \
           accession numbers to select by specifying the parameter 'col_ids'.\n \
           'col_ids' can be a str with the column name or column number (int - starting from 0)")
        # if it has one column: convert a df pandas series to list
        acc_ids = acc_ids[acc_ids.columns[0]].tolist()
    if type(acc_ids) != list: 
        sys.exit("'acc_ids' (" + str(acc_ids_name) + ") is of type " + str(type(acc_ids)) + ".\n \
            Please, provide a list of accession numbers (type list)\n \
            or a text file with accession numbers (type str of the file name)")  

    # get taxonomy
    lineage = [] # save taxonomy
    for acc_id in tqdm(acc_ids):
        try: # try to retrieve taxonomy with the biopython Entrez efetch method 
            efetch_acc_id = Entrez.efetch(db="nucleotide", id=acc_id, retmode="xml")
            seq_tax = Entrez.read(efetch_acc_id)
            lineage.append(seq_tax[0]["GBSeq_taxonomy"])
            efetch_acc_id.close()
        except: # raise exception if the accession no. provided is a bad request
            lineage.append("Bad request: " + acc_id)
            print("Bad request (accession number): " + acc_id)
            print("Please check the accession number provided!")
            continue
    
    # save data
    acc_id_tax = zip(acc_ids, lineage)
    data2save = pd.DataFrame(acc_id_tax, columns=['Accession_no', 'Taxonomy'])

    if type(file_name) == str:
        data2save.to_csv(file_name, sep = "\t", index=False) 
    else: 
        return(data2save)

#-------------------------------------------------------------------------------------------------------

#import os
#import sys
#import pandas as pd

def tax_summary(taxonomy_tbl, file_name = None, sep = '; '):
    
    '''
    'tax_summary()': summarizes a taxonomy table provided. It sums up the
    taxonomy to get the 'Abundance' of each unique full taxonomic path. 
    The taxonomy table provided needs to have a column named 'Taxonomy' 
    with the full taxonomic path, i.e., with all the taxonomic ranks from 
    'Kingdom' to 'Species' in one-line separated by default by '; ' ('sep'), 
    such as: 'Bacteria; Firmicutes; Clostridia; Clostridiales; Ruminococcaceae; 
    Faecalibacterium'. It returns (by default) or saves to a file a table
    with one column of abundance 'Abundance', with the abundance of unique 
    full taxonomic paths and the taxonomy split by taxonomic ranks possibly 
    named from: 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 
    'Species'.
    
    ---

    Parameters:

    'taxonomy_tbl' (mandatory): a pandas DataFrame obtained with the function 
    'get_taxonomy_from_NCBI()' or with one column named 'Taxonomy' with the 
    full taxonomic path. Also, it can be provided a string to the file path 
    of the taxonomic table to be imported with pandas.

    'file_name': if provided the file name (str) of the summarized taxonomic 
    table will be saved in tab-separated format instead of returned. Otherwise 
    by default the table/result is returned.  

    'sep': separator field of the full taxonomic path in the 'Taxonomy' 
    column to parse. By default is '; '.
    '''

    # import the table as pandas dataframe 
    if type(taxonomy_tbl)==str and os.path.isfile(taxonomy_tbl):
        taxonomy_tbl = pd.read_table(taxonomy_tbl)

    if not isinstance(taxonomy_tbl, pd.DataFrame): 
        sys.exit("The 'taxonomy_tbl' given - '" + str(taxonomy_tbl) + "' - is not a pandas Data Frame!\n \
           Please provide a pandas Data Frame or a file name (str) to\n \
           import the taxonomy table to summarize!")

    # count unique (strings) taxonomy: abundance
    tax_parse = taxonomy_tbl['Taxonomy'].value_counts().reset_index()
    tax_parse.columns = ['Taxonomy', 'Abundance'] # rename columns

    # split taxonomy by ranks 
    taxa_ranks = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'] # tax ranks 
    tax_parse_tax = tax_parse.Taxonomy.str.split(sep, expand = True) # split full taxonomy
    tax_parse[taxa_ranks[:tax_parse_tax.shape[1]]] = tax_parse_tax # add tax ranks by the col len of 'tax_parse_tax'
    tax_parse = tax_parse.drop(['Taxonomy'], axis=1) # delete 'Taxonomy' column

    # save the result or return 
    if type(file_name) == str:
        tax_parse.to_csv(file_name, sep = "\t", index=False) 
    else: 
        return(tax_parse)

#-------------------------------------------------------------------------------------------------------

""" 
## WARNING: remote blast with biopython (INCOMPLETE)

# import modules
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML

def blast_fasta(fasta, blast_tool = 'blastn'): 

    '''
    'blast_fasta' blasts a fasta file against the NCBI nt database.
    'fasta': fasta file to blast agains the NCBI nt database.
    'blast_tool' (one from):'blastn', 'blastp', 'blastx', 'tblast' and 'tblastx'.
    '''

    seq = SeqIO.read(fasta, format = "fasta") # import fasta
    result_handle = NCBIWWW.qblast(blast_tool, "nt", seq.format("fasta"))

    with open("blast_results.xml", "w") as out_handle: # save xml result
        out_handle.write(result_handle.read())
    result_handle.close()

    result_xml = open("blast_results.xml")
    blast_record = NCBIXML.read(result_xml)

    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...") """

#-------------------------------------------------------------------------------------------------------

#import re
#import os 
#import sys
#import gzip # for '.gz'
#import tarfile # for '.tar.gz'
#from zipfile import ZipFile # for '.zip' 

def decompress_file(file_name):

    '''
    'decompress_file()': decompress a given file 
    by identifying its file extension. Compression 
    formats supported are: '.zip', '.gz', '.tar.gz'.
    The decompressed file will have the same name 
    as the file name of the compressed file provided.

    Parameters

    ---

    'file_name' (mandatory): given a file name/path (str) 
    to decompress depending on the file extension format 
    supported. 
    '''

    # check file extension
    extensions = [".tar.gz", ".gz", ".zip"]
    file_extension = [ ext for ext in extensions if ext in file_name ]
    if len(file_extension) == 0: 
        sys.exit("The file to decompress ('" + file_name + "') is not supported!\n \
            The formats supported are: '.tar.gz', '.gz', '.zip'.\n \
            Please provide a compressed file supported!")
    file_extension = file_extension[0]

    # check if output file name exists in current wd; if so, abort!
    file_name_out = re.sub(file_extension + "$", "", file_name)
    for file_dir in os.listdir():
        if file_dir == file_name_out: 
            sys.exit("The file " + file_name_out + " exists in the current directory!\n\
            Decompressing the " + file_name + " would overwrite it.\n\
            Remove the file " + file_name_out + " from the current directory\n\
            or rename it. Aborting!")

    ## files
    # .zip
    if file_extension == ".zip": 
        with ZipFile(file_name, 'r') as zip:
            zip.extractall()

    # .gz
    if file_extension == ".gz":
        ungz_file_name = re.sub(".gz$", "", file_name)
        ungz_file = open(ungz_file_name, "w")
        with gzip.open(file_name, 'rb') as gz: # read compressed '.gz' file
            line = gz.read() # convert gz to bytes
            ungz_file.write(line.decode('UTF-8')) # write bytes converted into strings
        ungz_file.close()

    # .tar.gz
    if file_extension == ".tar.gz": 
        tar = tarfile.open(file_name, "r:gz")
        tar.extractall()
        tar.close()    

#-------------------------------------------------------------------------------------------------------

def check_tool(tool): 

    '''
    'check_tool': checks if command-line tool is in your PATH with 'which' command.
    It exits if 'which' command returns 1 meaning that the program is not in your 
    PATH.
    '''

    soft_call = subprocess.call(["which", tool], stdout = subprocess.PIPE)
    if soft_call != 0:
        print(tool + " is not in your PATH!")
        sys.exit("Install or export " + tool + " to your PATH! Exiting...")

#-------------------------------------------------------------------------------------------------------

def cmd_log(script_file, program, pattern_version, args, parse_args_dic):

    soft_call = subprocess.call(["which", program], stdout = subprocess.PIPE)
    program_install = subprocess.run(["which", program], stdout = subprocess.PIPE)
    program_version = subprocess.run([program, pattern_version], stdout = subprocess.PIPE)
    if soft_call != 0:
        print(program + " is not in your PATH!")
        sys.exit("Install or export " + program + " to your PATH! Exiting...")
    else: 
        log_file = open(script_file + ".log", "a")
        #if os.stat(script_file + ".log").st_size == 0: # if file is empty      
        log_file.write(script_file + " run: " + str(datetime.datetime.now()) + "\n")
        cmd_line = "./" + script_file
        for arg in vars(args):
            if getattr(args, arg) is not None and type(getattr(args, arg)) is not bool: 
                cmd_line = cmd_line + " " + str(parse_args_dic[arg]) + " " + str(getattr(args, arg))
            if getattr(args, arg) is True: 
                cmd_line = cmd_line + " " + str(parse_args_dic[arg])
        log_file.write(script_file + " command-line: " + cmd_line)
        log_file.write("\n\n---\n\n")
        log_file.write(program_version.stdout.decode())
        log_file.write(program + " installed: " + program_install.stdout.decode() )
        log_file.write("\n---\n\n")
        log_file.close()

#-------------------------------------------------------------------------------------------------------
