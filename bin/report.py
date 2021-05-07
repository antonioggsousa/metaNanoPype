#!/usr/bin/env python3

#-------------------------------------------------------------------------------------------------------
# author: António Sousa
# written: 27/01/2021
# updated: 27/01/2021
#
# programs: python3 packages (re, os, sys, time, markdown)
#
# TODO: in function 'parse_log_2_md()' add ' version:' tag to log files in order to avoid long 
# if statement below.
# 
#-------------------------------------------------------------------------------------------------------

# import modules
#import re
import os
import sys
import time
import subprocess

def grab_logs(dir_log = "./", log_list = None): 

    '''
    'grab_logs()': grab a list of log files from 
    the current dir_logectory by default ('dir_log = "./"')
    or from a 'log_list' if specified, and it returns
    a list of tuples with the file path name and time 
    of creation by order of creation.

     ---

    Parameters: 

    'dir_log': dir_logectory where the log files are located. 
    By default checks the current dir_logectory. It will 
    grab any text file with the extension '.log' under
    the dir_logectory specified. 

    'log_list': list of log files separated by comma (',') 
    to grab. If specified, the 'dir_log' option/parameter will
    be ignored. The absolute or relative path to the log
    files needs to be specified. 
    '''

    # list of log files to import
    if log_list is not None: 
        log_files = log_list.split(",") 
        log_files = [ file if os.path.isfile(file) else sys.exit("File " + file + " does not exists!\nExiting...") for file in log_files ]
    else: 
        files_2_check = os.listdir(dir_log)
        log_files = [ file for file in files_2_check if file.endswith(".log") ]
        if len(log_files) == 0: 
            sys.exit("It was not found log files under the dir_logectory provided:\n" +\
                dir_log + "\nPlease provide a dir_logectory containing the log files.")
    
    # list of tuples with log files and time of creation
    log_files_time = [ (file, os.path.getctime(file)) for file in log_files ]
    time_2_order = [ file[1] for file in log_files_time ]
    time_2_order.sort()
    log_files_time.sort(key = lambda i: time_2_order.index(i[1]))
    log_files_time = [ (file[0], time.ctime(file[1])) for file in log_files_time ]

    return(log_files_time)

#-------------------------------------------------------------------------------------------------------

from datetime import date
import re

def build_report(dir_log = "./", log_list = None, filename = None, 
                 report_name = "metaNanoPype reproducible report", 
                 author = "", references = None): 

    '''
    'build_report()': build a reproducible report in markdown and 
    html format.

     ---

    Parameters: 

    'dir_log': dir_logectory where the log files are located. 
    By default checks the current dir_logectory. It will 
    grab any text file with the extension '.log' under
    the dir_logectory specified. 

    'log_list': list of log files separated by comma (',') 
    to grab. If specified, the 'dir_log' option/parameter will
    be ignored. The absolute or relative path to the log
    files needs to be specified. 

    'filename': file name (str) to give to the markdown and html
    report built. By default is 'None'. 

    'report_name': title name of the report (str). By default is 
    "metaNanoPype reproducible report".

    'author': name of the author building the report (str). 
    By default is '""'.

    'references': reference file name (str) to import and add to the 
    report.
    '''

    ## check

    # date
    date_today = date.today()
    date_today = date_today.strftime("%d_%m_%Y")
    
    # file name
    if filename is None: 
        filename = "metaNanoPype_report_" + date_today + ".md"
    else: 
        filename = filename + ".md"

    md_file = open(filename, "w")

    # build header
    md_file.write("\n---\n")
    md_file.write("\n<br>\n\n")
    md_file.write("## Report: " + report_name + "\n")
    md_file.write("#### Author: " + author + "\n")
    md_file.write("#### Date: " + re.sub("_", "/",date_today) + "\n")
    md_file.write("\n<br>\n")
    md_file.write("\n---\n")
    md_file.write("\n<br>\n")

    # log files to build the report
    log_files = grab_logs(dir_log = dir_log, log_list = log_list)

    # loop over log files and structure them
    for file in log_files: 
        up_file = parse_log_2_md(file[0])
        for line in up_file: 
            md_file.write(line)
        
        # highlight differently different script runs
        md_file.write("\n<br>\n")
        md_file.write("\n---\n")
    
    # add references
    if references is None:
        try: 
            references = subprocess.run(["which", "report-py"], stdout = subprocess.PIPE)
            references = os.path.dirname(references.stdout.decode())
            references = references + "/ref/references.md"
        except: 
            references = None
    if references is not None and os.path.isfile(references): 

        md_file.write("\n<br>\n")
        md_file.write("\n<br>\n")
        md_file.write("\n---\n")
        md_file.write("\n<br>\n")
        ref_file = open(references, "r")
        for line in ref_file: 
            md_file.write(line)
        ref_file.close()
        md_file.write("\n<br>\n")
        md_file.write("\n---\n")

    # build footer
    md_file.write("\n<br>\n\n")
    md_file.write("\n<br>\n")
    md_file.write("\n---\n")    
    md_file.write("\n<br>\n\n")
    md_file.write("**metaNanoPype** project page: <a href='https://antonioggsousa.github.io/metaNanoPype' target='_blank'>https://antonioggsousa.github.io/metaNanoPype</a>\n")
    #md_file.write("GitHub **metaNanoPype** project page: [https://github.com/antonioggsousa/metaNanoPype](https://github.com/antonioggsousa/metaNanoPype)\n")
    md_file.write("\n<br>\n")
    #md_file.write("\n---\n")

    md_file.close()

    # convert md into html
    md_2_html(md_file = filename, html_file = None)

#-------------------------------------------------------------------------------------------------------

import re

def parse_log_2_md(log_file): 

    # TODO: add ' version:' tag to log files in order to avoid 
    # long if statement below.   

    '''
    'parse_log_2_md()': parse log file given to markdown. It 
    returns a list with a log file parsed to markdown. 

     ---

    Parameters: 

    'log_file': log file name to parse.
    '''

    # read log file
    file = open(log_file, "r")

    # list parsed
    parsed_file = []
    for line in file: 

        # once per run
        # run
        if " run:" in line: 
            parsed_file.append("\n")
            parsed_line = line.split(" ")
            parsed_file.append("<br>\n\n### " + parsed_line[0] + "\n" + "\n<br>\n\n")
            parsed_file.append("   + **" + parsed_line[1] + "** " + parsed_line[2] + " " + \
                               re.sub("\.\w+\n", "", parsed_line[3]) + "\n")

        # command-line
        if " command-line:" in line: 
            parsed_file.append("\n")
            parsed_line = line.split(": ")
            parsed_file.append("   + **" + parsed_line[0].split(" ")[1] + ":** `" + parsed_line[1].strip("\n") + "`\n")

        # multiple per run
        # version
        if "---\n" not in line and line != "\n" and " run:" not in line and " command-line:" not in line and " installed:" not in line:
            parsed_file.append("\n")
            parsed_file.append("   + **version:** " + line.strip("\n") + "\n")           

        # installed
        if " installed:" in line: 
            parsed_file.append("\n")
            parsed_line = line.split(": ")
            parsed_file.append("   + **" + parsed_line[0].split(" ")[1] + ":** `" + parsed_line[1].strip("\n") + "`\n")            

    file.close()
    return(parsed_file)

#-------------------------------------------------------------------------------------------------------

import markdown

def md_2_html(md_file, html_file = None): 

    '''
    'md_2_html()': converts markdown file into html file.

    ---

    Parameters: 

    'md_file' (mandatory): markdown file name (str) to be
    imported.    

    'html_file': html file name (str) to be exported. 
    By default 'None'. If 'html_file = None', the name of 
    markdown file is used instead.
    '''

    # import and read md file
    with open(md_file, "r") as line: 
        md = line.read()
        html = markdown.markdown(md)

    # export html file
    if html_file is None: 
        html_file = re.sub("\.\w+$", "", md_file)
        html_file = html_file + ".html"

    with open(html_file, 'w') as line:
        line.write(html)

#-------------------------------------------------------------------------------------------------------
