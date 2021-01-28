#-------------------------------------------------------------------------------------------------------
# author: António Sousa
# written: 19/12/2020
# updated: 19/12/2020
#
# programs: porechop (v.0.2.4) & NanoFilt (v.2.7.1)
#
# TODO: ...
#-------------------------------------------------------------------------------------------------------

import os
import subprocess
import re
import datetime

def filter_and_trim(fastq_path, out_folder, threads = 1, min_len = None, 
                    max_len = None, qs = None, min_gc = None, max_gc = None, 
                    trim_lead = None, trim_end = None): # filter/trim primers/adapters with Porechop
    fastq_name = os.path.basename(fastq_path)
    fastq_name = re.sub(".fastq.gz|.fastq", "", fastq_name)
    if not fastq_path.endswith(("fastq", "fastq.gz")): # exit if a file is not .fastq or .fastq.gz
        print("fastq file does not ends in fastq or fastq.gz!", "\nExiting...")
        exit()
    # porechop: trim primers + adapters
    subprocess.run(["porechop", "-i", fastq_path, "-o", out_folder + "/" + \
            fastq_name + "_adapter_clipped.fastq", "-v", "0", "-t", threads])
    # NanoFilt: trim length + Q score + GC%
    nanofilt_names = ['min_len', 'max_len', 'qs', 'min_gc', 'max_gc', 'trim_lead', 'trim_end']
    nanofilt_options = [min_len, max_len, qs, min_gc, max_gc, trim_lead, trim_end]
    trans_opts = {'min_len' : '-l', 'max_len' : '--maxlength', 'qs' : '-q', 
                  'min_gc' : '--minGC', 'max_gc' : '--maxGC', 'trim_lead' : '--headcrop', 
                  'trim_end' : '--tailcrop'}
    if nanofilt_options.count(None) == len(nanofilt_options): 
        print("Options to trim length, Q score, GC% not provided:", ', '.join(nanofilt_names), "\nExiting...")
        exit()
    parse_opts = ["NanoFilt", out_folder + "/" + fastq_name + "_adapter_clipped.fastq"]
    for opt in range(len(nanofilt_options)): 
            if nanofilt_options[opt] is not None:
                    parse_opts.append(trans_opts[nanofilt_names[opt]])
                    parse_opts.append(str(nanofilt_options[opt]))
    #print(' '.join(parse_opts))
    with open(out_folder + "/" + fastq_name + "_trimmed.fastq", "w") as out_file: 
        subprocess.run(parse_opts, stdout = out_file) 

def get_list(folder, pattern): # get path list of files under folder with a specific pattern
    files = subprocess.run(["ls " +  folder + "/" + pattern], \
            shell = True, stdout = subprocess.PIPE)
    files_list = files.stdout.decode().split('\n')
    files_list_final = list(filter(None, files_list))
    return(files_list_final)


def cmd_log(script_file, program, pattern_version, args, parse_args_dic):
        soft_call = subprocess.call(["which", program], stdout = subprocess.PIPE)
        program_install = subprocess.run(["which", program], stdout = subprocess.PIPE)
        program_version = subprocess.run([program, pattern_version], stdout = subprocess.PIPE)
        if soft_call != 0:
                print(program + " is not in your PATH!")
                sys.exit("Install or export " + program + " to your PATH! Exiting...")
        else: 
                log_file = open(script_file + ".log", "a")
                if os.stat(script_file + ".log").st_size == 0: # if file is empty      
                        log_file.write(script_file + " run: " + str(datetime.datetime.now()) + "\n")
                        cmd_line = "./" + script_file
                        for arg in vars(args):
                                if getattr(args, arg) is not None:
                                        cmd_line = cmd_line + " " + parse_args_dic[arg] + " " + getattr(args, arg)
                        log_file.write(script_file + " command-line: " + cmd_line)
                        log_file.write("\n\n---\n\n")
                log_file.write(program_version.stdout.decode())
                log_file.write(program + " installed: " + program_install.stdout.decode() )
                log_file.write("\n---\n\n")
                log_file.close()
