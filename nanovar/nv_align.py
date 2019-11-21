"""
Functions for index and blast alignment.

Copyright (C) 2019 Tham Cheng Yong

This file is part of NanoVar.

NanoVar is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

NanoVar is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with NanoVar.  If not, see <https://www.gnu.org/licenses/>.
"""


import os
import logging
import distutils.spawn
from subprocess import Popen, PIPE, STDOUT


# Check paths and executables
def check_exe(path, exe):
    if path is None:
        if distutils.spawn.find_executable(exe):
            return exe
        else:
            logging.critical("Error: %s executable is not in PATH" % exe)
            raise Exception("Error: %s executable is not in PATH" % exe)
    else:
        if distutils.spawn.find_executable(path):
            return path
        else:
            logging.critical("Error: %s path do not exist" % path)
            raise Exception("Error: %s path do not exist" % path)


def log_subprocess(out):
    for line in iter(out.readline, ''):
        if line != '\n':
            logging.debug(line.strip())


# Create blast index
def blast_index(ref, mdb):
    process = Popen([mdb, '-in', ref, '-input_type', 'fasta', '-dbtype', 'nucl'], universal_newlines=True,
                    stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    # Process finished
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: makeblastdb failed")
        raise Exception("Error: makeblastdb failed, see log")


# Create FMD index
def fmd_index(ref, hsb):
    process = Popen([hsb, 'index', ref], universal_newlines=True, stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: hs-blastn index failed")
        raise Exception("Error: hs-blastn index failed, see log")


# Create windowmasker counts
def window_counts(ref, counts_path, wmk):
    process = Popen([wmk, '-in', ref, '-infmt', 'blastdb', '-mk_counts', '-out', counts_path],
                    universal_newlines=True,
                    stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    # Process finished
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: windowmasker counts failed")
        raise Exception("Error: windowmasker counts failed, see log")


# Create counts binary
def window_binary(counts_path, obinary_path, wmk):
    process = Popen([wmk, '-in', counts_path, '-sformat', 'obinary', '-out', obinary_path, '-convert'],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT)
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: windowmasker obinary failed")
        raise Exception("Error: windowmasker obinary failed, see log")


# Run hs-blastn
def hsblastn(ref, read_para, tab_path, obinary_path, threads, hsb):
    # process = Popen(['hs-blastn', 'align', '-db', ref, '-window_masker_db', obinary_path, '-query', read_para, '-out',
    #                  tab_path, '-outfmt', '6', '-num_threads', str(threads), '-max_target_seqs', '3', '-gapopen', '0',
    #                  '-gapextend', '4', '-penalty', '-3', '-reward', '2'],
    #                 universal_newlines=True, stdout=PIPE, stderr=STDOUT)
    process = Popen([hsb + ' align -db ' + ref + ' -window_masker_db ' + obinary_path + ' -query ' + read_para + ' -out '
                     + tab_path + ' -outfmt 6 -num_threads ' + str(threads) + ' -max_target_seqs 3 -gapopen 0 '
                     '-gapextend 4 -penalty -3 -reward 2'],
                    universal_newlines=True, stdout=PIPE, stderr=STDOUT, shell=True, executable='/bin/bash')
    # cmd = process.args[0]
    with process.stdout:
        log_subprocess(process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        logging.critical("Error: hs-blastn alignment failed")
        raise Exception("Error: hs-blastn alignment failed, see log")


# Alignment pipeline
def master_align(ref, wk_dir, ref_name, read_para, read_name, threads, mode, mdb, wmk, hsb):
    counts_path = os.path.join(wk_dir, ref_name + '.counts')
    obinary_path = counts_path + '.obinary'
    tab_path = os.path.join(wk_dir, '%s-%s-blast.tab' % (read_name, ref_name))
    if mode:
        logging.info('Building blast database')
        blast_index(ref, mdb)
        logging.info('Windowmasker counts')
        window_counts(ref, counts_path, wmk)
        logging.info('Windowmasker binary')
        window_binary(counts_path, obinary_path, wmk)
        logging.info('Make FMD-index')
        fmd_index(ref, hsb)
        logging.info('hs-blastn alignment')
        hsblastn(ref, read_para, tab_path, obinary_path, threads, hsb)
        return ["hs-blastn align -db " + ref + " -window_masker_db obinary_path -query " + read_para + " -out " + tab_path +
                " -outfmt 6 -num_threads " + str(threads) + " -max_target_seqs 3 -gapopen 0 -gapextend 4 -penalty -3 -reward 2",
                tab_path]
    else:
        if os.path.isfile(ref + '.nsq') and os.path.isfile(ref + '.nhr') and os.path.isfile(ref + '.nin'):
            # Only check if index files exist, they might be empty
            logging.debug("Make blast index skipped")
        else:
            logging.info('Building blast database')
            blast_index(ref, mdb)
        if os.path.isfile(counts_path):  # Only check if files exist, they might be empty
            logging.debug("Windowmasker counts skipped")
        else:
            logging.info('Windowmasker counts')
            window_counts(ref, counts_path, wmk)
        if os.path.isfile(obinary_path):  # Only check if files exist, they might be empty
            logging.debug("Windowmasker obinary skipped")
        else:
            logging.info('Windowmasker binary')
            window_binary(counts_path, obinary_path, wmk)
        if os.path.isfile(ref + '.bwt') and os.path.isfile(ref + '.header') and os.path.isfile(ref + '.sa') and os.path.isfile(
                ref + '.sequence'):  # Only check if index files exist, they might be empty
            logging.debug("Make hs-blastn FMD-index skipped")
        else:
            logging.info('Make FMD-index')
            fmd_index(ref, hsb)
        if os.path.isfile(tab_path):  # Only check if index file exist, they might be empty
            logging.debug("hs-blastn alignment skipped")
        else:
            logging.info('hs-blastn alignment')
            hsblastn(ref, read_para, tab_path, obinary_path, threads, hsb)
        return ["hs-blastn align -db " + ref + " -window_masker_db obinary_path -query " + read_para + " -out " + tab_path +
                " -outfmt 6 -num_threads " + str(threads) + " -max_target_seqs 3 -gapopen 0 -gapextend 4 -penalty -3 -reward 2",
                tab_path]
