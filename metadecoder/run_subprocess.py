import gzip
import os
from multiprocessing import Pool
from shutil import move, rmtree
from stat import S_IRWXU, S_IRGRP, S_IXGRP, S_IROTH, S_IXOTH
from subprocess import DEVNULL, PIPE, STDOUT, run

from .fasta_utility import split_fasta
from .make_file import make_file


def set_permissions(file):
    if not os.access(file, os.R_OK | os.X_OK): # only for linux #
        os.chmod(file, S_IRWXU + S_IRGRP + S_IXGRP + S_IROTH + S_IXOTH) # 0o755 #
    return None


def worker(command, message):
    assert not run(command, stdout = DEVNULL, stderr = DEVNULL).returncode, message
    return None


def get_program_information(program):
    process = run([program, ], stdout = PIPE, stderr = STDOUT)
    if not process.stdout:
        process.stdout = b''
    return process.stdout


def run_fraggenescan(fraggenescan, input_fasta, output_fasta, threads):
    '''
    Run FragGeneScan to predict all protein sequences.
    '''
    set_permissions(fraggenescan)
    open_file = open(input_fasta, 'rb')
    magic_code = open_file.read(2)
    open_file.close()
    if magic_code == b'\x1f\x8b':
        flag = True
        open4r = gzip.open(input_fasta, mode = 'r')
        input_fasta = make_file()
        open4w = open(input_fasta, 'wb')
        while open4w.write(open4r.read(1048576)):
            pass
        open4w.close()
        open4r.close()
    else:
        flag = False
    worker(
        [fraggenescan, '-s', input_fasta, '-o', output_fasta, '-w', '0', '-t', 'complete', '-p', str(threads)],
        'An error has occured while running fraggenescan.'
    )
    if flag:
        os.remove(input_fasta)
    os.remove(output_fasta + '.ffn')
    os.remove(output_fasta + '.out')
    os.replace(output_fasta + '.faa', output_fasta)
    return None


def run_prodigal(prodigal, input_fasta, output_fasta, threads):
    '''
    Run prodigal to predict all protein sequences.
    '''
    input_fastas = list()
    output_fastas = list()
    process_pool = Pool(os.cpu_count())
    for INPUT_FASTA in split_fasta(input_fasta, threads):
        input_fastas.append(INPUT_FASTA)
        output_fastas.append(make_file())
        process_pool.apply_async(
            worker,
            (
                [prodigal, '-a', output_fastas[-1], '-i', INPUT_FASTA, '-p', 'meta'],
                'An error has occured while running prodigal.'
            )
        )
    process_pool.close()
    process_pool.join()
    for INPUT_FASTA in input_fastas:
        os.remove(INPUT_FASTA)
    open4w = open(output_fasta, 'wb')
    for OUTPUT_FASTA in output_fastas:
        open4r = open(OUTPUT_FASTA, 'rb')
        while open4w.write(open4r.read(1024)):
            pass
        open4r.close()
        os.remove(OUTPUT_FASTA)
    open4w.close()
    return None


def run_hmmsearch(hmmsearch, input_hmm, input_fasta, output_file, threads):
    '''
    Run Hmmsearch to map hmms to sequences.
    '''
    set_permissions(hmmsearch)
    input_fastas = list()
    output_files = list()
    process_pool = Pool(os.cpu_count())
    for INPUT_FASTA in split_fasta(input_fasta, threads):
        input_fastas.append(INPUT_FASTA)
        output_files.append(make_file())
        process_pool.apply_async(
            worker,
            (
                [hmmsearch, '--cpu', '1', '--noali', '--domtblout', output_files[-1], input_hmm, INPUT_FASTA],
                'An error has occured while running hmmsearch.'
            )
        )
    process_pool.close()
    process_pool.join()
    for INPUT_FASTA in input_fastas:
        os.remove(INPUT_FASTA)
    open4w = open(output_file, 'wb')
    for OUTPUT_FILE in output_files:
        open4r = open(OUTPUT_FILE, 'rb')
        while open4w.write(open4r.read(1024)):
            pass
        open4r.close()
        os.remove(OUTPUT_FILE)
    open4w.close()
    return None


def run_makeblastdb(makeblastdb, input_fasta, output_file):
    worker(
        [
            makeblastdb, '-in', input_fasta, '-input_type', 'fasta',
            '-dbtype', 'nucl', '-hash_index', '-out', output_file
        ],
        'An error has occured while running makeblastdb.'
    )
    return None


def run_blast(blast, input_database, input_fasta, output_file, threads):
    input_fastas = list(split_fasta(input_fasta, threads))
    output_files = list()
    process_pool = Pool(os.cpu_count())
    for INPUT_FASTA in input_fastas:
        output_files.append(make_file())
        process_pool.apply_async(
            worker,
            (
                [
                    blast, '-db', input_database, '-num_threads', '1',
                    '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send slen pident score',
                    '-strand', 'both', '-query', INPUT_FASTA, '-out', output_files[-1]
                ],
                'An error has occured while running blast.'
            )
        )
    process_pool.close()
    process_pool.join()
    open4w = open(output_file, 'wb')
    for INPUT_FASTA in input_fastas:
        os.remove(INPUT_FASTA)
    for OUTPUT_FILE in output_files:
        open4r = open(OUTPUT_FILE, 'rb')
        while open4w.write(open4r.read(1024)):
            pass
        open4r.close()
        os.remove(OUTPUT_FILE)
    open4w.close()
    return None


def run_idbaud(idbaud, input_fastqs_1, input_fastqs_2, output_fasta, threads):
    input_fasta = make_file()
    open4w = open(input_fasta, 'w')
    for input_fastq_1, input_fastq_2 in zip(input_fastqs_1, input_fastqs_2):
        open4r1 = open(input_fastq_1, 'r')
        open4r2 = open(input_fastq_2, 'r')
        while True:
            read1 = open4r1.readline().rstrip('\n')
            read2 = open4r2.readline().rstrip('\n')
            if not (read1 and read2):
                break
            open4w.write('>' + read1.split(' ', maxsplit = 1)[0][1:] + ' 1\n')
            open4w.write(open4r1.readline())
            open4w.write('>' + read2.split(' ', maxsplit = 1)[0][1:] + ' 2\n')
            open4w.write(open4r2.readline())
            open4r1.readline()
            open4r1.readline()
            open4r2.readline()
            open4r2.readline()
        open4r1.close()
        open4r2.close()
    open4w.close()
    worker(
        [
            idbaud, '--pre_correction', '--num_threads', str(threads),
            '--read', input_fasta, '--out', output_fasta + '.idbaud'
        ],
        'An error has occured while running idba_ud.'
    )
    os.remove(input_fasta)
    move(os.path.join(output_fasta + '.idbaud', 'contig.fa'), output_fasta)
    rmtree(output_fasta + '.idbaud')
    return None


def run_bowtie2_build(bowtie2_build, input_fasta, output_database, threads):
    worker(
        [bowtie2_build, '--threads', str(threads), input_fasta, output_database],
        'An error has occured while running bowtie2-build.'
    )
    return None


def run_bowtie2(bowtie2, input_database, input_fastq_1, input_fastq_2, output_sam, threads):
    worker(
        [
            bowtie2, '--threads', str(threads), '-x', input_database,
            '-1', input_fastq_1, '-2', input_fastq_2, '-S', output_sam
        ],
        'An error has occured while running bowtie2.'
    )
    return None


def run_spades(spades, input_fastqs_1, input_fastqs_2, output_fasta, threads):
    #create yami
    input_yami = make_file()
    open4w = open(input_yami, 'w')
    open4w.write('[\n')
    open4w.write('    {\n')
    open4w.write('        orientation: "fr",\n')
    open4w.write('        type: "paired-end",\n')
    open4w.write('        left reads: [{0}],\n'.format(','.join('"' + input_fastq_1 + '"' for input_fastq_1 in input_fastqs_1)))
    open4w.write('        right reads: [{0}]\n'.format(','.join('"' + input_fastq_2 + '"' for input_fastq_2 in input_fastqs_2)))
    open4w.write('    }\n')
    open4w.write(']\n')
    open4w.close()
    #run spades.py
    worker(
        [spades, '--meta', '--threads', str(threads), '--dataset', input_yami, '-o', output_fasta + '.metaspades'],
        'An error has occured while running metaspades.'
    )
    os.remove(input_yami)
    move(os.path.join(output_fasta + '.metaspades', 'contigs.fasta'), output_fasta)
    rmtree(output_fasta + '.metaspades')
    return None