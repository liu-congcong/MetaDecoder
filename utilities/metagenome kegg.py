#!/usr/bin/env python3
import os
import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from json import load
from math import ceil, inf
from multiprocessing import Process
from shutil import which
from subprocess import DEVNULL, run
from tempfile import mkstemp

import numpy


def __init__(parameters):
    parser = ArgumentParser(
        prog = sys.argv[0],
        formatter_class = RawTextHelpFormatter,
        description = 'KEGG for a metagenome.',
        epilog = 'Liucongcong congcong_liu@icloud.com.'
    )
    subparsers = parser.add_subparsers(
        title = 'program', dest = 'program', required = True, help = 'Run %(prog)s mapping|abundance|pathway --help for help.'
    )

    mapping_parser = subparsers.add_parser(
        'mapping', formatter_class = RawTextHelpFormatter, help = 'Map hmms to assembliy to generate the kmap file.\nUsage: %(prog)s mapping -assembly ASSEMBLY -output KMAP -threads 30.'
    )
    mapping_parser.add_argument(
        '-prodigal', '--prodigal', type = str, default = which('prodigal'), required = False, metavar = 'file',
        help = 'The path of "prodigal".\nDefault: /usr/bin/env prodigal.'
    )
    mapping_parser.add_argument(
        '-hmmsearch', '--hmmsearch', type = str, default = which('hmmsearch'), required = False, metavar = 'file',
        help = 'The path of "hmmsearch".\nDefault: /usr/bin/env hmmsearch.'
    )
    mapping_parser.add_argument(
        '-ko_list', '--ko_list', type = str, default = os.getenv('KO_LIST'), required = False, metavar = 'file',
        help = 'The path of "ko_list" file.\nwget https://www.genome.jp/ftp/db/kofam/ko_list.gz\ngunzip ko_list.gz\nDefault: /usr/bin/env KO_LIST.'
    )
    mapping_parser.add_argument(
        '-k_hmm', '--k_hmm', type = str, default = os.getenv('K_HMM'), required = False, metavar = 'file',
        help = 'The path of single combined hmm file.\nwget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz\ntar xvf profiles.tar.gz\ncd profiles\ncat K*.hmm > ../k.hmm\nDefault: /usr/bin/env K_HMM.'
    )
    mapping_parser.add_argument(
        '-assembly', '--assembly', type = str, required = True, metavar = 'file',
        help = 'The path of fasta format assembly file.'
    )
    mapping_parser.add_argument(
        '-output', '--output', type = str, required = True, metavar = 'file'
    )
    mapping_parser.add_argument(
        '-threads', '--threads', type = int, default = 1, required = False, metavar = 'threads'
    )

    abundance_parser = subparsers.add_parser('abundance', formatter_class = RawTextHelpFormatter, help = 'Calculate abundance in a kmap file on the basis of sorted and indexed bam file.\nUsage: %(prog)s abundance -bam SORTED.BAM -kmap KMAP -output AKMAP -threads 30.')
    abundance_parser.add_argument(
        '-samtools', '--samtools', type = str, default = which('samtools'), required = False, metavar = 'file',
        help = 'The path of "samtools".\nDefault: /usr/bin/env samtools.'
    )
    abundance_parser.add_argument(
        '-bam', '--bam', type = str, required = True, metavar = 'file',
        help = 'The path of sorted and indexed bam file.'
    )
    abundance_parser.add_argument(
        '-kmap', '--kmap', type = str, required = True, metavar = 'file',
        help = 'The path of kmap file.'
    )
    abundance_parser.add_argument(
        '-output', '--output', type = str, required = True, metavar = 'file'
    )
    abundance_parser.add_argument(
        '-threads', '--threads', type = int, default = 1, required = False, metavar = 'threads'
    )

    pathway_parser = subparsers.add_parser('pathway', formatter_class = RawTextHelpFormatter, help = 'Construct the KEGG pathways of clusters.\nUsage: %(prog)s pathway -clusters *CLUSTER -kmap KMAP -output PATHWAY.')
    pathway_parser.add_argument(
        '-ko_json', '--ko_json', type = str, default = os.getenv('KO_JSON'), required = False, metavar = 'file',
        help = 'The path of "ko*.json" file.\nwget -O ko00001.json "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001&format=json"\nDefault: /usr/bin/env KO_JSON.'
    )
    pathway_parser.add_argument(
        '-clusters', '--clusters', type = str, nargs = '+', required = True, metavar = 'file',
        help = 'The fasta format cluster files.'
    )
    pathway_parser.add_argument(
        '-kmap', '--kmap', type = str, required = True, metavar = 'file',
        help = 'The kmap file with abundance information.'
    )
    pathway_parser.add_argument(
        '-output', '--output', type = str, required = True, metavar = 'file'
    )
    return parser.parse_args(parameters)


def train_model(prodigal, input_fasta):
    output_model = make_file()
    run_process = run([prodigal, '-g', '11', '-i', input_fasta, '-m', '-p', 'single', '-t', output_model], stdout = DEVNULL, stderr = DEVNULL)
    assert not run_process.returncode, 'Run Prodigal field.'
    return output_model


def make_file():
    file_descriptor, file_name = mkstemp(dir = os.getcwd())
    os.close(file_descriptor)
    os.remove(file_name)
    return file_name


def split_fasta(input_fasta, output_files):
    total_size = os.path.getsize(input_fasta)
    block_size = ceil(total_size / output_files) # block_size <= total_size #
    file_position = 0
    file_position_ = 0
    open4r = open(input_fasta, 'rb')
    while file_position < total_size:
        line = open4r.readline()
        file_position += len(line)
        if line.startswith(b'>'):
            file_position -= len(line)
            if file_position > 0:
                open4r.seek(file_position_, os.SEEK_SET)
                output_file = make_file()
                open4w = open(output_file, 'wb')
                while file_position_ < file_position:
                    file_position_ += open4w.write(open4r.read(min(10485760, file_position - file_position_)))
                open4w.close()
                yield output_file
                # file_position_ will be equal to file_position, open4r.tell() will be equal to file_position #
            file_position = open4r.seek(min(file_position + block_size, total_size), os.SEEK_SET)
    open4r.seek(file_position_, os.SEEK_SET)
    output_file = make_file()
    open4w = open(output_file, 'wb')
    while file_position_ < file_position:
        file_position_ += open4w.write(open4r.read(min(10485760, file_position - file_position_)))
    open4w.close()
    yield output_file
    open4r.close()
    return None


def predict_protein(prodigal, model_file, input_fasta, output_fasta):
    run_process = run([prodigal, '-a', output_fasta, '-g', '11', '-i', input_fasta, '-m', '-p', 'single', '-t', model_file], stdout = DEVNULL, stderr = DEVNULL)
    assert not run_process.returncode, 'Run Prodigal field.'
    os.remove(input_fasta)
    return None


def combine_files(input_files, output_file):
    sequence_id2position = dict()
    open4w = open(output_file, 'w')
    for input_file in input_files:
        open4r = open(input_file, 'r')
        for line in open4r:
            line = line.rstrip('\n')
            if line.startswith('>'):
                sequence_id, start, end, _ = line[1 : ].split(' # ', maxsplit = 3)
                sequence_id2position[sequence_id] = (start, end)
            open4w.write(line + '\n')
        open4r.close()
        os.remove(input_file)
    open4w.close()
    return sequence_id2position


def run_hmmsearch(hmmersearch, hmm_file, protein_file, hit_file):
    run_process = run([hmmersearch, '--tblout', hit_file, '--noali', '--notextw', '-T', '0', '--cpu', '1', hmm_file, protein_file], stdout = DEVNULL)
    assert not run_process.returncode, 'Run hmmsearch field.'
    os.remove(protein_file)
    return None


def read_ko_list(input_file):
    k_number2score = dict()
    open4r = open(input_file, 'r')
    open4r.readline()
    # knum threshold score_type profile_type F-measure nseq nseq_used alen mlen eff_nseq re/pos definition #
    for line in open4r:
        k_number, score, score_type, _ = line.rstrip('\n').split('\t', maxsplit = 3)
        if score_type == 'full':
            k_number2score[k_number] = (5, float(score))
        elif score_type == 'domain':
            k_number2score[k_number] = (8, float(score))
        else:
            k_number2score[k_number] = (5, -inf)
    open4r.close()
    return k_number2score


def read_hit_file(input_file, score_hash):
    open_file = open(input_file, 'r')
    for line in open_file:
        if not line.startswith('#'):
            lines = line.rstrip('\n').split(maxsplit = 9)
            score_index, score_threshold = score_hash[lines[2]]
            if float(lines[score_index]) >= score_threshold:
                yield (lines[0], lines[2])
    open_file.close()
    return None


def output_kmap(hit_files, k_number2score, sequence_id2position, output_file):
    open_file = open(output_file, 'w')
    open_file.write('\t'.join(['Sequence ID', 'Start', 'End', 'K number']) + '\n')
    for hit_file in hit_files:
        for sequence_id, k_number in read_hit_file(hit_file, k_number2score):
            start, end = sequence_id2position[sequence_id]
            open_file.write('\t'.join([sequence_id.rsplit('_', maxsplit = 1)[0], start, end, k_number]) + '\n')
        os.remove(hit_file)
    open_file.close()
    return None


def split_kmap(input_kmap, output_files):
    open4r = open(input_kmap, 'r')
    open4r.readline()
    total_lines = 0
    for line in open4r:
        total_lines += 1
    open4r.seek(0, 0)
    open4r.readline()
    lines_per_file = ceil(total_lines / output_files)
    for output_file_index in range(0, total_lines, lines_per_file):
        output_file = make_file()
        open4w = open(output_file, 'w')
        for read_line in range(lines_per_file):
            line = open4r.readline().rstrip('\n')
            if not line:
                break
            sequence_id, start, end, k_number  = line.split('\t')
            open4w.write('\t'.join([sequence_id, str(int(start) - 1), end, k_number]) + '\n')
        open4w.close()
        yield output_file
    open4r.close()
    return None


def run_samtools(samtools, input_bam, input_bed, output_file):
    open_file = open(output_file, 'wb')
    run_process = run([samtools, 'bedcov', input_bed, input_bam], stdout = open_file)
    assert not run_process.returncode, 'Run samtools field.'
    os.remove(input_bed)
    return None


def output_kmap_with_abundance(input_files, output_file):
    open4w = open(output_file, 'w')
    open4w.write('\t'.join(['Sequence ID', 'Start', 'End', 'K number', 'Abundance']) + '\n')
    for input_file in input_files:
        open4r = open(input_file, 'r')
        for line in open4r:
            lines = line.rstrip('\n').split('\t')
            start = int(lines[1])
            end = int(lines[2])
            open4w.write('\t'.join([lines[0], str(start + 1), lines[2], lines[3], str(float(lines[4]) / (end - start))]) + '\n')
        open4r.close()
        os.remove(input_file)
    open4w.close()
    return None


def read_kmap(input_file):
    k_index2sequence_id2abundance = list()
    k_number2k_index = dict()
    k_index = 0
    open_file = open(input_file, 'r')
    open_file.readline()
    for line in open_file:
        lines = line.rstrip('\n').split('\t')
        sequence_id = lines[0]
        k_number = lines[3]
        if k_number not in k_number2k_index: # a new k_number #
            k_number2k_index[k_number] = k_index
            k_index2sequence_id2abundance.append(dict())
            k_index += 1
        sequence_id2abundance = k_index2sequence_id2abundance[k_number2k_index[k_number]]
        sequence_id2abundance[sequence_id] = sequence_id2abundance.get(sequence_id, 0.0) + float(lines[4])
    open_file.close()
    return (k_number2k_index, k_index2sequence_id2abundance)


def read_ko_json(input_file, k_number2k_index):
    levela2levelb2levelc2k_indices = dict()
    open_file = open(input_file, 'r')
    root_node = load(open_file)
    open_file.close()
    for node_a in root_node.get('children', list()):
        level_a = node_a['name'].split(maxsplit = 1)[1]
        if level_a in ('Metabolism', 'Genetic Information Processing', 'Environmental Information Processing', 'Cellular Processes', 'Organismal Systems', 'Human Diseases'):
            levela2levelb2levelc2k_indices[level_a] = dict()
            levelb2levelc2k_indices = levela2levelb2levelc2k_indices[level_a]
            for node_b in node_a.get('children', list()):
                level_b = node_b['name'].split(maxsplit = 1)[1]
                levelb2levelc2k_indices[level_b] = dict()
                levelc2k_indices = levelb2levelc2k_indices[level_b]
                for node_c in node_b.get('children', list()):
                    level_c = node_c['name'].split(maxsplit = 1)[1]
                    levelc2k_indices[level_c] = list()
                    k_indices = levelc2k_indices[level_c]
                    for node_d in node_c.get('children', list()):
                        k_number = node_d['name'].split(maxsplit = 1)[0]
                        if k_number in k_number2k_index:
                            k_indices.append(k_number2k_index[k_number])
    return levela2levelb2levelc2k_indices


def read_fastas(input_files):
    sequence_id2file_index = dict()
    for file_index, input_file in enumerate(input_files):
        open_file = open(input_file, 'r')
        for line in open_file:
            if line.startswith('>'):
                sequence_id2file_index[line.rstrip('\n').split(maxsplit = 1)[0][1 : ]] = file_index
        open_file.close()
    return sequence_id2file_index


def calculate_relative_abundance(abundance, k_index2sequence_id2abundance, sequence_id2cluster_index):
    for k_index, sequence_id2abundance in enumerate(k_index2sequence_id2abundance):
        for sequence_id, abundance_ in sequence_id2abundance.items():
            if sequence_id in sequence_id2cluster_index:
                cluster_index = sequence_id2cluster_index[sequence_id]
                abundance[k_index][cluster_index] += abundance_
    abundance /= numpy.sum(abundance, axis = 0) + numpy.finfo(numpy.float64).eps
    abundance *= 100
    return abundance


def output_pathway(abundance, levela2levelb2levelc2k_indices, cluster_names, output_file):
    levela_k_indices = list()
    levelb_k_indices = list()
    open_file = open(output_file, 'w')
    open_file.write('\t'.join(['KEGG level', 'KEGG level 1', 'KEGG level 2', 'KEGG level 3'] + cluster_names) + '\n')
    for levela, levelb2levelc2k_indices in levela2levelb2levelc2k_indices.items():
        for levelb, levelc2k_indices in levelb2levelc2k_indices.items():
            for levelc, k_indices in levelc2k_indices.items():
                open_file.write('\t'.join(['3', levela, levelb, levelc] + numpy.sum(abundance[k_indices], axis = 0).astype(numpy.str_).tolist())+ '\n')
                levelb_k_indices.extend(k_indices)
            open_file.write('\t'.join(['2', levela, levelb, '-'] + numpy.sum(abundance[list(set(levelb_k_indices))], axis = 0).astype(numpy.str_).tolist())+ '\n')
            levela_k_indices.extend(levelb_k_indices)
            levelb_k_indices.clear()
        open_file.write('\t'.join(['1', levela, '-', '-'] + numpy.sum(abundance[list(set(levela_k_indices))], axis = 0).astype(numpy.str_).tolist())+ '\n')
        levela_k_indices.clear()
    open_file.close()


if __name__ == '__main__':
    parameters = __init__(sys.argv[1 : ])

    if parameters.program == 'mapping':
        assert parameters.prodigal and os.access(parameters.prodigal, os.R_OK), 'You need to specify the path of "prodigal".'
        assert parameters.hmmsearch and os.access(parameters.prodigal, os.R_OK), 'You need to specify the path of "hmmsearch".'
        assert parameters.ko_list and os.access(parameters.ko_list, os.R_OK), 'You need to specify the path of "ko_list".'
        assert parameters.k_hmm and os.access(parameters.k_hmm, os.R_OK), 'You need to specify the path of "k.hmm".'
        assert os.access(parameters.assembly, os.R_OK), 'The {0} file is not readable.'.format(parameters.assembly)

        # Train model using Prodigal. #
        model_file = train_model(parameters.prodigal, parameters.assembly)

        processes = list()
        protein_files = list()
        for fasta_file in split_fasta(parameters.assembly, parameters.threads):
            protein_files.append(make_file())
            # Predict protein sequences using Prodigal. #
            processes.append(
                Process(
                    target = predict_protein,
                    args = (parameters.prodigal, model_file, fasta_file, protein_files[-1])
                )
            )
            processes[-1].start()
        for process in processes:
            process.join()
        os.remove(model_file)
        sequence_id2position = combine_files(protein_files, parameters.output + '.proteins.fasta')

        hit_files = list()
        for fasta_file in split_fasta(parameters.output + '.proteins.fasta', parameters.threads):
            hit_files.append(make_file())
            # Map proteins to KEGG hmms using hmmsearch. #
            processes.append(
                Process(
                    target = run_hmmsearch,
                    args = (parameters.hmmsearch, parameters.k_hmm, fasta_file, hit_files[-1])
                )
            )
            processes[-1].start()
        for process in processes:
            process.join()

        k_number2score = read_ko_list(parameters.ko_list)

        output_kmap(hit_files, k_number2score, sequence_id2position, parameters.output)

    elif parameters.program == 'abundance':
        assert parameters.samtools, 'You need to specify the path of "samtools".'
        assert os.access(parameters.kmap, os.R_OK), 'The {0} file is not readable.'.format(parameters.kmap)
        assert os.access(parameters.bam, os.R_OK), 'The {0} file is not readable.'.format(parameters.bam)

        processes = list()
        kmap_files_with_abundances = list()
        for bed_file in split_kmap(parameters.kmap, parameters.threads):
            kmap_files_with_abundances.append(make_file())
            # Calculate abundance using samtools. #
            processes.append(
                Process(
                    target = run_samtools,
                    args = (parameters.samtools, parameters.bam, bed_file, kmap_files_with_abundances[-1])
                )
            )
            processes[-1].start()
        for process in processes:
            process.join()

        output_kmap_with_abundance(kmap_files_with_abundances, parameters.output)

    else:
        assert parameters.ko_json and os.access(parameters.ko_json, os.R_OK), 'You need to specify the path of "ko*.json".'
        assert os.access(parameters.kmap, os.R_OK), 'The {0} file is not readable.'.format(parameters.kmap)
        for cluster_file in parameters.clusters:
            assert os.access(cluster_file, os.R_OK), 'The {0} file is not readable.'.format(cluster_file)

        k_number2k_index, k_index2sequence_id2abundance = read_kmap(parameters.kmap)
        levela2levelb2levelc2k_indices = read_ko_json(parameters.ko_json, k_number2k_index)
        sequence_id2cluster_index = read_fastas(parameters.clusters)
        abundance = numpy.zeros((len(k_number2k_index), len(parameters.clusters)), dtype = numpy.float64)
        calculate_relative_abundance(abundance, k_index2sequence_id2abundance, sequence_id2cluster_index)
        output_pathway(abundance, levela2levelb2levelc2k_indices, parameters.clusters, parameters.output)
