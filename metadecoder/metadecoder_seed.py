import os
from datetime import datetime
from math import inf

from metadecoder.fasta_utility import read_fasta_file
from metadecoder.make_file import make_file
from metadecoder.run_subprocess import get_program_information, run_fraggenescan, run_hmmsearch, run_prodigal


def parse_sequence_id(file, trim):
    '''
    fraggenescan: >id_start_end_strand
    prodigal: >id_index
    '''
    gene = 1
    output = os.path.basename(file) + '.proteins'
    open_file = open(output, 'w')
    for sequence_id, sequence in read_fasta_file(file):
        open_file.write('>' + str(gene) + '_' + sequence_id.rsplit('_', maxsplit = trim)[0] + '\n')
        open_file.write(sequence + '\n')
        gene += 1
    open_file.close()
    return output


def read_hmm_file(input_hmm):
    model2tc = dict()
    open_file = open(input_hmm, 'r')
    for line in open_file:
        if line.startswith('NAME'):
            model = line.rstrip('\n').split()[1]
        elif line.startswith('TC'):
            score1, score2 = line.rstrip(';\n').split()[1 : ]
            model2tc[model] = (float(score1), float(score2))
    open_file.close()
    return model2tc


def get_seeds(file, model2tc, coverage, accuracy, output):
    model2sequences = dict()
    open_file = open(file, 'r', encoding = 'utf-8')
    for line in open_file:
        if not line.startswith('#'):
            lines = line.rstrip('\n').split()
            tc = model2tc.get(lines[3], (-inf, -inf))
            if float(lines[7]) >= tc[0] and float(lines[13]) >= tc[1] and (int(lines[16]) - int(lines[15]) + 1) / int(lines[5]) >= coverage and float(lines[21]) >= accuracy:
                model2sequences.setdefault(lines[3], list()).append(lines[0].split('_', maxsplit = 1)[1])
    open_file.close()
    open_file = open(output, 'w')
    for model, sequences in model2sequences.items():
        open_file.write(model + '\t' + '\t'.join(list(set(sequences))) + '\n')
    open_file.close()


def main(parameters):
    PROTEIN = make_file()
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Identifying protein sequences.', flush = True)
    if b'fraggenescan' in get_program_information(parameters.protein_predictor).lower():
        # FragGeneScan #
        run_fraggenescan(parameters.protein_predictor, parameters.fasta, PROTEIN, parameters.threads)
        protein = parse_sequence_id(PROTEIN, 3)
    else:  # prodigal #
        run_prodigal(parameters.protein_predictor, parameters.fasta, PROTEIN, parameters.threads)
        protein = parse_sequence_id(PROTEIN, 1)
    os.remove(PROTEIN)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Mapping marker genes to protein sequences.', flush = True)
    hmmsearch_output = make_file()
    run_hmmsearch(parameters.hmmsearch, parameters.marker, protein, hmmsearch_output, parameters.threads)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Done.', flush = True)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Writing to file.', flush = True)
    os.remove(protein)
    model2tc = read_hmm_file(parameters.marker)
    get_seeds(hmmsearch_output, model2tc, parameters.coverage, parameters.accuracy, parameters.output)
    os.remove(hmmsearch_output)
    print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '->', 'Finished.', flush = True)