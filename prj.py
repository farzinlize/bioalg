from TrieFind import TrieNode

from getopt import getopt
import sys

DATA_FLU = ['./Influenza/AB470663.fasta']
DATA_MIT = ['./mitochondrial/AF010406_Sheep.fasta']

def read_sequence(sequence_type='FLU', data_index=0):
    sequence = ''

    if sequence_type == 'FLU':        
        dataset_filename = DATA_FLU[data_index]
    elif sequence_type == 'MIT':
        dataset_filename = DATA_MIT[data_index]

    with open(dataset_filename, 'r') as fasta:

        # ommit first line
        fasta.readline()

        for line in fasta:
            sequence += line.replace('\n', '')

    return sequence


def generates_blocks(sequence, block_size, overlap, fasta_filename='blocks.fasta'):
    report = '> block size\n%d\n> overlap\n%d\n'%(block_size, overlap)
    block_index = 0
    block = ''
    for nucleotide in sequence:
        if len(block) == block_size:
            report += '> block_%d\n%s\n'%(block_index, block)
            block_index += 1
            block = block[-overlap:]
            if overlap == 0:
                block = ''
        else:
            block = block + nucleotide
    if len(block) != 0:
        report += '> block_%d\n%s\n'%(block_index, block)
        block_index += 1

    report = '> block count\n%d\n'%block_index + report

    # writing on file
    fasta = open(fasta_filename, 'w')
    fasta.write(report)
    fasta.close


def read_blocks(filename):
    blocks = []

    with open(filename, 'r') as f:
        for _ in range(7):
            f.readline()

        for line in f:
            if '>' not in line:
                blocks += [line]
    
    return blocks


def count_kmer(blocks, k):
    tree = TrieNode()
    for block in blocks:
        kmer_start = 0
        kmer_end = k
        while kmer_end < len(block):
            kmer = block[kmer_start:kmer_end]
            tree.add_frame(kmer)
            kmer_start += 1
            kmer_end +=1
    return tree


def report_count(tree, fasta_filename='kmers.fasta'):
    with open(fasta_filename, 'w') as fasta:
        fasta.write(tree.report(tree.max_min_count()))
    

########################################
#            main call                 #
########################################

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise Exception('request command must be specified (read the description for supported commands)')

    shortopt = 'b:o:f:s:c:k:t:'
    longopts = ['block-size=', 'overlap=', 'filename=', 'sequence-index=', 'category=', 'block-target=']

    # default arguments
    arg_dics = {'block-size':100 , 'overlap':0, 'filename':'default.out', 'k':5,
            'sequence-index':0, 'category':'FLU', 'block-target':'default.out'}

    command = sys.argv[1]

    opts, _ = getopt(sys.argv[2:], shortopt, longopts)
    for o, a in opts:
        if o in ['-b', '--block-size']:
            arg_dics.update({'block-size':int(a)})
        elif o in ['-o', '--overlap']:
            arg_dics.update({'overlap':int(a)})
        elif o in ['-f', '--filename']:
            arg_dics.update({'filename':a})
        elif o in ['-s', '--sequence-index']:
            arg_dics.update({'sequence-index':int(a)})
        elif o in ['-c', '--category']:
            arg_dics.update({'category':a})
        elif o in ['-t', '--block-target']:
            arg_dics.update({'block-target':a})
        elif o == '-k':
            arg_dics.update({'k':int(a)})
    
    if command == 'BLK':
        generates_blocks(read_sequence(arg_dics['category'], arg_dics['sequence-index']), arg_dics['block-size'], arg_dics['overlap'], arg_dics['filename'])
    elif command == 'CNT':
        tree = count_kmer(read_blocks(arg_dics['block-target']), arg_dics['k'])
        report_count(tree, arg_dics['filename'])