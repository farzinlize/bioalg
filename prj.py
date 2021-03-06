from TrieFind import TrieNode

from deprecated import deprecated
from functools import reduce
from getopt import getopt
import sys

# dataset file addresses
DATASET = {
    'FLU':
    ['AB470663.fasta', 'AB546159.fasta', 'AB684161.fasta', 'AF509102.fasta', 'AM157358.fasta', 
    'AM914017.fasta', 'AY646080.fasta', 'CY005540.fasta', 'CY014788.fasta', 'CY039321.fasta', 
    'CY076231.fasta', 'CY129336.fasta', 'CY138562.fasta', 'CY140047.fasta', 'CY149630.fasta', 
    'CY186004.fasta', 'DQ017487.fasta', 'EF541464.fasta', 'EU026046.fasta', 'EU500854.fasta', 
    'EU635875.fasta', 'FJ357114.fasta', 'FM177121.fasta', 'GQ411894.fasta', 'GU186511.fasta', 
    'HM370969.fasta', 'HQ185381.fasta', 'HQ185383.fasta', 'HQ897966.fasta', 'JF699677.fasta', 
    'JX081142.fasta', 'KC608160.fasta', 'KC609801.fasta', 'KF259688.fasta', 'KF259734.fasta', 
    'KF572435.fasta', 'KF938945.fasta', 'KM244078.fasta'],

    'MIT':
    ['AF010406_Sheep.fasta', 'AF303110_Brown_Bear.fasta', 'AF303111_Polar_Bear.fasta', 
    'AF533441_Goat.fasta', 'AJ001562_Dormouse.fasta', 'AJ001588_Rabbit.fasta', 'AJ002189_Pig.fasta', 
    'AJ238588_Squirrel.fasta', 'AY488491_Buffalo.fasta', 'AY863426_Ver_Monkey.fasta', 'D38113_Com_Chim.fasta', 
    'D38114_Gorilla.fasta', 'D38115_Bor_Oran.fasta', 'D38116_Pig_Chim.fasta', 'DQ402478_Black_Bear.fasta', 
    'EF212882_Giant_panda.fasta', 'EF551002_Leopard.fasta', 'EF551003_Tiger.fasta', 'EU442884_Wolf.fasta', 
    'NC_001321_Fin_Whale.fasta', 'NC_001640_horse.fasta', 'NC_001788_donkey.fasta', 'NC_002083_Sum_Oran.fasta', 
    'NC_002764_Macaca_Thibet.fasta', 'NC_005268_Bowhead_Whale.fasta', 'NC_005270_Gray_Whale.fasta', 
    'NC_005275_Indus_river_dolphin.fasta', 'NC_006931_North_pacific_whale.fasta', 'NC_007441_Chiru.fasta', 
    'NC_008830_Common_warthog.fasta', 'NC_010640_Taiwan_serow.fasta', 'U20753_Cat.fasta', 'U96639_Dog.fasta', 
    'V00654_Cow.fasta', 'V00662_Human.fasta', 'X72204_Blue_Whale.fasta', 'X88898_Hedgehog.fasta', 
    'X97336_Indian_Rhino.fasta', 'X99256_Gibbon.fasta', 'Y07726_White_Rhino.fasta', 'Y18001_Baboon.fasta']
}

# extract sequence
def read_sequence(sequence_type='FLU', data_index=0):
    sequence = ''

    if sequence_type == 'FLU':        
        dataset_filename = './Influenza/' + DATASET['FLU'][data_index]
    elif sequence_type == 'MIT':
        dataset_filename = './mitochondrial/' + DATASET['MIT'][data_index]

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

@deprecated('related to phase #2 of project')
def count_kmer_blocks(blocks, k):
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


def count_kmer(sequence, k):
    tree = TrieNode()
    kmer_start = 0
    kmer_end = k
    while kmer_end < len(sequence):
        kmer = sequence[kmer_start:kmer_end]
        tree.add_frame(kmer)
        kmer_start += 1
        kmer_end +=1
    return tree


def report_count(tree, fasta_filename='kmers.fasta'):
    with open(fasta_filename, 'w') as fasta:
        fasta.write(tree.report(tree.max_min_count()))
    

def distance_manhattan(sequence_1, sequence_2, k):
    tree_1 = count_kmer(sequence_1, k)
    tree_2 = count_kmer(sequence_2, k)
    
    allkmers_1 = len(sequence_1)-k+1
    allkmers_2 = len(sequence_2)-k+1

    ma_diff = lambda count_1, count_2: abs((count_1/allkmers_1)-(count_2/allkmers_2))

    diffs = tree_1.measure_diff(tree_2, ma_diff)
    return reduce(lambda a,b:a+b, diffs + [count/allkmers_2 for count in tree_2.extract_counts()])


def block_distance_matrix(blocks_1, blocks_2, k, report_filename):
    # matrix = [[distance_manhattan(blocks_1, blocks_2) for block_2 in blocks_2] for block_1 in blocks_1]
    with open(report_filename, 'w') as report:
        for block_1 in blocks_1:
            for block_2 in blocks_2:
                d = distance_manhattan(block_1, block_2, k)
                report.write('%.6f   '%d)
            report.write('\n')


def dataset_distance_matrix(category, k, report_filename):
    with open(report_filename, 'w') as report:
        for i in range(len(DATASET[category])):
            for j in range(len(DATASET[category])):
                d = distance_manhattan(read_sequence(category, i), read_sequence(category, j), k)
                report.write('%.6f   '%d)
            report.write('\n')


########################################
#            main call                 #
########################################

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise Exception('request command must be specified (read the description for supported commands)')

    shortopt = 'b:o:f:s:c:k:t:n:h:'
    longopts = ['block-size=', 'overlap=', 'filename=', 'sequence-index=', 
            'category=', 'block-target=', 'other-index=', 'other-target=']

    # default arguments
    arg_dics = {'block-size':100 , 'overlap':0, 'filename':'default.out', 'k':5,
            'sequence-index':0, 'category':'FLU', 'block-target':'default.out', 
            'other-index':0, 'other-target':'default.out'}

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
        elif o in ['-n', '--other-index']:
            arg_dics.update({'other-index':int(a)})
        elif o in ['-h', '--other-target']:
            arg_dics.update({'other-target':a})
    
    if command == 'BLK':
        generates_blocks(
            read_sequence(arg_dics['category'], arg_dics['sequence-index']), 
            arg_dics['block-size'], 
            arg_dics['overlap'], 
            arg_dics['filename'])
    elif command == 'CNT':
        tree = count_kmer_blocks(read_blocks(arg_dics['block-target']), arg_dics['k'])
        report_count(tree, arg_dics['filename'])
    elif command == 'CMP':
        print('manhattan distance of sequences -> %f'%(distance_manhattan(
            read_sequence(arg_dics['category'], arg_dics['sequence-index']), 
            read_sequence(arg_dics['category'], arg_dics['other-index']), 
            arg_dics['k'])))
    elif command == 'CMP-B':
        block_distance_matrix(
            read_blocks(arg_dics['block-target']), 
            read_blocks(arg_dics['other-target']),
            arg_dics['k'], 
            arg_dics['filename']
        )
    elif command == 'CMP-Z':
        dataset_distance_matrix(arg_dics['category'], arg_dics['k'], arg_dics['filename'])
    else:
        print('ERROR - not supported command (check the description)')
