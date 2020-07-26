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


########################################
#            main call                 #
########################################

if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise Exception('request command must be specified (read the description for supported commands)')

    shortopt = 'b:o:f:s:c:'
    longopts = ['block-size=', 'overlap=', 'filename=', 'sequence-index=', 'category=']

    # default arguments
    arg_dics = {'block-size':100 , 'overlap':0, 'filename':'default.out', 'sequence-index':0, 'category':'FLU'}

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
            arg_dics.update({'filename':a})
    
    if command == 'BLK':
        generates_blocks(read_sequence(arg_dics['category'], arg_dics['sequence-index']), arg_dics['block-size'], arg_dics['overlap'], arg_dics['filename'])