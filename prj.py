from getopt import getopt
import sys


def generates_blocks(sequence, block_size, overlap, fasta_filename='blocks.fasta'):
    report = '> block size\n%d\n> overlap\n%d\n'%(block_size, overlap)
    block_index = 0
    block = ''
    for nucleotide in sequence:
        if len(block) == block_size:
            report += '> block_%d\n%s\n'%(block_index, block)
            block_index += 1
            block = block[-overlap:]
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
        

def test_one():
    s_file = open('dm01g.fasta', 'r')
    seq = s_file.readline()
    s_file.close()
    
    generates_blocks(seq, 200, 10)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise Exception('request command must be specified (read the description for supported commands)')

    shortopt = 'b:o:f:'
    longopts = ['block-size=', 'overlap=', 'filename=']

    arg_dics = {}

    command = sys.argv[1]

    opts, _ = getopt(sys.argv[2:], shortopt, longopts)
    for o, a in opts:
        if o in ['-b', '--block-size']:
            arg_dics.update({'block-size':int(a)})
        elif o in ['-o', '--overlap']:
            arg_dics.update({'overlap':int(a)})
        elif o in ['-f', '--filename']:
            arg_dics.update({'filename':a})
    
    if command == 'BLK':

        # test file TODO: remove this part after implementing standard reading sequences
        s_file = open('dm01g.fasta', 'r')
        seq = s_file.readline()
        s_file.close()

        generates_blocks(_, arg_dics['block-size'], arg_dics['overlap'], arg_dics['filename'])