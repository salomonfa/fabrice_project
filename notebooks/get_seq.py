import os, json
from Bio import SeqIO
import argparse

def get_seq(graph, paths, output, blockID):
    """
    

    Args:
        graph: path to a .json containing a partial pangraph.
        paths: path to a folder containing the .fasta files containing the full path sequences
        output: path to the location were the sequence of the the partial pangraph is stored
        blockID: string of the ID of the interesting block
    """

    with open(graph, 'r') as file:
        graph = json.load(file)
    
    strand_dict = dict()
    end_dict = dict()
    
    with open(output, 'w') as outfile:
        
        for path in os.listdir(paths):
        
            content = SeqIO.read(paths + '/' + path, 'fasta')

            for i, p in enumerate(graph['paths']):
                if p['name'] == content.id:
                    start = p['position'][0]
                    end = p['position'][-1]
                    for b in p['blocks']:
                        if p['name'] in strand_dict.keys():
                            strand_dict[content.id].append((b['id'], b['strand']))
                        else:
                            strand_dict[content.id] = list()
                            strand_dict[content.id].append((b['id'], b['strand']))

            end_dict[content.id] = (start, end)
        
        for p in strand_dict:
            ib = True
            for b in strand_dict[p]:
                if b[0] == blockID:
                    ib = b[1]
            start = end_dict[p][0]
            end = end_dict[p][1]
            if (strand_dict[p][0][1] == True and strand_dict[p][-1][1] == True):
                content = SeqIO.read(paths + '/' + p.split('.')[0] + '_all_na.fa' ,
                                      'fasta')[start:end]
            elif (strand_dict[p][0][1] == True and strand_dict[p][-1][1] == False):
                if ib == True:
                    content = SeqIO.read(paths + '/' + p.split('.')[0] + '_all_na.fa' ,
                                          'fasta')[start:end]
                else:
                    content = SeqIO.read(paths + '/' + p.split('.')[0] + '_all_na.fa' ,
                                          'fasta')[start:end]
                    content.seq = content.seq.reverse_complement()
            elif (strand_dict[p][0][1] == False and strand_dict[p][-1][1] == True):
                content = SeqIO.read(paths + '/' + p.split('.')[0] + '_all_na.fa' ,
                                      'fasta')[start:end]
                content.seq = content.seq.reverse_complement()
            else:
                content = SeqIO.read(paths + '/' + p.split('.')[0] + '_all_na.fa' ,
                                      'fasta')[start:end]
                content.seq = content.seq.reverse_complement()
            SeqIO.write(content, outfile, 'fasta')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="extract part of pangraph",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--graph", required=True, help="JSON containing pangraph")
    parser.add_argument("--paths", required=True, help="JSON containing path names and there species")
    parser.add_argument("--output", required=True, help="JSON containing pangraph")
    parser.add_argument("--blockID", required=False, help="block ID of interesting block")

    args = parser.parse_args()

    get_seq(args.graph, args.paths, args.output, args.blockID)


