import json
import pypangraph as pp
import argparse

def extract_json_graph(graph_path, strains_path):
    """
    Takes two json files extracts a pangraph object from the first and a dictionary containing 
    the beloning of strains to bacteria type.

    Args:
        graph_path: path to a json file containing a pangraph object.
        strains_path: path to a json file containing information about the bacterial strains.

    Returns:
    full_graph: a pangraph object.
    strains: a dictionary containing the bacteria as key and the paths as values.

    """
    
    full_graph = pp.Pangraph.load_json(graph_path)
    with open(strains_path, 'r') as file:
        strains = json.load(file)

    return full_graph, strains

def make_dictionary(full_graph, strains):
    """
    Takes a pangraph onject and a dictionary containing the belonging of paths. Creates a dictionary
    using Block IDs as keys and lists of the tuples as values. The tuples contain the path ID of the 
    path that the block appears in and an integer showing how many times said block appears in said 
    path. Also creates two sets, one containing the paths belonging to coli and one containing the 
    paths belonging to jejuni.

    Args:
        full_graph: a pangraph object.
        strains: a dictionary containing the bacteria as key and the paths as values.

    Returns:
        dictionary: dictionary using block ID as key and lists of tuples as values.
        set_coli: a set containing path IDs belonging to coli.
        set_jejuni: a set containing path IDs belonging to jejuni.
    """
    df_seq = full_graph.to_blockcount_df()
    blocks = df_seq.columns.values.tolist()
    paths = df_seq.index.tolist()
    
    dictionary = dict()
    set_coli = set() 
    set_jejuni = set()

    for path in paths:
        if path.split('.')[0] in strains['coli']:
            set_coli.add(path.split('.')[0])
        elif path.split('.')[0] in strains['jejuni']:
            set_jejuni.add(path.split('.')[0])
        else:
            set_coli.add(path.split('.')[0])


    for block in blocks:
        for path in paths:
            if df_seq[block][path] != 0.0:
                if block in dictionary:
                    dictionary[block].append((path.split(".")[0], df_seq[block][path]))
                else:
                    dictionary[block] = []
                    dictionary[block].append((path.split(".")[0], df_seq[block][path]))
    
    return dictionary, set_coli, set_jejuni


def extract_division_dicts(full_graph, dictionary, set_coli, set_jejuni):
    """
    Takes a pangraph object, dictionary and two sets. The dictionary and the sets can be created 
    using make_dictionary. Divides the the block IDs into categorys depending on there prevelance
    in coli and jejuni. 'core' blocks that are seen in all paths, 'core additional' essentialy core 
    but with a up to 2 paths missing, 'all coli' seen in all coli paths but no jejuni, 'all jejuni' 
    seen in all jejuni paths but no coli, 'part coli' seen in some coli paths but no coli, 
    'part jejuni' seen in some jejuni but no coli, 'accesory coli': seen in all coli paths and some 
    jejuni, 'accesory jejuni' seen in all jejuni paths and some coli, 'accesory unspecific' seen in
    some coli and some jejuni, 'singles total' all blocks that are only seen once, 'single coli' coli
    blocks that are only seen once, 'single jejuni' jejuni blocks that are only seen once.
    The name of the division is used as key for both dictionarys. The value associated with that key 
    is a list of tuples for the division_dict and a list of lists for the division_len_dict. Each 
    tuple in the division_dict contains the block ID and the length of a block. The list in the 
    division_len_dict contains the number of block in the division and the total length of all the
    blocks.

    Args:
        full_graph: a pangraph object.
        dictionary: dictionary using block ID as key and lists of tuples as values.
        set_coli: a set containing path IDs belonging to coli.
        set_jejuni: a set containing path IDs belonging to jejuni.

    Returns:
        division_dict: a dictionary using the division as key and lists of tupels as values
        division_len_dict: a dictionary using the division as key and lists as values
    """

    df_len = df_len = full_graph.to_blockstats_df()
    
    division_dict = {'core': [], 'core additional': [], 'all coli': [], 'all jejuni': [],
                    'part coli': [], 'part jejuni': [], 'accesory coli': [], 'accesory jejuni': [],
                    'accesory unspecific': [],'singles total': [], 'single coli': [], 
                    'single jejuni': []}

    division_len_dict = {'total with singlets': [0,0], 'total no singlets': [0,0], 'core': [0,0], 
                        'core additional': [0,0], 'all coli': [0,0], 'all jejuni': [0,0], 
                        'part coli': [0,0], 'part jejuni': [0,0], 'accesory coli': [0,0], 
                        'accesory jejuni': [0,0],'accesory unspecific': [0,0], 'singles total': [0,0],
                        'single coli': [0,0], 'single jejuni': [0,0]}

    single_list = list()

    for block in dictionary:
        tmp_set = set([path[0] for path in dictionary[block]])
        if len(tmp_set) > 1:
            if tmp_set.issubset(set_coli):
                if tmp_set == set_coli:
                    division_dict['all coli'].append((block, df_len['len'][block]))
                    division_len_dict['all coli'][0] += 1
                    division_len_dict['all coli'][1] += df_len['len'][block]
                else:
                    division_dict['part coli'].append((block, df_len['len'][block]))
                    division_len_dict['part coli'][0] += 1
                    division_len_dict['part coli'][1] += df_len['len'][block]
            elif tmp_set.issubset(set_jejuni):
                if tmp_set == set_jejuni:
                    division_dict['all jejuni'].append((block, df_len['len'][block]))
                    division_len_dict['all jejuni'][0] += 1
                    division_len_dict['all jejuni'][1] += df_len['len'][block]
                else:
                    division_dict['part jejuni'].append((block, df_len['len'][block]))
                    division_len_dict['part jejuni'][0] += 1
                    division_len_dict['part jejuni'][1] += df_len['len'][block]
            elif df_len['core'][block] and df_len['n. strains'][block] == (len(set_jejuni) + len(set_coli)):
                division_dict['core'].append((block, df_len['len'][block]))
                division_len_dict['core'][0] += 1
                division_len_dict['core'][1] += df_len['len'][block]
            else:
                if tmp_set.intersection(set_jejuni.union(set_coli)) >= set_jejuni.union(set_coli):
                    division_dict['core additional'].append((block, df_len['len'][block]))
                    division_len_dict['core additional'][0] += 1
                    division_len_dict['core additional'][1] += df_len['len'][block]
                elif set_coli.issubset(tmp_set):
                    division_dict['accesory coli'].append((block, df_len['len'][block]))
                    division_len_dict['accesory coli'][0] += 1
                    division_len_dict['accesory coli'][1] += df_len['len'][block]
                elif set_jejuni.issubset(tmp_set):
                    division_dict['accesory jejuni'].append((block, df_len['len'][block]))
                    division_len_dict['accesory jejuni'][0] += 1
                    division_len_dict['accesory jejuni'][1] += df_len['len'][block]
                else:
                    division_dict['accesory unspecific'].append((block, df_len['len'][block]))
                    division_len_dict['accesory unspecific'][0] += 1
                    division_len_dict['accesory unspecific'][1] += df_len['len'][block]
        else:
            single_list.extend(list(tmp_set))
            division_dict['singles total'].append((block, df_len['len'][block]))
            division_len_dict['singles total'][0] += 1
            division_len_dict['singles total'][1] += df_len['len'][block]
            if tmp_set.issubset(set_coli):
                division_dict['single coli'].append((block, df_len['len'][block]))
                division_len_dict['single coli'][0] += 1
                division_len_dict['single coli'][1] += df_len['len'][block]
            elif tmp_set.issubset(set_jejuni):
                division_dict['single jejuni'].append((block, df_len['len'][block]))
                division_len_dict['single jejuni'][0] += 1
                division_len_dict['single jejuni'][1] += df_len['len'][block]
                

    for key in division_len_dict:
        if not key.startswith('single'):
            division_len_dict['total no singlets'][0] += division_len_dict[key][0]
            division_len_dict['total no singlets'][1] += division_len_dict[key][1]
    division_len_dict['total with singlets'][0] = division_len_dict['total no singlets'][0] + \
        division_len_dict['singles total'][0]
    division_len_dict['total with singlets'][1] = division_len_dict['total no singlets'][1] + \
        division_len_dict['singles total'][1]

    return division_dict, division_len_dict




def extr_regular(path, block_id, index, content_set, division_dict, end_list, key, side):
    """
    extracts index of the next and previous core block starting from the block of interest.

    Args:
        path: dictionary containing the blocks of one path and there positions
        block_id: string of the ID of the block of intrest
        index: index of the block of intrest
        content_set: set containing all the blocks which will be includet in the partial pangraph
        division_dict: a dictionary using the division as key and lists of tupels as values
        end_list: dictionary using the species as key to a list contaiing the block at the start and end of the partial graph
        key: key for the end_list
        side: boolean clarifying if the end or start is to be found

    Returns:
        out_index: index of either the next or previous core block in the path
    """
    out_index = None
    
    if side == True:
        for j, b in enumerate(reversed(path['blocks'][:index])):
            content_set.add(b['id'])
            if b['id'] in [bl[0] for bl in division_dict[key]] and b['id'] != block_id:
                out_index = block_index(b['id'], index, j, end_list, key, True)
                break
    if side == False:    
        for j, b in enumerate(path['blocks'][index:]):
            content_set.add(b['id'])
            if b['id'] in [bl[0] for bl in division_dict[key]] and b['id'] != block_id:
                out_index = block_index(b['id'], index, j, end_list, key, False)
                break
     
    return out_index


                
def block_index(id, acc_index, core_index, end_list, key, side):
    """
    extracts index of the next and previous core block starting from the block of interest.

    Args:
        id: string of the ID of the block of intrest
        acc_index: index of the block of intrest
        core_index: index of either the next or previous core block starting from the block of interst
        end_list: dictionary using the species as key to a list contaiing the block at the start and end of the partial graph
        key: key for the end_list
        side: boolean clarifying if the end or start is to be found

    Returns:
        start: index of the previous core block in the path
        end: index of either the next core block in the path
    """
    if side == True:
        start = acc_index - core_index - 1
        end_list[key].add(id)
        return start
        
    elif side == False:
        end = acc_index + core_index + 1
        end_list[key].add(id)
        return end


def extr_start_end(path, block_id, division_dict, content_set, end_list, key, side):
    """
    extracts index of the next or previous core block starting from the block of interest. Is used if 
    the end has lower index then start or the sstart has a higher index then end. 


    Args:
        path: dictionary containing the blocks of one path and there positions
        block_id: string of the ID of the block of intrest
        division_dict: a dictionary using the division as key and lists of tupels as values
        content_set: set containing all the blocks which will be includet in the partial pangraph
        end_list: dictionary using the species as key to a list contaiing the block at the start and end of the partial graph
        key: key for the end_list
        side: boolean clarifying if the end or start is to be found

    Returns:
        out_index: index of either the next or previous core block in the path
    """
    out_index = None
    
    if side == True: 
        iter_list = reversed(path['blocks'])
    elif side == False:
        iter_list = path['blocks']

    for j, b in enumerate(iter_list):   
        content_set.add(b['id'])
        if b['id'] in [bl[0] for bl in division_dict[key]] and b['id'] != block_id:
            out_index = start_end_index(b['id'], j, path, key, side, end_list)
            break

    return out_index



def start_end_index(id, index, path, key, side, end_list):
    """
    extracts index of the next or previous core block starting from the block of interest. Is used if 
    the end has lower index then start or the sstart has a higher index then end. 

    Args:
        id: string of the ID of the block of intrest
        index: index of the block of intrest
        path: dictionary containing the blocks of one path and there positions
        key: key for the end_list
        side: boolean clarifying if the end or start is to be found
        end_list: dictionary using the species as key to a list contaiing the block at the start and end of the partial graph

    Returns:
        out_index: index of either the next or previous core block in the path
    """
    if side == True:
        out_index = len(path['blocks']) - index
        end_list[key].add(id)
        return out_index - 1

    elif side == False:
        out_index = index
        end_list[key].add(id)
        return out_index + 1

    


def extr_other(path, end_list, content_set):
    """
    extracts index of the two core blocks based from the end list. Is used if the path doesn't containe 
    the block of interest.

    Args:
        path: dictionary containing the blocks of one path and there positions
        end_list: dictionary using the species as key to a list contaiing the block at the start and end of the partial graph
        content_set: set containing all the blocks which will be includet in the partial pangraph

    Returns:
        This is a description of what is returned.

    Raises:
        start: index of the start block from the end_list
        end: index of the end block from the end_list
    """
    start = None
    end = None
    
    for i, b in enumerate(path['blocks']):

        if b['id'] in end_list and start == None:
            start = i
        elif b['id'] in end_list and start != None:
            end = i 
        
        if start != None and end != None:
            break

    if end > start and (len(path['blocks']) - end) + start < \
    (len(path['blocks']) - start) - (len(path['blocks']) - end):
        end, start = start, end

    if start > end and  (len(path['blocks']) - start) + end > \
    (len(path['blocks']) - end) -(len(path['blocks']) - start):
        end, start = start, end

    for i, b in enumerate(path['blocks']):
        if i >= start and i <= end:
            content_set.add(b['id'])

    return start, end + 1



                  


def delete_blocks(graph, content_set):
    """
    Deletes blocks from the partial pangraph which aren't in the content_set.

    Args:
        graph: full pangraph
        content_set: set containing all the blocks which will be includet in the partial pangraph
    """
    block_del_list = list()

    for block in graph['blocks']:
        if not block['id'] in content_set:

            block_del_list.append(block)
            continue

    graph['blocks'][:] = [x for x in graph['blocks'] if x not in block_del_list] 


def extr_partial_graph(graph_json, division_dict, block_id, set_coli, set_jejuni, output):
    """
    Extracts a partial pangraph from a full pangraph. the partial pangraph goes from a core block to another 
    core block and containes the block of interest if the full graph contained it.

    Args:
        graph_json: path to the full pangraph
        division_dict: a dictionary using the division as key and lists of tupels as values
        block_id: string of the ID of the block of intrest
        set_coli: a set containing path IDs belonging to coli.
        set_jejuni: a set containing path IDs belonging to jejuni.
        output: path to the location at which the partial pangraph is to be stored.
    """
    with open(graph_json, 'r') as file:
        graph = json.load(file)

    end_list = {'all coli': set(), 'all jejuni': set()}
    del_list = list()
    content_set = set()


    for n, path in enumerate(graph['paths']):
        start = None
        end = None

        for i, block in enumerate(path['blocks']):

            if block['id'] == block_id:
                if path['name'].split('.')[0] in set_coli:
                    start = extr_regular(path, block_id, i, content_set, division_dict, end_list,
                                                        'all coli', True)
                    end = extr_regular(path, block_id, i, content_set, division_dict, end_list,
                                                        'all coli', False)

                if path['name'].split('.')[0] in set_jejuni:
                    start = extr_regular(path, block_id, i, content_set, division_dict, end_list,
                                                        'all jejuni', True)
                    end = extr_regular(path, block_id, i, content_set, division_dict, end_list,
                                                        'all jejuni', False)
                    
                if start == None:
                    if path['name'].split('.')[0] in set_coli:
                        start = extr_start_end(path, block_id, division_dict, content_set, end_list,
                                                            'all coli', True)
    
                    if path['name'].split('.')[0] in set_jejuni:
                        start = extr_start_end(path, block_id, division_dict, content_set, end_list,
                                                            'all jejuni', True)
                    
                if end == None:
                    if path['name'].split('.')[0] in set_coli:
                        end = extr_start_end(path, block_id, division_dict, content_set, end_list,
                                                        'all coli', False)
                    if path['name'].split('.')[0] in set_jejuni:
                        end = extr_start_end(path, block_id, division_dict, content_set, end_list,
                                                        'all jejuni', False)
        
                break 

        
        if start == None or end == None:
            if path['name'].split('.')[0] in set_coli:
                if len(end_list['all coli']) != 0:
                    start, end = extr_other(path, end_list['all coli'], content_set)

            elif path['name'].split('.')[0] in set_jejuni:
                if len(end_list['all jejuni']) != 0:
                    start, end = extr_other(path, end_list['all jejuni'], content_set)

        if (path['name'].split('.')[0] in set_coli and len(end_list['all coli']) == 0) or \
            (path['name'].split('.')[0] in set_jejuni and len(end_list['all jejuni']) == 0): 
            graph.setdefault('paths').append(graph['paths'][n])
            del_list.append(path)
            continue
        
        if end > start and (len(path['blocks']) - end) + start < \
        (len(path['blocks']) - start) - (len(path['blocks']) - end):
            end, start = start, end

        if start > end and  (len(path['blocks']) - start) + end > \
        (len(path['blocks']) - end) -(len(path['blocks']) - start):
            end, start = start, end

        if start > end:
            graph['paths'][n]['circular'] = False
            graph['paths'][n].setdefault('blocks').extend(graph['paths'][n]['blocks'][:end])
            graph['paths'][n]['blocks'] = graph['paths'][n]['blocks'][start:]

            graph['paths'][n]['position'] = graph['paths'][n]['position'][:-1]
            graph['paths'][n].setdefault('position').extend(graph['paths'][n]['position'][:end + 1])
            graph['paths'][n]['position'] = graph['paths'][n]['position'][start:]

        else:
            graph['paths'][n]['circular'] = False
            graph['paths'][n]['blocks'] = graph['paths'][n]['blocks'][start:end]
            graph['paths'][n]['position'] = graph['paths'][n]['position'][start:end + 1]

    for path in del_list:
        graph.setdefault('paths').remove(path)

    delete_blocks(graph, content_set)

    with open(output, 'w') as file:
        file.write(json.dumps(graph))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="extract part of pangraph",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--pangraph", required=True, help="JSON containing pangraph")
    parser.add_argument("--strains", required=True, help="JSON containing path names and there species")
    parser.add_argument("--output", required=True, help="JSON of partial pangraph")
    parser.add_argument("--blockID", required=True, help="ID of block of interest")

    args = parser.parse_args()

    full_graph, strains = extract_json_graph(args.pangraph, args.strains)
    dictionary, set_coli, set_jejuni = make_dictionary(full_graph, strains)
    division_dict, division_len_dict = extract_division_dicts(full_graph, dictionary,
                                                            set_coli, set_jejuni)
    extr_partial_graph(args.pangraph, division_dict, args.blockID, set_coli, 
                   set_jejuni, args.output)

