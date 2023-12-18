#!/usr/bin/env python3
'''
This programm:
1) Compare probe id for 1 gene and find max_averege
2) Get a list with unike probe ID (only 1 probe with max_averege for 1 gene)
3)Find columns (In assignment for gene areas LHM and PHA, but programm can use different names)
4)Collect data for 1 column (or any columns for 1 name)
5)Compare data: a) Choose all above or equal 15   b) Choose all unique probe ID
6)Find shared probe between LHM and PHA (simular probe ID) 
7)Find unique probe ID both for LHM and PHA regions

input data:
#Please, use the correct directory paths in constants
write a names of Gene locations that we want compare (use only 2 structure acronims)
example: python3 assignment_3_vlad.py LHM PHA
 
output data:
information about this areas (shared and unique probes)
'''

__author__ = 'Vlad Pazenko'

# imports
import argparse  #used instead of sys

# constants
DIRECTORY1 = '../../New_Folder/normalized_microarray_donor9861/Probes.csv'
DIRECTORY2 = '../../New_Folder/normalized_microarray_donor9861/MicroarrayExpression.csv'
DIRECTORY3 = '../../New_Folder/normalized_microarray_donor9861/SampleAnnot.csv'

# functions


def open_file(directory):
    '''Function open directory and split each line to list 
    Input: directory 
    Output: data from file. Line by line, splited in "," 
    '''
    with open(directory, mode="r", encoding="utf-8") as data:
        return [line.strip('"').split(',') for line in data]


def create_dict_probe_id_averege(directory):
    '''Function that create and return a dictionary where key = probe_id, value = average
    Input: MicroarrayExpression.csv directory
    Output: dictionary[probe_id]=average
    '''
    dictionary_probe_aver = {}
    line_list = open_file(directory)
    for line in line_list:
        probe_id = line[0]
        line_no_id = list((map(float, line[1:])))
        dictionary_probe_aver[probe_id] = sum(line_no_id) / len(line_no_id)  # =average
    return dictionary_probe_aver


def create_dict_gene_probe_aver(dictionary_probe_aver, directory):
    '''Function that return a dictionary where key = gene name, value = list of lists contain 
    all pair [probe_id, average] connecting with this gene name
    Input: dictionary[probe_id] = average, Probes.csv directory
    Output: dictionary[gene_name] = [[probe1_id, average],[probe2_id, average]]
    '''
    dict_gene_and_list_id_aver = {}
    line_list = open_file(directory)
    for line in line_list[1:]:
        keys = dict_gene_and_list_id_aver.keys()
        if line[3] in keys: # if gene name exist in dictionary than add new list[ID, aver] in value
            dict_gene_and_list_id_aver[line[3]].append([line[0], dictionary_probe_aver[line[0]]])
        else:  # if gene name doesn't exist in dictionary than create new key-value pair, where key=gene name, value=list of lists
            dict_gene_and_list_id_aver[line[3]] = [[line[0], dictionary_probe_aver[line[0]]],]
    return dict_gene_and_list_id_aver


def compare_probe_for_gene(gene_dict):
    '''Function return a list of probes_id with max averege for each gene
    Input: dictionary[gene_name] = [[probe1_id, average],[probe2_id, average]]
    Output: list_probe_id
    '''
    list_probe_id = []
    aver = 0
    for i in gene_dict:
        max_aver = 0
        for j in range(0, len(gene_dict[i])):
            aver = gene_dict[i][j][1]
            if aver > max_aver:
                max_aver = aver
                probe_id = gene_dict[i][j][0]
        list_probe_id.append(probe_id)
    return list_probe_id


def find_columns(gene_zone_names, directory):
    '''Function return a dictionary where key = name of gene region (in our case 'LHM' and 'PHA') 
    and value = list of column numbers where ('LHM' and 'PHA') located
    Input: gene_zones_names #list of str, SampleAnnot.csv directory 
    Output: dictionary[gene_zone_name] = [column1, column2] #column in file MicroarrayExpression
    '''
    gene_zone_dict = {}
    gene_name1 = ''
    gene_name2 = ''
    for gene_zone in gene_zone_names:
        gene_zone_dict[gene_zone] = []
    line_list = open_file(directory)
    for count, line in enumerate(line_list[1:]):
        gene_abbrev = str(line[4]).replace('"', "")
        gene_name = str(line[5])
        for gene_zone in gene_zone_names:
            if gene_abbrev == gene_zone:
                # skip 1 line, and countdown start with 0
                gene_zone_dict[gene_zone].append(count+1)
        if gene_abbrev == gene_zone_names[0]:
            gene_name1 = gene_name.replace('"', "")
        elif gene_abbrev == gene_zone_names[1]:
            gene_name2 = gene_name.replace('"', "")
    return gene_zone_dict, gene_name1, gene_name2


def collect_probes_in_sets(list_probe_id, gene_zone_dict, directory):
    '''Function return two sets with unique probe_id for both names 
    in gene_zone_dict that greater or equal 15
    Input: List unique probes id, dictionary[gene_zone_name] = [column1, column2],
    MicroarrayExpression.csv directory
    Output: set1 for first arg (gene_zone), set2 for second arg (gene zone)
    '''
    set_arg1 = set()
    set_arg2 = set()
    list_gene_zones = list(gene_zone_dict.keys())
    line_list = open_file(directory)
    for line in line_list:
        if line[0] in list_probe_id:
            probe_id = line[0]
            line_no_id = list((map(float, line[1:])))
            gene_zone = gene_zone_dict[list_gene_zones[0]]
            for column in gene_zone:
                if line_no_id[column-1] >= 15:  # if expression in columns more then 15
                    set_arg1.add(probe_id)  # add probe ID to set
            gene_zone = gene_zone_dict[list_gene_zones[1]]
            for column in gene_zone:
                if line_no_id[column-1] >= 15:
                    set_arg2.add(probe_id)
    return set_arg1, set_arg2


def compare_two_sets(set_arg1, set_arg2, name_arg1, name_arg2):
    '''Function compare two sets and print result information
    Input: two sets of unique probes_id that greater 15 for gene zones that user input
    two names of these gene_zones 
    Output: none (during the function the data will be displayed on the screen)
    '''

    '''
    Input:  1. set_arg1: ....
            2.  set_arg2: ....
            3. name_arg1: ....
            4. name_arg2: ....
        
    Description: Function compare two sets and print result information

    Output: none (during the function the data will be displayed on the screen)

    
    
    
    '''
    shared_set = [element for element in set_arg1 if element not in (set_arg1 - set_arg2)]
    print(f'Shared probes:\n{shared_set}\n')
    print(len(shared_set))
    uniq_arg1 = set_arg1 - set_arg2  #Find list with unique probe_id for arg1
    print(f'Probes unique in {name_arg1}:\n{uniq_arg1}\n')
    print(len(uniq_arg1))
    uniq_arg2 = set_arg2 - set_arg1  #Find list with unique probe_id for arg2
    print(f'Probes unique in {name_arg2}:\n{uniq_arg2}')
    print(len(uniq_arg2))
    #return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process 2 gene zones abbreviation.',
                                    usage='e.g. "assignment_3_vlad.py arg1 arg2" use only 2 gene regions for these program')
    parser.add_argument('gene_region', type=str, nargs=2,
                        help='Please input only 2 gene zones to compare')
    gene_zones = parser.parse_args().gene_region
    if len(gene_zones[0]) > 11 or len(gene_zones[1]) > 11:
        print('Please use gene zones abbreviations')
    else:
        dict_probe_averege = create_dict_probe_id_averege(DIRECTORY2)
        dict_gene = create_dict_gene_probe_aver(dict_probe_averege, DIRECTORY1)
        list_probe = compare_probe_for_gene(dict_gene)
        gene_column_dict, name1, name2 = find_columns(gene_zones, DIRECTORY3)
        arg1, arg2 = collect_probes_in_sets(list_probe, gene_column_dict, DIRECTORY2)
        compare_two_sets(arg1, arg2, name1, name2)

