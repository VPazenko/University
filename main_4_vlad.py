#!/usr/bin/env python3
'''
Change the constants (path) if nedeed (in both files)
print in the terminal "main_4_vlad.py gene_zone1 gene_zone2 cut_off" 

MicroArray class (imported)
should store the Structure id, Structure acronym, Structure name
has a list of probe objects for a given experiment
should store if the probe is above background (see PACall.csv file)
above information is mandatory to pass when creating a microarray object
should keep track of the number of created microarrays
when printing a microarray object it should show the Structure id, Structure acronym \
    and Structure name nicelly formatted
has a method to get the probes that are above a certain expression value (cut_off) \
    (give a result as list of probe_id)

Probe class(imported)
should store expression value (in my case it stored only for two input regions)
should store probe_id, gene id and gene name, chromosome
These item use in MicroArray class because of minimizing runtime. 
(should store if the probe is above background (see PACall.csv file))
it should be mandatory to pass all information when the probe is created
when a probe is printed it should show a nicelly formatted string containing the value of the attributes
should keep track of the number of created probes
'''

__author__ = 'Vlad Pazenko'

# imports
import argparse # used instead of sys
import yaml #enable to process configuration.yaml file.

# import classes
# import all that we have in file classes_assignment4_vlad (2 classes 'Probe', "MicroArray" and 1 function 'open_file')
from classes_assignment4_vlad import *

"""
also we can use:
import classes_assignment4_vlad
but in the future we have to use command like 'classes_assignment4_vlad.Probe'
or line 'Probe = classes_assignment4_vlad.Probe'

also we can import only needed function/classes:
from classes_assignment4_vlad import Probe, MicroArray
"""

# constants
def open_config_yaml(directory):
    '''
    Input: 1. directory - config.yaml file

    Function return a dictionary from directory file

    Output: dictionary from (config file)
    '''
    with open(directory, 'r') as config_info:
        config_dict = yaml.safe_load(config_info)
    return config_dict

#fuctions

def find_columns(gene_zone_names, directory):
    '''
    Input:  1. gene_zones_names - list of str (user's input) 
            2. directory - SampleAnnot.csv file
    
    Function return a dictionary where key = name of gene region (in our case 'LHM' and 'PHA')
        and value = list of column numbers where ('LHM' and 'PHA') located
    
    Output: dictionary[gene_zone_name] = [column1, column2] #column in file MicroarrayExpression
    '''
    gene_zone_dict = {}
    for gene_zone in gene_zone_names:
        gene_zone_dict[gene_zone] = []
    line_list = open_file(directory)

    for count, line in enumerate(line_list[1:]):
        gene_abbrev = str(line[4]).replace('"', "")
        for gene_zone in gene_zone_names:
            if gene_abbrev == gene_zone:
                # skip 1 line, and countdown start with 0
                gene_zone_dict[gene_zone].append(count + 1)
    return gene_zone_dict


def find_probe_id_and_expressions(gene_zone_dict, directory):
    '''
    Input:  1. gene_zone_dict (key = name of gene region and value = list of column numbers where ('LHM' and 'PHA')
            located), 
            2. directory - MicroarrayExpression.csv file

    Function return a list of probes where each probe is a member of class Probe and contain 3 attributes
        probe_id, expression_values_region1, expression_values_region2 (the other 3 attributes are '').
        For region1 and region2 we realise regions from command line (in our case 'LHM' and 'PHA').
    
    Output: probes = list of probes where each element is a member of class Probe
    '''
    probes = []
    line_list = open_file(directory)
    for line in line_list:
        expression_values_reg1 = []
        expression_values_reg2 = []
        probe_id = line[0]
        line_no_id = list(map(float, line[1:]))
        average = sum(line_no_id) / len(line_no_id)
        list_gene_zones = list(gene_zone_dict.keys())
        for column in gene_zone_dict[list_gene_zones[0]]: #first elem in dict
            expression_values_reg1.append(line_no_id[column-1])
        for column in gene_zone_dict[list_gene_zones[1]]: #second elem in dict
            expression_values_reg2.append(line_no_id[column-1])     #we have to use index column-1 !!!
                                                                    #because we search expr in line_no_id list !
        #create an object probe in class Probe and give 3 attr (other attr are empty on this step)
        probe = Probe(probe_id, expression_values_reg1, expression_values_reg2, '', '', '', average)
        probes.append(probe)
    return probes


def find_probes_description(gene_zone_dict, directory1, directory2):
    '''
    Input:  1. gene_zone_dict (key = name of gene region and value = list of column numbers where ('LHM' and 'PHA')
            located),  #use a function "find_probe_id_and_expressions" inside this one.
            2. directory1 - MicroarrayExpression.csv file
            3. directory2 - probe.csv file

    Function return a list of probes where each probe is a member of class Probe and contain all attributes.
        First 3 attr we find through another function "find_probe_id_and_expressions"
    
    Output: probes = list of probes where each element is a member of class Probe
    '''
    gene_dict = {}
    probes = find_probe_id_and_expressions(gene_zone_dict, directory1)
    line_list = open_file(directory2)
    for count, line in enumerate(line_list[1:]):
        gene_id = line[2]
        chromosome = line[-1]
        gene_name = ','.join(line[4:-2])
        probe = probes[count]
        probe.gene_id = gene_id
        probe.gene_name = gene_name
        probe.chromosome = chromosome
        if gene_name in gene_dict:
            gene_dict[gene_name].append(probe)
        else:
            gene_dict[gene_name] = [probe,]
    return probes


def find_unique_probes(gene_dict):
    '''
    Input: probes = list of probes where each element is a member of class Probe

    Function return a list of probe (objects of class Probe) with max averege for each gene
    
    Output: unique_probes = list of probes (objects) with max averege for each gene
    '''
    unique_probes = []
    for value in gene_dict.values():
        value.sort(reverse=True, key=lambda elem: elem.average.max()) 
        unique_probes.append(value[0])
    return unique_probes


def find_micro_array_description(probes, gene_zone_dict, cut_off_value, directory):
    '''
    Input:  1. probes - list of objects class probe (we need attr: probe_id, expression_values_region1, expression_values_region2)
            2. gene_zone_dict (key = name of gene region and value = list of column numbers where ('LHM' and 'PHA') located
            3. cut_off_value - certain express value (users input), default = 1 
            4. directory - SampleAnnot.csv directory

    Function create 2 objects in class MicroArray (gene zones, that was printed in command line)
        (self, structure_id, structure_acronim, structure_name, probes_info, gene_zone_dict):
    
    Output: 2 objects (class MicroArray)    
    '''
    gene_zone1 = None
    gene_zone2 = None
    one_time1 = True
    probes_info = [[],[]]
    probes_info2 = []    
    gene_zone_list = []
    line_list = open_file(directory)
# =================================================
    for probe in probes:
        probes_info[0].append([probe.expression_values_region1, probe.probe_id])
        probes_info[1].append([probe.expression_values_region2, probe.probe_id])
    for count, key, value in enumerate(gene_zone_dict):
        chosen_line = line_list[value[0]]
        gene_abbrev = key
        structure_id = chosen_line[0]
        structure_name = ','.join(chosen_line[5:-7])
        gene_zone = MicroArray(structure_id, gene_abbrev, structure_name, probes_info[count], 
                                gene_zone_dict, cut_off_value)
        gene_zone_list.append(gene_zone)
# =================================================
    return gene_zone_list


def compare_two_gene_zones(zone_object1, zone_object2):
    '''
    Input:  1. zone_object1: object class MicroArray with all attributes for input gene zone 1
            2. zone_object2: object class MicroArray with all attributes for input gene zone 2
        
    Description: Function compare two sets(probe_id) from 2 objects and print result information

    Output: none (during the function the data will be displayed on the screen)
    '''
    set_arg1 = set(zone_object1.list_probe_id)
    set_arg2 = set(zone_object2.list_probe_id)
    shared_set = [element for element in set_arg1 if element not in (set_arg1 - set_arg2)]
    print(f'Shared probes:\n{shared_set}\n')
    #print(len(shared_set))
    uniq_arg1 = set_arg1 - set_arg2  #Find list with unique probe_id for arg1
    print(f'Probes unique in {zone_object1.structure_name}:\n{uniq_arg1}\n')
    #print(len(uniq_arg1))
    uniq_arg2 = set_arg2 - set_arg1  #Find list with unique probe_id for arg2
    print(f'Probes unique in {zone_object2.structure_name}:\n{uniq_arg2}')
    #print(len(uniq_arg2))
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process 2 gene zones abbreviation and cut_off value for compare.',
                                    usage='e.g. "main_4_vlad.py arg1 arg2 cut_off" use only 2 gene regions and 1 value for these program')
    parser.add_argument('gene_region', type=str, nargs=2,
                        help='Please input only 2 gene zones to compare')
    parser.add_argument('cut_off', type=int, nargs='?', default=1,
                        help='Please input only 1 cut_off value')    
    gene_zones = parser.parse_args().gene_region
    cut_off_value = parser.parse_args().cut_off
    if len(gene_zones[0]) > 11 or len(gene_zones[1]) > 11:
        print('Please use gene zones abbreviations')
    else:
        gene_zone_dict = find_columns(gene_zones, DIRECTORY1)
        probes = find_probes_description(gene_zone_dict, DIRECTORY2, DIRECTORY3)
        unique_probes = find_unique_probes(probes)
        zone1, zone2 = find_micro_array_description(unique_probes, gene_zone_dict, cut_off_value, DIRECTORY1)
        compare_two_gene_zones(zone1, zone2)



 