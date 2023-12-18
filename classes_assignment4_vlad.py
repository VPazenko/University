#!/usr/bin/env python3
'''
Change the constants (path) if nedeed

Creation of 2 classes:

MicroArray class
Probe class
'''

__author__ = 'Vlad Pazenko'

# imports
import yaml #enable to process configuration.yaml file.

# constants
#DIRECTORY4 = '../../New_Folder/normalized_microarray_donor9861/PACall.csv'

class Probe():
    '''
    Probe class
    should store expression value (in my case it stored only for two input regions)
    should store probe_id, gene id and gene name, chromosome
    These item use in MicroArray class because of minimizing runtime. 
    (should store if the probe is above background (see PACall.csv file))
    it should be mandatory to pass all information when the probe is created
    when a probe is printed it should show a nicelly formatted string containing the \
    value of the attributes
    should keep track of the number of created probes
    '''
    counter = 0

    def __init__(self, probe_id, expression_values_region1, expression_values_region2,
                 gene_id, gene_name, chromosome, average):
        self.probe_id = probe_id
        self.expression_values_region1 = expression_values_region1
        self.expression_values_region2 = expression_values_region2
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chromosome = chromosome
        self.average = average
        Probe.counter += 1

    def __str__(self):
        return f'''probe ID: {self.probe_id}, contain information about expression \
                in the first region {self.expression_values_region1} and second region \
                {self.expression_values_region2}. Also, this sample belongs to the gene \
                {self.gene_name} with gene id {self.gene_id} in chromosome {self.chromosome}'''


class MicroArray:
    '''
    MicroArray class
    should store the Structure id, Structure acronym, Structure name
    has a list of probe objects for a given experiment
    should store if the probe is above background (see PACall.csv file)
    above information is mandatory to pass when creating a microarray object
    should keep track of the number of created microarrays
    when printing a microarray object it should show the Structure id, Structure acronym and \
    Structure name nicelly formatted
    has a method to get the probes that are above a certain expression value (cut_off) \
    (give a result as list of probe_id)
    '''
    counter = 0

    def __init__(self, structure_id, structure_acronim, structure_name, probes_info, \
                  gene_zone_dict, cut_off):
        self.structure_id = structure_id
        self.structure_acronim = structure_acronim
        self.structure_name = structure_name
        self.cut_off = cut_off
        # comment line below if you want shut down compare with background
        try:
            self.compare_probe_with_background(probes_info, gene_zone_dict, structure_acronim)
        except FileNotFoundError:
            print('''File PACall.csv not found, please, check the directory.
                Programm will be process without compare_probe_with_background''')
        self.list_of_probes = self.find_list_of_probes(probes_info)
        self.list_probe_id = self.find_list_probe_id_cutoff(self.list_of_probes,
                                                            self.cut_off)

    def compare_probe_with_background(self, probes_info, gene_zone_dict, structure_acronim, \
                                      directory=DIRECTORY4):
        '''
        Input:  1. probes_info - list of expressions in genes area + probe_id (looks like \
                [[[expr1, expr2], probe1_id], [...]])
                2. gene_zone_dict (key=name of gene region and value = list of column numbers where
                (in our case 'LHM' and 'PHA', but it may be any zones) located)
                3. structure_acronim - short name of gene region (same as keys in gene_zone_dict)
                4. directory for file which contain info about fone signal (PACall.csv)

        Function change expression value to 0 if signal equal fone noise

        Output: none (we change expression values in probes_info list)
        '''
        line_list = open_file(directory)
        # each element looks like [[expr1, expr2],[expr1, expr2, expr3]] expr for 2 regions
        for count1, element in enumerate(probes_info):
            line = line_list[count1]
            # first elem in dict
            for count2, column in enumerate(gene_zone_dict[structure_acronim]):
                if int(line[column]) == 0:
                    # expr = 0, because signal equal fone
                    element[0][count2] = 0

    def find_list_of_probes(self, probes_info):
        '''
        Input:  1. probes_info - list of expressions in gene area + probe_id (looks like \
        [[[expr1, expr2], probe1_id], [...]])

        Function give a list of probes where at least 1 expression value more than 0.

        Output: list_of_probes - list of expr + probe_id for this region (only if this probe \
        have significant (>fone) expression value )
        '''
        list_of_probes = []
        for element in probes_info:
            if not all([e == 0 for e in element[0]]):
                list_of_probes.append(element)
        return list_of_probes

    def find_list_probe_id_cutoff(self, list_of_probes, cut_off):
        '''
        Input: 1. list_of_probes - list of expr + probe_id for this region \
        (only if this probe have significant (>fone) expression value )

        Function modifies list_of_probes and return list of probe id, where at least \
        one expression more than cut off value

        Output: list_probe_id (list of probe id, where at least one expression \
        (for this gene zone) more than cut off value)
        '''
        list_probe_id = []
        for element in list_of_probes:
            # (element[0] >= cut_off):
            if any([e >= cut_off for e in element[0]]):
                list_probe_id.append(element[1])
        return list_probe_id

    def __str__(self):
        return f'{self.structure_acronim} collect probes with expr more than \
                    {self.cut_off}: {self.list_probe_id}'


# fuctions

def open_file(directory):
    '''
    Input: 1. directory - file location ''

    Function open directory and split each line to list

    Output: data from file. Line by line, splited in ","
    '''
    with open(directory, mode="r", encoding="utf-8") as data:
        if 'probe' in directory.lower():
            pass
        return [line.strip('"\n').split(',') for line in data]


def open_config_yaml(directory):
    '''
    Input: 1. directory - config.yaml file

    Function return a dictionary from directory file

    Output: dictionary from (config file)
    '''
    with open(directory, 'r') as config_info:
        config_dict = yaml.safe_load(config_info)
    return config_dict



# main


if __name__ == '__main__':
    print(0)
