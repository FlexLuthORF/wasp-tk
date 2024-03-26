# read all the bed files in a given directory

import os
import csv

# This specifies the window outside the region delimited by the franken co-ordinates which will be
# provided in the 'genes' field and can therefore be searched if the allele is not aligned with the franken
GENE_WINDOW_SIZE_3 = 0     # 3' end
GENE_WINDOW_SIZE_5 = 0     # 5' end


def get_gene_type(label):
    gene = label.split('*')[0]

    if gene[3] == 'J' or gene[3] == 'V':
        return gene[3]
    elif gene[3] == 'D' and '-' in gene and '_' not in gene:
        return 'D'
    else:
        return 'C'


def read_beds(assembly_name, sense, dir, assemblies=None):
    """
     Read co-ordinates for an assembly from all BED files in a given directory

     :param sense: The sense of the assembl(ies) listed in the BED files, where + is 5' to 3', - is 3' to 5' and +- is both (eg IGK) or unknown/mixed.
     :type sense: str
     :param dir: The directory containing the BED files
     :type dir: str
     :param assemblies: Sequences of the assemblies listed in the BED files. If not provided, the output will not contain motif sequences
     :type assemblies: dict

     :return: Data structure containing, for each assembly, and for each gene in the assembly, the co-ordinates of the features denoted by the BED files.
     :rtype: dict of dicts
    """
    beds = {}
    found_D = False
    print('read beds sense is:' + sense)
    if sense not in ('+', '-', '+-'):
        raise ValueError(f"Sense must be one of +, -, +-")

    if assemblies is None:
        assemblies = {}

    for entry in os.scandir(dir):
        if entry.is_file() and '.bed' in entry.name:
            with open(os.path.join(dir, entry.name), 'r') as fi:
                el_type = entry.name.replace('.bed', '').upper()
                reader = csv.DictReader(fi, delimiter='\t', fieldnames=['assembly_name', 'start', 'end', 'gene', 'sense'])
                for row in reader:
                    if not row['gene'] or row['assembly_name'] != assembly_name:
                        continue

                    if row['start'] and row['end']:
                        refname = row['assembly_name']
                        row['el_type'] = el_type

                        row['start'] = int(row['start'])
                        row['end'] = int(row['end'])

                        row['seq'] = assemblies[refname][row['start']:row['end']] if refname in assemblies else ''

                        if refname not in beds:
                            beds[refname] = {}

                        if row['gene'] not in beds[refname]:
                            beds[refname][row['gene']] = {}

                        if row['el_type'] in beds[refname][row['gene']]:
                            if row['el_type'] in ('NONAMER', 'SPACER', 'HEPTAMER') and get_gene_type(row['gene']) == 'D':
                                found_D = True
                                if sense == '-':
                                    beds[refname][row['gene']]['3_' + row['el_type']] = beds[refname][row['gene']][row['el_type']]
                                    del beds[refname][row['gene']][row['el_type']]
                                    beds[refname][row['gene']]['5_' + row['el_type']] = row
                                else:
                                    beds[refname][row['gene']]['5_' + row['el_type']] = beds[refname][row['gene']][row['el_type']]
                                    del beds[refname][row['gene']][row['el_type']]
                                    beds[refname][row['gene']]['3_' + row['el_type']] = row
                            else:
                                print(f"Error: {refname} {row['gene']} {row['el_type']} is multiply defined")
                        elif el_type == 'CONSTANT':
                            beds[refname][row['gene']]['GENE'] = row

                        else:
                            beds[refname][row['gene']][row['el_type']] = row

    # In the above we assumed that any Ds were in the + orientation if sense was +-. If sense was specified as +-, check the sense
    # now, based on the first J we encountered in each ref (this assumes that Ds, where present, are always in consistent 
    # sense with the Js - if this isn't true the sense will have to be specified explicitly

    if sense == '+-' and found_D:
        for refname in beds:
            j_sense = None
            for gene in beds[refname]:
                gene_type = get_gene_type(gene)
                if gene_type == 'J':
                    row = beds[refname][gene]
                    if 'NONAMER' in row and 'HEPTAMER' in row:
                        j_sense = '+' if row['NONAMER']['start'] < row['HEPTAMER']['start'] else '-'
                        break

            if j_sense and j_sense == '-':
                for gene in beds[refname]:
                    gene_type = get_gene_type(gene)
                    if gene_type == 'D':
                        for el in ['HEPTAMER', 'SPACER', 'NONAMER']:
                            if '3_' + el in beds[refname][gene] and '5_' + el in beds[refname][gene]:
                                (beds[refname][gene]['3_' + el], beds[refname][gene]['5_' + el]) = (beds[refname][gene]['5_' + el], beds[refname][gene]['3_' + el])
                            else:
                                print(f"Warning: {refname} {gene} missing 3_ or 5_ {el}")


    # Sanity checks

    if sense == '-':
        g_s = 'start'
        g_e = 'end'
    else:
        g_s = 'end'
        g_e = 'start'

    for refname in beds:
        for gene in list(beds[refname]):
            if not gene:
                continue

            gene_type = get_gene_type(gene)
            row = beds[refname][gene]
            if gene_type == 'V':
                complete = True
                for el in ['EXON_1', 'INTRON', 'EXON_2', 'HEPTAMER', 'SPACER', 'NONAMER']:
                    if el not in row:
                        print(f'element {el} missing from gene {gene}')
                        complete = False

                if complete:
                    if sense == '+-':
                        if row['NONAMER']['start'] > row['HEPTAMER']['start']:
                            g_s = 'end'         # + sense gene
                            g_e = 'start'
                        else:
                            g_s = 'start'       # - sense gene
                            g_e = 'end'

                    if row['NONAMER'][g_e] != row['SPACER'][g_s]:
                        print(f"maths problem in {gene}: row['NONAMER'][{g_e}] != row['SPACER'][{g_s}]")
                    if row['SPACER'][g_e] != row['HEPTAMER'][g_s]:
                        print(f"maths problem in {gene}: row['SPACER'][{g_e}] != row['HEPTAMER'][{g_s}]")
                    if row['HEPTAMER'][g_e] != row['EXON_2'][g_s]:
                        print(f"maths problem in {gene}: row['HEPTAMER'][{g_e}] != row['EXON_2'][{g_s}]")
                    if row['EXON_2'][g_e] != row['INTRON'][g_s]:
                        print(f"maths problem in {gene}: row['EXON_2'][{g_e}] != row['INTRON'][{g_s}]")
                    if row['INTRON'][g_e] != row['EXON_1'][g_s]:
                        print(f"maths problem in {gene}: row['INTRON'][{g_e}] != row['EXON_1'][{g_s}]")

                    # add additional co-ords
                    row_sense = sense
                    if sense == '+-':
                        if row['NONAMER']['start'] > row['HEPTAMER']['start']:
                            row_sense = '+'
                        else:
                            row_sense = '-'

                    if row_sense == '-':
                        if 'V-REGION' not in row:
                            row['V-REGION'] = {}
                            row['V-REGION']['start'] = row['EXON_2']['start']
                            row['V-REGION']['end'] = row['EXON_2']['end'] - 11
                        if 'L-PART2' not in row:
                            row['L-PART2'] = {}
                            row['L-PART2']['start'] = row['EXON_2']['end'] - 11
                            row['L-PART2']['end'] = row['EXON_2']['end']
                    else:
                        if 'V-REGION' not in row:
                            row['V-REGION'] = {}
                            row['V-REGION']['start'] = row['EXON_2']['start'] + 11
                            row['V-REGION']['end'] = row['EXON_2']['end']
                        if 'L-PART2' not in row:
                            row['L-PART2'] = {}
                            row['L-PART2']['start'] = row['EXON_2']['start']
                            row['L-PART2']['end'] = row['EXON_2']['start'] + 11

                # provide a window around the gene in the allele sequence if one is specified
                if 'GENE' in row:
                    row['GENE']['start'] = row['GENE']['start'] - GENE_WINDOW_SIZE_5
                    row['GENE']['end'] = row['GENE']['end'] + GENE_WINDOW_SIZE_3

            elif gene_type == 'J':
                for el in ['HEPTAMER', 'SPACER', 'NONAMER']:
                    complete = True
                    if el not in row:
                        print(f'element {el} missing from gene {gene}')
                        complete = False

                if complete:
                    if sense == '+-':
                        if row['NONAMER']['start'] < row['HEPTAMER']['start']:
                            g_s = 'end'
                            g_e = 'start'
                        else:
                            g_s = 'start'
                            g_e = 'end'

                    if row['HEPTAMER'][g_e] != row['SPACER'][g_s]:
                        print(f"maths problem in {gene}: row['HEPTAMER'][{g_e}] != row['SPACER'][{g_s}]")
                    if row['SPACER'][g_e] != row['NONAMER'][g_s]:
                        print(f"maths problem in {gene}: row['SPACER'][{g_e}] != row['NONAMER'][{g_s}]")

                if 'GENE' in row:
                    row['GENE'][g_e] -= GENE_WINDOW_SIZE_5
                    row['GENE'][g_s] += GENE_WINDOW_SIZE_3

            elif gene_type == 'D':
                for el in ['3_HEPTAMER', '3_SPACER', '3_NONAMER', '5_HEPTAMER', '5_SPACER', '5_NONAMER']:
                    complete = True
                    if el not in row:
                        print(f'element {el} missing from gene {gene}')
                        complete = False

                if complete:
                    if row['3_NONAMER'][g_e] != row['3_SPACER'][g_s]:
                        print(f"maths problem in {gene}: row['3_NONAMER'][{g_e}] != row['3_SPACER'][{g_s}]")
                    if row['3_SPACER'][g_e] != row['3_HEPTAMER'][g_s]:
                        print(f"maths problem in {gene}: row['3_SPACER'][{g_e}] != row['3_HEPTAMER'][{g_s}]")
                    if row['5_HEPTAMER'][g_e] != row['5_SPACER'][g_s]:
                        print(f"maths problem in {gene}: row['5_HEPTAMER'][{g_e}] != row['5_SPACER'][{g_s}]")
                    if row['5_SPACER'][g_e] != row['5_NONAMER'][g_s]:
                        print(f"maths problem in {gene}: row['5_SPACER'][{g_e}] != row['5_NONAMER'][{g_s}]")

                if 'GENE' in row:
                    row['GENE'][g_s] -= GENE_WINDOW_SIZE_5
                    row['GENE'][g_e] += GENE_WINDOW_SIZE_3

    return beds


def write_beds(beds, dir):
    bed_files = {}
    for refname in beds:
        for gene in beds[refname]:
            for el_type in beds[refname][gene]:
                row = beds[refname][gene][el_type]
                if el_type not in bed_files:
                    bed_files[el_type] = {}
                if refname not in bed_files[el_type]:
                    bed_files[el_type][refname] = []
                if 'el_type' in row:
                    del row['el_type']
                if 'seq' in row:
                    del row['seq']
                bed_files[el_type][refname].append(row)
                
    for el_type in bed_files:
        with open(os.path.join(dir, f'{el_type}.bed'), 'w', newline='') as fo:
            writer = csv.DictWriter(fo, delimiter='\t', fieldnames=['assembly_name', 'start', 'end', 'gene', 'sense'])
            for refname in bed_files[el_type]:
                rows = bed_files[el_type][refname]
                rows.sort(key=lambda x: x['start'])
                writer.writerows(rows)
