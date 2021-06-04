from __future__ import print_function

import numpy as np
import gzip
import re
import time
import argparse

VERSION = 1.01

##Program to read in multiple admixfrog SNP files and then write them out as VCF

parser = argparse.ArgumentParser(description="This program will attempt to combine admixfrog snp output files into a VCF.")


parser.add_argument('--target_snpfiles', '-targs', dest='target_snpfiles', type=str, action = 'store', nargs='*', default = [],
                    help='space-separated list of locations of all snp files that will be included in the VCF. The full directory structure is required.')
parser.add_argument('--target_names', '-names', dest='target_names', type=str, action = 'store', nargs='*', default = [],
                    help='space-separated list of sample names for the VCF file. If these are not given then the file names from target_snpfiles are used.')
parser.add_argument('--admixfrog_reffile', '-ref', dest='admixfrog_reffile', type=str, action = 'store', nargs='?', default = '',
                    help='the location of the admixfrog referene file. This is used for the Ref and Alt values for each SNP. The full directory structure is required.')
parser.add_argument('--admixfrog_reffile_type', '-reftype', dest='admixfrog_reffile_type', type=str, action = 'store', nargs='?', default = 'reference',
                    help='the type of reference file. Usually best to leave as is, ie the refernece file you actually used when calling admixfrog. However, you can use a VCF.')

parser.add_argument('--save_snps_mode', '-mode', dest='save_snps_mode', type=str, action = 'store', nargs='?', default = 'only_called',
                    help='should the VCF: ["all_in_ref": save all SNPs in the reference file (including missing data as ./.)  location of the admixfrog referene file; "only_called" : Onlt the SNPs that are called in at least one individual].')
parser.add_argument('--outfile', '-out', dest='outfile', type=str, action = 'store', nargs='?', default = '/home/guy/outtmp.vcf.gz',
                    help='the location of the output VCF file.')

args = parser.parse_args()

#args.target_snpfiles = ['/media/guy/sf_Genetic_Data/aDNA_archaic_Helen/Data/mes1_snp_eg2.csv.gz', '/media/guy/sf_Genetic_Data/aDNA_archaic_Helen/Data/deni8_snp_eg2.csv.gz']

run_program = True


#Input error checking

if len(args.target_snpfiles) == 0:
    print("No target SNP files given, cannot run.")
    run_program = False
elif args.admixfrog_reffile == '':
    print("No reference file given, cannot run.")
    run_program = False
elif len(args.target_snpfiles) != len(args.target_names):
    print("No names given. Defaulting to file names.")
    args.target_names = [target_snpfile.split('/')[-1].split('.')[0] for target_snpfile in args.target_snpfiles]


#Functions to read in the SNP files and the relevant bits of the Reference file

def readin_admixfrogsnp_csvtsv(infile):
    #Read in a csv (comma separated) or tsv (tab separated) SNP file as output by admixfrog and return it as a numpy float array
    open_in = gzip.open if infile[-3:] == '.gz' else open
    read_data = []
    with open_in(infile, 'rb') as f_in:
        headers = f_in.readline()[0:-1].split(','.encode())
        for line in f_in:
            line = line.decode().replace('NA', 'nan')
            split_line = re.split('[\t,]', line[0:-1])
            chrom = int(split_line[1])
            pos = int(split_line[6])
            g0 = split_line[7]
            g1 = split_line[8]
            g2 = split_line[9]
            p = float(split_line[10])
            random_read = int(split_line[11])
            read_data.append([chrom, pos, g0, g1, g2, p, random_read])
    read_data = np.array(read_data, dtype = float)
    return read_data

def readin_reference_extract_chrposrefalt(infile, reference_file_type = 'reference'):
    #Read in reference file to create a 'chrom_pos' : [ref, alt] dictionary for each position called.
    open_in = gzip.open if infile[-3:] == '.gz' else open
    refalt_dict = dict()
    with open_in(infile, 'rb') as f_in:
        if reference_file_type == 'reference':
            #We expect one header line
            f_in.readline()
            for line in f_in:
                line_decoded = line.decode()
                split_line = line_decoded.split(',')
                refalt_dict['_'.join(split_line[0:2])] = [split_line[2],split_line[3]]
        elif reference_file_type == 'vcf':
            for line in f_in:
                line_decoded = line.decode()
                if line_decoded[0] == '#':
                    #headers
                    pass
                else:
                    #Try to read the line in
                    split_line = line_decoded.split('\t')
                    refalt_dict['_'.join(split_line[0:2])] = [split_line[3],split_line[4]]
    return refalt_dict


#Functions to work on the SNP data, by calling it based on the maximum likelihood genotype, and retrieving and merging the chr_pos lists

def call_csv_in(csv_array):
    #This calls the CSV, converting it to [chr, pos, geno] where geno is 0/1/2
    called_array = []
    for snp_idx in range(len(csv_array)):
        call = np.array([0,1,2])[np.argmax(csv_array[snp_idx][2:5])]
        called_array.append([csv_array[snp_idx][0], csv_array[snp_idx][1], call])
    called_array = np.array(called_array, dtype = int)
    return called_array

def convert_called_array_to_dict(called_array):
    #This converts a called array into a dict, which will be helpful when filling in the VCFs ie knowing if we have missing data
    keys = np.array(list(map('_'.join, np.array(called_array[::,0:2], dtype = str))))
    called_dict = dict(zip(keys, called_array[::,2]))
    return called_dict

def retreive_and_sort_chrpos(called_dict_list = []):
    #To write a VCF I need to know what positions were called, over all files.
    #Read in all the chr_pos dict keys and add to a set. This avoids duplicates.
    keys_set = set()
    for called_dict in called_dict_list:
        keys_set.update(called_dict.keys())
    keys_list = list(keys_set)
    #Split the key strings on _ and convert to int aray
    keys_list = np.array(list(map(lambda x : x.split('_'), keys_list)), dtype = int)
    #Sort on position then chromosome
    keys_list_sorted = keys_list[np.lexsort((keys_list[::,1], keys_list[::,0]))]
    keys_list_sorted = np.array(list(map('_'.join, np.array(keys_list_sorted, dtype = str))))
    return keys_list_sorted

"""
def convert_haplodip(called_array):
    #This converts 1s randomly to 0s or 2s in a called array
    hets = called_array[::,2] == 1
    num_hets = np.sum(hets)
    called_array[::,2][hets] = np.random.randint(0,2, num_hets) * 2
    return called_array
"""

#Function to write out the VCF

def write_vcf_from_chrposlist_dicts(outfile, chrpos_sorted, refalt_dict, called_dict_list = [], sample_name_list = [], save_snps_mode = 'only_called'):
    open_out = gzip.open if outfile[-3:] == '.gz' else open
    geno_map = ['0/0', '1/0', '1/1', './.']
    with open_out(outfile, 'w') as f_out:
        #Write the header
        f_out.write('##fileformat=VCFv4.1\n')
        currtime = time.localtime()
        f_out.write('##fileDate=%02d%02d%d_%02dh%02dm%02ds\n' %(currtime[2],currtime[1],currtime[0],currtime[3],currtime[4],currtime[5]))
        f_out.write('##source=admixfrog_snp_conversion.v%.1f\n' %(VERSION))
        f_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Unphased Genotype">\n')
        f_out.write(('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + sample_name_list) + '\n'))
        snps_to_write = chrpos_sorted if save_snps_mode == 'only_called' else refalt_dict.keys() if save_snps_mode == 'all_in_ref' else []
        for chrpos in snps_to_write:
            chr_pos_vals = chrpos.split('_')
            calls = [geno_map[called_dict.get(chrpos, 3)] for called_dict in called_dict_list]
            try:
                line_to_write = [chr_pos_vals[0], chr_pos_vals[1], '.', refalt_dict[chrpos][0], refalt_dict[chrpos][1], '.', 'PASS', '', 'GT']
            except KeyError:
                raise KeyError("Ref/Alt not known for SNP %s %s, make sure reference includes this SNP?" %(chr_pos_vals[0], chr_pos_vals[1]))
            line_to_write = line_to_write + calls
            f_out.write('\t'.join(line_to_write) + '\n')
    return outfile


if run_program == True:
    
    called_dict_list = []
    #Read in each SNP file in the list
    for target_infile in args.target_snpfiles:
        print("Reading in %s" %(target_infile))
        snp_data = readin_admixfrogsnp_csvtsv(target_infile)
        print("Got to readin+_admixfrogsnp_csvtsv")
        snp_data_called = call_csv_in(snp_data)
        print("Got to call_csv_in(snp_data)")
        snp_data_called_dict = convert_called_array_to_dict(snp_data_called)
        print("Got to convert_called_array")
        called_dict_list.append(snp_data_called_dict)
        print("Got to snp_data_called_dict")

    chrpos_sorted = retreive_and_sort_chrpos(called_dict_list)
    print("Got to retrieve_and_sort_chrpos")
    refalt_dict = readin_reference_extract_chrposrefalt(args.admixfrog_reffile, reference_file_type = args.admixfrog_reffile_type)
    print("Got to readin_reference_extract_chrpos")
    
    write_vcf_from_chrposlist_dicts(outfile = args.outfile, chrpos_sorted = chrpos_sorted, refalt_dict = refalt_dict, called_dict_list = called_dict_list, sample_name_list = args.target_names, save_snps_mode = args.save_snps_mode)
    print("Got to write_vcf_from_chrposlist_dicts")

else:
    raise RuntimeError()
