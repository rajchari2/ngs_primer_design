# script to design primers using primer3

# input: list of targets with chromsomal position and strand from sgRNA scorer
# input: genome file
# input: primer3 location
# input: last PSID used
# output: tab file with everything for analysis
# output: list of oligos to order
# output: primer info file upload for db


import sys
import argparse
import re
import glob
import subprocess
import os
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq

def design_ngs_primers(target_file, genome_file, design_set, primer3_location, primer3_settings_file_location, psid, ivt_num):
	
	# make directories if needed
	if os.path.exists('./primer3_input/')==False:
		mkdir_command1 = 'mkdir primer3_input'
		p = subprocess.Popen(mkdir_command1,shell=True)
		p.communicate()
	if os.path.exists('./primer3_output/')==False:
		mkdir_command2 = 'mkdir primer3_output'
		p = subprocess.Popen(mkdir_command2,shell=True)
		p.communicate()

	# primer3 input and output file
	primer3_input_file = './primer3_input/' + design_set + '_primer3_input.txt'
	primer3_input = open(primer3_input_file,'w')
	ps_to_region = defaultdict(str)

	# variable to hold absolute region coordinates
	region_db = defaultdict(str)

	# identify which chromosomes to store
	chr_list = []
	chr_db = defaultdict(str)

	# increment psid
	psid = psid.replace('PS','')
	num_psid = int(psid)
	num_psid += 1

	# increment ivt number
	ivt_number = int(ivt_num) + 1

	# find longest set of guides by position
	guide_dict = defaultdict(dict)
	guide_to_psid = defaultdict(str)
	name_to_target = defaultdict(str)
	guide_list = []
	ivt_to_name = defaultdict(str)

	print('Going through target file')

	# go through targets
	for line in target_file:
		line = line.rstrip('\r\n')
		name,target,scores,hits,position = line.split('\t')
		chrom,start_end,strand = position.split(':')
		start,end = start_end.split('-')
		if chrom not in chr_list:
			chr_list.append(chrom)
		guide_dict[chrom][int(end)] = line
		guide_list.append(line)
		name_to_target[name] = target
	target_file.close()

	print('Storing relevant genome files')

	# store relevant chrs
	chr_list_copy = chr_list
	for record in SeqIO.parse(genome_file,'fasta'):
		if str(record.id) in chr_list:
			chr_db[str(record.id)] = record.seq
			chr_list_copy.remove(str(record.id))
		if len(chr_list_copy)==0:
			break

	print('Building clusters')

	# identify longest window of guides within 130 bp
	for chrom in guide_dict:
		# guide clusters
		guide_cluster = defaultdict(list)
		positions = sorted(guide_dict[chrom].keys())
		index1 = 0
		index2 = 0
		chain = defaultdict(list)
		while index1 < len(positions):
			index2 = index1 + 1
			current_node = positions[index1]
			chain[current_node].append(current_node)
			while index2 < len(positions):
				if positions[index2] - positions[index1] <= 150:
					chain[current_node].append(positions[index2])
				index2 += 1
			index1 += 1
		
		# now go through and make clusters
		while len(positions) > 0:
			largest_node_size = 0
			largest_node = ''
			for node in chain:
				if len(chain[node]) > largest_node_size:
					largest_node = node
					largest_node_size = len(chain[node])

			# add largest node to cluster
			guide_cluster[largest_node] = []
			for neighbor in chain[largest_node]:
				if neighbor in positions and neighbor not in guide_cluster[largest_node]:
					guide_cluster[largest_node].append(neighbor)
					positions.remove(neighbor)

			# if it's an empty cluster, delete
			if len(guide_cluster[largest_node])==0:
				del guide_cluster[largest_node]

			# remove key from chain
			del chain[largest_node]

		# go through and design primers
		for gc in sorted(guide_cluster.keys()):
			# find start and end
			first_pos = 0
			last_pos = 0
			if len(guide_cluster[gc]) > 1:
				# grab the info
				first_pos = guide_cluster[gc][0] - 23
				last_pos = guide_cluster[gc][-1]
			elif len(guide_cluster[gc])==1:
				first_pos = guide_cluster[gc][0] - 23
				last_pos = guide_cluster[gc][0]

			# grab the minimal spanning region and the 
			position = chrom + ':' + str(first_pos) + '-' + str(last_pos)
			region = chrom + ':' + str(first_pos-200) + '-' + str(last_pos+200)
			relevant_dna = str(chr_db[chrom][first_pos-200-1:last_pos+200]).upper()
			size = last_pos - first_pos + 30

			# write to primer3 input
			psid_name = 'PS' + str(num_psid)
			primer3_input.write('SEQUENCE_ID=' + psid_name + '\n')
			primer3_input.write('SEQUENCE_TEMPLATE=' + relevant_dna + '\n')
			primer3_input.write('SEQUENCE_TARGET=170,' + str(size) + '\n')
			primer3_input.write('=\n')

			# increment psid #
			num_psid += 1

			# store region
			ps_to_region[psid_name] = relevant_dna
			region_db[psid_name] = region
			
			# go through and populate PSID
			for guide in guide_cluster[gc]:
				if guide in guide_dict[chrom]:
					entry = guide_dict[chrom][guide]
					parts = entry.split('\t')
					guide_to_psid[parts[0]] = psid_name
					ivt_to_name[ivt_number] = parts[0]
					ivt_number += 1


	# close primer3_input
	primer3_input.close()

	# write to output file
	primer3_output_file =  './primer3_output/' + design_set + '_primer3_output.txt'

	# run primer3
	print('Running primer3')
	primer3_command = primer3_location + ' --p3_settings_file=' + primer3_settings_file_location + ' --output=' + primer3_output_file + ' ' + primer3_input_file
	p = subprocess.Popen(primer3_command,shell=True)
	p.communicate()

	print('Parsing primer3 output')

	# open primer3 output and parse and store information
	primer3_output = open(primer3_output_file,'r')
	primer_db = defaultdict(dict)
	primer_left = ''
	primer_right = ''
	for line in primer3_output:
		line = line.rstrip('\r\n')
		if 'SEQUENCE_ID' in line:
			seq_name = line.replace('SEQUENCE_ID=','')
		elif 'PRIMER_LEFT_0_SEQUENCE' in line:
			primer_left = line.replace('PRIMER_LEFT_0_SEQUENCE=','')
		elif 'PRIMER_RIGHT_0_SEQUENCE' in line:
			primer_right = line.replace('PRIMER_RIGHT_0_SEQUENCE=','')
		elif line=='=':
			# record ended
			seq_name = seq_name.replace('PS','')
			ps = int(seq_name)
			primer_db[ps] = defaultdict(str)
			if primer_left != '':
				primer_db[ps]['F'] = primer_left
				primer_db[ps]['R'] = primer_right
			else:
				primer_db[ps]['F'] = 'TO-DESIGN'
				primer_db[ps]['R'] = 'TO-DESIGN'
			# reset variables
			seq_name = ''
			primer_left = ''
			primer_right = ''

	# close output
	primer3_output.close()

	# tab file and oligo file
	oligo_file = design_set + '_oligo_order.tab'
	analysis_file = design_set + '_analysis_file.tab'
	oligo_output = open(oligo_file,'w')
	analysis_output = open(analysis_file,'w')

	# adapter sequences
	left_ivt = 'GGATCCTAATACGACTCACTATAG '
	right_ivt = ' GTTTTAGAGCTAGAA'
	left_ngs = 'TCCCTACACGACGCTCTTCCGATCT '
	right_ngs = 'GTTCAGACGTGTGCTCTTCCGATCT '

	# cell type
	if 'hg38' in genome_file.name:
		cell_type = '293T'
	else:
		cell_type = 'P19'

	# write the primers
	amplicon_db = defaultdict(str)
	for ps in sorted(primer_db.keys()):
		# get the name
		ps_name = 'PS' + str(ps)

		# get the region
		region = ps_to_region[ps_name]
		region_coords = region_db[ps_name]

		# forward and reverse primers
		forward_primer = primer_db[ps]['F']
		reverse_primer = primer_db[ps]['R']
		# find positions
		if forward_primer != 'TO-DESIGN' and reverse_primer != 'TO-DESIGN':
			# get the indices
			start_index = region.find(forward_primer)
			end_index = region.find(str(Seq(reverse_primer).reverse_complement())) + len(reverse_primer)

			# amend coordinates
			chrom,start_end = region_coords.split(':')
			start,end = start_end.split('-')
			amp_start = int(start) + start_index
			amp_end = int(start) + end_index - 1

			# store in amplicon db
			amplicon_db[ps_name] = chrom + ':' + str(amp_start) + '-' + str(amp_end)

			# write to idt file
			oligo_output.write(ps_name + '-NGS-F\t' + left_ngs + forward_primer + '\n')
			oligo_output.write(ps_name + '-NGS-R\t' + right_ngs + reverse_primer + '\n')
		else:
			oligo_output.write(ps_name + '-NGS-F\t' + left_ngs + 'TO-DESIGN' + '\n')
			oligo_output.write(ps_name + '-NGS-R\t' + right_ngs + 'TO-DESIGN' + '\n')
			amplicon_db[ps_name] = 'TBD'
	# make output files, first write oligo 
	for ivtnum in sorted(ivt_to_name.keys()):
		# name
		guide_name = ivt_to_name[ivtnum]
		psname = guide_to_psid[guide_name]

		# get the target
		target = name_to_target[guide_name]
		spacer = target[0:20]

		# write to oligo output file
		oligo_output.write(guide_name + '_IVT-' + str(ivtnum) + '\t' + left_ivt + spacer + right_ivt + '\n')

		# write to analysis file
		analysis_output.write('IVT-' + str(ivtnum) + '-R1\t' + genome_file.name + '\t' + amplicon_db[psname] + '\t' + target + '\t' + cell_type + '-' + psname + '-Control\n')
		analysis_output.write('IVT-' + str(ivtnum) + '-R2\t' + genome_file.name + '\t' + amplicon_db[psname] + '\t' + target + '\t' + cell_type + '-' + psname + '-Control\n')

	# close file handles
	oligo_output.close()
	analysis_output.close()


def make_long_psid(primer_pair_count):
	primer_pair_id = ''
	if primer_pair_count < 10:
		primer_pair_id = '0000000'
	elif primer_pair_count < 100:
		primer_pair_id = '000000'
	elif primer_pair_count < 1000:
		primer_pair_id = '00000'
	elif primer_pair_count < 10000:
		primer_pair_id = '0000'
	elif primer_pair_count < 100000:
		primer_pair_id = '000'
	elif primer_pair_count < 1000000:
		primer_pair_id = '00'
	elif primer_pair_count < 10000000:
		primer_pair_id = '0'
	primer_pair_id = 'PS' + primer_pair_id + str(primer_pair_count)
	return primer_pair_id

def main(argv):
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-t','--target_file',type=argparse.FileType('r'),required=True)
	parser.add_argument('-g','--genome_file',type=argparse.FileType('r'),required=True)	
	parser.add_argument('-d','--design_set',required=True)
	parser.add_argument('-e','--primer3_location',required=True)
	parser.add_argument('-s','--primer3_settings_file',required=True)
	parser.add_argument('-p','--psid',required=True)
	parser.add_argument('-i','--ivt_num',required=True)
	opts = parser.parse_args(argv)
	design_ngs_primers(opts.target_file, opts.genome_file, opts.design_set, opts.primer3_location, opts.primer3_settings_file, opts.psid, opts.ivt_num)

if __name__ == '__main__':
	main(sys.argv[1:])


