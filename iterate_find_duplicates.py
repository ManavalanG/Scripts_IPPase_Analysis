# This script identifies %x identical sequences among FASTA sequences using BLAST and writes them into a CSV file.
# It first gathers all sequences by genera, runs BLAST for them and then identifies %x identical seqs based on certain criteria for seq length and %identity.
# In output CSV file, certain sequences are marked 'unique' based on criteria mentioned in this script.
# Obviously, %x identical seqs are not found by this script if they do not share 'Genus' name.
# This script needs BLAST installed locally.

# Author: MANAvalan Gajapathy

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import csv
import StringIO
import operator
import os

# read COMMENT ranking file and get it into dictionary
ranked_file = "refseq_ids_ranked.csv"
dict_rank_ids = {}
with open(ranked_file, "Ur") as rank_handle:
	for line in rank_handle:
		line = line.replace("\n", "")
		line_list = line.split(',')
		dict_rank_ids[ line_list[0] ] = int( line_list[2] )


input_file_fasta = "22Feb_21755seqs.fasta"

iterate_or_not = True
iteration_round = 1
	
while iterate_or_not:
	print "\nIteration round: %i\n" %iteration_round
	
	if iteration_round != 1:
		input_file_fasta = unique_seqs_filename
	dict_seq_100perc_matches = {}
	count = 0

	# gathers sequences into a dict by their genus name
	dict_grouping_by_genus = {}
	dict_seqlength = {}
	for seq_record in SeqIO.parse(input_file_fasta, "fasta"):
		genus = ((seq_record.id).split('|') )[4]
		genus = genus.lower()
		# one genus in test data had slash symbol and produced error when it was part of file path. So replace it.
		genus = genus.replace('/', '_')
		if genus in dict_grouping_by_genus:
			dict_grouping_by_genus[genus].append(seq_record)
		else:
			dict_grouping_by_genus[genus]= [seq_record]

		dscrptn = seq_record.description
		dscrptn = dscrptn.split('|')
		dict_seqlength[ dscrptn[0] ] = len(seq_record.seq)

	total_genus = len(dict_grouping_by_genus)
	print "\tTotal no. of genus in file:  ", total_genus


	# For each genera of seqs, BLAST is run first and then their results are read to identify sequences that match in length and content (eg. 100% identical considering seq length)
	progress = 1
	for genus in dict_grouping_by_genus:
		print "\t%i/%i Genus being processed: %s  " %(progress, total_genus, genus)
				
		try:
			os.makedirs('temp')
		except OSError as exception:
			pass	
				
		out_filename = 'temp/%s.fasta' % genus
		Output_handle = open(out_filename, 'w')
		SeqIO.write(dict_grouping_by_genus[genus], Output_handle, 'fasta')
		Output_handle.close()

		xml_filename = 'temp/' + genus + '_Blast_results.xml'
		blastp_cline = NcbiblastpCommandline(query= out_filename, subject= out_filename, outfmt= 5, out= xml_filename)
		stdout, stderr = blastp_cline()
		if stdout or stderr:
			print stdout, stderr

		result_handle = open(xml_filename)
		blast_records = NCBIXML.parse(result_handle)

		identical_match_IDs = []		# once a set of matching seqs is found, this list is used to avoid repeating the same set of matching in different combination
		for record in NCBIXML.parse(result_handle):
			# defines range of acceptable seq length on basis of query seq's length		
			temp =  (record.query_length)/100
			diff = (temp + 1) * 2
			low = record.query_length - diff
			high = record.query_length + diff
			
			if record.query not in dict_seq_100perc_matches:
				dict_seq_100perc_matches[record.query] = []
			
			for alignment in record.alignments:
				if alignment.hit_def not in dict_seq_100perc_matches and \
						alignment.hit_def not in identical_match_IDs and \
						record.query not in identical_match_IDs:		# to avoid duplicates in different combinations in results. However duplicates will be written in results in some cases

		# 				if record.query != alignment.hit_def and record.query_length == alignment.length:	# if IDs dont match and query's length is same as subject's length
					if record.query != alignment.hit_def and (low <= alignment.length <= high):	# if IDs dont match and query's length meets the threshold size
						for hsp in alignment.hsps:
		# 					if hsp.identities == alignment.length:	# if no. of identities is same as subject seq's length ie. 100% identity
							if hsp.identities >= (0.98 * alignment.length):	# if no. of identities match the threshold set
								identical_match_IDs.append(alignment.hit_def)
								title_list = alignment.hit_def.split(' ')		# to separate seq ID to write in output file
								dict_seq_100perc_matches[record.query].append( [title_list[0], alignment.length, hsp.identities] )\
		
		
			# deletes key if no duplicates were found for a query sequence
			if not len(dict_seq_100perc_matches[record.query]):
				dict_seq_100perc_matches.pop(record.query, None)
		
		result_handle.close()
		progress += 1

			
				
	 
	# this dictionary is used just to have the output file written in descending order
	dict_4_sorting_order = {}
	for key in dict_seq_100perc_matches:
		dict_4_sorting_order[key] = len(dict_seq_100perc_matches[key])

	out_csv_data = 'Following are shown for each sequence in this order: \nSeq length, #of identities, %identities.\n\n'
	duplicates_size = []
	all_ids_in_csv_data = []	
	all_unique_ids_in_csv_data = []
	for key, value in sorted(dict_4_sorting_order.iteritems(), key = lambda (k,v): (v,k), reverse = True):	# helps to sort
		key_list = key.split(' ')
		key_list[0] = key_list[0].replace(',', '_')
		if key_list[0][-1] == '|':	# for cases where Seq ID ends with pipe symbol
			key_list[0] = key_list[0][:-1]
		key_list = key_list[0].split('|')
		query_seq = ''
		for ele in key_list:
			query_seq += "%s," %ele
		query_seq += "%s\n" %dict_seqlength[key_list[0]] 
			
		
		dict_species = {}		# to gather IDs corresponding to each species
		species = "%s %s" %(key_list[4], key_list[5])
		species = str(species)
		dict_species[species] = [ str(key_list[0]) ]

		data_similar = ''
		for item in dict_seq_100perc_matches[key]:	
			count += 1
								
			item[0] = item[0].replace(',', '_')
			if item[0][-1] == '|':	# for cases where Seq ID ends with pipe symbol
				item[0] = item[0][:-1]
			seq_id_list = item[0].split('|')
			data_similar += ','.join(seq_id_list)
			data_similar += ",%s,%s,%0.1f\n" %(item[1], item[2], ( float(item[2])/float(item[1]) ) *100 )
			
			species = '%s %s' %(seq_id_list[4], seq_id_list[5])
			species = str(species)
			if species in dict_species:
				dict_species[species].append( str(seq_id_list[0]) )
			else:
				dict_species[species] = [ str(seq_id_list[0]) ]

		
		# sort results before writing into output file
		csv_data = csv.reader(StringIO.StringIO(data_similar))
		low = len(seq_id_list) - 1
		high = low + 3		
		sorted_data = sorted(csv_data, key=operator.itemgetter(*range(low, high)))	# sorts by species name and then by the numbers in following columns
		

		# identifies sequences that fit criteria in order to be saved/removed later on from results
		# Sequences are marked 'unique' using following criteria:
			#	1. Sequence record with highest ranking on basis of COMMENT status is chosen
			#	2. If multiple sequences are present for a species, seq with largest seq length is chosen.
			#	3. If multiple seqs pass above criteria, a seq is randomly chosen
		unique_ids_list = []
		for species in dict_species:	
			rank_dict = {k:[] for k in range(1,6)}
			seqlength_dict = {k:[] for k in range(1,6)}
			
			for identifier in dict_species[species]:
				if identifier not in all_ids_in_csv_data:
					all_ids_in_csv_data.append(identifier)
					
				rank_id = dict_rank_ids[identifier]
				rank_dict[rank_id].append( identifier )
				seqlength_dict[rank_id].append( dict_seqlength[identifier] )

			for no in range(1,6):
				if rank_dict[no]:
					max_seq_length = max( seqlength_dict[no] )
					best_id_index = seqlength_dict[no].index( max_seq_length)
					best_id = rank_dict[no][best_id_index]
					break
												
			# If only genus part is known but species part is not, such seqs are declared unique only when such genus doesn't have any species with known complete species name part
			if species.split(' ')[1] == 'NA':
				if len(dict_species) == 1:		# if only species in list has 'NA' as its species name
					unique_ids_list.append(best_id)
			else:
				unique_ids_list.append(best_id)
		all_unique_ids_in_csv_data += unique_ids_list
		
							
		# Sequences marked as unique are written into a output file
		ref_id = str(query_seq).split(',')
		if ref_id[0] in unique_ids_list:
			out_csv_data += 'Unique,'
		else:
			out_csv_data += ','
		out_csv_data += query_seq
		
		for line in sorted_data:
			if line[0] in unique_ids_list:
				out_csv_data += 'Unique,'
			else:
				out_csv_data += ','
			out_csv_data += ( ','.join(line) + '\n' )
		out_csv_data += '\n\n'
			
		duplicates_size.append(len(dict_seq_100perc_matches[key]))


	if len(duplicates_size):
		big_no = max(duplicates_size)
	else:
		big_no = 0
		
	summary = '\n\tTotal no. of duplicate seqs:  %i\n' %count + \
	 '\tTotal no. of Queries for which there is/are duplicate seqs:  %i\n' %len(dict_seq_100perc_matches) + \
	 '\tHighest no. of duplicates found for a query:  %i' %big_no

	out_csv_data += '\n%s' %summary

	duplicate_ids_all_data = set(all_ids_in_csv_data) - set(all_unique_ids_in_csv_data)

	dict_unique_seqs = {}
	extracted_seqs_by_id = []
	for seq_record in SeqIO.parse(input_file_fasta, "fasta"):
		identifier = (seq_record.id).split('|')
		if identifier[0] not in duplicate_ids_all_data:
			extracted_seqs_by_id.append(seq_record)
	
	folder_name = "iter_%i" %iteration_round
	try:
		os.makedirs(folder_name)
	except OSError as exception:
		pass	

	output_File = '%s/seqs_grouped_round_%i.csv' %(folder_name, iteration_round)
	output_handle_csv = open(output_File, 'w')
	output_handle_csv.write(out_csv_data)
	output_handle_csv.close()
	
	unique_seqs_filename = "%s/Unique_seqs_extracted_round_%i.fasta" %(folder_name, iteration_round)
	output_seqs_handle = open(unique_seqs_filename, "w")
	SeqIO.write(extracted_seqs_by_id, output_seqs_handle, 'fasta')
	output_seqs_handle.close()

	print '\n\tOutput is written into folder: %s' %folder_name
	print "\n\tNo. of seqs present in input fasta file:   %i" %len(dict_seqlength)
	print "\tNo. of non-duplicate sequences found and written into output fasta file:   %i" %len(extracted_seqs_by_id)
	print summary
	
	if len(extracted_seqs_by_id) != len(dict_seqlength):
		iteration_round += 1
	else:	
		iterate_or_not = False
	
	

print "\nSequences in Iteration# %i does not have anymore duplicate sequences.\n" %(iteration_round)
