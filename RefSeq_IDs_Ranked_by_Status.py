# Reads GenBank/GenPept file and ranks accession IDs in them on the basis of status in 'COMMENT' section
# Author: 'Mana'valan Gajapathy

from Bio import SeqIO   # BioPython module

status_codes = ["REVIEWED REFSEQ", "VALIDATED REFSEQ", "PROVISIONAL REFSEQ", "PREDICTED REFSEQ",
					"MODEL REFSEQ", "INFERRED REFSEQ", "REFSEQ"]
status_ranks= [1, 2, 3, 4, 4, 4, 5]
dict_status = dict( zip(status_codes, status_ranks) )

filename = "22Feb_RefSeq_22173seqs.gp"      # Input GenPept file
out_handle = open("refseq_ids_ranked.csv", "w")
for record in SeqIO.parse(filename, "genbank"):
	comment = record.annotations["comment"]
	comment = comment.split(":")
	rank_status = dict_status[comment[0]]
	
	out_handle.write("%s,%s,%i\n" %(record.id, comment[0], rank_status) )

out_handle.close()