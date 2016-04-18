# read two fasta files and write a output file without any duplicate seqs that may have been present in two files combined
# Author: 'Mana'valan Gajapathy

from Bio import SeqIO

file_1 = 'Ecoli_FamI_BlastFiltered.fasta'
file_2 = 'Yeast_FamI_BlastFiltered.fasta'

ids_file1 = []
all_sequence_records = []
unique_seqs = 0
file1_seqs = 0
for seq_record in SeqIO.parse(file_1, "fasta"):
	if seq_record.id not in ids_file1:
		ids_file1.append(seq_record.id)
		all_sequence_records.append(seq_record)
		unique_seqs += 1
	file1_seqs += 1

file2_seqs = 0
for seq_record in SeqIO.parse(file_2, "fasta"):
	if seq_record.id not in ids_file1:
		all_sequence_records.append(seq_record)
		unique_seqs += 1
	file2_seqs += 1

out_handle= open('Unique_seqs_combined.fasta', 'w')
SeqIO.write(all_sequence_records, out_handle, 'fasta')
out_handle.close()

print 'Number of seqs in 1st file: %i' % file1_seqs
print 'Number of seqs in 2nd file: %i' % file2_seqs
print 'Number of unique seqs found when two files combined: %i' % unique_seqs