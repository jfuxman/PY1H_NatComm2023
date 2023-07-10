with open('Gene_Transcript_Seq.txt') as f1, open ('Transcript_Seq.fasta', 'w') as o1:
	f1.readline()
	for line in f1:
		line = line.rstrip().split('\t')
		if line[2] == 'Sequence unavailable':
			continue
		o1.write('>' + line[0] + '_' + line[1] + '\n' + line[2] + '\n')
