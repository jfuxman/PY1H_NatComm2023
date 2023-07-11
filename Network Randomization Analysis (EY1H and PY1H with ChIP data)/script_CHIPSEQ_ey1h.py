import pandas as pd
import subprocess
import os
import bbi

# wget http://gtrd.biouml.org:8888/downloads/current/metadata/ChIP-seq.metadata.txt
metadata = pd.read_table('ChIP-seq.metadata.txt')
metadata = metadata[metadata.specie == 'Homo sapiens']

### Events
df = pd.read_excel('py1h ey1h chip-seq overlap list for Luis.xlsx')

tfs = df['TF'].tolist()
tfs = list(set(tfs))

cytokines = list(set(df['Bait gene'].tolist()))
print(cytokines)
cytokines = [str(i) for i in cytokines]
print(cytokines)

# Do blast Cytokine sequence vs human genome to get coordinates
info = dict(list(zip(df['Bait gene'], df['bait sequence'])))
info = {str(i):j for i,j in info.items()}

o1 = open('tmp2.fasta', 'w')
for k,v in info.items():
	#print(k,v)
	o1.write('>' + str(k) + '\n' + v + '\n')
o1.close()

subprocess.call('blastn -db GRCh38.fa -query tmp2.fasta -out blast.results.out -outfmt "6 qgi qacc sacc qlen qstart qend sstart send evalue length mismatch"', shell = True)

files = os.listdir('gtrd/gtrd2/egrid/bigBeds/hg38/ChIP-seq/Peaks/MACS2')

results = pd.DataFrame()

for cytokine in cytokines:

	print('cytokine', cytokine)
	sequence = info[cytokine]

	seq_len = len(sequence)

	tmp_df = pd.read_table('blast.results.out', names = ['index', 'cytokine', 'chromosome', 'lenx', 'qstart', 'qend', 'start', 'end', 'evalue', 'qlen', 'x'])
	tmp_df = tmp_df[(tmp_df['cytokine'] == cytokine) & (tmp_df['qend'] == seq_len) & (tmp_df['qstart'] <= 2)]

	print(tmp_df)
	if tmp_df.shape[0] > 1:
		print('More than 1 row')
		tmp_df = tmp_df.head(1)

	chromosome = str(tmp_df['chromosome'].iloc[0])
	start = tmp_df['start'].iloc[0]
	end = tmp_df['end'].iloc[0]

	if end <= start:
		real_end = start
		start = end
		end = real_end

	subdf = df[df['Bait gene'] == cytokine]

	for tf in tfs:

		tf_files = [i for i in files if i.split('_')[1] == tf]
		print('tf', tf)
		for file in tf_files:
			peakID = file.split('_')[0]
			file = 'gtrd/gtrd2/egrid/bigBeds/hg38/ChIP-seq/Peaks/MACS2/' + file
			BBIFile = bbi.open(file)
			all_chromosomes = list(dict(BBIFile.chromsizes).keys())
			if chromosome not in all_chromosomes:
				continue

			peaks = BBIFile.fetch_intervals(chromosome, start-200000, end+200000)

			if peaks.shape[0] == 0:
				continue

			tmp_info = metadata[metadata['input'] == peakID]

			if tmp_info.shape[0] == 0:
				print(tf, 'N peaks:', peaks.shape[0], 'No info metadata')
				continue

			peaks = peaks[(peaks['abs_summit'] >= start) & (peaks['abs_summit']<= end)]

			peaks['cytokine'] = [cytokine] * peaks.shape[0]
			peaks['id'] = [ tmp_info.iloc[0, 0] ] * peaks.shape[0]
			peaks['treatment'] = [ tmp_info.iloc[0, 3] ] * peaks.shape[0]
			peaks['antibody'] = [ tmp_info.iloc[0, 1] ] * peaks.shape[0]
			peaks['cell_line'] = [tmp_info.iloc[0, 8]] * peaks.shape[0]

			results = pd.concat([results, peaks])

print(results)

results.to_excel('MACS2_Peaks_TFs_final_26_june.xlsx')

