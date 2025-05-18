#!/usr/bin/env python3

import optparse
from Bio import SeqIO
from collections import defaultdict

################################# Command line options

desc='Calculate the best hit for each query based on bitscore and return the coverage in percentage on the best hit. Good to check quantitative view of how well your assembled transcripts or queries align to a reference'

parser = optparse.OptionParser(description=desc, version='%prog version 0.1 - 16-05-2025 - Author: FCicconardi')

parser.add_option('-b', '--Blast-hits', dest='bhits', help='Blast/Diamond Protein hits in  output-6 format. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-q', '--Query-Fasta', dest='qfasta', help='Query fasta file. Mandatory opt.', action='store', metavar='FILE')
parser.add_option('-t', '--Target-Fasta', dest='tfasta', help='Target fasta file. Mandatory opt.', action='store', metavar='FILE')

(opts, args) = parser.parse_args()

mandatories = ['bhits','qfasta','tfasta']
for m in mandatories:
        if not opts.__dict__[m]:
                print("\nWARNING! One or more options not specified\n")
                parser.print_help()
                exit(-1)

############################## Reading files and parametersfrom sys import argv

def ParseFastas(fasta_file):
	DictDB=defaultdict(list)
	for record in SeqIO.parse(fasta_file, "fasta"):
		DictDB[record.id].append((' '.join(record.description.split(' ')[1:]),len(record.seq)))

	return(DictDB)

#QueryDB=defaultdict(list)
#Target=defaultdict(list)


with open(opts.qfasta) as fasta_file:
	QueryDB=ParseFastas(fasta_file)
	
with open(opts.tfasta) as fasta_file:
	TargetDB=ParseFastas(fasta_file)


HitsDB=defaultdict(list)

with open(opts.bhits) as hits:
	for el in hits:
		query_id, db_id, pident, _, _, _, qstart, qend, sstart, send, evalue, bitscore = el.strip().split('\t')[:12]
		HitsDB[query_id].append((db_id, pident, qstart, qend, sstart, send, evalue, bitscore))


print('#query_id\ttarget_id\tpident\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tpct_aln_ratio_q_t\tpct_query_len_aligned\tpct_hit_len_aligned\tdesc')

for query in HitsDB:
	bestHit=defaultdict(list)
	for i,hit in enumerate(HitsDB[query]):
		bestHit[HitsDB[query][i][7]].append(hit)

	BestScore=max(bestHit.keys())
	BestHit=bestHit[BestScore][0]

	db_id,qstart,qend,tstart,tend=BestHit[0],int(BestHit[2]),int(BestHit[3]),int(BestHit[4]),int(BestHit[5])

	if tend > tstart: TargetCov=tend-tstart
	else: TargetCov=tstart-tend

	if qend > qstart: QueryCov=qend-qstart
	else: QueryCov=tstart-tend

	Tlen=int(TargetDB[db_id][0][1])
	TDesc=TargetDB[db_id][0][0]

	Qlen=int(QueryDB[query][0][1])


	PctHitLenAli=TargetCov/Tlen*100

	PctQueryLenAli=QueryCov/Qlen*100

	PctQTratio=PctHitLenAli/PctQueryLenAli

	Best='\t'.join(BestHit)	
	print(f"{query}\t{Best}\t{PctQTratio}\t{PctQueryLenAli}\t{PctHitLenAli}\t{TDesc}")
#	print(query,'\t'.join(BestHit))
