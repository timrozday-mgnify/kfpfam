import json,csv
from common import *
frag=dict(seq_fasta('truth/frag50.faa'))
rows=list(csv.DictReader(open('truth/frag50.tsv'),delimiter='\t'))
det=set(l.strip() for l in open('truth/frag50.detected'))
out=[dict(id=r['id'],orfs=[frag[r['id']]],true=r['acc'],hmm=r['id'] in det) for r in rows[2000:2300]]
json.dump(out,open('truth/true300.json','w'))
b=json.load(open('truth/reads10k.json'))[:200]
json.dump(b,open('truth/bg200.json','w'))
print(len(out),len(b))
