from collections import defaultdict
from common import *
d=np.load('pfam.npz'); ga={a.split('.')[0]:g for a,g in zip(d['accs'],d['ga'])}
H=defaultdict(dict)
for h in domtbl('truth/reads10k.domtbl'):
    r=h['t'].split('_')[0]; H[r].setdefault(h['acc'],[]).append((h['domscore'],h['evalue'],ga[h['acc']]))
n=10000
print('reads with any model E<=10:',len(H)/n)
for name,f in (('E<=1e-3',lambda x:x[1]<=1e-3),('GA',lambda x:x[0]>=x[2]),('E<=1e-5',lambda x:x[1]<=1e-5)):
    c=[sum(1 for v in hs.values() if any(f(x) for x in v)) for hs in H.values()]
    c=np.array(c);print(name,'reads with a hit %.3f; models per hit read mean %.2f p90 %.0f max %d'%((c>0).sum()/n,c[c>0].mean(),np.percentile(c[c>0],90),c.max()))
