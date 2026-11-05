import csv
from collections import defaultdict
from common import *
d=np.load('pfam.npz'); ga={a.split('.')[0]:g for a,g in zip(d['accs'],d['ga'])}
hits=defaultdict(dict)
for h in domtbl('truth/fragall.domtbl'):
    q=h['t']
    if h['domscore']>hits[q].get(h['acc'],(-1e9,))[0]: hits[q][h['acc']]=(h['domscore'],h['domi'],h['evalue'])
for pre,tag in (('A','frag50'),('B','frag33')):
    rows=list(csv.DictReader(open(f'truth/{tag}.tsv'),delimiter='\t'))
    n=len(rows);e10=ge3=gaok=top=topga=0;nh=[];scores=[];ng=[]
    det=open(f'truth/{tag}.detected','w');detga=open(f'truth/{tag}.detga','w')
    for r in rows:
        q=pre+r['id'][1:];h=hits.get(q,{});a=r['acc']
        nh.append(len(h)); ng.append(sum(1 for k,v in h.items() if v[0]>=ga[k]))
        if a in h:
            e10+=1;s=h[a][0];scores.append(s)
            ge3+=h[a][2]<=1e-3
            if s>=ga[a]: gaok+=1;detga.write(r['id']+'\n')
            if h[a][2]<=1e-3: det.write(r['id']+'\n')
            if s>=max(v[0] for v in h.values()): top+=1
    print(tag,'n',n,'true model reported E<=10 %.3f; E<=1e-3 (Z=40k) %.3f; >=GA %.3f; true model is best %.3f; mean models hit E<=10 %.2f, >=GA %.3f; median true score %.1f'%(e10/n,ge3/n,gaok/n,top/n,np.mean(nh),np.mean(ng),np.median(scores)))
    gas=[ga[r['acc']] for r in rows];print(' GA of the fragments families: median',np.median(gas),'p10',np.percentile(gas,10),'p90',np.percentile(gas,90))
