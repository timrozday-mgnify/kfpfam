"""Recall against the true model on large samples: seed, two-hit, ungapped diagonal score.
Groups: what hmmsearch itself reports on that fragment (E<=1e-3; >=GA), and all domain fragments."""
import csv,json
from collections import defaultdict
from common import *
d,P,LO=load(); offs=d['offs']; acc2i={a.split('.')[0]:i for i,a in enumerate(d['accs'])}; ga=d['ga']
SEEDS=json.load(open('seeds.json'))
XS=[14,16,18,20,22,25,30]
def score(S,shape,Ts):
    M,L=S.shape;span=shape[-1]+1;res={}
    if M<span or L<span:
        for T in Ts: res[T]=(False,False,0.)
        return res
    W=np.zeros((M-span+1,L-span+1),np.float32)
    for o in shape: W+=S[o:o+M-span+1,o:o+L-span+1]
    for T in Ts:
        ii,jj=np.nonzero(W>=T)
        if len(ii)==0: res[T]=(False,False,0.);continue
        dg=ii-jj;p=False
        if len(ii)>1:
            for a in range(len(ii)):
                for b in range(a+1,len(ii)):
                    if abs(int(jj[a])-int(jj[b]))>=span and abs(int(dg[a])-int(dg[b]))<=3: p=True;break
                if p: break
        best=0.
        for dd in np.unique(dg):
            x=np.arange(max(0,-dd),min(L,M-dd))
            v=S[dd+x,x];c=0.;bb=0.
            for s in v:
                c=max(0.,c+float(s));bb=max(bb,c)
            best=max(best,bb)
        res[T]=(True,p,best)
    return res
stat=defaultdict(lambda:defaultdict(float))
def run(sets,name):
    for gname,seqs in sets.items():
        for m,seq in seqs:
            f=enc(seq);S=LO[offs[m]:offs[m+1]][:,f]
            for sn,(shape,Ts) in SEEDS.items():
                r=score(S,shape,Ts)
                for T,(a,b,c) in r.items():
                    st=stat[(name,gname,sn,T)]
                    st['n']+=1;st['s']+=a;st['p']+=b
                    for x in XS: st[x]+=c>=x
for tag,N in (('frag50',4000),('frag33',4000)):
    frag=dict(seq_fasta(f'truth/{tag}.faa'));rows=list(csv.DictReader(open(f'truth/{tag}.tsv'),delimiter='\t'))[:N]
    det=set(l.strip() for l in open(f'truth/{tag}.detected'));dg=set(l.strip() for l in open(f'truth/{tag}.detga'))
    sets={'all':[],'E3':[],'GA':[]}
    for r in rows:
        x=(acc2i[r['acc']],frag[r['id']]);sets['all'].append(x)
        if r['id'] in det: sets['E3'].append(x)
        if r['id'] in dg: sets['GA'].append(x)
    run(sets,tag)
# real reads: hmmsearch hits on the read's own ORFs
orfs=dict(seq_fasta('truth/reads10k.faa'))
E3={};GA={}
for h in domtbl('truth/reads10k.domtbl'):
    m=acc2i[h['acc']]
    if h['evalue']<=1e-3: E3[(h['t'],m)]=orfs[h['t']]
    if h['domscore']>=ga[m]: GA[(h['t'],m)]=orfs[h['t']]
import random
random.seed(1)
e=list(E3.items());random.shuffle(e);g=list(GA.items());random.shuffle(g)
run({'E3':[(k[1],v) for k,v in e[:3000]],'GA':[(k[1],v) for k,v in g[:3000]]},'reads')
with open('recall3.tsv','w') as o:
    o.write('set\tgroup\tseed\tT\tn\tsingle\tpair\t'+'\t'.join(f'd{x}' for x in XS)+'\n')
    for k,st in sorted(stat.items(),key=str):
        n=st['n'];o.write('\t'.join(map(str,k))+f"\t{int(n)}\t{st['s']/n:.3f}\t{st['p']/n:.3f}\t"+'\t'.join(f'{st[x]/n:.3f}' for x in XS)+'\n')
