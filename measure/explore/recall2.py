"""Recall by family traits, per read and per family; real data vs the HMM's own emissions;
sampled k-mers vs enumerated ones at the same index size."""
import csv,sys
from collections import defaultdict
from common import *
d,P,LO=load(); offs=d['offs']; acc2i={a.split('.')[0]:i for i,a in enumerate(d['accs'])}
rel=(P*np.log2(P/Q)).sum(1)
effn=d['effn'];leng=d['leng']
frag=dict(seq_fasta('truth/frag50.faa'));rows=list(csv.DictReader(open('truth/frag50.tsv'),delimiter='\t'))
det=set(l.strip() for l in open('truth/frag50.detected'))
detga=set(l.strip() for l in open('truth/frag50.detga'))
def seeds(S,k):
    M,L=S.shape
    if M<k or L<k: return None
    W=np.zeros((M-k+1,L-k+1),np.float32)
    for o in range(k): W+=S[o:o+M-k+1,o:o+L-k+1]
    return W
def diag(S,W):
    """Best ungapped local score on any diagonal the seed set W>=T touches; returns (maxseed, best diag score per seed bar)"""
    return None
res=defaultdict(list)
fam=defaultdict(list)
for r in rows:
    m=acc2i.get(r['acc'])
    f=enc(frag[r['id']]); S=LO[offs[m]:offs[m+1]][:,f]
    W5=seeds(S,5)
    if W5 is None: continue
    best=W5.max()
    # pair at T: any two non-overlapping seeds >=T on diagonal +-3
    def pair(T):
        ii,jj=np.nonzero(W5>=T)
        for a in range(len(ii)):
            for b in range(a+1,len(ii)):
                if abs(int(jj[a])-int(jj[b]))>=5 and abs((int(ii[a])-int(jj[a]))-(int(ii[b])-int(jj[b])))<=3: return True
        return False
    p10=pair(10)
    i2=(r['id'] in det)
    rr=dict(b=best,p10=p10,hmm=i2,ga=r['id'] in detga,relmean=rel[offs[m]:offs[m+1]].mean(),effn=effn[m],leng=leng[m],maxw=np.convolve(rel[offs[m]:offs[m+1]],np.ones(5),'valid').max() if leng[m]>=5 else 0,acc=r['acc'],hf=int(r['hf']),ht=int(r['ht']))
    res['all'].append(rr);fam[r['acc']].append(rr)
R=res['all'];import statistics as st
def report(rs,name):
    n=len(rs);h=[x for x in rs if x['hmm']]
    print(f'{name:28s} n={n:6d}  hmmE-3={len(h)/max(n,1):.2f} | seed k5>=10 all={np.mean([x["b"]>=10 for x in rs]):.2f} hmm={np.mean([x["b"]>=10 for x in h]) if h else 0:.2f} | k5>=8 hmm={np.mean([x["b"]>=8 for x in h]) if h else 0:.2f} | pair10 hmm={np.mean([x["p10"] for x in h]) if h else 0:.2f}')
report(R,'all')
for lo,hi in ((0,6),(6,8),(8,11),(11,16),(16,99)):
    report([x for x in R if lo<=x['maxw']<hi],f'max 5-window relent {lo}-{hi}')
for lo,hi in ((0,1),(1,3),(3,10),(10,100),(100,1e9)):
    report([x for x in R if lo<=x['effn']<hi],f'EFFN {lo}-{hi}')
for lo,hi in ((0,60),(60,150),(150,300),(300,1e9)):
    report([x for x in R if lo<=x['leng']<hi],f'model length {lo}-{hi}')
# family-level
A=[];B=[];C=[];D=[]
for a,rs in fam.items():
    h=[x for x in rs if x['hmm']]
    if len(h)<3: continue
    A.append(np.mean([x['b']>=10 for x in h])); B.append(np.mean([x['b']>=8 for x in h])); C.append(np.mean([x['p10'] for x in h]));D.append(len(h))
A=np.array(A);B=np.array(B);C=np.array(C)
print('families with >=3 hmm-detected fragments:',len(A))
for nm,X in (('k5>=10',A),('k5>=8',B),('pair10',C)):
    print(f' per-family recall {nm}: mean {X.mean():.2f}, p10 {np.percentile(X,10):.2f}, med {np.median(X):.2f}; zero-recall families {np.mean(X==0):.3f}; <0.5 {np.mean(X<0.5):.3f}')
print('distinct families sampled',len(fam))
# Index-size budget: enumerated top-LR vs random hmmemit-style sampling at one model sample
rng=np.random.default_rng(3)
Tgrid=[8,10,12]
ms=rng.choice(len(leng),200,replace=False)
print('\nsampling vs enumerating, k=5, per window, windows pooled over 200 models')
tot=defaultdict(float);nw=0
from itertools import product
for m in ms:
    p=P[offs[m]:offs[m+1]]
    if len(p)<5: continue
    for st_ in rng.choice(len(p)-4,min(3,len(p)-4),replace=False):
        w=p[st_:st_+5]
        # whole distribution over 3.2M 5-mers
        pk=w[0]
        for c in w[1:]: pk=np.outer(pk,c).ravel()
        qk=Q.copy()
        for c in range(4): qk=np.outer(qk,Q).ravel()
        lr=np.log2(pk/qk)
        nw+=1
        for T in Tgrid:
            sel=lr>=T;n=sel.sum()
            tot[('enum_LR',T,'n')]+=n;tot[('enum_LR',T,'cover')]+=pk[sel].sum();tot[('enum_LR',T,'bg')]+=qk[sel].sum()
        for n in (100,1000,10000):
            # top-n by probability
            idx=np.argpartition(-pk,n)[:n]
            tot[('top_P',n,'cover')]+=pk[idx].sum();tot[('top_P',n,'bg')]+=qk[idx].sum()
            idx=np.argpartition(-lr,n)[:n]
            tot[('top_LR',n,'cover')]+=pk[idx].sum();tot[('top_LR',n,'bg')]+=qk[idx].sum()
            # n emission samples, deduplicated: coverage = sum p(1-(1-p)^n), unique= sum 1-(1-p)^n
            inc=1-np.power(1-pk,n)
            tot[('sampled',n,'cover')]+=(pk*inc).sum();tot[('sampled',n,'bg')]+=(qk*inc).sum();tot[('sampled',n,'uniq')]+=inc.sum()
for k,v in sorted(tot.items(),key=str): print(k,'%.4g'%(v/nw))
print('windows',nw)
