"""Fragment recall against the true model for seeds and thresholds, single and two-hit."""
import sys,csv,json
from common import *
from collections import defaultdict
d,P,LO=load()
offs=d['offs']; acc2i={a.split('.')[0]:i for i,a in enumerate(d['accs'])}
frag=dict(seq_fasta(sys.argv[1]+'.faa'))
rows=list(csv.DictReader(open(sys.argv[1]+'.tsv'),delimiter='\t'))
detected=set(l.strip() for l in open(sys.argv[2])) if len(sys.argv)>2 and sys.argv[2]!='-' else None
L=lambda n:None
SEEDS={'k4':(0,1,2,3),'k5':(0,1,2,3,4),'k6':(0,1,2,3,4,5),'k7':(0,1,2,3,4,5,6),
       's5a':(0,1,3,4,6),'s5b':(0,2,3,5,6),'s5c':(0,1,2,4,5),'s6a':(0,1,2,4,6,7),'s6b':(0,2,3,4,6,8)}
GRID={'k4':[6,8,10],'k5':[8,10,12,14],'k6':[10,12,14,16],'k7':[12,14,16,18],'s5a':[8,10,12],'s5b':[8,10,12],'s5c':[8,10,12],'s6a':[10,12,14],'s6b':[10,12,14]}
SETS={'k4':['k4'],'k5':['k5'],'k6':['k6'],'k7':['k7'],'s5x3':['s5a','s5b','s5c'],'k5+s5x3':['k5','s5a','s5b','s5c'],'s6x2':['s6a','s6b']}
stat=defaultdict(lambda:[0,0,0]) 
# keys: (set,T,mode,group)
def windows(S,seed):
    span=seed[-1]+1
    if S.shape[0]<span or S.shape[1]<span: return None,span
    M,Ln=S.shape
    W=np.zeros((M-span+1,Ln-span+1),dtype=np.float32)
    for o in seed: W+=S[o:o+M-span+1,o:o+Ln-span+1]
    return W,span
out={}
for r in rows:
    m=acc2i.get(r['acc'])
    if m is None: continue
    f=enc(frag[r['id']]); S=LO[offs[m]:offs[m+1]][:,f]
    hits={}
    for sd,shape in SEEDS.items():
        W,span=windows(S,shape)
        hits[sd]=(W,span)
    grp=['all']
    if detected is not None: grp.append('hmm' if r['id'] in detected else 'nohmm')
    for sn,members in SETS.items():
        Ts=GRID[members[0]]
        for T in Ts:
            pts=[]
            for sd in members:
                W,span=hits[sd]
                if W is None: continue
                ii,jj=np.nonzero(W>=T)
                pts+= [(int(i),int(j),span) for i,j in zip(ii,jj)]
            single=len(pts)>0
            pair=False
            if len(pts)>1:
                for a in range(len(pts)):
                    for b in range(a+1,len(pts)):
                        i1,j1,sp1=pts[a];i2,j2,sp2=pts[b]
                        if abs(j2-j1)>=max(sp1 if j2>j1 else sp2,1) and abs((i2-j2)-(i1-j1))<=3:
                            pair=True;break
                    if pair: break
            for g in grp:
                st=stat[(sn,T,g)]
                st[0]+=1;st[1]+=single;st[2]+=pair
print('set\tT\tgroup\tn\tsingle\tpair')
for k in sorted(stat): 
    n,s,p=stat[k]
    print(f'{k[0]}\t{k[1]}\t{k[2]}\t{n}\t{s/n:.3f}\t{p/n:.3f}')
