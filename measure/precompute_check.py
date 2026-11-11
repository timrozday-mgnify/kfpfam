import numpy as np, itertools, random
A=b"ACDEFGHIKLMNPQRSTVWY"; T=np.full(256,20,np.uint8); T[list(A)]=np.arange(20)
def orfs(p,limit=None):
    with open(p,'rb') as f:
        for i,l in enumerate(f):
            if limit and i>=limit*2: return
            if l[0]!=62: yield T[np.frombuffer(l.strip(),np.uint8)]
K=4; n=160000
idx=open('data/index/k4T7.idx','rb'); h=np.fromfile(idx,np.uint64,2); k=int(h[0]); tot=int(h[1]); starts=np.fromfile(idx,np.uint64,n+1).astype(np.int64)
entries=np.memmap('data/index/k4T7.idx',dtype=np.uint32,mode='r',offset=16+8*(n+1),shape=(tot,))
lens=np.diff(starts)
lo=np.memmap('data/prep/profile.f32',dtype=np.float32,mode='r').reshape(-1,21); offs=np.fromfile('data/prep/offs.u32',np.uint32).astype(np.int64)
top=lo[:,:20].max(1); print('columns',len(top),'mean top bits %.2f p10 %.2f p90 %.2f'%(top.mean(),*np.quantile(top,[.1,.9])))
print('entries %d; list length mean %.0f p50 %.0f p99 %.0f max %d'%(tot,lens.mean(),np.median(lens),np.quantile(lens,.99),lens.max()))
# skew and window repeats, 50k reads
freq=np.zeros(n,np.int64); win={w:[] for w in (8,12,16,20)}; res=0
pw=20**np.arange(K-1,-1,-1)
R='data/real/ERR10889147/'
for s in orfs(R+'orfs.subset.faa'):
    res+=len(s)
    if len(s)<K: continue
    v=np.lib.stride_tricks.sliding_window_view(s,K); ok=(v<20).all(1); np.add.at(freq,(v[ok]*pw).sum(1),1)
    for w,l in win.items():
        if len(s)>=w: l.append(np.lib.stride_tricks.sliding_window_view(s,w).astype(np.uint8))
hitk=freq*lens; o=np.argsort(-hitk); c=np.cumsum(hitk[o])/hitk.sum()
print('hits per read %.0f; hits from top 0.1%%/1%%/10%% of 4-mers: %.3f %.3f %.3f'%(hitk.sum()/50000,c[159],c[1599],c[15999]))
print('top 4-mers',[''.join(chr(A[(i//20**(3-p))%20]) for p in range(4)) for i in o[:8]], 'residues',res)
for w,l in win.items():
    a=np.concatenate(l); u=np.unique(a.view(np.dtype((np.void,w))))
    print(f'windows {w}: {len(a)} total, distinct {len(u)/len(a):.3f}')
del win
# bounds on hits from 200 ORFs of a random sample
random.seed(1); all_=[s for s in orfs('data/index/first300.faa') if len(s)>=K]; print('orfs',len(all_),'residues',sum(map(len,all_)))
cum=np.concatenate([[0],np.cumsum(top.astype(np.float64))])
mm=np.stack([lo[offs[m]:offs[m+1],:20].max(0) for m in range(len(offs)-1)]) # model max per residue
mm=np.hstack([mm,np.zeros((len(mm),1))])
H=0;A_=0;B=0;C=0; samp=[]
for s in all_:
    for j in range(len(s)-K+1):
        q=s[j:j+K]
        if (q>=20).any(): continue
        code=int((q*pw).sum()); e=np.asarray(entries[starts[code]:starts[code+1]]); 
        if not len(e): continue
        m=(e>>16).astype(np.int64); col=(e&0xffff).astype(np.int64); a=offs[m]; cols=offs[m+1]-a
        back=np.minimum(j,col); ln=np.minimum(len(s)-j+back,cols-col+back); st=a+col-back
        ubA=cum[st+ln]-cum[st]; H+=len(e); A_+=(ubA<22).sum()
        # bound B: per-model best per residue, over the read's residues on the diagonal
        tab=np.cumsum(np.hstack([np.zeros((len(m),1)),mm[m][:,s]]),1); i0=j-back; ub=tab[np.arange(len(m)),i0+ln]-tab[np.arange(len(m)),i0]
        B+=(ub<22).sum(); C+=(np.minimum(ubA,ub)<22).sum()
        if random.random()<0.02: samp.append((s,j,back,ln,st,ub,ubA))
print('hits',H,'per read',H/300,'reject column-top bound %.3f model-residue bound %.3f either %.3f'%(A_/H,B/H,C/H))
print('ratio bound/22 median',np.median(np.concatenate([x[5] for x in samp])), 'bound p1',np.quantile(np.concatenate([x[5] for x in samp]),.01))
