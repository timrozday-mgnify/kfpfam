"""Exact per-window k-mer statistics: for each threshold T on window log-odds bits,
N = kmers in set, C = share of the model's emissions covered, F = background hit rate."""
import numpy as np, sys
from itertools import product
A='ACDEFGHIKLMNPQRSTVWY'
Q=np.array([.0787945,.0151600,.0535222,.0668298,.0397062,.0695071,.0229198,.0590092,.0594422,.0963728,.0237718,.0414386,.0482904,.0395639,.0540978,.0683364,.0540687,.0673417,.0114135,.0304133],dtype=np.float64)
Q/=Q.sum()
d=np.load(sys.argv[1]); P=d['p'].astype(np.float64); offs=d['offs']; leng=d['leng']; effn=d['effn']
LO=np.log2(P/Q)
rel=(P*LO).sum(1)
print('per-column relative entropy (bits): mean %.3f  p10 %.3f med %.3f p90 %.3f p99 %.3f'%(rel.mean(),*np.percentile(rel,[10,50,90,99])))
ent=-(P*np.log2(P)).sum(1)
print('per-column entropy (bits): mean %.3f (20**k eff size: 2**(k*H)) => eff letters %.1f'%(ent.mean(),2**ent.mean()))
print('max per-column LO: mean %.2f, p90 %.2f'%(LO.max(1).mean(),np.percentile(LO.max(1),90)))
mrel=np.array([rel[offs[i]:offs[i+1]].mean() for i in range(len(leng))])
print('per-model mean relent: min %.2f p10 %.2f med %.2f p90 %.2f max %.2f'%(mrel.min(),*np.percentile(mrel,[10,50,90]),mrel.max()))
print('max k=5 window relent per model (bits): ',end='')
mw=[]
for i in range(len(leng)):
    r=rel[offs[i]:offs[i+1]]
    if len(r)<5: mw.append(r.sum());continue
    mw.append(np.convolve(r,np.ones(5),'valid').max())
mw=np.array(mw);print('p10 %.1f med %.1f p90 %.1f'%tuple(np.percentile(mw,[10,50,90])))
BIN=0.25; NB=int(80/BIN); OFF=int(40/BIN)
def window_tables(p,lo):
    # three-way DP on binned scores: count, model mass, bg mass
    c=np.zeros(NB);m=np.zeros(NB);q=np.zeros(NB);c[OFF]=1;m[OFF]=1;q[OFF]=1
    for col in range(p.shape[0]):
        cc=np.zeros(NB);mm=np.zeros(NB);qq=np.zeros(NB)
        sh=np.round(lo[col]/BIN).astype(int)
        for a in range(20):
            s=sh[a]; 
            if s>=0:
                cc[s:]+=c[:NB-s]; mm[s:]+=m[:NB-s]*p[col,a]; qq[s:]+=q[:NB-s]*Q[a]
            else:
                cc[:s]+=c[-s:]; mm[:s]+=m[-s:]*p[col,a]; qq[:s]+=q[-s:]*Q[a]
        c,m,q=cc,mm,qq
    return c,m,q
rng=np.random.default_rng(1)
models=rng.choice(len(leng),400,replace=False)
Ts=np.arange(0,32,1.0)
for k in (4,5,6,7,8):
    Nt=np.zeros(len(Ts));Ct=np.zeros(len(Ts));Ft=np.zeros(len(Ts));Fmax=[];w=0;Nwin=0
    bigC=np.zeros(len(Ts))
    for mi in models:
        a,b=offs[mi],offs[mi+1]
        if b-a<k: continue
        idx=rng.choice(b-a-k+1,min(10,b-a-k+1),replace=False)  # 10 windows a model
        wt=(b-a-k+1)/len(idx)
        for st in idx:
            c,m,q=window_tables(P[a+st:a+st+k],LO[a+st:a+st+k])
            ccum=np.cumsum(c[::-1])[::-1];mcum=np.cumsum(m[::-1])[::-1];qcum=np.cumsum(q[::-1])[::-1]
            I=np.minimum((OFF+np.round(Ts/BIN)).astype(int),NB-1)
            Nt+=wt*ccum[I];Ct+=wt*mcum[I];Ft+=wt*qcum[I];Nwin+=wt
    scale=leng.sum()/Nwin   # to all of Pfam (windows ≈ states)
    print(f'\n k={k}  windows sampled-weight {Nwin:.0f}; 20^k={20**k:.3g}')
    print(' T(bits)  N/window  all-Pfam-entries  meanCover  bgHit/window  coverage/bg')
    for j,T in enumerate(Ts):
        if T%2: continue
        nw=Nt[j]/Nwin;cw=Ct[j]/Nwin;fw=Ft[j]/Nwin
        if nw<1e-3: continue
        print(f' {T:5.0f} {nw:12.3g} {nw*leng.sum():14.3g} {cw:9.3f} {fw:12.3g} {cw/max(fw,1e-30):10.3g}')
