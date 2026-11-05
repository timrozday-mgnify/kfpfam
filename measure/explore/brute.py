"""Whole-Pfam seed scan by brute force. For each read (one or several ORFs) lists the models
passing: one seed >=T; two seeds on one diagonal (+-3), non-overlapping. Also the same after a
per-read rank rule: only models within D bits of the read's best seed score.
A k-mer is 'in the index' exactly when its score at some model window is >=T, so no enumeration."""
import sys,json,time
from common import *
d,P,LO=load()
offs=d['offs']; nmodel=len(offs)-1
model_of=np.repeat(np.arange(nmodel),np.diff(offs))
SEEDS=json.loads(sys.argv[3])          # {"k5":[[0,1,2,3,4],[10,12,14]],...}
reads=json.load(open(sys.argv[1]))     # [{"id":..,"orfs":[...],"true":acc or null}]
acc2i={a.split('.')[0]:i for i,a in enumerate(d['accs'])}
CH=500000
BW=4
res={}
def scan(seq,shape,Tmin):
    f=enc(seq);span=shape[-1]+1;L=len(f)
    if L<span: return np.zeros(0,int),np.zeros(0,int),np.zeros(0,np.float32)
    II=[];JJ=[];SS=[]
    M=LO.shape[0]
    for a in range(0,M,CH):
        b=min(M,a+CH+span-1)
        S=LO[a:b][:,f]
        n=b-a-span+1
        W=np.zeros((n,L-span+1),dtype=np.float32)
        for o in shape: W+=S[o:o+n,o:o+L-span+1]
        i,j=np.nonzero(W>=Tmin)
        II.append(i+a);JJ.append(j);SS.append(W[i,j])
    I=np.concatenate(II);J=np.concatenate(JJ);Sc=np.concatenate(SS)
    ok=model_of[I]==model_of[I+span-1]
    return I[ok],J[ok],Sc[ok]
def group_stats(I,J,Sc,span):
    """models with one seed; models with a non-overlapping pair on the diagonal +-3, best pair score."""
    m=model_of[I]
    single=np.zeros(nmodel,np.float32); np.maximum.at(single,m,Sc)
    pair={}
    if len(I)<2: return single,pair
    diag=I-offs[m]-J+1000
    for shift in (0,BW//2):
        key=(m.astype(np.int64)<<16)|((diag+shift)//BW)
        o=np.argsort(key,kind='stable');k=key[o]
        st=np.flatnonzero(np.r_[True,k[1:]!=k[:-1]]);en=np.r_[st[1:],len(k)]
        for s,e in zip(st[en-st>1],en[en-st>1]):
            idx=o[s:e];jj=J[idx]
            if jj.max()-jj.min()<span: continue
            ss=Sc[idx];best=0.
            for x in range(len(idx)):
                for y in range(x+1,len(idx)):
                    if abs(jj[x]-jj[y])>=span and ss[x]+ss[y]>best: best=float(ss[x]+ss[y])
            mi=int(m[idx[0]])
            if best>pair.get(mi,0): pair[mi]=best
    return single,pair
XS=[14,16,18,20,22,25,30]
def diagscores(I,J,Sc,seq):
    """Ungapped local score along each seeded diagonal (one row per diagonal): the best seed on it, its model, its score."""
    if len(I)==0: return np.zeros(0,int),np.zeros(0,np.float32),np.zeros(0,np.float32)
    f=enc(seq);L=len(f);dg=I-J
    u,inv=np.unique(dg,return_inverse=True)
    best=np.zeros(len(u),np.float32);np.maximum.at(best,inv,Sc)
    mods=model_of[I[np.unique(inv,return_index=True)[1]]]
    st=u[:,None]+np.arange(L)[None,:]
    bad=(st<0)|(st>=len(model_of))
    st2=np.clip(st,0,len(model_of)-1)
    M=LO[st2,f[None,:].repeat(len(u),0)].astype(np.float32)
    M[bad|(model_of[st2]!=mods[:,None])]=-100.0
    c=np.cumsum(M,axis=1)
    rmin=np.minimum.accumulate(np.hstack([np.zeros((len(u),1),np.float32),c[:,:-1]]),axis=1)
    return mods,best,(c-rmin).max(axis=1)
out=[];t0=time.time()
for ri,r in enumerate(reads):
    row={'id':r['id'],'true':r.get('true')}
    tm=acc2i.get(r.get('true'))
    for sn,(shape,Ts) in SEEDS.items():
        Tmin=min(Ts);span=shape[-1]+1
        hits=[scan(o,shape,Tmin) for o in r['orfs']]
        for T in Ts:
            single=np.zeros(nmodel,np.float32);pair={}
            for I,J,Sc in hits:
                ok=Sc>=T
                if not ok.any(): continue
                s,p=group_stats(I[ok],J[ok],Sc[ok],span)
                single=np.maximum(single,s)
                for k,v in p.items():
                    if v>pair.get(k,0): pair[k]=v
            ds=np.zeros(nmodel,np.float32)
            for (I,J,Sc),o in zip(hits,r['orfs']):
                ok=Sc>=T
                mods,bs,sc=diagscores(I[ok],J[ok],Sc[ok],o)
                if len(mods): np.maximum.at(ds,mods,sc)
            dmax=float(ds.max())
            hitm=np.flatnonzero(single>0)
            row[f'{sn}.{T}']=dict(
              n1=int(len(hitm)),n2=len(pair),
              dmax=dmax,dtrue=float(ds[tm]) if tm is not None else None,
              dn={x:int((ds>=x).sum()) for x in XS},
              dr10=int((ds>=max(dmax-10,0.001)).sum()) if dmax>0 else 0,
              dr5=int((ds>=max(dmax-5,0.001)).sum()) if dmax>0 else 0,
              # rank rules: within 10 bits of the read's best
              n1r=int((single>=single.max()-10).sum()) if len(hitm) else 0,
              n1r5=int((single>=single.max()-5).sum()) if len(hitm) else 0,
              n2r=sum(1 for v in pair.values() if v>=max(pair.values())-10) if pair else 0,
              t1=bool(tm is not None and single[tm]>0),
              t2=bool(tm is not None and tm in pair),
              t1r=bool(tm is not None and single[tm]>0 and single[tm]>=single.max()-10),
              t1r5=bool(tm is not None and single[tm]>0 and single[tm]>=single.max()-5),
              t2r=bool(tm is not None and tm in pair and pair[tm]>=max(pair.values())-10))
    out.append(row)
    if ri%20==0: print(ri,round(time.time()-t0),flush=True); json.dump(out,open(sys.argv[2],'w'))
json.dump(out,open(sys.argv[2],'w'))
print('done',round(time.time()-t0))
