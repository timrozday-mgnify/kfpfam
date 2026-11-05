import numpy as np
A='ACDEFGHIKLMNPQRSTVWY'
AIDX=np.full(256,20,dtype=np.int8)
for i,c in enumerate(A): AIDX[ord(c)]=i
Q=np.array([.0787945,.0151600,.0535222,.0668298,.0397062,.0695071,.0229198,.0590092,.0594422,.0963728,.0237718,.0414386,.0482904,.0395639,.0540978,.0683364,.0540687,.0673417,.0114135,.0304133])
Q/=Q.sum()
def load(path='pfam.npz'):
    d=np.load(path)
    P=d['p'].astype(np.float64)
    LO=np.log2(P/Q)
    LO=np.hstack([LO,np.full((LO.shape[0],1),-100.0)]).astype(np.float32)
    return d,P,LO
def enc(s): return AIDX[np.frombuffer(s.encode(),dtype=np.uint8)]
def seq_fasta(p):
    n=None;s=[]
    for l in open(p):
        if l[0]=='>':
            if n: yield n,''.join(s)
            n=l[1:].split()[0];s=[]
        else: s.append(l.strip())
    yield n,''.join(s)
def domtbl(p):
    for l in open(p):
        if l[0]=='#': continue
        f=l.split()
        yield dict(t=f[0],tlen=int(f[2]),q=f[3],acc=f[4].split('.')[0],qlen=int(f[5]),evalue=float(f[6]),score=float(f[7]),domi=float(f[12]),domscore=float(f[13]),hf=int(f[15]),ht=int(f[16]),af=int(f[17]),at=int(f[18]),ef=int(f[19]),et=int(f[20]))
