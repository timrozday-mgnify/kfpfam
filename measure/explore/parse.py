"""Parse Pfam-A.hmm.gz match emissions into one npz: scores in bits vs HMMER's background."""
import gzip, sys, numpy as np
BG=None
ms=[];names=[];accs=[];leng=[];ga=[];effn=[];nseq=[];offs=[0]
ins_open=[];
with gzip.open(sys.argv[1],'rt') as f:
    state=0
    for l in f:
        if l.startswith('NAME '): names.append(l.split()[1])
        elif l.startswith('ACC '): accs.append(l.split()[1])
        elif l.startswith('LENG '): L=int(l.split()[1]);leng.append(L)
        elif l.startswith('GA '): ga.append(float(l.split()[1]))
        elif l.startswith('EFFN '): effn.append(float(l.split()[1]))
        elif l.startswith('NSEQ '): nseq.append(int(l.split()[1]))
        elif l.startswith('HMM  '):
            next(f); # transitions header
            comp=next(f); ins=next(f); tr=next(f)
            for i in range(L):
                m=next(f).split(); next(f); t=next(f).split()
                ms.append(np.array(m[1:21],dtype=np.float32))
                ins_open.append(float(t[1]) if t[1]!='*' else 99.)
            offs.append(len(ms))
ms=np.stack(ms); P=np.exp(-ms)
print(ms.shape, P.sum(1).min(), P.sum(1).max())
np.savez(sys.argv[2],p=P.astype(np.float32),offs=np.array(offs),names=np.array(names),accs=np.array(accs),leng=np.array(leng),ga=np.array(ga),effn=np.array(effn),nseq=np.array(nseq),ins=np.array(ins_open,dtype=np.float32))
