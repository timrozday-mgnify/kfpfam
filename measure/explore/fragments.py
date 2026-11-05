"""Fragments from GA-passing domains of MGnify proteins: abundance-weighted (each domain equally)."""
import sys, random
from common import *
domfile,fafile,out,L,N,minov=sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6])
random.seed(11+L)
seqs=dict(seq_fasta(fafile))
doms=[d for d in domtbl(domfile) if d['t'] in seqs]
print(len(doms),'domains',len({d['acc'] for d in doms}),'families',len({d['t'] for d in doms}),'proteins with a hit of',len(seqs))
n=0
with open(out+'.faa','w') as fa,open(out+'.tsv','w') as tv:
    tv.write('id\tprotein\tacc\tdom_af\tdom_at\thf\tht\tstart\tend\tdomscore\n')
    while n<N:
        d=random.choice(doms); s=seqs[d['t']]
        if len(s)<L: continue
        lo=max(0,d['af']-1-(L-minov)); hi=min(len(s)-L,d['at']-minov)
        if lo>hi: 
            if d['at']-d['af']+1<minov: continue
            lo=max(0,d['af']-1);hi=min(len(s)-L,d['at']-1)
            if lo>hi: continue
        st=random.randint(lo,hi)
        f=s[st:st+L]
        fa.write(f'>F{n}\n{f}\n')
        tv.write(f"F{n}\t{d['t']}\t{d['acc']}\t{d['af']}\t{d['at']}\t{d['hf']}\t{d['ht']}\t{st}\t{st+L}\t{d['domscore']}\n")
        n+=1
print('wrote',n)
