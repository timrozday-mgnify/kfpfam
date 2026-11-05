"""Six-frame ORFs of reads (split at stops, keep >=minlen), as JSON for brute.py and FASTA for hmmsearch."""
import sys,json,random,itertools
bases='TCAG';aa='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
T={a+b+c:aa[i] for i,(a,b,c) in enumerate(itertools.product(bases,repeat=3))}
rc=str.maketrans('ACGTN','TGCAN')
fq,N,minlen,oj,of=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],sys.argv[5]
lines=open(fq).read().split('\n')
seqs=[(lines[i][1:].split()[0],lines[i+1]) for i in range(0,len(lines)-3,4)]
random.seed(5);random.shuffle(seqs);seqs=seqs[:N]
out=[]
with open(of,'w') as fa:
    for ri,(name,s) in enumerate(seqs):
        orfs=[]
        for st in (s,s[::-1].translate(rc)):
            for fr in range(3):
                p=''.join(T.get(st[i:i+3],'X') for i in range(fr,len(st)-2,3))
                orfs+=[o for o in p.split('*') if len(o)>=minlen]
        out.append(dict(id=f'R{ri}',orfs=orfs,true=None))
        for j,o in enumerate(orfs): fa.write(f'>R{ri}_{j}\n{o}\n')
json.dump(out,open(oj,'w'))
import statistics as st
print(len(out),'reads; orfs per read',st.mean(len(r['orfs']) for r in out),'residues per read',st.mean(sum(map(len,r['orfs'])) for r in out))
