import random,sys
random.seed(7)
recs=[];n=None;s=[]
for l in open(sys.argv[1]):
    if l[0]=='>':
        if n: recs.append((n,''.join(s)))
        n=l[1:].split()[0];s=[]
    else: s.append(l.strip())
recs.append((n,''.join(s)))
random.shuffle(recs)
sel=recs[:int(sys.argv[2])]
with open(sys.argv[3],'w') as o:
    for n,s in sel: o.write(f'>{n}\n{s}\n')
print(len(sel),sum(len(s) for _,s in sel))
