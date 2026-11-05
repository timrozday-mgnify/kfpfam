import json,numpy as np
bg=json.load(open('truth/bg60.out.json')); tr=json.load(open('truth/true150.out.json'))
reads=json.load(open('truth/true150.json'));hm={r['id']:r['hmm'] for r in reads}
hm_tr=[r for r in tr if hm[r['id']]]
print('reads',len(bg),'true',len(tr),'hmm-detected true',len(hm_tr))
keys=[k for k in bg[0] if '.' in k]
print('seed.T | bg candidates per read: one seed mean/med | two-hit mean/med | d>=20 d>=22 d>=25 mean/med | d>=22 & within10 of top mean/med || true (E3): seed two-hit d22 d25 top10')
for k in keys:
    def mm(f,src=bg): 
        x=[f(r[k]) for r in src];return f'{np.mean(x):8.1f}/{np.median(x):7.0f}'
    within=lambda v,src=None:None
    def rate(f): return np.mean([f(r[k]) for r in hm_tr])
    print(f'{k:7s} {mm(lambda v:v["n1"])} {mm(lambda v:v["n2"])} {mm(lambda v:v["dn"]["20"])} {mm(lambda v:v["dn"]["22"])} {mm(lambda v:v["dn"]["25"])} dr10 {mm(lambda v:v["dr10"])} || {rate(lambda v:v["t1"]):.2f} {rate(lambda v:v["t2"]):.2f} {rate(lambda v:(v["dtrue"] or 0)>=22):.2f} {rate(lambda v:(v["dtrue"] or 0)>=25):.2f} top10 {rate(lambda v:(v["dtrue"] or 0)>0 and v["dtrue"]>=v["dmax"]-10):.2f} top5 {rate(lambda v:(v["dtrue"] or 0)>0 and v["dtrue"]>=v["dmax"]-5):.2f}')
# rule: d>=X and within 10 of best: need per-read scores; approximated by min of dn[X], dr10
k='k5.10'
for k in ('k5.10','k5.8','k4.8','k6.10'):
    x=[min(r[k]['dn']['22'],r[k]['dr10']) for r in bg];print(k,'min(d22,top10) mean/med',np.mean(x),np.median(x))
    x=[r[k]['dmax'] for r in bg];print('  bg read best diag score: pct',np.percentile(x,[10,50,90]).round(1))
    x=[r[k]['dtrue'] for r in hm_tr];print('  true read true diag: pct',np.percentile(x,[10,50,90]).round(1))
    x=[r[k]['dmax']-r[k]['dtrue'] for r in hm_tr];print('  true read best minus true: pct',np.percentile(x,[50,90,95,99]).round(1))
