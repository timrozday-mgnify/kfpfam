from common import *
d,P,LO=load(); offs=d['offs'];L=d['leng']
mx=LO[:,:20].max(1)
for k,Ts in ((4,[8,10]),(5,[8,10,12]),(6,[10,12]),(7,[12,14])):
    c=np.cumsum(np.r_[0,mx]);best=[]
    for i in range(len(L)):
        a,b=offs[i],offs[i+1]
        best.append((c[a+k:b+1]-c[a:b-k+1]).max() if b-a>=k else -1)
    best=np.array(best)
    print(k,{T:f'{(best<T).mean():.4f} of models have no seed' for T in Ts},'median best possible %.1f'%np.median(best))
# which amino acid columns: share of windows that can reach T=10 at k=5
c=np.cumsum(np.r_[0,mx]);w=c[5:]-c[:-5]
print('share of 5-windows able to reach 10 bits: %.3f; 12 bits %.3f; 8 bits %.3f'%((w>=10).mean(),(w>=12).mean(),(w>=8).mean()))
