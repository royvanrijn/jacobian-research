import pickle,re,numpy as np
p=71
S=pickle.load(open('hne0_polred.pkl','rb'))
pat=re.compile(r'^C\|(\d+)\|(\d+)\|(\d+)\|(\d+)\|')
cols=[]
for l in open('hcert71.out'):
    m=pat.match(l)
    if m: cols.append(tuple(map(int,m.groups())))
assert len(cols)==385
rows={(1,0,0)}
for i,a,b,c in cols:
    for e in S[i-1]: rows.add((a+e[0],b+e[1],c+e[2]))
rows=sorted(rows,key=lambda m:(sum(m),m)); ri={m:i for i,m in enumerate(rows)}
assert len(rows)==402
# reduce coefficient dict to F_p^5
def coeff_vec(c):
    out=np.zeros(5,dtype=np.int64)
    for j,(n,d) in c.items(): out[j]=(n%p)*pow(d%p,-1,p)%p
    return out
# minpoly: w^5 - w^4 +3w^3+3w^2+26 =0
# w^k reduction, return 5 coeffs
poww=[]
for k in range(9):
    v=[0]*(max(k+1,5));v[k]=1
    for d in range(len(v)-1,4,-1):
        q=v[d]%p
        if q:
            v[d]=0
            # w^d = w^(d-5)*(w^4 -3w^3 -3w^2 -26)
            s=d-5
            for off,co in [(4,1),(3,-3),(2,-3),(0,-26)]:
                if s+off>=len(v):v.extend([0]*(s+off-len(v)+1))
                v[s+off]=(v[s+off]+q*co)%p
    poww.append(np.array(v[:5],dtype=np.int64)%p)
def mulmat(a):
    # column j = a*w^j
    M=np.zeros((5,5),dtype=np.int64)
    for j in range(5):
        v=np.zeros(5,dtype=np.int64)
        for i,ai in enumerate(a):
            if ai: v=(v+ai*poww[i+j])%p
        M[:,j]=v
    return M
A=np.zeros((len(rows)*5,len(cols)*5),dtype=np.int16)
for cj,(i,a,b,c) in enumerate(cols):
    for e,co in S[i-1].items():
        rr=ri[(a+e[0],b+e[1],c+e[2])]
        A[rr*5:(rr+1)*5,cj*5:(cj+1)*5]=mulmat(coeff_vec(co)).astype(np.int16)
# target b for h monomial, coefficient 1
bvec=np.zeros(len(rows)*5,dtype=np.int16); bvec[ri[(1,0,0)]*5]=1
np.savez_compressed('fixed_matrix_p71.npz',A=A,b=bvec,rows=np.array(rows,dtype=np.int16),cols=np.array(cols,dtype=np.int16))
print('shape',A.shape,'nnz',np.count_nonzero(A),'targetrow',ri[(1,0,0)])
