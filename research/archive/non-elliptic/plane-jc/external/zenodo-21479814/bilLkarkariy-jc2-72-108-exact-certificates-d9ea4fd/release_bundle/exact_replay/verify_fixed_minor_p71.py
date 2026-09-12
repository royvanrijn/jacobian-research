import numpy as np,time
p=71
z=np.load('fixed_matrix_p71.npz');Aall=z['A'].astype(np.int64);ball=z['b'].astype(np.int64)
sel=np.load('pivot_scalar_rows_p71.npy')
A=Aall[sel].copy();b=ball[sel].copy();n=A.shape[0]
perm=np.arange(n);det=1;sign=1;t=time.time()
# forward elimination
for c in range(n):
    nz=np.flatnonzero(A[c:,c])
    if not nz.size: raise SystemExit(f'singular at {c}')
    q=c+int(nz[0])
    if q!=c:
        A[[c,q]]=A[[q,c]];b[[c,q]]=b[[q,c]];sign=-sign
    piv=int(A[c,c]);det=det*piv%p
    inv=pow(piv,-1,p)
    A[c,c:]=A[c,c:]*inv%p;b[c]=b[c]*inv%p
    idx=np.flatnonzero(A[c+1:,c])+c+1
    if idx.size:
        fac=A[idx,c].copy()
        for s in range(0,len(idx),128):
            ii=idx[s:s+128];ff=fac[s:s+128]
            A[ii,c:]=(A[ii,c:]-ff[:,None]*A[c,c:])%p
            b[ii]=(b[ii]-ff*b[c])%p
    if (c+1)%400==0:print('elim',c+1,round(time.time()-t,1),flush=True)
if sign<0:det=(-det)%p
# back substitute (diag 1)
x=np.zeros(n,dtype=np.int64)
for i in range(n-1,-1,-1):
    x[i]=(b[i]-int(np.dot(A[i,i+1:],x[i+1:])%p))%p
res=(Aall@x-ball)%p
print('det_mod_71',det)
print('selected_residual_nonzero',np.count_nonzero((Aall[sel]@x-ball[sel])%p))
print('all_residual_nonzero',np.count_nonzero(res),'max',res.max())
np.save('fixed_solution_p71.npy',x.astype(np.int16))
open('fixed_minor_p71_meta.txt','w').write(f'prime=71\nrows={Aall.shape[0]}\ncols={Aall.shape[1]}\nminor_size={n}\ndet_mod_71={det}\nall_residual_nonzero={np.count_nonzero(res)}\n')
