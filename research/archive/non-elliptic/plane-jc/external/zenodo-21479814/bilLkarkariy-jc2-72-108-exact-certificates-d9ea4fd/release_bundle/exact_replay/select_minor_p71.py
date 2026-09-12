import numpy as np,time,sys
p=71
z=np.load('fixed_matrix_p71.npz')
A=z['A'].astype(np.int64)
orig=np.arange(A.shape[0],dtype=np.int32)
r=0;t=time.time();nrows,ncols=A.shape
for c in range(ncols):
    nz=np.flatnonzero(A[r:,c])
    if nz.size==0:
        print('NO PIVOT col',c,'rank',r);break
    q=r+int(nz[0])
    if q!=r:
        A[[r,q]]=A[[q,r]];orig[[r,q]]=orig[[q,r]]
    inv=pow(int(A[r,c]),-1,p)
    A[r,c:]=(A[r,c:]*inv)%p
    idx=np.flatnonzero(A[r+1:,c])+r+1
    if idx.size:
        fac=A[idx,c].copy()
        # block rows to cap temp memory
        for s in range(0,len(idx),128):
            ii=idx[s:s+128];ff=fac[s:s+128]
            A[ii,c:]=(A[ii,c:]-ff[:,None]*A[r,c:])%p
    r+=1
    if (c+1)%100==0:
        print('pivot',c+1,'elapsed',round(time.time()-t,1),flush=True)
else:
    print('FULL',r,'elapsed',round(time.time()-t,1))
sel=orig[:r]
np.save('pivot_scalar_rows_p71.npy',sel)
# also save corresponding extension-row and component
np.savetxt('pivot_scalar_rows_p71.txt',np.c_[sel,sel//5,sel%5],fmt='%d',header='scalar_row extension_row component')
print('saved',len(sel),'unique ext rows',len(set((sel//5).tolist())))
