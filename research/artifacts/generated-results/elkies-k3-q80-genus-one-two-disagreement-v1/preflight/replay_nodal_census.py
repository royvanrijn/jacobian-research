"""Independent vectorized census; reconstruct square roots from constant terms."""
import json
from pathlib import Path
import resource
import time
import numpy as np

resource.setrlimit(resource.RLIMIT_CPU,(40,45))
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
P=131
def census(A,B,root):
    z0=np.tile(np.arange(P,dtype=np.int64),P)
    z1=np.repeat(np.arange(P,dtype=np.int64),P)
    inv=np.array([0]+[pow(i,-1,P) for i in range(1,P)],dtype=np.int64)
    found=[]
    for z2 in range(P):
        x=np.zeros((5,P*P),dtype=np.int64)
        x[0]=(root[0]+88*z0)%P;x[1]=(root[1]+62*z0+88*z1)%P
        x[2]=(z0+62*z1+88*z2)%P;x[3]=(z1+62*z2)%P;x[4]=z2
        square=np.zeros((9,P*P),dtype=np.int64)
        for i in range(5):
            for j in range(5):square[i+j]+=x[i]*x[j]
        square%=P
        f=np.zeros((13,P*P),dtype=np.int64)
        for i in range(9):
            for j in range(5):f[i+j]+=(square[i]+A[i])*x[j]
        for i in range(13):f[i]+=B[i]
        f%=P
        s=np.zeros((11,P*P),dtype=np.int64)
        for i in range(12,1,-1):
            s[i-2]=f[i];f[i-1]=(f[i-1]-62*f[i])%P;f[i-2]=(f[i-2]-88*f[i])%P
        if np.any(f[:2]):raise ValueError('nodal-root divisibility')
        # The actual fibre at t=0 has no rational2-torsion, so this is a unit.
        if np.any(s[0]==0):raise ValueError('nonzero constant term required')
        normalized=s*inv[s[0]]%P
        y=np.zeros((6,P*P),dtype=np.int64);y[0]=1
        for i in range(1,6):
            coefficient=np.zeros(P*P,dtype=np.int64)
            for j in range(1,i):coefficient+=y[j]*y[i-j]
            y[i]=(normalized[i]-coefficient)*66%P
        good=np.ones(P*P,dtype=bool)
        for i in range(6,11):
            coefficient=np.zeros(P*P,dtype=np.int64)
            for j in range(max(0,i-5),6):coefficient+=y[j]*y[i-j]
            good &= coefficient%P==normalized[i]
        for k in np.flatnonzero(good):
            found.append({'z':[int(z0[k]),int(z1[k]),z2],
                          'x':[int(v) for v in x[:,k]],
                          'constant_scalar':int(s[0,k]),
                          'constant_normalized_root':[int(v) for v in y[:,k]]})
    return found
if __name__=='__main__':
    out=Path(__file__).resolve().parent.parent
    data=json.loads((out/'integral-preview.json').read_text());start=time.process_time()
    rows=census(data['A'],data['B'],data['simple_root'])
    result={'status':'PREVIEW_ONLY','tried':P**3,'records':rows,'cpu_seconds':time.process_time()-start,
            'method':'constant-term square-root reconstruction, independent of leading-term C++ producer'}
    with (out/'nodal-independent-preview.json').open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps(result,sort_keys=True),flush=True)
