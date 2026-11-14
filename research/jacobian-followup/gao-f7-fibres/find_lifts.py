"""Bounded F13 preflight for F7 multiple-root configurations, not exclusions.

Positive points are only candidate characteristic-zero certificates until
their exact multiplicities and Jacobian ranks are separately checked.
"""
import itertools
import json
from pathlib import Path
import time
from flint import nmod_poly

HERE=Path(__file__).resolve().parent
start=time.monotonic()
found={}
rounds=[]


def pivots(matrix,p):
    mat=[r[:] for r in matrix];out=[];row=0
    for col in range(len(mat[0])):
        pivot=next((i for i in range(row,len(mat)) if mat[i][col]%p),None)
        if pivot is None:continue
        mat[row],mat[pivot]=mat[pivot],mat[row]
        inv=pow(mat[row][col]%p,-1,p)
        mat[row]=[(v*inv)%p for v in mat[row]]
        for i in range(len(mat)):
            if i!=row:
                d=mat[i][col]%p
                mat[i]=[(a-d*b)%p for a,b in zip(mat[i],mat[row])]
        out.append(col);row+=1
        if row==len(mat):break
    return out


for p in [13,17,19,23]:
 hist={};tested=0
 for Y1,Y2,Y3,Y4 in itertools.product(range(p),repeat=4):
    a=nmod_poly([Y3,-2*Y2,Y1-Y4,-1,3*Y2,-2*Y1,1],p)
    b=nmod_poly([Y4,0,-3*Y2,2*Y1,-1],p)
    R=a*a+a+b*nmod_poly([0,1],p)+nmod_poly([-Y2,Y1,0,1],p)
    gcd=R.gcd(R.derivative())
    dg=gcd.degree()
    tested+=1
    dd=gcd.gcd(R.derivative().derivative()).degree()
    count=12-2*dg+dd
    hist[str(count)]=hist.get(str(count),0)+1
    if count not in [4,5,6] or str(count) in found:
        continue
    groups={}
    for f,m in R.factor()[1]:groups[m]=groups.get(m,nmod_poly([1],p))*f
    if sum(f.degree()*(m-1) for m,f in groups.items())>4:continue
    tt=nmod_poly([0,1],p)
    da=[tt**2-2*tt**5,-2*tt+3*tt**4,nmod_poly([1],p),-tt**2]
    db=[2*tt**3,-3*tt**2,nmod_poly([0],p),nmod_poly([1],p)]
    cols=[2*a*da[i]+da[i]+tt*db[i]+[tt,nmod_poly([-1],p),nmod_poly([0],p),nmod_poly([0],p)][i] for i in range(4)]
    for m,f in sorted(groups.items()):
        other=nmod_poly([1],p)
        for n,h in groups.items():other*=h**(n-1 if n==m else n)
        for j in range(f.degree()):cols.append(-m*tt**j*other)
    matrix=[[int(col[i]) for col in cols] for i in range(12)]
    pivot=pivots(matrix,p)
    if len(pivot)!=12:continue
    found[str(count)]={'modulus':p,'Y':[Y1,Y2,Y3,Y4],
                      'multiplicity_factors':{str(m):[int(f[i]) for i in range(f.degree()+1)] for m,f in sorted(groups.items())},
                      'jacobian_mod_p':matrix,'pivot_columns':pivot}
    print('candidate',count,found[str(count)],flush=True)
    if len(found)==3:break
 rounds.append({'modulus':p,'targets_tested':tested,'total_target_box':p**4,
                'exhaustive':tested==p**4,'simple_root_histogram':hist})
 if len(found)==3:break
receipt={'rounds':rounds,'wall_seconds':time.monotonic()-start,
         'positive_candidates':found,
         'supersedes':'lift-preflight-v1.json has an incorrect histogram for gcd degree2; its raw positive points remain finite-field candidates only. The corrected preflight always computes the second gcd degree.',
         'boundary':'Full-rank finite-field candidates. A separate checker and written Hensel argument are required for characteristic-zero claims. Failed finite boxes are not nonexistence certificates.'}
(HERE/'lift-preflight.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
