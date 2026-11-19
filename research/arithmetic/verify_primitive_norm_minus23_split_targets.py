#!/usr/bin/env python3
"""Exact finite support for the -23 split-target coordinate-gcd operator lemma."""
if not __debug__:
    raise RuntimeError('This certificate requires Python assertions enabled.')
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import gcd,isqrt
import json
import argparse
import hashlib
from pathlib import Path


def mmul(a,b):
    return (a[0]*b[0]+a[1]*b[2],a[0]*b[1]+a[1]*b[3],
            a[2]*b[0]+a[3]*b[2],a[2]*b[1]+a[3]*b[3])


def local_rank_counts(p):
    bydet=Counter();bad=Counter();rank1=axes=0
    for a,b,c,d in product(range(p),repeat=4):
        det=(a*d-b*c)%p
        if det:
            bydet[det]+=1
            if ((a==d==0) or (b==c==0)):bad[det]+=1
        elif a or b or c or d:
            rank1+=1
            if sum(bool(v) for v in(a,b,c,d))==1:axes+=1
    assert rank1==(p-1)*(p+1)**2
    assert axes==4*(p-1)
    assert F(axes,rank1)==F(4,(p+1)**2)
    for det in range(1,p):
        assert bydet[det]==p*(p*p-1)
        assert F(bad[det],bydet[det])==F(2,p*(p+1))
    return {'p':p,'rank1':rank1,'rank1_bad':axes,
            'rank1_bad_density':str(F(axes,rank1)),
            'unit_bad_density':str(F(2,p*(p+1)))}


def primitive_path_cosets(p,maximum):
    result=[]
    for t in range(1,maximum+1):
        full=chosen=0
        for a in range(t+1):
            c=t-a
            for s in range(p**a):
                M=(p**a,s,0,p**c)
                if all(v%p==0 for v in M):continue
                full+=1
                for z in list(range(p))+[None]:
                    h=(1,0,z,1) if z is not None else(0,1,1,0)
                    g=mmul(M,h)
                    # Image=kernel=e1: precisely the upper-right matrix axis.
                    if g[0]%p==g[2]%p==g[3]%p==0 and g[1]%p:
                        chosen+=1
        assert full==(p+1)*p**(t-1)
        assert chosen==p**(t-1)
        result.append({'p':p,'t':t,'full_primitive_degree':full,
                       'chosen_endpoint_degree':chosen})
    return result


def norm_vectors(bound):
    rows=[[]for _ in range(bound+1)]
    ymax=isqrt(4*bound//23)
    for y in range(-ymax,ymax+1):
        root=isqrt(4*bound-23*y*y)
        lo=-((root+y)//2);hi=(root-y)//2
        for x in range(lo,hi+1):
            n=x*x+x*y+6*y*y
            if 0<=n<=bound:rows[n].append((x,y))
    return rows


def signed_global_controls(bound=256):
    vectors=norm_vectors(bound)
    result=[]
    for p in(2,3,13,59):
        roots=[r for r in range(p) if (r*r-r+6)%p==0]
        assert len(roots)==2
        r,s=roots
        cases=0;positive=0
        for N in range(p,bound+1,p):
            count=chosen=0
            primitive=0
            for n in range(N+1):
                for x,y in vectors[n]:
                    for z,w in vectors[N-n]:
                        if x%p==y%p==z%p==w%p==0:continue
                        primitive+=1
                        matrix=((x+r*y)%p,(z+r*w)%p,
                                -(z+s*w)%p,(x+s*y)%p)
                        bad=((x%p==y%p==0)or(z%p==w%p==0))
                        if bad:
                            count+=1
                            assert sum(bool(v)for v in matrix)==1
                        if matrix[0]==matrix[2]==matrix[3]==0 and matrix[1]:
                            chosen+=1
            assert count==4*chosen
            assert chosen%2==0 and count%8==0
            # The primitive whole-four-vector restriction is independently a
            # subtraction of the all-divisible shell; no split anisotropy used.
            whole=sum(len(vectors[n])*len(vectors[N-n])for n in range(N+1))
            removed=0
            if N%(p*p)==0:
                M=N//(p*p)
                removed=sum(len(vectors[n])*len(vectors[M-n])for n in range(M+1))
            assert primitive==whole-removed
            cases+=1;positive+=bool(count)
        result.append({'p':p,'target_cases':cases,'positive_rank1_bad_cases':positive,
                       'axis_case_factor':4,'edge_stabilizer_units':2,
                       'total_return_factor':8})
    assert len(vectors[4])==2
    assert all(gcd(x,y)>1 for x,y in vectors[4])
    return result


def higher_digit_lifts(p,h):
    modulus=p**h
    targets=[p**t for t in range(1,h)]
    data={N:[0,0]for N in targets}
    for a,b,c,d in product(range(modulus),repeat=4):
        if a%p==b%p==c%p==d%p==0:continue
        N=(a*d-b*c)%modulus
        if N not in data:continue
        data[N][0]+=1
        if ((a%p==d%p==0)or(b%p==c%p==0)):
            data[N][1]+=1
    for N,(total,bad)in data.items():
        assert total==(p-1)*(p+1)**2*p**(3*(h-1))
        assert bad==4*(p-1)*p**(3*(h-1))
        assert F(bad,total)==F(4,(p+1)**2)
    return {'p':p,'h':h,'targets':data}


def fmul(a,b):
    return(a[0]*b[0]-6*a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1])


def fconj(a):return(a[0]+a[1],-a[1])


def qmul(a,b):
    x=fmul(a[:2],b[:2]);y=fmul(a[2:],fconj(b[2:]))
    z=fmul(a[:2],b[2:]);w=fmul(a[2:],fconj(b[:2]))
    return(x[0]-y[0],x[1]-y[1],z[0]+w[0],z[1]+w[1])


def in_ideal(q,p,s):
    x,y,z,w=q
    return(-s*x-z-w)%p==0 and(-s*z+x+y)%p==0


def ideal_basis(p,s):
    return[(1,-1,0,-s),(0,s,1,-1),(0,p,0,0),(0,0,0,p)]


def left_order_norm_elements(p,s,k):
    basis=ideal_basis(p,s);bound=k*p*p
    rows=norm_vectors(bound);found=[]
    for n in range(bound+1):
        for a in rows[n]:
            for b in rows[bound-n]:
                q=a+b
                for v in basis:
                    product=qmul(q,v)
                    if any(x%p for x in product):break
                    if not in_ideal(tuple(x//p for x in product),p,s):break
                else:found.append(q)
    return found


def multiply(A,B):
    return[[sum(a*b for a,b in zip(row,col))for col in zip(*B)]for row in A]


def identity(n):return[[int(i==j)for j in range(n)]for i in range(n)]


def complete_bad_count_two(N,rows):
    count=0
    for n in range(N+1):
        for x,y in rows[n]:
            for z,w in rows[N-n]:
                if x%2==y%2==z%2==w%2==0:continue
                if x%2==y%2==0 or z%2==w%2==0:count+=1
    return count


def independent_brandt_controls():
    I2units=left_order_norm_elements(2,0,1)
    I2norm2=left_order_norm_elements(2,0,2)
    I3units=left_order_norm_elements(3,1,1)
    I3norm2=left_order_norm_elements(3,1,2)
    assert len(I2units)==2 and len(I3units)==6
    assert len(I2norm2)==2 and len(I3norm2)==0
    assert F(1,2)+F(1,3)+1==F(11,6)
    # Class order is principal(units4), I3(units6), I2(units2).
    weights=[F(3,11),F(2,11),F(6,11)]
    coarse=[[1,0,2],[0,0,3],[1,1,1]]
    assert all(sum(row)==3 for row in coarse)
    assert all(weights[i]*coarse[i][j]==weights[j]*coarse[j][i]
               for i in range(3)for j in range(3))
    assert [coarse[i][i]for i in range(3)]==[1,len(I3norm2)//6,len(I2norm2)//2]
    # All integer possibilities with these diagonals/degrees/masses: unique.
    possible=[]
    from itertools import product
    for x,y,z in product(range(4),repeat=3):
        A=[[1,x,2-x],[y,0,3-y],[z,2-z,1]]
        if all(weights[i]*A[i][j]==weights[j]*A[j][i]
               for i in range(3)for j in range(3)):
            possible.append(A)
    assert possible==[coarse]
    edges=[(0,0),(0,2),(1,2),(2,0),(2,1),(2,2)]
    reverse=[0,3,4,1,2,5]
    edgeweights=[F(x,11)for x in(1,2,2,2,2,2)]
    flip=[[int(j==reverse[i])for j in range(6)]for i in range(6)]
    B=[]
    for i,(u,v)in enumerate(edges):
        row=[]
        for j,(s,t)in enumerate(edges):
            row.append((coarse[v][t]if s==v else 0)-int(j==reverse[i]))
        assert min(row)>=0 and sum(row)==2
        B.append(row)
    FB_F=multiply(multiply(flip,B),flip)
    assert all(edgeweights[i]*FB_F[i][j]==edgeweights[j]*B[j][i]
               for i in range(6)for j in range(6))
    rows=norm_vectors(4096)
    power=identity(6);controls=[]
    for t in range(1,13):
        V=multiply(power,flip)
        assert all(sum(row)==2**(t-1)for row in V)
        expected=8*V[1][1]
        actual=complete_bad_count_two(2**t,rows)
        assert expected==actual,(t,expected,actual)
        controls.append({'t':t,'N':2**t,'brandt_return':V[1][1],
                         'predicted_bad_count':expected,'coordinate_bad_count':actual})
        power=multiply(power,B)
    result={'I2_unit_numerators_over2':I2units,
            'I2_norm2_numerators_over2':I2norm2,
            'I3_unit_numerators_over3':I3units,
            'I3_norm2_numerators_over3':I3norm2,
            'projective_mass':'11/6','complete_class_unit_counts':[4,6,2],
            'coarse_brandt_T2':coarse,'oriented_edge_classes':edges,
            'edge_weights':[str(x)for x in edgeweights],
            'flip_matrix':flip,'nonbacktracking_matrix':B,
            'rank1_bad_return_controls':controls,
            'scope':'Exact left-order units and norm2 elements, mass-complete maximal class list, unique coarse T2, induced six-state Iwahori flip/nonbacktracking matrices and actual rank1 coordinate count comparison at N=2^t through4096.'}


    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("primitive_norm_minus23_split_target_checks.json"))
    args=parser.parse_args()
    result={'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'independent_brandt_two':independent_brandt_controls(),
            'local_rank_counts':[local_rank_counts(p)for p in(2,3,13)],
            'path_cosets':[r for p,m in[(2,6),(3,5),(13,3)]
                           for r in primitive_path_cosets(p,m)],
            'higher_lifts':[higher_digit_lifts(2,3),higher_digit_lifts(3,2)],
            'signed_global_controls':signed_global_controls(),
            'scope':'Exact finite local densities and lifts, primitive HNF path degrees, signed shell controls, and a mass-complete three-class Brandt computation with actual Iwahori returns at p=2 through N=4096. These computations do not prove the spectral norm bound or unrestricted global correspondence.'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
