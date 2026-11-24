#!/usr/bin/env python3
"""Exact controls for a positive-measure merger obstruction, not an SFC counterexample."""
if not __debug__:
    raise RuntimeError("Assertions are required; do not run with -O.")

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, factorial, prod
from pathlib import Path
import json
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
EXPECTED_SOURCE_SHA256 = {
    "SFC_ROW_MERGING_OBSTRUCTIONS.md": "731e5dcd27f8cfe6f688c351b55888facffe0756c797c192f3f1970a7113edee",
    "SFC_BALANCED_PRIVATE_AMPLIFICATION.md": "f4ebaa9028cb04258cf65be178685650d53e2852644ed0678f1a655d134c065e",
    "SFC_CYCLIC_POLYNOMIAL_TILTS.md": "ba7cc5fe27558993b93eb87d4a08e2cb29adb1c30d86ab8868b20cbfc5cfde0d",
    "UNRESTRICTED_SFC_MECHANISMS.tex": "f8b3dfc5947601fb002aee0ba942b0fc20b7ef46ca8e0d4f69b5a3e5e1361785"
}


def trim(p):
    p=list(p)
    while len(p)>1 and not p[-1]:p.pop()
    return p


def add(a,b):
    out=[F(0)]*max(len(a),len(b))
    for i,c in enumerate(a):out[i]+=c
    for i,c in enumerate(b):out[i]+=c
    return trim(out)


def mul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,c in enumerate(a):
        for j,d in enumerate(b):out[i+j]+=c*d
    return trim(out)


def divide(a,b):
    a,b=trim(a),trim(b)
    assert any(b)
    q=[F(0)]*max(1,len(a)-len(b)+1)
    while any(a) and len(a)>=len(b):
        k,c=len(a)-len(b),a[-1]/b[-1]
        q[k]+=c
        for j,d in enumerate(b):a[k+j]-=c*d
        a=trim(a)
    return trim(q),a


def det(A):
    A=[list(map(F,row)) for row in A]
    value=F(1)
    for i in range(len(A)):
        pivot=next((j for j in range(i,len(A)) if A[j][i]),None)
        if pivot is None:return F(0)
        if pivot != i:A[i],A[pivot]=A[pivot],A[i];value=-value
        value*=A[i][i]
        for j in range(i+1,len(A)):
            factor=A[j][i]/A[i][i]
            for k in range(i,len(A)):A[j][k]-=factor*A[i][k]
    return value


def resultant(a,b):
    a,b=list(reversed(trim(a))),list(reversed(trim(b)))
    m,n=len(a)-1,len(b)-1
    rows=[]
    for i in range(n):rows.append([0]*i+a+[0]*(n-1-i))
    for i in range(m):rows.append([0]*i+b+[0]*(m-1-i))
    return det(rows)


def power_linear(a,b,n):
    return [F(comb(n,j))*a**(n-j)*b**j for j in range(n+1)]


def jacobi_from_rodrigues(p,q):
    assert p<=q
    raw=mul(power_linear(1,-1,p),power_linear(1,1,q))
    derivative=[raw[k]*factorial(k)//factorial(k-p) for k in range(p,len(raw))]
    quotient,remainder=divide(derivative,power_linear(1,1,q-p))
    assert not any(remainder)
    J=[x*F((-1)**p,2**p*factorial(p)) for x in quotient]
    transformed=[F(0)]
    for j,c in enumerate(J):
        term=mul(power_linear(1,1,j),power_linear(1,-1,p-j))
        transformed=add(transformed,[c*x for x in term])
    expected=[F(comb(p,j)*comb(q,j)) for j in range(p+1)]
    assert transformed==expected
    assert sum(expected)==comb(p+q,p)
    return transformed


def scalar_poly_mul(a,b):
    out={}
    for ea,ca in a.items():
        for eb,cb in b.items():
            e=tuple(x+y for x,y in zip(ea,eb))
            out[e]=out.get(e,F(0))+ca*cb
    return {e:c for e,c in out.items() if c}


ATOMS=((F(1),F(2),F(121,296)),(F(1),F(-2),F(121,296)),
       (F(11,3),F(0),F(27,148)))
DELTA=F(1,100)


def atom_moment(a,b):
    return sum(weight*A**a*B**b for A,B,weight in ATOMS)


def disk_moment(a,b):
    if a%2 or b%2:return F(0)
    r,s=a//2,b//2
    return DELTA**(a+b)*F(factorial(a)*factorial(b),
                         4**(r+s)*factorial(r)*factorial(s)*factorial(r+s+1))


def moment_value(poly):
    return sum(c*atom_moment(e[0],e[1])*disk_moment(e[2],e[3])
               for e,c in poly.items())


U={(1,0,0,0):F(1,2),(0,1,0,0):F(1,6),
   (0,0,1,0):F(1,2),(0,0,0,1):F(1,6)}
V={e:(-c if e[1] or e[3] else c) for e,c in U.items()}


@lru_cache(None)
def uv_moment(a,b):
    P={(0,0,0,0):F(1)}
    for _ in range(a):P=scalar_poly_mul(P,U)
    for _ in range(b):P=scalar_poly_mul(P,V)
    return moment_value(P)


def gi_add(a,b):return (a[0]+b[0],a[1]+b[1])


def gi_mul(a,b):return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])


def gi_power(a,n):
    out=(F(1),F(0))
    for _ in range(n):out=gi_mul(out,a)
    return out


def binary_gi(coeffs,x,y):
    m=len(coeffs)-1
    out=(F(0),F(0))
    for k,c in enumerate(coeffs):
        term=gi_mul(gi_power(x,k),gi_power(y,m-k))
        out=gi_add(out,(c*term[0],c*term[1]))
    return out


def verify_sources():
    found = {}
    for name, expected in EXPECTED_SOURCE_SHA256.items():
        actual = sha256((HERE/name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError("Source SHA256 mismatch: "+name)
        found[name] = actual
    return found


def verify_rejection_controls():
    optimized = subprocess.run([sys.executable, "-O", str(Path(__file__).resolve())],
                               capture_output=True, text=True)
    assert optimized.returncode != 0 and "do not run with -O" in optimized.stderr
    with tempfile.TemporaryDirectory(prefix="sfc_merger_obstruction_guard_") as directory:
        temporary = Path(directory)
        script = temporary/Path(__file__).name
        script.write_bytes(Path(__file__).read_bytes())
        for name in EXPECTED_SOURCE_SHA256:
            (temporary/name).write_bytes((HERE/name).read_bytes())
        note = temporary/"SFC_ROW_MERGING_OBSTRUCTIONS.md"
        note.write_bytes(note.read_bytes()+b"\nAltered-note rejection control.\n")
        altered = subprocess.run([sys.executable,str(script)], capture_output=True, text=True)
        assert altered.returncode != 0
        assert "Source SHA256 mismatch: SFC_ROW_MERGING_OBSTRUCTIONS.md" in altered.stderr
    return dict(optimized_mode_rejected=True, altered_note_rejected=True)


def main():
    sources=verify_sources()
    rejection_controls=verify_rejection_controls()
    beta_count,jacobi_count=0,0
    for P in range(9):
        for Q in range(9):
            beta=sum(F((-1)**j*comb(Q,j),P+j+1) for j in range(Q+1))
            assert beta==F(factorial(P)*factorial(Q),factorial(P+Q+1))
            assert beta*(P+Q+1)*comb(P+Q,P)==1
            beta_count+=1
            jacobi_from_rodrigues(min(P,Q),max(P,Q))
            jacobi_count+=1
    forward_square=comb(2,2)+comb(2,0)-2*comb(2,1)
    assert forward_square==-2
    assert det([[1,1,1],[1,1,2],[1,2,1]])==-1
    assert sum(weight for _,_,weight in ATOMS)==1
    assert atom_moment(1,0)==F(55,37)
    assert atom_moment(2,0)==atom_moment(0,2)==F(121,37)
    assert atom_moment(3,0)==3*atom_moment(1,2)==F(363,37)
    assert min(A-abs(B)/3 for A,B,_ in ATOMS)-F(4,3)*DELTA==F(8,25)
    v=DELTA**2/4
    H2=F(121,37)+v
    H3=F(363,37)+F(165,37)*v
    source2=[2*uv_moment(k,2-k) for k in range(3)]
    source3=[6*uv_moment(k,3-k) for k in range(4)]
    assert source2==[H2*c/9 for c in [5,4,5]]
    assert source3==[H3*c/18 for c in [15,13,13,15]]
    merged2=[comb(2,k)*c for k,c in enumerate(source2)]
    merged3=[comb(3,k)*c for k,c in enumerate(source3)]
    assert merged2==[H2*c/9 for c in [5,8,5]]
    assert merged3==[H3*c/18 for c in [15,39,39,15]]
    x,y=(F(1),F(3)),(F(1),F(-3))
    assert binary_gi(merged2,x,y)==binary_gi(merged3,x,y)==(0,0)
    assert resultant([5,4,5],[15,13,13,15])==29400
    rho_records=[]
    for rho in map(F,[-1,0,F(1,2),1,F(3,2),2]):
        S=[F(5),4*(1+rho),F(5)]
        T=[F(15),13*(1+2*rho),13*(1+2*rho),F(15)]
        value=resultant(S,T)
        assert value==9800*(3-2*rho)*(rho-1)**2
        rho_records.append(dict(rho=str(rho),normalized_resultant=str(value)))
    holder_count,suppression_count,gcd_count=0,0,0
    polynomials=[]
    for m in range(1,7):
        Hm=uv_moment(m,0)
        assert Hm==uv_moment(0,m)>0
        coefficient=[]
        for k in range(m+1):
            h=uv_moment(k,m-k)
            assert 0<h<=Hm
            holder_count+=1
            C=comb(m,k)
            raw=C*factorial(12*k)*factorial(12*(m-k))*h
            normalized=raw/(Hm*factorial(12*m))
            assert normalized<=F(C)**(-11)
            suppression_count+=1
            coefficient.append(normalized)
        assert coefficient[0]==coefficient[-1]==1
        polynomials.append(coefficient)
        # Both merged moments vanish after every common private power.
        if m in (2,3):
            common=[F(comb(m,k))*uv_moment(k,m-k)*factorial(12*m)
                    for k in range(m+1)]
            assert binary_gi(common,x,y)==(0,0)
    for a,b in zip(polynomials,polynomials[1:]):
        while any(b):a,b=b,divide(a,b)[1]
        assert len(trim(a))==1 and a[0]
        gcd_count+=1
    assert F(72*3**3,2**11)==F(243,256)<1
    assert F(288*3,2**11)==F(27,64)<1
    for D in (1,12,13):
        source=[(D,0,1,0),(0,D,0,1)]
        target=[(D,1,0),(D,0,1)]
        assert len(set(source))==len(set(target))==2
    receipt=dict(status="PASS",scope="positive-measure row-merger obstruction, not ordinary factorial SFC",
                 note_sha256=sha256((HERE/'SFC_ROW_MERGING_OBSTRUCTIONS.md').read_bytes()).hexdigest(),
                 checker_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                 source_sha256=sources,rejection_controls=rejection_controls,
                 beta_identities=beta_count,rodrigues_jacobi_identities=jacobi_count,
                 forward_mellin_kernel_square=-2,
                 atoms=[dict(A=str(A),B=str(B),probability=str(p)) for A,B,p in ATOMS],
                 noise_radius=str(DELTA),strict_positive_feature_lower_bound="4/25",
                 source_normalized_forms=[[5,4,5],[15,13,13,15]],
                 target_normalized_forms=[[5,8,5],[15,39,39,15]],
                 source_normalized_resultant=29400,
                 common_projective_zero=["1+3i","1-3i"],
                 physical_correlation_resultant_controls=rho_records,
                 holder_controls=holder_count,private_D12_coefficient_controls=suppression_count,
                 private_D12_consecutive_gcd_controls=gcd_count,
                 all_order_amplification_input="written extension of SFC_BALANCED_PRIVATE_AMPLIFICATION.md; finite checks alone do not prove it",
                 ordinary_factorial_counterexample=False)
    output=HERE/'sfc_row_merging_obstructions_verification.json'
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status="PASS",beta_identities=beta_count,jacobi_identities=jacobi_count,
                          source_resultant=29400,holder_controls=holder_count,
                          D12_gcd_controls=gcd_count,receipt=str(output)),indent=2))


if __name__=='__main__':main()
