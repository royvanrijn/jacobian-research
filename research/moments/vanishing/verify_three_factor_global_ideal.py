#!/usr/bin/env python3
"""Exact invariant-ring proof data for global three-factor Gaussian moment ideals."""
from pathlib import Path
from itertools import product
from math import gcd
from functools import reduce
import hashlib,json,time
import sympy as S
from sympy.polys.matrices import DomainMatrix

HERE=Path(__file__).resolve().parent
INPUT='verify_nodal_product_reduction.py'
EXPECTED_INPUT_SHA256='7254eff7dad198b2e15fa727da0703bab5f8e3b8704697b27a0ed212527821d6'
assert hashlib.sha256((HERE/INPUT).read_bytes()).hexdigest()==EXPECTED_INPUT_SHA256
from verify_nodal_product_reduction import moments,pair_count_moment


def weight_monomials(weights):
    if min(weights)<0:return []
    out=[]
    for d in range(min(weights[0],weights[1])+1):
        for e in range(min(weights[0]-d,weights[2])+1):
            for f in range(min(weights[1]-d,weights[2]-e)+1):
                left=[weights[0]-d-e,weights[1]-d-f,weights[2]-e-f]
                if min(left)<0 or any(v%2 for v in left):continue
                out.append(tuple(v//2 for v in left)+(d,e,f))
    return out


def standard_basis(G):
    lm=[tuple(p.LM(order=G.order).exponents) for p in G.polys]
    bounds=[min(v[i] for v in lm if v[i] and sum(v)==v[i]) for i in range(len(G.gens))]
    base=[v for v in product(*(range(k) for k in bounds)) if not any(all(x>=y for x,y in zip(v,w)) for w in lm)]
    hf=[sum(sum(v)==k for v in base) for k in range(max(map(sum,base))+1)]
    return base,hf


def symmetric_moments(a,b,c,A,E,s,q):
    normalized=moments([[a,1,1],[1,b,1],[1,1,c]],10);answer=[]
    for m in range(1,6):
        phi,rem,mapping=S.symmetrize(normalized[2*m-1],[a,b,c],formal=True);assert rem==0
        phi=phi.subs({v:k for (v,_),k in zip(mapping,[s,q,A])})
        terms=S.Poly(phi,s,q,A).terms()
        assert all(m-i-2*j-k>=0 for (i,j,k),coef in terms)
        answer.append(S.expand(sum(coef*s**i*q**j*A**k*E**(m-i-2*j-k) for (i,j,k),coef in terms)))
    return answer


def semigroup_checks():
    matrix=S.Matrix([[1,1,0,0,0],[1,0,1,0,0],[1,0,0,1,0],[0,0,0,2,1],[0,0,2,0,1],[0,2,0,0,1]])
    assert matrix.rank()==4 and matrix*S.Matrix([-1,1,1,1,-2])==S.zeros(6,1)
    checked=0
    for m in range(7):
        for exp in weight_monomials((2*m,)*3):
            a,b,c,d,e,f=exp;odd=d%2;assert e%2==f%2==odd
            n=m-odd;p,q,r=(d-odd)//2,(e-odd)//2,(f-odd)//2
            k=max(0,p+q+r-n);assert k<=min(p,q,r)
            powers=(n-p-q-r+k,r-k,q-k,p-k,odd+2*k)
            assert min(powers)>=0 and sum(powers)==m
            assert matrix*S.Matrix(powers)==S.Matrix(exp)
            checked+=1
    return {'generator_exponent_matrix':[[int(v) for v in row] for row in matrix.tolist()],'integer_kernel':[-1,1,1,1,-2],
            'relation':'B*C*D=A*E^2','balanced_monomials_replayed_through_degree_6':checked,
            'proof_boundary':'The nonnegative factorization formula in the note proves every degree; finite replay checks its transcription.'}


def dual_nonmembership(coords,polys,M10,label,multiplier,delta):
    target=tuple(10+k for k in delta);base=weight_monomials(target);lookup={v:i for i,v in enumerate(base)};cols=[];labels=[]
    for j,poly in enumerate(polys,1):
        for exponent in weight_monomials(tuple(v-2*j for v in target)):
            col=[0]*len(base)
            for exp,coefficient in poly.terms():col[lookup[tuple(x+y for x,y in zip(exp,exponent))]]+=int(coefficient)
            cols.append(col);labels.append([2*j,list(exponent)])
    matrix=S.Matrix.hstack(*(S.Matrix(v) for v in cols))
    rhs=S.Poly(multiplier*M10.as_expr(),*coords);vec=S.Matrix([rhs.coeff_monomial(v) for v in base])
    null=DomainMatrix.from_Matrix(matrix.T).nullspace().to_Matrix()
    witness=next(null.row(i) for i in range(null.rows) if (null.row(i)*vec)[0]!=0)
    den=S.ilcm(*[S.denom(v) for v in witness]);integer=[int(v*den) for v in witness]
    common=reduce(gcd,(abs(v) for v in integer if v));integer=[v//common for v in integer]
    if next(v for v in integer if v)<0:integer=[-v for v in integer]
    row=S.Matrix([integer]);assert row*matrix==S.zeros(1,matrix.cols)
    value=(row*vec)[0];assert value!=0
    return {'target':label,'weight':list(target),'ambient_monomials':len(base),'all_possible_cofactor_columns':len(cols),
            'witness':[[list(exp),str(v)] for exp,v in zip(base,integer) if v],
            'evaluation_on_target':str(value),'annihilates_every_I8_cofactor_column':True}


def main():
    started=time.monotonic();a,b,c,d,e,f=S.symbols('a b c d e f');coords=[a,b,c,d,e,f]
    A,B,C,D,E,s,q=S.symbols('A B C D E s q');fm=symmetric_moments(a,b,c,A,E,s,q)
    covariance=[[a,d,e],[d,b,f],[e,f,c]];raw=moments(covariance,10)
    assert S.expand(raw[9]-pair_count_moment(10,covariance))==0
    simple={A:a*b*c,B:a*f*f,C:b*e*e,D:c*d*d,E:d*e*f}
    sb=B+C+D;qb=B*C+B*D+C*D
    for m,phi in enumerate(fm,1):
        pulled=S.expand(phi.subs({s:sb,q:qb}).subs(simple,simultaneous=True))
        assert S.expand(pulled-raw[2*m-1])==0
    AA=-2*s-8*E;QQ=(7*s*s+88*s*E+288*E*E)/20
    F=107*s**3+2196*s*s*E+11280*s*E*E+17152*E**3
    G=337*s**4+12736*s**3*E+109968*s*s*E*E+359424*s*E**3+406272*E**4
    substituted=[S.expand(v.subs(A,AA).subs(q,QQ)) for v in fm]
    assert substituted[:2]==[0,0]
    assert S.expand(substituted[2]+1620*F)==0 and S.expand(substituted[3]-381024*G)==0
    IG=S.groebner([F,G],s,E,domain=S.QQ,order='grevlex');base8,hf8=standard_basis(IG)
    remainder=S.expand(IG.reduce(substituted[4])[1])
    assert remainder==-S.Rational(225834448896000,92791)*E**5
    assert len(base8)==12 and hf8==[1,2,3,3,2,1]
    assert all(IG.reduce(s**i*E**(6-i))[1]==0 for i in range(7))
    I10G=S.groebner([F,G,substituted[4]],s,E,domain=S.QQ,order='grevlex');base10,hf10=standard_basis(I10G)
    assert len(base10)==11 and hf10==[1,2,3,3,2]
    # A global monic S-module presentation proves the six free basis elements.
    c_relation=C**2+B*C+B**2-s*C-s*B+q
    b_relation=B**3-s*B**2+q*B-A*E**2
    assert S.expand(c_relation.subs({s:sb,q:qb}))==0
    assert S.expand(b_relation.subs({s:sb,q:qb})-(B*C*D-A*E**2))==0
    freegb=S.groebner([c_relation,b_relation],C,B,domain=S.QQ.poly_ring(A,E,s,q),order='lex')
    assert {tuple(g.LM(order=freegb.order).exponents) for g in freegb.polys}=={(2,0),(0,3)}
    assert all(g.LC()==1 for g in freegb.polys)
    balanced=[]
    for cutoff,expected in [(4,[1,4,9,14,16,14,9,4,1]),(5,[1,4,9,14,16,13,7,2])]:
        ideal=[B*C*D-A*E**2]+[phi.subs({s:sb,q:qb}) for phi in fm[:cutoff]]
        bg=S.groebner(ideal,A,B,C,D,E,domain=S.QQ,order='grevlex');base,hf=standard_basis(bg)
        assert hf==expected
        balanced.append({'last_moment_order':2*cutoff,'length':len(base),'Hilbert_function':hf,
                         'irrelevant_ideal_nilpotence_index':len(hf),'standard_monomials':base,'Groebner_basis':list(map(str,bg))})
    # Exact global radical-power upper bounds by monomial extraction.
    extraction=[]
    for power,k in [(17,8),(19,9)]:
        count=0
        for i in range(power+1):
            for j in range(power+1-i):
                for l in range(power+1-i-j):
                    for aa in range(power+1-i-j-l):
                        ee=power-i-j-l-aa;odd=sum(v%2 for v in [i,j,l])
                        extracted=sum(v//2 for v in [i,j,l])+aa+ee+(2 if odd==3 else 0)
                        assert extracted>=k;count+=1
        extraction.append({'J_power':power,'balanced_power':k,'all_generator_monomials':count})
    # A genuine eighth radical power survives, even after M10 is added.
    z=S.symbols('z');ar=-2*(d*d+e*e+4*z*d*e)/(1+2*z*z)
    row=[S.cancel(raw[n-1].subs({b:1,c:1,f:z}).subs(a,ar)) for n in [4,6,8,10]]
    rowgb=S.groebner(row[:2],d,e,domain=S.QQ.frac_field(z),order='grevlex')
    assert rowgb.reduce(row[2])[1]==rowgb.reduce(row[3])[1]==0
    assert rowgb.reduce(e**8)[1]==e**8
    polys=[S.Poly(raw[n-1],*coords) for n in [2,4,6,8]]
    duals=[dual_nonmembership(coords,polys,S.Poly(raw[9],*coords),label,mult,delta)
           for label,mult,delta in [('M10',S.Integer(1),(0,0,0)),('a*M10',a,(2,0,0)),('d*M10',d,(1,1,0))]]
    cert={'field':'QQ; all claims characteristic zero','semigroup':semigroup_checks(),
          'symmetric_moments_F1_to_F5':list(map(str,fm)),'F':str(F),'G':str(G),
          'invariant_elimination':{'A':str(AA),'q':str(QQ),'M10_remainder':str(remainder),'I8_basis':list(map(str,IG)),
                                   'I8_standard_monomials':base8,'I8_Hilbert':hf8,'I10_basis':list(map(str,I10G)),
                                   'I10_standard_monomials':base10,'I10_Hilbert':hf10},
          'balanced_rings':balanced,'free_S_module_relations':[str(c_relation),str(b_relation)],'free_S_module_basis':['1','B','B^2','C','B*C','B^2*C'],'radical_power_extraction':extraction,
          'radical_power_lower_witness':'(b*e)^8 pulls back to e^8, nonzero over QQ(z) after b=c=1,f=z and M2 elimination',
          'full_Gram_nonmembership_duals':duals,
          'all_order_bridge':'Every M_(2m) is weighted homogeneous of degree m in the symmetric invariant ring. I8 annihilates all residual degrees>=6, so every M_(2m) for m>=6 is in I8; M10 is the sole new generator.',
          'boundary':'No full embedded-primary decomposition, exact global nilpotence index, positive-characteristic extension or claim for affine factors.'}
    artifact=HERE/'three_factor_global_ideal_certificate.json';artifact.write_text(json.dumps(cert,indent=2)+'\n')
    receipt={'all_moment_ideal':'I_infinity=I10=(M2,M4,M6,M8,M10), strictly larger than I8',
             'sharp_initial_ideal_cutoff':10,'set_theoretic_cutoff_retained':8,
             'all_moments_from_order12':'Already in I8 by the written invariant grading proof',
             'symmetric_invariant_lengths':[12,11],'balanced_lengths':[72,66],
             'global_nilradical_index_bounds':{'I8':[9,19],'I10_equals_I_infinity':[9,17]},
             'embedded_structure':'Nonzero I10/I8 vanishes at all six generic components by the retained generic algebra proof, so I8 has embedded associated structure.',
             'M10_remainder':str(remainder),'full_Gram_dual_certificates':len(duals),
             'certificate_sha256':hashlib.sha256(artifact.read_bytes()).hexdigest(),
             'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'input_hashes':{INPUT:EXPECTED_INPUT_SHA256},'sympy_version':S.__version__,'elapsed_seconds':time.monotonic()-started,
             'assurance':'Internal written all-order invariant proof and exact algebra. No independent replay, external review or formal verification; new characteristic claims are not made.'}
    (HERE/'three_factor_global_ideal_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
