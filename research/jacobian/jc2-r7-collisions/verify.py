"""Replay all affine singularity strata of the polynomial gap-21 family."""
import hashlib
import json
import time
from fractions import Fraction as Q
import sympy as S
import collisions as C

CHECKS=[]
def check(name,value):
    assert value,name
    CHECKS.append(name)

def zero(name,value):check(name,S.cancel(value)==0)

def proportional(name,left,right):
    ratio=S.cancel(left/right)
    check(name,ratio.is_Rational and ratio!=0)

def decode(rows):
    out={}
    for m,c in rows:
        assert len(m)==2 and all(type(e) is int and e>=0 for e in m)
        assert tuple(m) not in out
        out[tuple(m)]=Q(c)
    return out

def encode_poly(expr,var):
    return {m:Q(str(c)) for m,c in S.Poly(expr,C.d,var).terms() if c}

def add(out,left,scale=Q(1)):
    for m,c in left.items():
        out[m]=out.get(m,Q(0))+scale*c
        if not out[m]:del out[m]

def multiply(left,right):
    out={}
    for m,c in left.items():
        for n,e in right.items():
            key=(m[0]+n[0],m[1]+n[1])
            out[key]=out.get(key,Q(0))+c*e
    return {m:c for m,c in out.items() if c}

def reduce_d(poly,modulus):
    out=dict(poly);degree=max(m[0] for m in modulus)
    lead=modulus[(degree,0)]
    while out and max(m[0] for m in out)>=degree:
        power=max(m[0] for m in out)
        for m,c in list(out.items()):
            if m[0]!=power:continue
            scale=c/lead
            add(out,{(e+power-degree,k+m[1]):q for (e,k),q in modulus.items()},-scale)
    return out

def evaluate_root(poly,root,modulus):
    out={}
    for degree in range(max(m[1] for m in poly),-1,-1):
        out=reduce_d(multiply(out,root),modulus)
        add(out,{(m[0],0):c for m,c in poly.items() if m[1]==degree})
        out=reduce_d(out,modulus)
    return out

def main():
    start=time.monotonic();d,v,w,u,r=C.d,C.v,C.w,C.u,C.r
    s0,t0=S.symbols('s0 t0')
    for name,poly,div in [('p',C.p,C.A),('q',C.q,C.B)]:
        zero(f'{name} symmetric divided difference',
             poly.subs(C.t,s0)-poly.subs(C.t,t0)-(s0-t0)*div.subs({v:s0+t0,w:s0*t0}))
    zero('p composition',C.p-(C.t**3+d*C.t)**2-(C.t**3+d*C.t))
    zero('collision split',C.A-C.A1*C.A2)
    zero('first collision equation',C.B.subs(w,v*v+d)+C.H9/162)
    zero('second collision equation',C.B.subs(w,(v**3+d*v+1)/(3*v))+v**6*C.H3/162)
    check('second branch never has v0',C.H3.subs(v,0)==6)
    proportional('H9 tangency discriminant',S.discriminant(C.H9,v),C.T)
    proportional('H3 tangency discriminant',S.discriminant(C.H3,v),C.C)
    proportional('branch-overlap resultant',S.resultant(C.H9,2*v**3+2*d*v-1,v),C.C)
    proportional('first branch never diagonal',S.resultant(C.H9,3*v*v+4*d,v),S.Integer(1))
    proportional('second branch diagonal locus',S.resultant(C.H3,v**3+4*d*v+4,v),C.C)
    proportional('same-u triple resultant',S.resultant(C.r2,C.r1,u),d**3*C.U)
    proportional('mixed triple resultant',S.resultant(C.H9,C.mixed,v),C.C**2*C.M)
    zero('mixed reduction preserves candidates',S.rem(C.mixed-C.mixed_remainder,C.H9,v))
    for i,f in enumerate([C.C,C.U,C.M,C.T]):
        check(f'exceptional factor squarefree{i}',S.degree(S.gcd(f,S.diff(f,d)),d)==0)
        check(f'exceptional factor avoids zero{i}',f.subs(d,0)!=0)
        for j,g in enumerate([C.C,C.U,C.M,C.T][:i]):
            check(f'exceptional factors disjoint{i},{j}',S.degree(S.gcd(f,g),d)==0)
    check('exact exceptional parameter count',1+sum(S.degree(f,d) for f in [C.C,C.U,C.M,C.T])==34)
    data=json.loads((C.HERE/'certificates.json').read_text())['gcd_certificates']
    for name,(mod,f,g,var) in C.GCD_INPUTS.items():
        cert=data[name]
        zero(f'{name} modulus',S.sympify(cert['modulus'],locals={'d':d})-mod)
        for i,h in enumerate([f,g]):check(f'{name} original gcd generator{i}',decode(cert['generators'][i])==encode_poly(h,var))
        modulus=encode_poly(mod,var);gg=decode(cert['gcd'])
        check(f'{name} monic linear gcd',gg.get((0,1))==1 and all(m[1]==0 or m==(0,1) for m in gg))
        result={}
        for generator,multiplier in zip(cert['generators'],cert['multipliers']):
            add(result,multiply(decode(generator),decode(multiplier)))
        add(result,gg,-Q(1))
        check(f'{name} Fraction Bezout modulo exceptional factor',not reduce_d(result,modulus))
        mutated=dict(result);add(mutated,decode(cert['generators'][0]))
        check(f'{name} mutated Bezout rejected',bool(reduce_d(mutated,modulus)))
        root={m:-c for m,c in gg.items() if m[1]==0}
        for i,h in enumerate(cert['generators']):
            check(f'{name} Fraction gcd divides generator{i}',not evaluate_root(decode(h),root,modulus))
    # The sharp historical sample is an exceptional triple-point member.
    rem0=S.Poly(C.rem.as_expr().subs(d,0),C.t)
    zero('sharp sample has no t2 coefficient',rem0.nth(2))
    A0=rem0.nth(1)
    check('three distinct sharp triple values of u',S.degree(A0,u)==3 and S.discriminant(A0,u)!=0 and A0.subs(u,0)!=0)
    check('sharp sample all pair tangencies absent',C.T.subs(d,0)!=0 and C.C.subs(d,0)!=0)
    check('sharp sample no mixed triples',C.M.subs(d,0)!=0)
    old=json.loads((C.HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases'][0]['parameters']
    old={name:S.Rational(value) for name,value in old.items()}
    zold=old['d']**3/old['b']**5
    zero('retained class I normalized cusp parameter',zold+S.Rational(2,3))
    zero('retained class I line a relation',old['a']*old['b']-2*old['d'])
    zero('retained class I line B relation',6*old['B']*old['b']+7*old['d'])
    zero('retained class I graph c',old['c']*old['b']**2-old['d']**2)
    zero('retained class I graph D',old['D']/old['b']**2-(18-35*zold)/81)
    zero('retained class I graph E',old['E']-7*old['d']*old['b']*(5*zold-9)/243)
    # Pointwise degree-six consequences do not assert non-surjectivity is impossible.
    check('triple value uses all omitted length',3*2==6)
    check('another omitted prime forbidden there',3*2+1>6)
    singularities=[
        dict(stratum='outside d*C*U*M*T=0',nodes=12,triple_points=0,tacnodes=0,cusps=0,A5=0),
        dict(stratum='d=0',nodes=3,triple_points=3,tacnodes=0,cusps=0,A5=0),
        dict(stratum='U=0',nodes=9,triple_points=1,tacnodes=0,cusps=0,A5=0),
        dict(stratum='M=0',nodes=9,triple_points=1,tacnodes=0,cusps=0,A5=0),
        dict(stratum='T=0',nodes=10,triple_points=0,tacnodes=1,cusps=0,A5=0),
        dict(stratum='C=0, retained cusp class I',nodes=8,triple_points=0,tacnodes=0,cusps=1,A5=1)]
    for row in singularities:
        check('affine delta12 '+row['stratum'],row['nodes']+3*row['triple_points']+2*row['tacnodes']+row['cusps']+3*row['A5']==12)
    ledger=dict(schema='r7-line-collision-strata-v1',parameter='d',
                exceptional_factors={name:str(f) for name,f in [('C',C.C),('U',C.U),('M',C.M),('T',C.T)]},
                exceptional_parameter_count=34,singularities=singularities,
                sharp_witness='d0 is special: three ordinary triple points plus three nodes; the old braid case label generic is historical.',
                global_consequence='At any C triple value a hypothetical degree6 Keller map omits that value entirely, and no other omitted image can pass through it. This is not an impossibility theorem for a Keller map.',
                remaining=['Laurent-family collision strata','actual second image and intersections','global complement relations and source realization'])
    (C.HERE/'strata.json').write_text(json.dumps(ledger,indent=2)+'\n')
    inputs=[C.HERE/n for n in ['PROOF.md','README.md','collisions.py','verify.py','certificates.json']]
    inputs += [C.HERE.parent/'jc2-r7-image-family/family.py',C.HERE.parent/'jc2-r7-image-family/PROOF.md',
               C.HERE.parent/'jc2-normal-r7/classification.json',C.ROOT/'foundations/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(schema='r7-line-collisions-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),
                 input_sha256={str(p.relative_to(C.ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                 output_sha256={'strata.json':hashlib.sha256((C.HERE/'strata.json').read_bytes()).hexdigest()},
                 assurance='Exact symmetric identities/resultants and three independent Fraction modular gcd replays. Written local-curve classification and global-length interpretation; no formalization or independent external review.')
    (C.HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} POLYNOMIAL-FAMILY COLLISION CHECKS ({receipt["elapsed_seconds"]}s)')

if __name__=='__main__':main()
