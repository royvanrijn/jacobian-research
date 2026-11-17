"""Exact reciprocal-contact certificates; no elimination is run by default.

The base-map/contact bridge and inherited geometry are written proof.
All endpoint memberships/units/powers replay with Fraction arithmetic.
"""
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
import hashlib,json,time

import sympy as S
from sympy.polys.rings import ring
from contact import (P,a,b,c,d,B,D,E,h,v,w,p,q,rho,coefficient,
                     weights,degrees,columns,remainders)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
CHECKS=[]


def check(name,value):
    if not value:raise AssertionError(name)
    CHECKS.append(name)


def zero(name,value):check(name,S.expand(value)==0)


def decode(rows,n):
    result={}
    for exponent,coefficient in rows:
        key=tuple(exponent);value=Q(coefficient)
        assert len(key)==n and all(type(i) is int and i>=0 for i in key)
        assert key not in result and value
        result[key]=value
    return result


def add(a,b):
    out=dict(a)
    for m,c0 in b.items():out[m]=out.get(m,Q(0))+c0
    return {m:c0 for m,c0 in out.items() if c0}


def multiply(a,b):
    out={}
    for m,c0 in a.items():
        for n,d0 in b.items():
            key=tuple(x+y for x,y in zip(m,n));out[key]=out.get(key,Q(0))+c0*d0
    return {m:c0 for m,c0 in out.items() if c0}


def combination(generators,multipliers,n):
    result={}
    for f,g in zip(generators,multipliers):result=add(result,multiply(decode(f,n),decode(g,n)))
    return result


def dense(coefficients):return {(i,):Q(c0) for i,c0 in enumerate(coefficients) if Q(c0)}


def dense_combination(fs,ms):
    result={}
    for f,g in zip(fs,ms):result=add(result,multiply(dense(f),dense(g)))
    return result


def array_expr(rows,variables):
    return sum(S.Rational(c0)*S.prod(x**e for x,e in zip(variables,m)) for m,c0 in rows)


def dense_expr(rows,variable):return sum(S.Rational(c0)*variable**i for i,c0 in enumerate(rows))


def primitive(expression,variables):
    return S.Poly(S.expand(expression),*variables).clear_denoms()[1].primitive()[1].as_expr()


def determinant(matrix,rows,R):
    out=R.zero
    for perm in permutations(range(4)):
        term=R.ground_new((-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4)))
        for i,j in enumerate(perm):term*=matrix[rows[i]][j]
        out+=term
    return out


def augmented(equations,unknowns,R):
    matrix=[]
    for f in equations:
        coefficients=[S.diff(f,x) for x in unknowns]
        constant=f.subs(dict.fromkeys(unknowns,0))
        assert all(not g.has(*unknowns) for g in coefficients)
        assert S.expand(f-constant-sum(x*g for x,g in zip(unknowns,coefficients)))==0
        matrix.append([R.from_expr(g) for g in coefficients+[constant]])
    return matrix


def only_denominator_root(denominator,variable,value):
    f=S.Poly(denominator,variable,domain=S.QQ)
    while f.degree()>0:
        qq,rr=f.div(S.Poly(variable-value,variable))
        if rr:return False
        f=qq
    return f.as_expr()!=0


def main():
    start=time.monotonic();cert=json.loads((HERE/'certificates.json').read_text())
    sa,sb,sc,sd,sB,sD,sE,sh,sv,sw,sZ,st=S.symbols('a b c d B D E h v w Z t')
    check('universal column grading',all(
        all(sum(x*y for x,y in zip(weights,m))==degree-n for m in poly)
        for degree,column in zip(degrees,columns) for n,poly in enumerate(column)))
    check('cubic squared-remainder grading',all(
        sum(x*y for x,y in zip((1,1,3,2),m))+r==n
        for n,rows in enumerate(remainders) for r,poly in enumerate(rows) for m in poly))
    zero('birational parity gate',
         (rho[5].as_expr()+2*sa*rho[3].as_expr()/3+sc*rho[1].as_expr()/3).subs(sb,0)+5*sd**3/81)
    ch=-15*h**4-6*a*h*h-3*b*h;dh=24*h**5+8*a*h**3+3*b*h*h
    zero('marked first derivative',S.diff(p,st).subs(st,sh).subs(sc,ch.as_expr()).subs(sd,dh.as_expr()))
    zero('marked second derivative',S.diff(p,st,2).subs(st,sh).subs(sc,ch.as_expr()).subs(sd,dh.as_expr()))
    zero('universal v4 identity',coefficient((2,3,4,5),13).as_expr()+60*sv**4)
    marked={v:P.zero,c:ch,d:dh}
    zero('nonzero marked root forces w',coefficient((1,3,4,5),9,marked).as_expr()-S.Rational(40,27)*sh*(3*sw-sb)**4)

    # Centered chart: exact nonreduced scheme, not its radical ideal.
    scheme=cert['zero_scheme'];variables=S.symbols(' '.join(scheme['variables']))
    check('centered coefficient variables are explicit symbols',tuple(map(str,variables))==('E','D','B','a'))
    zero_sub={sb:1,sc:0,sd:0}
    fs=[primitive(rho[j].as_expr().subs(zero_sub),variables) for j in (1,3,5,7,9)]
    for i,(f,encoded) in enumerate(zip(fs,scheme['generators'])):zero(f'zero scheme generator{i}',f-array_expr(encoded,variables))
    check('zero scheme certificate dimensions',len(fs)==5 and len(scheme['targets'])==len(scheme['multipliers'])==4)
    expected=[3*sB-5*sa,9*sD-2,27*sE-4*sa,sa*sa]
    for i,(target,multipliers) in enumerate(zip(scheme['targets'],scheme['multipliers'])):
        zero(f'zero scheme target{i}',array_expr(target,variables)-expected[i])
        check(f'zero scheme Fraction membership{i}',combination(scheme['generators'],multipliers,4)==decode(target,4))
    for i,f in enumerate(fs):
        reduced=S.Poly(S.expand(f.subs({sB:5*sa/3,sD:S.Rational(2,9),sE:4*sa/27})),sa)
        check(f'zero scheme reverse containment{i}',all(m[0]>=2 for m,c0 in reduced.terms() if c0))

    zero_chart={h:P.zero,v:P.zero,a:P.zero,b:P.one,c:P.zero,d:P.zero,B:P.zero,D:P.ground_new(S.Rational(2,9)),E:P.zero}
    zcert=cert['zero_contact']
    for i,power in enumerate(zcert['Z_powers']):
        zero(f'centered contact coefficient Z{power}',coefficient(tuple(zcert['rows']),power,zero_chart).as_expr()-dense_expr(zcert['polynomials'][i],sw))
    check('centered Fraction unit',dense_combination(zcert['polynomials'],zcert['multipliers'])=={(0,):Q(1)})
    tampered=json.loads(json.dumps(zcert['multipliers']));tampered[0][0]=str(Q(tampered[0][0])+1)
    check('centered mutated unit rejected',dense_combination(zcert['polynomials'],tampered)!={(0,):Q(1)})

    # Nonzero chart h1, after its forced w=b/3.
    sub1={h:P.one,v:P.zero,w:b/3,c:-15-6*a-3*b,d:24+8*a+3*b}
    marked1={sc:-15-6*sa-3*sb,sd:24+8*sa+3*sb}
    gaps=[rho[j].as_expr().subs(marked1) for j in (1,3,5,7,9)]
    F10=coefficient((1,2,3,4),10,sub1).as_expr()
    V7=coefficient((1,3,4,5),7,sub1).as_expr()
    F9=coefficient((1,2,3,4),9,sub1).as_expr()
    T=sa*sb-12*sa-3*sb-36
    zero('T0 exceptional b12 is impossible',T.subs(sb,12)+72)
    # T0 rank minors. Clearing only known nonzero denominator factors.
    av=3*(sb+12)/(sb-12);eq0=[]
    for f in gaps+[F9]:
        numerator,denominator=S.fraction(S.cancel(f.subs(sa,av)))
        check('T0 denominator has only excluded root12',only_denominator_root(denominator,sb,12))
        eq0.append(primitive(numerator,(sB,sD,sE,sb)))
    R0,bb=ring('b',S.QQ);matrix0=augmented(eq0,(sB,sD,sE),R0)
    tcert=cert['T_zero'];factor=bb**tcert['common_b_power']*(bb-12)**tcert['common_b_minus_12_power']
    check('T0 certificate branch dimensions',len(tcert['rows'])==len(tcert['quotients'])==len(tcert['multipliers'])==4)
    for i,rows in enumerate(tcert['rows']):
        det=determinant(matrix0,rows,R0);quotient,remainder=det.div(factor)
        check(f'T0 exact determinant quotient{i}',not remainder and S.expand(quotient.as_expr()-dense_expr(tcert['quotients'][i],sb))==0)
    check('T0 independent Fraction unit',dense_combination(tcert['quotients'],tcert['multipliers'])=={(0,):Q(1)})
    tampered=json.loads(json.dumps(tcert['multipliers']));tampered[0][0]=str(Q(tampered[0][0])+1)
    check('T0 mutated unit rejected',dense_combination(tcert['quotients'],tampered)!={(0,):Q(1)})

    # T nonzero rank minors. All coefficient variables remain symbolic.
    extras=[S.cancel(F10/T),S.cancel(V7/T)]
    check('contact factors divide as polynomials',all(S.fraction(f)[1].is_Rational for f in extras))
    eq1=[primitive(f,(sB,sD,sE,sa,sb)) for f in gaps+extras]
    R1,aa,bbb=ring('a,b',S.QQ);matrix1=augmented(eq1,(sB,sD,sE),R1)
    check('all necessary systems affine-linear',len(matrix0)==6 and len(matrix1)==7)
    pc=cert['T_nonzero'];sr,sbb=S.symbols('r b');RR,rr,bb2=ring('r,b',S.QQ)
    check('power certificate dimensions',len(pc['generators'])==len(pc['multipliers'])==len(pc['minor_metadata'])==12)
    for i,item in enumerate(pc['minor_metadata']):
        det=determinant(matrix1,item['rows'],R1)
        unit_factor=bbb**item['b_power']*R1.from_expr(T)**item['T_power']
        quotient,remainder=det.div(unit_factor);check(f'offT removed factor exact{i}',not remainder)
        primitive_ab=quotient.primitive()[1].as_expr()
        actual=primitive(primitive_ab.subs(sa,(sr-sb-6)/2),(sr,sb))
        zero(f'offT exact normalized minor{i}',actual-array_expr(pc['generators'][i],(sr,sb)))
    target={(0,pc['target_power']):Q(1)}
    check('offT independent Fraction b-power identity',combination(pc['generators'],pc['multipliers'],2)==target)
    tampered=json.loads(json.dumps(pc['multipliers']));first=decode(tampered[0],2);first[(0,0)]=first.get((0,0),Q(0))+1
    tampered[0]=[[list(m),str(q0)] for m,q0 in sorted(first.items()) if q0]
    check('offT mutated power rejected',combination(pc['generators'],tampered,2)!=target)

    # Independent full-matrix control for the v4 coefficient.
    RF,tt,vv,ww,ZZ=ring('t,v,w,Z',S.QQ);divisor=(tt**3+vv*tt**2-ZZ*tt+ww)**2
    RM,vf,wf,Zf=ring('v,w,Z',S.QQ)
    pp=st**6+st**3;qq=st**10+S.Rational(5,3)*st**7+S.Rational(7,9)*st**4+S.Rational(7,81)*st
    raw=[]
    for f in [pp,qq,pp*pp,pp*qq]:
        _,rem=RF.from_expr(f).div(divisor)
        raw.append([RM.from_dict({m[1:]:co for m,co in rem.items() if m[0]==j}) for j in range(6)])
    full=[[column[j] for column in raw] for j in range(1,6)]
    det=determinant(full,[1,2,3,4],RM)
    coeff=RM.from_dict({m[:-1]+(0,):co for m,co in det.items() if m[-1]==13})
    zero('full-matrix v4 control',coeff.as_expr()+60*sv**4)
    good_p=-st**2;good_q=st**4+2*st
    zero('genuine lower-degree three-contact positive control',
         good_p*good_q+sZ*good_q+sZ*good_p**2+sZ**2*good_p-1+(st**3-sZ*st+1)**2)
    check('positive control is outside excluded degrees',(S.degree(good_p,st),S.degree(good_q,st))==(2,4))

    # Exact necessary global ledger, conditional on the written row exclusion.
    profiles=[dict(other_index=e,excess=X,unramified_primes=t0) for e in range(2,6)
              for X in range(6) for t0 in range(6) if 2+e+X+t0==5]
    check('three two-ramified-prime Euler ledgers',profiles==[
        dict(other_index=2,excess=0,unramified_primes=1),dict(other_index=2,excess=1,unramified_primes=0),dict(other_index=3,excess=0,unramified_primes=0)])
    check('a third ramified prime exceeds Euler budget',3*2>5)
    check('two ramified primes cannot share a two-branch image',4*2>6)

    inputs=[HERE/'PROOF.md',HERE/'README.md',HERE/'REVIEW.md',HERE/'contact.py',HERE/'certificates.json',Path(__file__).resolve(),
            ROOT/'foundations/scripts/verify_f2_degree_6_10_gap.py',
            HERE.parent/'jc2-r7-cubic-cover/PROOF.md',HERE.parent/'jc2-r7-universal-quartic/PROOF.md',
            HERE.parent/'jc2-omitted-divisors/GLOBAL_LEDGER.md']
    hashes={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    receipt=dict(schema='r7-reciprocal-contact-v1',passed=True,check_count=len(CHECKS),checks=CHECKS,
                 elapsed_seconds=round(time.monotonic()-start,3),input_sha256=hashes,
                 centered_scheme='3B-5a,9D-2,27E-4a,a^2; one complex point, length2 scheme',
                 nonzero_branches='T0 univariate unit; offT twelve-minor power certificate b^13',
                 single_ramified_prime_gap21_stratum='EXCLUDED under the declared degree6,e2f1 hypotheses',
                 necessary_two_ramified_prime_ledgers=profiles,
                 assurance='Exact algebraic identities and Fraction replays. Geometric inheritance, cubic pole/contact bridge and global second-image argument are written proofs; no formal verification or independent external review.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(f'PASS {len(CHECKS)} RECIPROCAL CONTACT CHECKS ({receipt["elapsed_seconds"]}s)')
    print('Both marked-root charts excluded; the declared entire single-ramified-prime stratum is impossible.')


if __name__=='__main__':main()
