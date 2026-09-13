"""First two local, level-preserving neighbours in the frozen first-chart pilot.

This is a small degree-two model experiment, not a full minimal-model enumerator.
The integral neighbour formula and final quartic transport are replayed exactly.
"""
import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import time

from sage.all import GF, PolynomialRing, QQ, ZZ
from cypari2 import Pari
from pointed_height_bounds import coefficients, substitute
from pointed_box_equivalence import box_key
from search_observability import multiply

R=PolynomialRing(QQ,'x');x=R.gen()


def invariants(P,Q):
    d=coefficients(4*P+Q*Q)
    a,b,c,e,f=d[4],d[3],d[2],d[1],d[0]
    I=12*a*f-3*b*e+c*c
    J=72*a*c*f+9*b*c*e-27*a*e*e-27*b*b*f-2*c**3
    return I,J/2


def build(data):
    start=time.process_time();old=data['mapping']
    P=R(list(map(QQ,old['reduced_P'])));Q=R(list(map(QQ,old['reduced_Q'])))
    D=4*P+Q*Q;inv=invariants(P,Q)
    # The retained Curve302 short equation is X=36*x_min+15.
    A,B=map(QQ,data['curve'][3:])
    if inv!=(-48*A/6**4,-864*B/6**6):
        raise ArithmeticError('first-chart minimal-model invariant binding differs')
    pari=Pari();pari.allocatemem(256000000,silent=True)
    seen={tuple(map(str,box_key(old['matrix'])))}
    models=[];attempts=[]
    for p in map(ZZ,data['certified_primes']):
        if p<3 or p>167:continue
        F=PolynomialRing(GF(p),'z');g=F(D).gcd(F(D.derivative()))
        roots=sorted(map(ZZ,g.roots(multiplicities=False))) if g.degree()>0 else []
        for r in roots:
            s=(-ZZ(Q(r))*pow(2,-1,int(p)))%p
            newP=(P(r+p*x)-s*Q(r+p*x)-s*s)/p**2
            newQ=(Q(r+p*x)+2*s)/p
            if any(c.denominator()!=1 for c in list(newP)+list(newQ)):
                attempts.append({'prime':int(p),'residue':int(r),'status':'NONINTEGRAL_NEIGHBOUR'});continue
            if invariants(newP,newQ)!=inv:
                raise ArithmeticError('neighbour changed c4/c6')
            polynomial=lambda f:'+'.join(f'({c})*x^{i}' for i,c in enumerate(f))
            ret=pari('my(m,C);C=hyperellred(['+polynomial(newP)+','+polynomial(newQ)+'],&m);[C,m]')
            reducedP=R([QQ(str(ret[0][0].polcoef(i))) for i in range(5)])
            reducedQ=R([QQ(str(ret[0][1].polcoef(i))) for i in range(3)])
            change=tuple(QQ(str(ret[1][1][i,j])) for i in range(2) for j in range(2))
            second=multiply((int(p),int(r),0,1),tuple(Fraction(str(c)) for c in change))
            matrix=multiply(tuple(map(Fraction,old['matrix'])),second)
            key=tuple(map(str,box_key(matrix)))
            if key in seen:
                attempts.append({'prime':int(p),'residue':int(r),'status':'DUPLICATE_BOX'});continue
            disc=4*reducedP+reducedQ*reducedQ
            transformed=R(substitute(old['raw_coefficients'],matrix).list())
            j=next(i for i in range(5) if disc[i]);ratio=transformed[j]/disc[j]
            if (invariants(reducedP,reducedQ)!=inv or transformed!=ratio*disc or
                ratio<=0 or not ratio.is_square()):
                raise ArithmeticError('reduction or composed chart changed')
            seen.add(key)
            result=deepcopy(data);result['policy']=f'minimal-neighbour-p{p}-r{r}'
            m=deepcopy(old)
            m.update(matrix=list(map(str,matrix)),first_matrix=old['matrix'],
                second_matrix=list(map(str,second)),
                reduced_P=list(map(str,coefficients(reducedP))),
                reduced_Q=list(map(str,coefficients(reducedQ,2))),
                discriminant_quartic=list(map(str,coefficients(disc))),square_ratio=str(ratio),
                coordinate_policy={'kind':'raw','matrix':list(map(str,matrix))},
                minimization_policy='integral level-preserving local neighbour, then hyperellred')
            result['mapping']=m
            witness={'prime':str(p),'abscissa_residue':str(r),'ordinate_residue':str(s),
                'unreduced_P':list(map(str,coefficients(newP))),
                'unreduced_Q':list(map(str,coefficients(newQ,2))),
                'reduction_matrix':list(map(str,change)),'c4':str(inv[0]),'c6':str(inv[1])}
            result['neighbour_witness']=witness;models.append(result)
            attempts.append({'prime':int(p),'residue':int(r),'status':'RETAINED_DISTINCT_NEIGHBOUR'})
            if len(models)==2:
                return models,{'attempts':attempts,'cpu_seconds':time.process_time()-start}
    return models,{'attempts':attempts,'cpu_seconds':time.process_time()-start}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if args.output_dir.exists():raise FileExistsError('preserve previous neighbour attempt')
    models,record=build(json.loads(args.input.read_text()))
    args.output_dir.mkdir()
    for i,model in enumerate(models):
        (args.output_dir/f'neighbour-{i}-input.json').write_text(json.dumps(model,indent=2)+'\n')
    record.update(input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),models=len(models))
    (args.output_dir/'construction.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record),flush=True)
