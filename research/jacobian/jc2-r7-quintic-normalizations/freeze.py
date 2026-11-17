"""Discovery exporter: turn Singular cofactors into literal integer identities."""
from pathlib import Path
import sys,json,gzip,time
import sympy as S
HERE=Path(__file__).resolve().parent
A,M,d,z=S.symbols('A M d z')

def freeze(label):
    family='L' if label=='L' else 'K';parameter=d if family=='L' else z
    row=json.loads((HERE/('reconstruction-'+family+'.json')).read_text())
    generators=[]
    for degree in range(9,3,-1):
        value=S.sympify(row['residual'][str(degree)],locals=dict(A=A,M=M,d=d,z=z))
        generators.append(S.Poly(value,A,M,parameter).clear_denoms()[1].primitive()[1].as_expr())
    cofactors=[S.sympify((HERE/f'cofactor-{label}-{i}.txt').read_text().replace('^','**'),locals=dict(A=A,M=M,d=d,z=z)) for i in range(1,7)]
    domain=S.QQ.frac_field(parameter)
    polynomials=[S.Poly(c,A,M,domain=domain) for c in cofactors]
    result=sum((polynomials[i]*S.Poly(generators[i],A,M,domain=domain) for i in range(6)),S.Poly(0,A,M,domain=domain))
    assert result==S.Poly(1,A,M,domain=domain)
    denominator=S.Poly(1,parameter,domain=S.QQ)
    for poly in polynomials:
        for co in poly.coeffs():
            denominator=S.lcm(denominator,S.Poly(S.fraction(S.cancel(co))[1],parameter,domain=S.QQ)).monic()
    integral=[S.Poly(S.expand(sum(S.cancel(co*denominator.as_expr())*A**i*M**j for (i,j),co in poly.terms())),A,M,parameter,domain=S.QQ) for poly in polynomials]
    numeric=S.ilcm(*[S.denom(co) for poly in integral+[S.Poly(denominator.as_expr(),A,M,parameter)] for co in poly.coeffs()])
    integral=[S.Poly(poly.as_expr()*numeric,A,M,parameter,domain=S.ZZ) for poly in integral]
    rhs=S.Poly(denominator.as_expr()*numeric,parameter,domain=S.ZZ)
    terms=lambda poly:[[list(m),str(co)] for m,co in poly.terms()]
    payload=dict(schema='quintic-normalization-necessary-unit-v1',family=family,label=label,
        variables=['A','M',str(parameter)],degrees=list(range(9,3,-1)),
        generators=[terms(S.Poly(v,A,M,parameter)) for v in generators],
        cofactors=[terms(poly) for poly in integral],rhs=terms(rhs),
        assurance='Discovery identity checked in QQ(parameter)[A,M]; canonical verifier independently convolves frozen integer coefficients.')
    (HERE/('unit-'+label+'.json.gz')).write_bytes(gzip.compress((json.dumps(payload,indent=2)+'\n').encode(),mtime=0))
    print(label,'terms',[len(v.terms()) for v in integral],'RHS degree',rhs.degree(),flush=True)
    return rhs

if __name__=='__main__':
    start=time.monotonic()
    for label in sys.argv[1:]:freeze(label)
    print('elapsed',time.monotonic()-start,flush=True)
