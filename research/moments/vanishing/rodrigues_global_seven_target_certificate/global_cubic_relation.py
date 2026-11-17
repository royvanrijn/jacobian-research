"""One bounded exact relation test among five character-six local vectors."""
import json
from pathlib import Path
import sympy as sp
root=Path(__file__).parent
data=json.loads((root/'rodrigues_local_length17_certificate/certificate_data.json').read_text())
targets=json.loads((root/'global_target_nilindices.json').read_text())
Ma=sp.Matrix([[sp.Rational(x) for x in row] for row in data['multiplication_a']])
Mb=sp.Matrix([[sp.Rational(x) for x in row] for row in data['multiplication_b']])
MM=[Ma**i*Mb**j for i,j in data['quotient_basis']]
normalized={name:sp.Matrix([sp.Rational(x) for x in vector]) for name,vector in targets['normalized_target_vectors'].items()}
def multiplication(v):return sum((x*MM[i] for i,x in enumerate(v) if x),sp.zeros(17))
def mul(v,w):return multiplication(v)*w
one=sp.eye(17)[:,0]
def power(v,n):return multiplication(v)**n*one
s0=normalized['s0'];s1=normalized['s1'];p=normalized['p']
s2=sp.Matrix([sp.Rational(x) for x in targets['normalized_s2']])
K1=4*p+5*s0
K6=4*normalized['U022']-13*normalized['U112']
vectors=[power(s0,3),power(s1,2),mul(s2,K1),mul(K1,power(s0,2)),K6]
M=sp.Matrix.hstack(*vectors)
null=M.nullspace()
print('Five-vector rank:',M.rank(),'nullity:',len(null),flush=True)
solutions=[v/v[0] for v in null if v[0]]
if solutions:
    v=solutions[0]
    assert v[0]==1 and M*v==sp.zeros(17,1)
    print('MONIC RELATION:',list(v),flush=True)
    payload={'terms':['s0^3','s1^2','s2*(4p+5s0)','(4p+5s0)*s0^2','4U022-13U112'],
             'coefficients':list(map(str,v)),'rank':M.rank(),'verified':True,
             'normalized_s2':list(map(str,s2))}
    common=sp.ilcm(*[x.q for x in v])
    ints=[int(common*x) for x in v]
    gcd=sp.igcd(*ints)
    ints=[int(x//gcd) for x in ints]
    assert M*sp.Matrix(ints)==sp.zeros(17,1)
    assert all(sp.Rational(coefficient,ints[0])==value for coefficient,value in zip(ints,v))
    payload['primitive_integer_coefficients']=ints
    payload['ordinary_polynomial_degree']=3
    payload['seven_target_maximum_degree']=5
    print('PRIMITIVE INTEGER RELATION:',ints,flush=True)
else:
    print('NO MONIC RELATION',flush=True)
    payload={'rank':M.rank(),'monic_relation':False,'nullspace':[list(map(str,v)) for v in null]}
Path(__file__).with_suffix('.json').write_text(json.dumps(payload,indent=2)+'\n')
