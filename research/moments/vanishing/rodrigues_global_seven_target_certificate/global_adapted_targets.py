"""Read the exact local target vectors and check simple linear adaptations."""
from pathlib import Path
import json
import sympy as s
root=Path(__file__).parent
d=json.loads((root/'rodrigues_local_length17_certificate/certificate_data.json').read_text())
g=json.loads((root/'global_target_nilindices.json').read_text())
Ma=s.Matrix([[s.Rational(x) for x in row] for row in d['multiplication_a']])
Mb=s.Matrix([[s.Rational(x) for x in row] for row in d['multiplication_b']])
MM=[Ma**i*Mb**j for i,j in d['quotient_basis']]
vectors={name:s.Matrix([s.Rational(x) for x in data['vector']]) for name,data in g['target_profile'].items()}
normalized={name:s.Matrix([s.Rational(x) for x in v]) for name,v in g['normalized_target_vectors'].items()}
def ni(v):
    M=sum((x*MM[i] for i,x in enumerate(v) if x),s.zeros(17))
    power=s.eye(17)
    for k in range(1,8):
        power=power*M
        if power==s.zeros(17):return k
    raise AssertionError
def vector_product(v,w):
    M=sum((x*MM[i] for i,x in enumerate(v) if x),s.zeros(17))
    return M*w
combos={
 '4p+5s0':{'p':4,'s0':5},
 '35s1-34r':{'s1':35,'r':-34},
 'U002-6U011':{'U002':1,'U011':-6},
 'U111+2U012':{'U111':1,'U012':2},
 '4U022-13U112':{'U022':4,'U112':-13},
}
out={}
for name,coeffs in combos.items():
    raw=sum((coef*vectors[key] for key,coef in coeffs.items()),s.zeros(17,1))
    norm=sum((coef*normalized[key] for key,coef in coeffs.items()),s.zeros(17,1))
    e=ni(raw)
    assert ni(norm)==e
    out[name]={'nilpotence_index':e,'raw_vector':list(map(str,raw)),'normalized_vector':list(map(str,norm))}
    print(name,e,flush=True)
Tprime=normalized['U222']+s.Rational(896,75)*vector_product(normalized['p'],normalized['p'])
assert ni(Tprime)==2
out['U222+128s2+(896/75)p^2']={'nilpotence_index':2,'normalized_vector':list(map(str,Tprime))}
print('U222+128s2+(896/75)p^2',2,flush=True)
(root/'global_adapted_targets.json').write_text(json.dumps(out,indent=2)+'\n')
