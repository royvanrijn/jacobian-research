#!/usr/bin/env python3
"""Exact semilinear ideal control and finite arithmetic/volume checks; not an all-degree proof replay."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form,smith_normal_form
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def integers(A):
    assert all(v.q==1 for v in A),'nonintegral asserted lattice map'
    return s.Matrix([[int(v) for v in row] for row in A.tolist()])
def serial(A):return [[int(v) for v in row] for row in A.tolist()]
def main():
    started=time.monotonic()
    # A cyclic-shift base implements an actual nontrivial holonomy automorphism.
    shifts=[1,2,4]
    ck('alpha2 cyclic-shift base intertwines C7',[shifts[(j+1)%3]%7 for j in range(3)]==[2*a%7 for a in shifts])
    # Full conductor23 model on a nonprincipal ideal, including its index and source lattice.
    z=s.symbols('z');d=22;Phi=s.cyclotomic_poly(23,z);qr=sorted({a*a%23 for a in range(1,23)});tau=sum(z**j for j in qr)
    poly=s.gcd(s.Poly(Phi,z,modulus=2),s.Poly(tau,z,modulus=2)).monic().as_expr()
    expected=z**11+z**9+z**7+z**6+z**5+z+1
    ck('prime ideal exact degree11 residue factor',poly==expected and s.Poly(poly,z,modulus=2).is_irreducible)
    def multiplication(f):
        return s.Matrix(d,d,lambda i,j:s.Poly(s.rem(f*z**j,Phi,z),z).coeff_monomial(z**i))
    R=multiplication(z);Tau=multiplication(tau);I=s.eye(d);F=multiplication(poly)
    Sigma=s.Matrix(d,d,lambda i,j:s.Poly(s.rem(z**(5*j),Phi,z),z).coeff_monomial(z**i))
    B=hermite_normal_form((2*I).row_join(F));A=I+Tau/2
    M=integers(B.inv()*A*Sigma*B);RI=integers(B.inv()*R*B)
    ck('prime ideal basis index2048',B.det()==2048)
    ck('semilinear endomorphism determinant minus2048',M.det()==-2048)
    ck('integral holonomy lattice order23',RI**23==I and RI!=I)
    ck('semilinear fifth-power intertwining',M*RI==RI**5*M)
    ck('primitive Galois generator5 has order22',s.n_order(5,23)==22)
    ck('quadratic period exact equation',Tau*Tau+Tau+6*I==s.zeros(d))
    ck('Galois generator exchanges quadratic periods',Sigma*Tau*Sigma.inv()==-I-Tau)
    ck('quadratic norm-eight element has full norm2^33',(2*I+Tau).det()==2**33)
    P2=hermite_normal_form((4*I).row_join(2*F).row_join(F*F))
    image=hermite_normal_form(integers(A*Sigma*B))
    ck('semilinear image is exactly P squared',image==P2 and P2.det()==2**22)
    ck('residue cyclic module generator from scalar2',hermite_normal_form((2*I).row_join(P2))==B)
    smith=[abs(int(v)) for v in smith_normal_form(M,domain=s.ZZ).diagonal()]
    ck('elementary degree2048 kernel exact Smith factors',smith==[1]*11+[2]*11)
    ck('residue order gives irreducible C23 module',s.n_order(2,23)==11)
    ck('base and fibre orientation give total Jacobian2048',Sigma.det()*M.det()==2048)
    # A deliberately small primitive-divisor control grid; the theorem uses the classical input.
    prime_controls=[]
    for p,e in [(2,2),(2,3),(2,5),(2,6),(2,9),(3,2),(3,3),(5,2),(7,2),(7,3)]:
        primes=sorted(int(q) for q in s.factorint(p**e-1) if q!=2 and s.n_order(p,int(q))==e)
        exceptional=(p,e)==(2,6) or (e==2 and s.isprime(p) and (p+1)&p==0)
        ck(f'primitive-divisor control p{p} e{e}',bool(primes)!=exceptional)
        prime_controls.append({'p':p,'e':e,'primitive_primes':primes,'classical_exception':bool(exceptional)})
    ck('Gaussian exception p3 and p7 irreducible',all(s.Poly(z*z+1,z,modulus=p).is_irreducible for p in [3,7]))
    ck('ninth-cyclotomic binary exception irreducible',s.Poly(s.cyclotomic_poly(9,z),z,modulus=2).is_irreducible and s.n_order(2,9)==6)
    # Exact punctured canonical volume and the fixed-mark determinant distinction.
    x,y,w=s.symbols('x y w',nonzero=True);coords=[x,y,w]
    lifted=s.Matrix([x**3,y**2,w/6]);density=1/(x*y)
    ck('punctured canonical map preserves actual top form',s.factor(lifted.jacobian(coords).det()*density.subs({x:x**3,y:y**2},simultaneous=True)/density)==1)
    swap=s.Matrix([y,x,-w]);ck('canonical lifted reflection preserves volume',s.factor(swap.jacobian(coords).det()*density.subs({x:y,y:x},simultaneous=True)/density)==1)
    ck('nontrivial determinant has no nonzero fixed covolume',s.solve(s.Eq(s.Symbol('q'),-s.Symbol('q')),s.Symbol('q'))==[0])
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'degree2048_model':{'base_ring':'Z[1/46]','conductor':23,'galois_automorphism_exponent':5,'source_dimension':44,'volume_source_dimension':45,'ideal_hermite_basis':serial(B),'holonomy_matrix':serial(RI),'semilinear_matrix':serial(M),'degree':2048,'smith_factors':smith,'residue_factor_mod2':str(poly),'geometric_group':'F2^11 semidirect C23','kernel_dimension':11,'source_ideal_nonprincipal':'Written quadratic norm/order-three argument; no class-number table.'},'primitive_divisor_controls':prime_controls,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact degree2048 ideal-lattice and semilinear map, image/Smith/residue controls, a small primitive-divisor grid, and actual punctured-volume identities. No formalization of the ambiguous class formula, all-degree primitive-divisor theorem, group cohomology, stable atomicity, spectrum or source geometry. No independent whole-proof replay, external review or publication-priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
