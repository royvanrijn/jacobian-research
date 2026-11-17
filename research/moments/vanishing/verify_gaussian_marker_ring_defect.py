"""One-use exact controls for written ring-defect proof. Never replay success.

Finite symbolic checks transcribe torsion/projector/guard formulas;
they do not establish the universal quotient or nonmembership theorem.
"""
from pathlib import Path
import hashlib,json,time
import sympy as S

root=Path(__file__).resolve().parent
output=root/'gaussian_marker_ring_defect_verification.json'
assert not output.exists(),'Successful checker sealed;never rerun.'
started=time.monotonic()
t,b,y,s,r,z=S.symbols('t b y s r z')
checks=[];packets=[]
def check(name,value):
    assert S.cancel(S.expand(value))==0,name
    checks.append(name)

for L in [1,2,5]:
    nu=(L+1)//2
    H=sum((-t*b*y)**k/S.factorial(k) for k in range(L))
    tail=(-1)**(L-1)*t**L*b**L*y**(L-1)/S.factorial(L-1)
    DYH=S.diff(H,y)-(1-2*t)*y*H+t*b*H
    check(f'torsion-remainder-L{L}',DYH+(1-2*t)*y*H-tail)
    A=sum((b*b*(2-s)/8)**k/S.factorial(k) for k in range(nu))
    F=sum((b*b*(1-s)**2/(8*s))**k/S.factorial(k) for k in range(nu))
    G=sum((b*b/(8*s))**k/S.factorial(k) for k in range(nu))
    product=S.Poly(S.expand(A*F),b)
    truncated=sum(c*b**powers[0] for powers,c in product.terms() if powers[0]<L)
    check(f'nilpotent-Gaussian-gauge-L{L}',truncated-G)
    projected=G
    for j in range(1,nu):
        # On h*f, E acts as h*(s*f'-f/2).
        projected=S.cancel((s*S.diff(projected,s)+j*projected)/j)
    check(f'Euler-projector-L{L}',projected-1)
    rel=b**L*(1-b*r);zz=b*r;e=zz**L
    geo=sum(zz**k for k in range(L))
    check(f'idempotent-witness-L{L}',e*e-e+r**L*rel*geo)
    assert S.expand(b*r*e-e+r**L*rel)==0
    assert S.expand(b**L*(1-e)-rel*geo)==0
    packets.append(dict(nilpotence_threshold=L,nu=nu,operator_order_bound=nu-1,
                        torsion_amplitude=str(H),torsion_remainder=str(tail),
                        Gaussian_gauge=str(A),Euler_projector=f'product_(j=1,...,{nu-1}) (E+j+1/2)/j',
                        idempotent=str(e),stabilization_guard=str(rel)+'=0'))

# Point-module normal form, with D multiplication by z and s acting -d/dz.
coeff=S.symbols('c0:7')
v=sum(coeff[k]*z**k for k in range(7))
def project(w):
    return S.expand(sum(z**j*(-1)**j*S.diff(w,z,j)/S.factorial(j) for j in range(7)))
check('point-module-projection',project(v)-coeff[0])
rebuilt=sum((-1)**j*z**j*project((-1)**j*S.diff(v,z,j))/S.factorial(j) for j in range(7))
check('point-module-reconstruction',rebuilt-v)

H=packets[-1]['torsion_amplitude'];HH=S.sympify(H)
assert S.expand((S.diff(HH,y)-(1-2*t)*y*HH)+(1-2*t)*y*HH-S.sympify(packets[-1]['torsion_remainder']))!=0
nu=3;bad=G
for j in range(1,nu):bad=S.cancel((s*S.diff(bad,s)+(S.Rational(-1,2)+j)*bad)/j)
assert S.cancel(bad-1)!=0
data=dict(date='2026-11-17',status='pass',scope='Bounded exact symbolic transcriptions,not independent verification of universal all-ring defect/nonmembership proof.',checks=checks,rejected_mutations=['Removed tb from source derivative','Removed half-residue correction from Euler projector'],packets=packets,point_model='D=z,s=-d/dz;finite polynomial vectors only,no false truncated Weyl representation.',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),elapsed_seconds=time.monotonic()-started,independent_replay=False,formal_verification=False,external_review=False,rerun='Successful checker sealed;never rerun.')
output.write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps(dict(status='pass',checks=len(checks),mutations=2,seconds=data['elapsed_seconds'])))
