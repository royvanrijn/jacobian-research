"""New exact identity controls; none certifies the written frame theorem."""
from pathlib import Path
import hashlib,json
import sympy as S

r,x,y,t,T,b,e=S.symbols('r x y t T b e', nonzero=True)
phi=r*(x+1/x);Y=x+1/x;Z=x-1/x
twisted=lambda q:S.expand(x*S.diff(q,x)+x*S.diff(phi,x)*q)
assert S.expand(r*Y+r*r*Y*Y-4*r*r-r*twisted(Z))==0
# The all-order logarithmic coefficient equation uses only H_n-H_(n-1)=1/n.
n,an,Hn=S.symbols('n an Hn',nonzero=True)
anm1=n*n*an;bn=-Hn*an;bnm1=-(Hn-1/n)*anm1
assert S.expand(n*n*bn-bnm1+n*an)==0
# Exact removed-component transport,including its scalar time term.
X=T-(1/t-b)/e;D=1-t*(b+e*T);kappa=-e
assert S.cancel(S.diff(X,T)-t*kappa*X/D)==0
w=1/D;R=S.diff(X,t)/S.diff(X,T)
assert S.cancel(S.diff(w,t)+w/t-S.diff(w*R,T))==0
# Central log volume and outside relative volume have different scaling.
central_phi=x*x/2-r*x
assert S.expand(central_phi.subs(x,y)-(y*y/2-r*y))==0
escaping_phi=r*r*x*x+x;center=-1/(2*r*r);delta=1/r
assert S.expand(escaping_phi.subs(x,center+delta*y)+1/(4*r*r)-y*y)==0
assert S.cancel(delta/center+2*r)==0
# A finite ramification multiplies,not kills,the nonzero logarithmic shear.
N=S.symbols('N',integer=True,positive=True)
assert S.Matrix([[1,N],[0,1]])!=S.eye(2)

p=Path(__file__);out=p.with_name('punctured_controls.json');assert not out.exists()
result=dict(date='2026-11-17',version='v004',status='passed',
    controls=['Universal two-ended seed ODE modulo spatial exactness',
      'All-order harmonic-number coefficient recurrence',
      'Removed-component spatial and parameter transport identities',
      'Central puncture and escaping-center scaling controls',
      'Nonzero logarithmic shear after symbolic positive ramification'],
    script_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
    not_proved=['Full punctured formal frame','Polynomial separation','Analytic logarithmic cut estimate','Full Gaussian boundary proof'],
    independent_replay=False,formal_verification=False,external_review=False)
out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(status='passed',universal_identity_groups=5,old_checks_replayed=False,full_frame_proved_by_computation=False)))
