"""Exact field and global-unit controls for KMFI77--80.

Universal geometric arguments are in LINEAR_PENCIL_OBSTRUCTION.md.
A finite group specialization is not used as a geometric-group proof.
"""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
import hashlib,json,platform,time
import sympy as sp
from linear_pencil import pencil_data,graph_pencil,homogenize,moving_field,monic_cofactor_family
from marked_wreath import polynomial_map
from coefficient_graph import passes_morse_gate
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
T,B,Q,R,Z,X,Y,A=sp.symbols("T B Q R Z X Y A")
pencils=[]
for n in range(4,11):
    # Exact T=B*Z normalization removes B from these homogeneous controls.
    for a,d,label,base in [(1,1,'neither_vanishes',0),(-1,0,'H_vanishes',0),
                           (0,-1,'J_vanishes',0),(-1,-1,'both_vanish',1)]:
        H=T**n+a*B**(n-2)*T*T
        J=T**(n-1)+d*B**(n-3)*T*T;L=T-B
        h=Z**n+a*Z*Z;j=Z**(n-1)+d*Z*Z;l=Z-1
        check(f"exact homogeneous field normalization H n{n},{label}",sp.cancel(H.subs(T,B*Z)/B**n-h)==0)
        check(f"exact homogeneous field normalization J n{n},{label}",sp.cancel(J.subs(T,B*Z)/B**(n-1)-j)==0)
        data=pencil_data(h,j,l,Z)
        check(f"common factor retained n{n},{label}",data['base_degree']==base)
        check(f"exact moving cover degree n{n},{label}",data['cover_degree']==1)
        check(f"exact moving image degree n{n},{label}",data['image_degree']==n-base)
        check(f"exact root count including fixed root n{n},{label}",data['moving_degree']+data['base_degree']==n)
        pencils.append({'n':n,'case':label,'base_degree':base,'cover_degree':data['cover_degree'],
          'moving_root_group':data['moving_root_group'],'degree_certificate':data['degree_certificate']})
for n in [5,7,11]:
    graphs={2:B*Q+B**3}
    H,J,L=graph_pencil(n,graphs);data=pencil_data(H,J,L,T,(B,))
    check(f"new prime pencil outside old strict weight gate n{n}",not passes_morse_gate(n,graphs))
    check(f"prime pencil with arbitrary offset has no base factor n{n}",data['base_degree']==0)
    check(f"prime degree certificate e1 with independent sections n{n}",
      data['cover_degree']==1 and data['image_degree']==n and
      data['degree_certificate']=='prime degree and three independent sections')
    pencils.append({'n':n,'case':'q_dependent_with_offset_outside_weight_gate','graphs':{2:str(graphs[2])},
      'base_degree':0,'cover_degree':1,'moving_root_group':data['moving_root_group'],
      'degree_certificate':data['degree_certificate']})
# A genuine imprimitive composite-degree control, not a missing gcd exception.
H=T**4;J=T**3+2*T*T;L=T*T/4+T+1
data=pencil_data(H,J,L,T)
check("D4 example no common factor",data['base_degree']==0)
check("D4 example degree2 map to conic",data['cover_degree']==2 and data['image_degree']==2)
phi=T*T/(T+2)
check("D4 exact composed pencil",sp.cancel((Q*H+R*J+L)/(T+2)**2-(Q*phi**2+R*phi+sp.Rational(1,4)))==0)
tau=-2*T/(T+2)
check("D4 block involution order2",sp.cancel(tau.subs(T,tau)-T)==0)
check("D4 block involution fixes phi",sp.cancel(phi.subs(T,tau)-phi)==0)
check("D4 block involution nonidentity",sp.cancel(tau-T)!=0)
u,v=sp.symbols("u v")
qrec=sp.Rational(1,4)/(u*v);rrec=-qrec*(u+v)
check("D4 top quadratic reconstructed",sp.cancel(qrec*u*u+rrec*u+sp.Rational(1,4))==0)
disc_u=u*(u+8);disc_v=v*(v+8)
check("D4 first Kummer valuation separates roots",sp.Poly(disc_u,u).terms()[-1][0][0]==1 and disc_v.subs(u,0)!=0)
check("D4 second Kummer valuation separates roots",sp.Poly(disc_v,v).terms()[-1][0][0]==1 and disc_u.subs(v,0)!=0)
from sympy.polys.numberfields import galois_group
arith,alt=galois_group(sp.Poly((H+L)*4,T))
check("D4 rational point group exact order8",arith.order()==8)
check("D4 rational point imprimitive transitive",arith.is_transitive() and not arith.is_primitive())
# Full homogeneous marking gradients on the original polynomial affine source.
for n in range(3,9):
    variables,TT,phi0,Fout=polynomial_map(n);x,y=variables[:2];p=1+x*y
    Ph=homogenize(phi0.as_expr(),TT,X,Z,n)
    PX=sp.diff(Ph,X).subs({X:x,Z:p},simultaneous=True)
    PZ=sp.diff(Ph,Z).subs({X:x,Z:p},simultaneous=True)
    check(f"entire homogeneous derivative X is p n{n}",sp.expand(PX-p)==0)
    check(f"entire homogeneous derivative Z is minus x n{n}",sp.expand(PZ+x)==0)
# A polynomial identity gives the cofactor inverse once the marked gradients
# and the chosen component equation are imposed.
x,y=sp.symbols("x y");p=1+x*y
for n in range(4,10):
    H=T**n-B**(n-2)*T*T;J=T**(n-1)-B**(n-3)*T*T;L=T-B;D=T-B
    moving=sp.cancel((Q*H+R*J+L)/D)
    Gh=homogenize(moving,T,X,Z,n-1);Dh=X-B*Z;Ph=sp.expand(Gh*Dh)
    sub={X:x,Z:p}
    G=Gh.subs(sub,simultaneous=True)
    gx=sp.diff(Gh,X).subs(sub,simultaneous=True);gz=sp.diff(Gh,Z).subs(sub,simultaneous=True)
    hx=sp.diff(Dh,X).subs(sub,simultaneous=True);hz=sp.diff(Dh,Z).subs(sub,simultaneous=True)
    px=sp.diff(Ph,X).subs(sub,simultaneous=True);pz=sp.diff(Ph,Z).subs(sub,simultaneous=True)
    unit=Dh.subs(sub,simultaneous=True)
    check(f"global polynomial cofactor inverse certificate n{n}",
      sp.expand(unit*(gx+y*gz)-1-((px-p)+y*(pz+x)-G*(hx+y*hz)))==0)
    mf=moving_field(H,J,L,T,D)
    check(f"moving marking has simple r-prime n{n}",mf['W']!=0 and sp.degree(mf['derivative'],R)==1)
    qval=mf['q']
    derivative=sp.diff(Q*H+R*J+L,T).subs(Q,qval)
    check(f"moving derivative field identity n{n}",sp.cancel(derivative-mf['derivative'])==0)
    m=n-2
    rval=sp.cancel((A**m/D-mf['V'])/mf['W'])
    check(f"moving marked rational field reconstruction n{n}",sp.cancel(mf['derivative'].subs(R,rval)-A**m)==0)
# Explicit polynomial p inverse modulo a monic root equation.
relative=[]
for d in range(2,12):
    f=T**d+T;P,D,C=monic_cofactor_family(f,T);n=d+3;m=d+1
    pp=sp.Poly(P,T)
    check(f"monic-factor graph linear coefficient1 d{d}",pp.nth(1)==1)
    check(f"monic-factor graph constant minus b d{d}",pp.nth(0)==-B)
    check(f"monic-factor graph top coefficients q,r d{d}",pp.nth(n)==Q and pp.nth(n-1)==R)
    check(f"monic static factor irreducible over k(b) d{d}",sp.degree(D,B)==1 and sp.diff(D,B)==-1)
    Dh=homogenize(D,T,X,Z,d).subs({X:x,Z:p},simultaneous=True)
    rr=sp.cancel((Dh-x**d)/p)
    ss=sp.Poly((Z-1)**d-(-1)**d,Z).exquo(sp.Poly(Z,Z)).as_expr().subs(Z,p)
    inv=(-1)**(d+1)*(rr*y**d+ss)
    check(f"explicit monic-component unit p certificate d{d}",
      sp.cancel(p*inv-1-(-1)**(d+1)*y**d*Dh)==0)
    fprime=sp.diff(f,T);c=sp.Poly(f,T).nth(d-1)
    qval=(A**m/fprime-1-R*T*T)/(T**3-c*T*T)
    actual=sp.diff(P,T).subs(B,f).subs(Q,qval)
    check(f"relative marked component rational field reconstruction d{d}",sp.cancel(actual-A**m)==0)
    check(f"relative marking simple q-prime d{d}",sp.diff(fprime*C,Q)!=0 and sp.degree(fprime*C,Q)==1)
    check(f"nonlinear monic curve has a nonunit derivative d{d}",sp.degree(fprime,T)==d-1 and d-1>0)
    relative.append({'root_degree':d,'ambient_n':n,'marking_order':m,'relative_degree':d*m,
      'curve_control':str(f),'relative_group':'C_m wr G_f',
      'field':'k(t,a,r), rational3','actual_source':'nonconstant unit p; not stably affine space'})
# Cofactor linear forms at distinct nonzero roots give independent Kummer
# divisors. Their coefficient pairs coincide only on the diagonal.
s,c=sp.symbols("s c")
g=sp.Poly(s*s-T*T,s,domain=sp.QQ.frac_field(T,c)).gcd(
  sp.Poly(s**3-c*s*s-T**3+c*T*T,s,domain=sp.QQ.frac_field(T,c))).monic()
check("relative cofactor linear forms distinct at generic distinct roots",sp.expand(g.as_expr()-(s-T))==0)
# Polynomial root frames preserve the marking gradient up to determinant.
# This guards the monic-factor obstruction against regular coordinate changes.
M=sp.Matrix([[1,B],[B,1+B*B]]);Minv=M.inv()
xp=M*sp.Matrix([x,p]);xp0,xp1=xp
check("parameter-dependent root frame determinant unit",sp.expand(M.det())==1)
check("root frame transforms homogeneous gradient",
  all(sp.expand(z)==0 for z in Minv.T*sp.Matrix([p,-x])-sp.Matrix([xp1,-xp0])))
bezout=(sp.Matrix([[-y,1]])*Minv)
check("root frame pair remains globally unimodular",
  sp.expand(bezout[0]*xp0+bezout[1]*xp1-1)==0)
U=sp.Symbol("U")
check("integral transformed root gives explicit denominator unit",
  sp.expand(xp1*(bezout[0]*U+bezout[1])-1-
            bezout[0]*(U*xp1-xp0))==0)
manifest=json.loads((ROOT/'sources/linear-pencil-manifest.json').read_text())
for src in manifest['sources']:
    path=ROOT/src['path']
    check("sectional primary source bytes",path.stat().st_size==src['bytes'])
    check("sectional primary source SHA256",hashlib.sha256(path.read_bytes()).hexdigest()==src['sha256'])
model={'claims':['KMFI77','KMFI78','KMFI79','KMFI80'],'pencil_controls':pencils,
 'composite_imprimitive_control':{'degree':4,'plane_cover_degree':2,'plane_image_degree':2,
   'geometric_root_group':'D4 of order8','rational_point_group_order':int(arith.order()),
   'coefficient_scope':'rational in b; tested polynomial q,r slice b=-1'},
 'relative_monic_factor_controls':relative,
 'universal_statements':{
 'primitive_moving_root_group_of_affine_linear_top_coefficient_graph':'always symmetric',
 'irreducible_prime_degree11_affine_linear_graph_PSL2_11_root_group':'excluded within this graph class',
 'no_offset_root_classification':'Sn, or fixed b plus S_(n-1) when H(b)=J(b)=0',
 'base_factor_moving_marked_component':'rational3, but a nonconstant cofactor unit forbids stable affine space',
 'monic_nonlinear_factor_actual_component':'nonconstant p unit forbids stable affine space',
 'polynomial_unimodular_root_frame_removes_monic_obstruction':False,
 'constant_nonzero_leading_coefficient_actual_affine_component':'only degree1 Keller maps can occur',
 'generic_imprimitive_pencil_full_wreath_equality_asserted':False,
 'relative_Cm_wr_Gf_group_claimed_as_absolute_Keller_occurrence':False},
 'probe_audit':{'dense_symbolic_n11_fiber_gcd':'interrupted explicitly after a confirmed-live session; no result or passed receipt inferred',
   'input_n':11,
   'input_graphs':{'2':'Q*(B+1)+R*(B**2-1)+B**3+2','3':'Q*(B-2)+R*B+1-B'},
   'method':'exact fiber gcd over Q(B,T), before adding the prime-degree shortcut',
   'replacement':'exact prime-degree and independent-section certificate, proved universally'},
 'exceptional_PSL2_11_degree11_absolute_occurrence_decided':False,
 'natural_PSL2_11_degree12_absolute_occurrence':'excluded by KMFI67',
 'full_action_classification_complete':False}
(ROOT/'linear_pencil_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['LINEAR_PENCIL_OBSTRUCTION.md','linear_pencil.py','verify_linear_pencil.py',
 'linear_pencil_model.json','sources/linear-pencil-manifest.json','sources/kadets-1809.07293v4.html']
receipt={'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),
 'sympy':sp.__version__,'runtime_seconds':round(time.monotonic()-started,3),
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written universal plane-map and source-unit proofs with attributed sectional-monodromy input and exact field, polynomial, Kummer and unit controls. No formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'linear_pencil_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'prime11_linear_pencil_root_route_excluded':True,'constant_leading_nontrivial_affine_source_excluded':True,
 'exceptional_PSL11_absolute_decided':False,'full_goal_complete':False}))
