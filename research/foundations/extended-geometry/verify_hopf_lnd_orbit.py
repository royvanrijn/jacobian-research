#!/usr/bin/env python3
"""Bounded exact controls/source integrity for the written Hopf LND theorem."""
import argparse,copy,hashlib,json,sys
from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
if not __debug__:
    raise SystemExit("Optimized Python is unsupported: exact assertions must execute.")
HERE=Path(__file__).resolve().parent
PN="HOPF_LND_ORBIT_DEGREE_PROPAGATION_AND_POWER_DESCENT.md"
CN="verify_hopf_lnd_orbit.py"
IN="hopf_lnd_orbit_inputs.json"
LN="hopf_lnd_orbit_source_ledger.json"
RN="hopf_lnd_orbit_verification.json"
PROOF_SHA="4c5aff7869bb933e4a3ef3b77d0321ab20525b2b0b375824ce400a31275a6032"
EXPECTED={'degree': {'all_k_formula_claim': False, 'formula': 'nu*k+c', 'meaning': 'least ordinary sphere filtration degree, not derivation height', 'only_nonzero_f': True, 'only_nu_positive': True, 'sufficient_threshold': 'k>=d-c+1'}, 'derivation': 'J=x partial_t-t partial_y', 'field': 'every characteristic-zero field and every coefficient-field extension', 'low_height': {'arbitrary_nonlinear_conjugates_claim': False, 'linear_orthogonal_conjugates_only': True, 'odd_height_uses_retained_full_quadratic_classification': True, 'odd_unsafe_nu_min': 3, 'unsafe_nu_min': 2}, 'mixed_box': {'all_odd_phases_have_a_height3_direction_claim': False, 'degree5_resolved': False, 'degree7_resolved': False, 'degree_bounds': ['degC<=d+1', 'degB<=d', 'degD<=d-1', 'degE<=d-2'], 'dimensions': {'5': 20, '7': 28}, 'nu_at_most': 3, 'odd_only_class_exhausts_mixed_box_claim': False, 'regularity': ['C0+D0=0', 'B0+E0=0'], 'shape': 'x^-1(C+tB+t^2D+t^3E)'}, 'moments': {'Kx_marker_unchanged': True, 'arbitrary_marker_preserves_every_even_order_claim': False, 'arbitrary_original_marker_first_balanced': True, 'marker_may_depend_on_k': True, 'one_fixed_marker_per_k': True, 'pure_all_positive_orders': True}, 'power_descent': {'all_extensions': True, 'both_bigrading_extremes_required': True, 'primitive_homogenization_required_for_iff': True, 'ring': 'FULL original four-variable polynomial ring', 'scalar_factors_allowed': True}, 'ring': 'K[x,y,t]/(t^2+2xy-1)', 'schema': 1, 'source_boundary': {'classical_LND_context': 'inspected context, not an imported logical premise', 'general_degree_and_power_proof': 'self-contained', 'old_even_or_prime_coverage': 'not a logical premise', 'quadratic_classification': 'logical only for odd-unsafe nu<=2 exclusion', 'source_authentication_is_not_a_finite_prefix_proof': True}, 'verification': {'bounded_exact_controls_only': True, 'external_review': False, 'finite_controls_prove_universal_quantifiers': False, 'formal_verification': False, 'independent_replay': False, 'worldwide_priority': False}, 'all_automorphisms': {'all_polynomial_A_automorphisms_preserve_Lambda': True, 'pure_safe_unsafe_transport': True, 'marker_transported_by_sigma': True, 'exact_form_vector_space_codimension': 1, 'exact_form_space_is_ideal_claim': False, 'pointwise_volume_form_preservation_claim': False, 'arbitrary_automorphism_preserves_ordinary_degree_claim': False, 'nu_k_plus_c_for_all_automorphisms_claim': False}}
RAW=[
 ("quadratic_classification","QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md","cc68360e2ae7d09d025e484433e0207df0498131ba2f2ee2bfd518f9bf9f86e4","LOGICAL_NU2_UNSAFE_EVEN_CLASSIFICATION_ONLY"),
 ("original_normalization","TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md","5dba378cfcde66b1ca7f862533c1ae3277cc2acafd9157c80bec2dd0063384a2","LOGICAL_RAW_IMAGE_KERNEL_AND_ORIGINAL_NORMALIZATION_ONLY"),
 ("formal_sphere_antecedent","TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md","3e3afcf39fbbf3ec1641843cc81b550248d15db531020f4635a2ea0f5759e825","AUTHENTICATED_RAW_MOMENTS_AND_MARKER_BRIDGE_ANTECEDENT"),
 ("graded_hopf_antecedent","HOPF_LIFT_CLASSIFICATION.md","1688857a80e9c09db5727bf66d8e59f85126d19a410db2dbef5c7929bd4f8a60","AUTHENTICATED_GRADED_DEGREE_CONTEXT_ONLY"),
 ("example_seed","HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md","e03807c520cd766507be7e12a225821390447524ac02099ca15ac2c106879056","AUTHENTICATED_EXAMPLE_SEED_ONLY_NOT_GENERAL_CORE"),
 ("odd_only_scope","ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md","51e8ccec1e80166c1c70f87e5cc83b0725a8d3bb04eb697931da18b11a06a826","AUTHENTICATED_SCOPE_BOUNDARY_ONLY"),
 ("odd_septic_scope","ODD_CUBIC_DEGREE_SEVEN_SAFETY.md","af8611c27bd832dedcbb77587fedd604d811c93f2ff79a210937bbf6ab03a958","AUTHENTICATED_SCOPE_BOUNDARY_ONLY"),
 ("odd_three_height_scope","THREE_HEIGHT_QUINTIC_IMAGE_SAFETY.md","ee811ec80c6d253e2fc3130fb3b74526f4b5d3bf7b7c99068dc3abcce962d3cb","AUTHENTICATED_SCOPE_BOUNDARY_ONLY"),
]
CLASSICAL=[{'author': 'Daniel Daigle', 'title': 'Polynomials f(X,Y,Z) of low LND-degree', 'url': 'https://aix1.uottawa.ca/~ddaigle/articles/lndeg.pdf', 'locator': 'Introduction printed p1; Section3.1(b) printed p4', 'scope': 'Derivation height additivity on characteristic-zero domains and local-slice localization. Context only; not ordinary sphere degree and not a new logical premise.', 'local_pdf_download_or_hash_claim': False}, {'author': 'Gene Freudenburg', 'title': 'Local slice constructions in k[X,Y,Z]', 'url': 'https://doi.org/10.18910/4752', 'pdf_url': 'https://ir.library.osaka-u.ac.jp/repo/ouka/all/4752/ojm34_04_02.pdf', 'locator': 'Osaka J. Math.34(1997), Proposition1 printed p757', 'scope': 'Integral characteristic-zero domain, LND D, Ds=1 imply R=(kerD)[s]. On A_x the slice for J is t/x. Context only; specialized argument is proved here.', 'local_pdf_download_or_hash_claim': False}, {'author': 'Frank Kutzschebauch and Matthias Leuenberger', 'title': 'Lie algebra generated by locally nilpotent derivations on Danielewski surfaces', 'url': 'https://arxiv.org/pdf/1311.1075v2', 'locator': 'Version2(24Nov2014), Section2 p3; Proposition2.1 pp3-4; Theorem3.9 proof p11', 'scope': 'Complex standard volume form and rank-one quadratic cohomology are geometric context only. Theorem2.7 requires degree>=3 and is NOT imported for the quadric. Algebraic all-field quotient argument is proved here.', 'local_pdf_download_or_hash_claim': False}]
FLAGS={"external_review":False,"independent_replay":False,"formal_verification":False,"worldwide_priority":False}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def strict_equal(a,b):
    return json.dumps(a,sort_keys=True,separators=(",",":"))==json.dumps(b,sort_keys=True,separators=(",",":"))
def require(ok,msg):
    if not ok:raise ValueError(msg)
def specs():
    out=[("main_proof","here",PN,"WRITTEN_UNIVERSAL_PROOF"),
         ("checker","here",CN,"BOUNDED_EXACT_CONTROLS_AND_SOURCE_AUTHENTICATION"),
         ("typed_inputs","here",IN,"EXACT_TYPED_SCOPE")]
    out += [(r,"repo","research/foundations/extended-geometry/"+p,b) for r,p,h,b in RAW]
    return out
def validate_input(d):
    require(strict_equal(d,EXPECTED),"Exact typed theorem/input scope mismatch")
def validate_ledger(d):
    require(set(d)=={"schema","sources","classical_primary_reading","assurance","logical_boundary"},"Source-ledger keys")
    require(type(d["schema"]) is int and d["schema"]==1,"Ledger schema")
    require(strict_equal(d["classical_primary_reading"],CLASSICAL),"Classical source scope changed")
    require(strict_equal(d["assurance"],FLAGS),"Assurance strengthening")
    require(strict_equal(d["logical_boundary"],EXPECTED["source_boundary"]),"Logical/evidence boundary")
    rows=d["sources"];sp=specs()
    require(type(rows) is list and len(rows)==len(sp),"Exact source cardinality")
    require(len({r.get("role") for r in rows if isinstance(r,dict)})==len(sp),"Duplicate/missing role")
    by={r["role"]:r for r in rows}
    for role,root,path,boundary in sp:
        r=by.get(role)
        require(type(r) is dict and set(r)=={"role","root","path","sha256","boundary"},"Source row schema")
        require(r["root"]==root and r["path"]==path and r["boundary"]==boundary,"Source path/root/scope")
        h=r["sha256"]
        require(type(h) is str and len(h)==64 and all(c in "0123456789abcdef" for c in h),"SHA format")
    return rows
def check_sources(repo):
    validate_input(load(HERE/IN));rows=validate_ledger(load(HERE/LN))
    pins={}
    hard={"main_proof":PROOF_SHA}
    hard.update({r:h for r,p,h,b in RAW})
    for r in rows:
        base=HERE if r["root"]=="here" else repo
        p=(base/r["path"]).resolve()
        try:p.relative_to(base)
        except ValueError:raise ValueError("Source escaped intended root")
        require(p.is_file(),"Missing mandatory source "+r["role"])
        actual=sha(p)
        require(actual==r["sha256"],"Bad current source hash "+r["role"])
        if r["role"] in hard:require(actual==hard[r["role"]],"Reviewed mathematical source changed "+r["role"])
        pins[r["role"]]=actual
    return rows,pins
# A small independent rational sparse-polynomial implementation.
def clean(p):return {m:Q(c) for m,c in p.items() if c}
def const(c,n):return {} if not c else {(0,)*n:Q(c)}
def var(i,n):
    a=[0]*n;a[i]=1;return {tuple(a):Q(1)}
def add(*ps):
    out={}
    for p in ps:
        for m,c in p.items():out[m]=out.get(m,Q(0))+c
    return clean(out)
def scale(p,c):return clean({m:v*c for m,v in p.items()})
def mul(a,b):
    out={}
    for m,v in a.items():
        for n,w in b.items():
            k=tuple(i+j for i,j in zip(m,n));out[k]=out.get(k,Q(0))+v*w
    return clean(out)
def power(p,n,nv=None):
    if nv is None:nv=len(next(iter(p))) if p else 3
    out=const(1,nv)
    for j in range(n):out=mul(out,p)
    return out
def deriv(p,i):
    out={}
    for m,v in p.items():
        if m[i]:
            k=list(m);k[i]-=1;out[tuple(k)]=v*m[i]
    return clean(out)
def sub(p,images):
    nv=len(next(iter(images[0]))) if images[0] else len(next(iter(images[1])))
    out={}
    for m,c in p.items():
        term=const(c,nv)
        for i,e in enumerate(m):term=mul(term,power(images[i],e,nv))
        out=add(out,term)
    return out
def nf(p):
    # Independently use xy=(1-t^2)/2; canonical monomials contain no xy.
    out={}
    for m,c in p.items():
        a,b,l=m;u=min(a,b)
        for r in range(u+1):
            k=(a-u,b-u,l+2*r)
            out[k]=out.get(k,Q(0))+c*Q(comb(u,r)*(-1)**r,2**u)
    return clean(out)
X,Y,T=[var(i,3) for i in range(3)]
ONE=const(1,3)
def J(p):return nf(add(mul(X,deriv(p,2)),scale(mul(T,deriv(p,1)),-1)))
def chart(p):
    out={}
    for (a,b,l),c in p.items():
        for r in range(b+1):
            k=(a-b,l+2*r);out[k]=out.get(k,Q(0))+c*Q(comb(b,r)*(-1)**r,2**b)
    return clean(out)
def lam(p):
    return sum((v/Q(l+1) for (a,l),v in chart(p).items() if a==0 and l%2==0),Q(0))
def degree(p):return max((sum(m) for m in nf(p)),default=-1)
def height(p):
    p=nf(p);nu=0
    while J(p):p=J(p);nu+=1
    return nu,p
def flow(p,k):
    xx=power(X,k);hh=power(X,k+1)
    yy=add(Y,scale(mul(T,xx),-1),scale(power(X,2*k+1),Q(-1,2)))
    return nf(sub(p,[X,yy,add(T,hh)]))
def exponential(p,k):
    out={};v=nf(p);j=0
    while v:
        out=add(out,scale(mul(power(X,k*j),v),Q(1,factorial(j))))
        v=J(v);j+=1
    return nf(out)
def digest(rows):return hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(",",":")).encode()).hexdigest()
def literal_lift(p,D):
    # Direct ORIGINAL four-variable expansion, no intermediate quotient.
    # xi1 exponent i, z1 exponent j encodes the balanced degree-D monomial.
    out={}
    for (a,b,l),c in nf(p).items():
        rr=D-a-b-l;require(rr>=0,"Lift degree too small")
        for i in range(l+1):
            for j in range(rr+1):
                ij=(a+i+j,b+i+j)
                out[ij]=out.get(ij,Q(0))+c*2**b*comb(l,i)*(-1)**(l-i)*comb(rr,j)
    return clean(out)
def balmul(a,b):
    out={}
    for (i,j),v in a.items():
        for (k,l),w in b.items():
            ij=(i+k,j+l);out[ij]=out.get(ij,Q(0))+v*w
    return clean(out)
def balpow(p,n):
    out={(0,0):Q(1)}
    for j in range(n):out=balmul(out,p)
    return out
def contract(p,D):return sum((v*factorial(i)*factorial(D-i) for (i,j),v in p.items() if i==j),Q(0))
def marker_z(p):return {(i+1,j):v for (i,j),v in p.items()}
def formal_identities():
    n=5;x,y,t,h,g=[var(i,n) for i in range(n)]
    q=add(power(t,2),scale(mul(x,y),2));one=const(1,n)
    def phi(a):return [x,add(y,scale(mul(t,a),-1),scale(mul(x,power(a,2)),Q(-1,2))),add(t,mul(x,a)),h,g]
    require(sub(q,phi(h))==q,"Universal quadric invariance")
    ph=phi(h);neg=phi(scale(h,-1))
    for i in range(3):require(sub(ph[i],neg)==[x,y,t][i],"Universal inverse flow")
    both=[sub(ph[i],phi(g)) for i in range(3)]
    for i in range(3):require(both[i]==phi(add(h,g))[i],"Universal additive flow")
    # Generic cubic regular representative t*B-2*t*y*D: J kills B,D.
    xx,yy,tt,B,D=[var(i,5) for i in range(5)]
    generic=add(mul(tt,B),scale(mul(mul(tt,yy),D),-2))
    def jraw(a):return add(mul(xx,deriv(a,2)),scale(mul(tt,deriv(a,1)),-1))
    require(jraw(jraw(jraw(generic)))==scale(mul(power(xx,2),D),6),"Universal cubic last derivative")
    # Symbolic exponent-gap identity for arbitrary c,nu,k,d,j.
    c,nu,k,d,j=[var(i,5) for i in range(5)]
    gap=add(c,mul(nu,k),scale(d,-1),scale(mul(j,k),-1))
    rhs=add(c,scale(d,-1),mul(k,add(nu,scale(j,-1))))
    require(gap==rhs,"Universal top-degree gap")
    # Universal first two height-one moments, with formal a,b.
    a,b,t=[var(i,3) for i in range(3)]
    ft=add(mul(a,t),b)
    def avg_t(v):
        o={}
        for (aa,bb,ll),cc in v.items():
            if ll%2==0:o[(aa,bb,0)]=o.get((aa,bb,0),Q(0))+cc/Q(ll+1)
        return clean(o)
    require(avg_t(ft)==b,"First height-one moment")
    require(avg_t(power(ft,2))==add(power(b,2),scale(power(a,2),Q(1,3))),"Second height-one moment")
    return ["quadric preservation with formal shear parameters","inverse shear on all generators","additive shear law on all generators","symbolic top-degree gap identity","formal first height-one moment","formal second height-one moment","generic regular cubic last derivative with invariant coefficient symbols"]

def top_form_controls():
    # Vectors use (dx^dy, dx^dt, dy^dt); dq wedge generators are explicit.
    zero={}
    omega=[scale(T,-1),scale(Y,2),zero]
    rx=[scale(X,-2),scale(T,-2),zero]
    ry=[scale(Y,2),zero,scale(T,-2)]
    rt=[zero,scale(Y,2),scale(X,2)]
    def vmul(v,p):return [mul(c,p) for c in v]
    def vadd(*vs):return [add(*(v[i] for v in vs)) for i in range(3)]
    def vanish(v):return all(not nf(c) for c in v)
    require(vanish(vadd(vmul(omega,X),[zero,scale(ONE,-1),zero],vmul(rx,scale(T,Q(-1,2))))),"Global x-omega identity")
    require(vanish(vadd(vmul(omega,T),[ONE,zero,zero],vmul(rx,Y))),"Global t-omega identity")
    require(vanish(vadd(vmul(omega,Y),[zero,zero,ONE],vmul(ry,scale(T,Q(1,2))),vmul(rt,scale(Y,-1)))),"Global y-omega identity")
    def coeff(v):return nf(add(scale(mul(T,v[0]),-1),mul(X,v[1]),scale(mul(Y,v[2]),-1)))
    def Jy(p):return nf(add(mul(Y,deriv(p,2)),scale(mul(T,deriv(p,0)),-1)))
    def H(p):return nf(add(mul(X,deriv(p,0)),scale(mul(Y,deriv(p,1)),-1)))
    rows=[]
    for a in range(5):
        for b in range(5):
            for l in range(5):
                if a+b+l>4:continue
                f={(a,b,l):Q(1)};fx,fy,ft=[deriv(f,i) for i in range(3)]
                require(coeff([scale(fy,-1),scale(ft,-1),zero])==scale(J(f),-1),"Exact d(f dx) coefficient")
                require(coeff([fx,zero,scale(ft,-1)])==Jy(f),"Exact d(f dy) coefficient")
                require(coeff([zero,fx,fy])==H(f),"Exact d(f dt) coefficient")
                require(lam(J(f))==lam(Jy(f))==lam(H(f))==0,"Exact coefficients not in ker Lambda")
                rows.append([a,b,l])
    # The d(f dx/dy/dt) coefficient identities with INDEPENDENT formal partials.
    xx,yy,tt,fx,fy,ft=[var(i,6) for i in range(6)]
    def wedge_form(a,b):
        return [add(mul(a[0],b[1]),scale(mul(a[1],b[0]),-1)),
                add(mul(a[0],b[2]),scale(mul(a[2],b[0]),-1)),
                add(mul(a[1],b[2]),scale(mul(a[2],b[1]),-1))]
    def to_omega(v):
        return add(scale(mul(tt,v[0]),-1),mul(xx,v[1]),scale(mul(yy,v[2]),-1))
    formal_partials=[fx,fy,ft]
    wants=[scale(add(mul(xx,ft),scale(mul(tt,fy),-1)),-1),
           add(mul(yy,ft),scale(mul(tt,fx),-1)),
           add(mul(xx,fx),scale(mul(yy,fy),-1))]
    for i in range(3):
        basis=[const(int(i==j),6) for j in range(3)]
        require(to_omega(wedge_form(formal_partials,basis))==wants[i],"Universal exact differential sign/coefficient")
    # A single polynomial identity in a FORMAL j parameter gives (O33).
    x,y,t,j=[var(i,4) for i in range(4)]
    q=add(power(t,2),scale(mul(x,y),2))
    left=add(scale(add(mul(j,mul(x,y)),scale(power(t,2),-1)),2),
             scale(add(j,scale(mul(add(j,scale(const(1,4),2)),power(t,2)),-1)),-1))
    require(left==mul(j,add(q,scale(const(1,4),-1))),"Universal J(yt^j) reduction")
    ev=Q(1);odd=Q(0)
    for r in range(1,21):
        ev*=Q(2*r-1,2*r+1);require(ev==Q(1,2*r+1),"Quotient moment recurrence")
    # V really is only a VECTOR space: t is exact but t^2 is not.
    require(J(Y)==scale(T,-1) and lam(T)==0 and lam(power(T,2))==Q(1,3),"V-is-an-ideal countercontrol")
    def wedge(a,b):
        return [add(mul(a[0],b[1]),scale(mul(a[1],b[0]),-1)),
                add(mul(a[0],b[2]),scale(mul(a[2],b[0]),-1)),
                add(mul(a[1],b[2]),scale(mul(a[2],b[1]),-1))]
    autos=[([X,Y,scale(T,-1)],Q(-1)),([Y,X,T],Q(-1))]
    for k in [0,1,2]:
        autos.append(([X,add(Y,scale(mul(T,power(X,k)),-1),scale(power(X,2*k+1),Q(-1,2))),add(T,power(X,k+1))],Q(1)))
    for images,want in autos:
        xx,yy,tt=images
        dx=[deriv(xx,i) for i in range(3)];dy=[deriv(yy,i) for i in range(3)];dt=[deriv(tt,i) for i in range(3)]
        pull=vadd(vmul(wedge(dx,dt),scale(yy,2)),vmul(wedge(dx,dy),scale(tt,-1)))
        require(coeff(pull)==const(want,3),"Automorphism volume multiplier")
        for fexp in [(0,0,0),(1,1,2),(2,0,3),(0,2,2)]:
            f={fexp:Q(1)}
            require(lam(sub(f,images))==lam(f),"Automorphism Lambda control")
    return {"global_form_identity_families":7,"exact_differential_and_kernel_controls":len(rows),
            "quotient_even_moment_controls":20,"literal_automorphism_volume_controls":len(autos),
            "volume_pointwise_preservation_counterexample":"t-reflection multiplies omega by -1 while preserving Lambda",
            "nonideal_counterexample":"t belongs to V, but Lambda(t^2)=1/3",
            "scope":"These controls support the universal written quotient argument; finitely many automorphisms do not establish the all-Aut theorem."}

def exact_controls():
    basis=[]
    for a in range(8):
        for b in range(8):
            for l in range(8):
                if a*b==0 and a+b+l<=7:basis.append({(a,b,l):Q(1)})
    flow_rows=[];pole_rows=[];degree_rows=[];box_rows=[];moment_rows=[];proper_rows=[]
    for p in basis:
        nu,last=height(p);ch=chart(p)
        require(nu==max(l for (a,l) in ch),"J height differs from Laurent t degree")
        require(all(b==l==0 for a,b,l in last),"Last derivative not in K[x]")
        c=degree(last);d=degree(p);K=d-c+1
        N=max(0,-min(a for a,l in ch))
        require(2*N<=nu,"General regularity pole bound")
        pole_rows.append([d,nu,N])
        if nu:
            for kk in sorted({K,K+1}):
                ff=flow(p,kk);ee=exponential(p,kk)
                require(ff==ee,"Literal/exponential shear mismatch")
                D=nu*kk+c;require(degree(ff)==D,"Eventual degree formula")
                top={m:v for m,v in ff.items() if sum(m)==D}
                lc=last.get((c,0,0));require(top=={(D,0,0):lc/Q(factorial(nu))},"Unique pure-x top")
                flow_rows.append([d,nu,c,kk,D,str(lc)])
        if nu<=3:
            num={(a+1,l):v for (a,l),v in ch.items()}
            require(all(a>=0 and l<=3 for a,l in num),"kerJ4 not in regular winding-one box")
            coeff=[[num.get((0,l),Q(0)) for l in range(4)]][0]
            require(coeff[0]+coeff[2]==0 and coeff[1]+coeff[3]==0,"Mixed regularity constants")
            bounds=[max((a for a,ll in num if ll==l),default=-10) for l in range(4)]
            D=max(bounds[l]+l-1 for l in range(4))
            require(D==d,"Exact mixed-box degree")
            box_rows.append([d,nu,bounds])
    # Degree additivity on the full nonconstant monomial normal-form basis through4.
    small=[p for p in basis if degree(p)<=4]
    for a in small:
        for b in small:
            require(degree(mul(a,b))==degree(a)+degree(b),"Ordinary degree additivity")
            degree_rows.append([degree(a),degree(b)])
    # Literal Delta_q expression of the sphere functional.
    def lap(v):return add(scale(deriv(deriv(v,0),1),2),deriv(deriv(v,2),2))
    delta_count=0
    for a in range(11):
        for b in range(11-a):
            for l in range(11-a-b):
                total=a+b+l;p={(a,b,l):Q(1)}
                if total%2:want=Q(0)
                else:
                    r=total//2;u=p
                    for jj in range(r):u=lap(u)
                    denom=2**r*factorial(r)
                    for jj in range(1,2*r+2,2):denom*=jj
                    want=u.get((0,0,0),Q(0))/denom
                require(lam(p)==want,"Formal Delta_q normalization")
                delta_count+=1
    # Full mixed-box dimensions by independent normal-form counting.
    dimensions={}
    for d in [3,4,5,7,9]:
        count=sum(1 for p in basis if degree(p)<=d and height(p)[0]<=3) if d<=7 else None
        direct=sum(d-l+1 for l in range(4))+2
        require(direct==4*d,"Mixed box parameter count")
        if count is not None:require(count==direct,"Mixed-box NF dimension")
        dimensions[str(d)]=direct
    # Source-independent polynomial examples, not a purity test for general f.
    quartic=nf(add(mul(add(ONE,X),Y),scale(mul(power(T,2),mul(add(ONE,X),add(scale(ONE,2),X))),Q(-1,2))))
    qnu,qlast=height(quartic);require(qnu==2 and degree(qlast)==4,"Quartic nu,c")
    for kk,maxm in [(1,3),(2,2)]:
        fk=flow(quartic,kk);D=degree(fk);lift=literal_lift(fk,D)
        require(lift.get((D,0))==Q(-1,2),"Original radial point")
        p=ONE;raw={(0,0):Q(1)}
        for m in range(1,maxm+1):
            p=nf(mul(p,fk));raw=balmul(raw,lift)
            pure=contract(raw,m*D);marked=contract(marker_z(raw),m*D+1)
            require(pure==factorial(m*D+1)*lam(p)==0,"Original pure contraction")
            require(marked==factorial(m*D+2)*lam(mul(X,p)),"Original marked normalization")
            require(lam(mul(X,p))==lam(mul(X,power(quartic,m))),"Fixed-x moment transport")
            require(marked!=0,"Quartic marker disappeared")
            moment_rows.append([kk,m,D,str(marked)])
    # Power ascent/descent controls in the FULL original ring.
    roots=[add(X,T),add(Y,T,ONE),add(power(X,2),mul(Y,T)),quartic]
    for g0 in roots:
        e=degree(g0);G=literal_lift(g0,e)
        for r in [2,3]:
            f=nf(scale(power(g0,r),Q(3,2)));D=degree(f)
            require(D==r*e,"Power degree")
            F=literal_lift(f,D)
            require(F==scale(balpow(G,r),Q(3,2)),"Original scalar-power ascent")
            require(all(0<=i<=e and 0<=j<=e for i,j in G),"Balanced root grading")
            proper_rows.append([e,r,D,len(G),len(F)])
    # Arbitrary-marker balancing checked directly monomialwise.
    marker_count=0
    for u in range(4):
        for v in range(u,5):
            L=v-u
            for aa in range(u+1):
                for bb in range(v+1):
                    # H=xi1^aa xi2^(u-aa) z1^bb z2^(v-bb).
                    if aa>bb or u-aa>v-bb:continue
                    coef=Q(factorial(bb)*factorial(v-bb),factorial(bb-aa)*factorial(v-bb-u+aa))
                    beta=(bb-aa,v-bb-u+aa)
                    require(sum(beta)==L,"Balanced marker degree")
                    full=factorial(bb)*factorial(v-bb)
                    require(coef*factorial(beta[0])*factorial(beta[1])==full,"Fixed coefficient derivative")
                    marker_count+=1
    # Explicit mathematical countercontrols to overstrong statements.
    bad=nf(add(power(X,8),Y));nu,last=height(bad)
    require(nu==2 and degree(last)==1 and degree(flow(bad,1))==8!=3,"Missing-threshold countercontrol")
    threshold_example=nf(add(power(X,3),T));nn,pp=height(threshold_example)
    require(nn==1 and degree(pp)==1 and degree(threshold_example)==3 and degree(flow(threshold_example,1))==3!=2,"Simple missing-threshold countercontrol")
    # f=x^2 is a proper power, but its padded original lift R*Z^2 is not:
    # on xi1=xi2=z2=1,z1=z the literal slice is z+1, with a simple zero.
    padded=balmul({(1,1):Q(1),(0,0):Q(1)},{(2,0):Q(1)})
    slice_coeff={}
    for (i,j),v in padded.items():slice_coeff[j]=slice_coeff.get(j,Q(0))+v
    require(slice_coeff=={1:Q(1),0:Q(1)},"Original padded-power counterexample slice")
    require(sum(v*(-1)**j for j,v in slice_coeff.items())==0 and sum(j*v*(-1)**(j-1) for j,v in slice_coeff.items() if j)==1,"Original simple zero")
    require(degree(Y)==1 and degree(flow(Y,1))==3,"Nonlinear degree-invariance countercontrol")
    # Actually expand H=xi1*(1+z1) in the FULL four-variable original ring.
    xx1,xx2,zz1,zz2=[var(i,4) for i in range(4)]
    original_H=mul(xx1,add(const(1,4),zz1))
    original_square=power(original_H,2,4)
    xi_degrees={a+b for a,b,c,d in original_square}
    z_degrees={c+d for a,b,c,d in original_square}
    require(xi_degrees=={2} and z_degrees=={0,1,2},"Both-grading countercontrol")
    return {"normal_form_basis":len(basis),"formal_degree_additivity_controls":len(degree_rows),
      "literal_vs_exponential_and_eventual_degree_controls":len(flow_rows),
      "general_pole_bound_controls":len(pole_rows),"mixed_box_controls":len(box_rows),
      "laplacian_normalization_controls":delta_count,"mixed_box_dimensions":dimensions,
      "original_four_variable_factorial_controls":moment_rows,"original_scalar_power_controls":proper_rows,
      "fixed_original_marker_grading_controls":marker_count,
      "control_stream_sha256":digest([flow_rows,pole_rows,degree_rows,box_rows]),
      "mathematical_negative_controls":["eventual degree formula fails without its k threshold","primitive hypothesis is necessary for power iff: RZ^2 has original simple slice zero","nonlinear shear does not preserve ordinary degree","one bihomogeneity check does not replace BOTH grading arguments","f=x^3+t at k1 has degree3 instead of nu*k+c=2 below its threshold"]}
def source_scope_controls(rows):
    count=0
    def rejects(fn):
        nonlocal count
        try:fn()
        except (ValueError,TypeError,KeyError,IndexError):count+=1;return
        raise AssertionError("A corrupt scope/source schema was accepted")
    for path,value in [
      (("degree","all_k_formula_claim"),True),
      (("power_descent","primitive_homogenization_required_for_iff"),False),
      (("power_descent","ring"),"selected invariant subring only"),
      (("power_descent","all_extensions"),False),
      (("moments","one_fixed_marker_per_k"),False),
      (("moments","arbitrary_marker_preserves_every_even_order_claim"),True),
      (("low_height","unsafe_nu_min"),3),
      (("low_height","odd_unsafe_nu_min"),2),
      (("low_height","arbitrary_nonlinear_conjugates_claim"),True),
      (("mixed_box","all_odd_phases_have_a_height3_direction_claim"),True),
      (("mixed_box","degree5_resolved"),True),
      (("mixed_box","degree7_resolved"),True),
      (("mixed_box","odd_only_class_exhausts_mixed_box_claim"),True),
      (("verification","finite_controls_prove_universal_quantifiers"),True),
      (("verification","formal_verification"),True),
      (("verification","external_review"),True),
      (("all_automorphisms","pointwise_volume_form_preservation_claim"),True),
      (("all_automorphisms","exact_form_space_is_ideal_claim"),True),
      (("all_automorphisms","nu_k_plus_c_for_all_automorphisms_claim"),True),
      (("all_automorphisms","arbitrary_automorphism_preserves_ordinary_degree_claim"),True),
      (("schema",),True),
    ]:
        a=copy.deepcopy(EXPECTED);q=a
        for k in path[:-1]:q=q[k]
        q[path[-1]]=value;rejects(lambda a=a:validate_input(a))
    template={"schema":1,"sources":rows,"classical_primary_reading":CLASSICAL,"assurance":FLAGS,"logical_boundary":EXPECTED["source_boundary"]}
    for role,root,path,b in specs():
        a=copy.deepcopy(template);a["sources"]=[r for r in a["sources"] if r["role"]!=role]
        rejects(lambda a=a:validate_ledger(a))
    for mutate in [
      lambda a:a.update(schema=2),
      lambda a:a.update(schema=True),
      lambda a:a.update(sources=[]),
      lambda a:a["sources"].append(a["sources"][0]),
      lambda a:a["sources"][0].update(root="outside"),
      lambda a:a["sources"][0].update(path="../bad"),
      lambda a:a["sources"][0].update(boundary="FINITE_TESTS_PROVE_UNIVERSAL_THEOREM"),
      lambda a:a["assurance"].update(formal_verification=True),
      lambda a:a["logical_boundary"].update(source_authentication_is_not_a_finite_prefix_proof=1),
    ]:
        a=copy.deepcopy(template);mutate(a);rejects(lambda a=a:validate_ledger(a))
    return count
def run(repo):
    rows,pins=check_sources(repo)
    ids=formal_identities();controls=exact_controls();forms=top_form_controls()
    report={"schema":1,"status":"PASS_SOURCE_GUARDED_HOPF_LND_ORBIT_AND_POWER_DESCENT",
      "scope":EXPECTED,"proof_sha256":pins["main_proof"],"checker_sha256":pins["checker"],
      "inputs_sha256":pins["typed_inputs"],"source_ledger_sha256":sha(HERE/LN),
      "mandatory_source_roles":len(rows),"pins":rows,"classical_primary_reading":CLASSICAL,
      "universal_polynomial_identity_families":ids,"exact_controls":controls,"all_automorphism_top_form_controls":forms,
      "source_scope_negative_controls":source_scope_controls(rows),"assurance":FLAGS,
      "verification_boundary":"Written proof supplies all fields/parameters/orders and power iff. Finite exact controls do not prove those quantifiers. Raw quadratic theorem is logical ONLY for odd-unsafe low-height exclusion; primary context is not a logical premise. No child checker or old finite-prefix result is imported."}
    return report
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo-root",required=True,type=Path)
    p.add_argument("--output",type=Path,default=HERE/RN)
    a=p.parse_args();repo=a.repo_root.resolve()
    require(repo.is_dir(),"Repository root does not exist")
    report=run(repo)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"mandatory_source_roles":report["mandatory_source_roles"],
      "universal_polynomial_identity_families":len(report["universal_polynomial_identity_families"]),
      "source_scope_negative_controls":report["source_scope_negative_controls"],
      "receipt_sha256":sha(a.output)},sort_keys=True))
if __name__=="__main__":
    try:main()
    except (ValueError,KeyError,TypeError,AssertionError) as exc:
        raise SystemExit("FAILED: "+str(exc))

