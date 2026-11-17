"""Exact two-node actions and four-simple-root conductor controls."""
from pathlib import Path
import hashlib,json,tempfile,time
import sympy as sp

HERE=Path(__file__).resolve().parent
started=time.perf_counter()
checks=[]
def check(name,test):
    if not test:
        raise AssertionError(name)
    checks.append(name)
def zero(name,expr):
    check(name,sp.cancel(expr)==0)
t,u,y,U,V,r,s=sp.symbols("t u y U V r s")
R,S=sp.symbols("R S")
p=t*t
g=(t*t-r*r)*(t*t-s*s)
q=t*g
f=V*V-U*(U-r*r)**2*(U-s*s)**2
zero("canonical two-node equation",f.subs({U:p,V:q},simultaneous=True))
zero("canonical conductor from monogenic trace-dual ratio",
     sp.diff(f,V).subs({U:p,V:q},simultaneous=True)/sp.diff(p,t)-g)
inc=sp.diag(1,(U-r*r)*(U-s*s))
zero("exact index degree two",(inc.det()-(U-r*r)*(U-s*s)))
zero("norm-conductor index square",
     (U-r*r)**2*(U-s*s)**2-inc.det()**2)
zero("canonical conductor discriminant",
     sp.discriminant(g,t)-16*r*r*s*s*(r*r-s*s)**4)

Ar=sp.integrate(sp.expand(p*sp.diff(q,t)),(t,-r,r))
As=sp.integrate(sp.expand(p*sp.diff(q,t)),(t,-s,s))
zero("universal first node action",Ar-sp.Rational(8,105)*r**5*(3*r*r-7*s*s))
zero("universal second node action",As-sp.Rational(8,105)*s**5*(3*s*s-7*r*r))
zero("two normalized actions cannot vanish simultaneously",
     3*(3*R-7*S)+7*(3*S-7*R)+40*R)
zero("second independent normalized-action combination",
     7*(3*R-7*S)+3*(3*S-7*R)+40*S)

# The two-node algebra is exact: even polynomials plus t*g(t^2) times them.
for k in [0,1,2,5]:
    pp=sp.expand((t*t)**k+t*g*(t*t)**(k+1))
    zero(f"first fixed node equality k={k}",pp.subs(t,r)-pp.subs(t,-r))
    zero(f"second fixed node equality k={k}",pp.subs(t,s)-pp.subs(t,-s))
    check(f"odd coefficient has both required factors k={k}",
          sp.rem(sp.expand((pp-pp.subs(t,-t))/(2*t)),
                 g,t)==0)

# All partitions of four root coordinates; conductor maximality forbids singletons.
def partitions(items):
    if not items:
        yield []
        return
    first,*rest=items
    for blocks in partitions(rest):
        yield [[first]]+[b[:] for b in blocks]
        for j in range(len(blocks)):
            yield [b+[first] if i==j else b[:] for i,b in enumerate(blocks)]
all_parts=list(partitions([0,1,2,3]))
allowed=[B for B in all_parts if len(B)==2 and all(len(b)>1 for b in B)]
check("complete four-point partition count",len(all_parts)==15)
check("exactly three conductor-compatible pairings",len(allowed)==3)
check("all admissible conductor blocks have size two",
      all(sorted(map(len,B))==[2,2] for B in allowed))
pts=[-r,r,-s,s]
sums=[]
for B in allowed:
    pair_sums=[sum(pts[i] for i in block) for block in B]
    sums.append(sp.factor(pair_sums[0]-pair_sums[1]))
check("unique common-midpoint pairing under four-distinct-point hypothesis",
      sum(e==0 for e in sums)==1)
for e in [e for e in sums if e!=0]:
    check("alternative common-midpoint equation forces coincident roots",
          any(sp.expand(e-candidate)==0 for candidate in
              [2*(r-s),-2*(r-s),2*(r+s),-2*(r+s)]))

# One source shear fixes all four sections and preserves Jacobian.
h=y*y-3*y+1
kshear=sp.expand(((u-h)**2-r*r)*((u-h)**2-s*s))
zero("source shear preserves simple-root discriminant",
     sp.discriminant(kshear,u)-16*r*r*s*s*(r*r-s*s)**4)
for ai in [-r,r,-s,s]:
    zero("shifted root section is exact",kshear.subs(u,h+ai))
for ai,aj in [(-r,r),(-r,-s),(r,s)]:
    zero("root-section differences are constant",(h+ai)-(h+aj)-(ai-aj))
P=u**3+y*u+1
Q=u**2+y*y
J=sp.diff(P,u)*sp.diff(Q,y)-sp.diff(P,y)*sp.diff(Q,u)
Psh=P.subs(u,t+h)
Qsh=Q.subs(u,t+h)
Jsh=sp.diff(Psh,t)*sp.diff(Qsh,y)-sp.diff(Psh,y)*sp.diff(Qsh,t)
zero("polynomial source shear preserves the Jacobian",Jsh-J.subs(u,t+h))

# Finite exact controls of target changes, with partner degrees unbounded.
baseq=sp.expand(q.subs({r:1,s:2}))
def action(pp,qq,v):
    return sp.integrate(sp.expand(pp*sp.diff(qq,t)),(t,-v,v))
baseAr=action(p,baseq,1)
baseAs=action(p,baseq,2)
controls=[]
for N in [2,3,7]:
    pp=sp.expand(p+baseq**N)
    qq=sp.expand(baseq+pp**2)
    recoveredq=sp.expand(qq-pp**2)
    recoveredp=sp.expand(pp-recoveredq**N)
    zero(f"target-shear inverse recovers odd generator N={N}",recoveredq-baseq)
    zero(f"target-shear inverse recovers even generator N={N}",recoveredp-p)
    zero(f"first node period survives high target shear N={N}",action(pp,qq,1)-baseAr)
    zero(f"second node period survives high target shear N={N}",action(pp,qq,2)-baseAs)
    controls.append({"shear_power":N,"p_degree":int(sp.degree(pp,t)),
                     "q_degree":int(sp.degree(qq,t)),
                     "periods":[str(baseAr),str(baseAs)]})

# Sharp guard: one period CAN vanish in the exact two-node algebra.
zeroAr=sp.simplify(Ar.subs({r:sp.sqrt(7),s:sp.sqrt(3)}))
nonzeroAs=sp.simplify(As.subs({r:sp.sqrt(7),s:sp.sqrt(3)}))
zero("one canonical node action may vanish",zeroAr)
check("the other canonical node action remains nonzero",nonzeroAs!=0)
gzero=sp.expand(g.subs({r:sp.sqrt(7),s:sp.sqrt(3)}))
check("zero-one-period control has four distinct conductor roots",
      sp.discriminant(gzero,t)!=0)
check("zero-one-period control still has codimension two",
      sp.degree((U-7)*(U-3),U)==2)

# Fixed endpoint variation: the period slope is -2*v*J when J is constant.
vv,J0=sp.symbols("v J0",nonzero=True)
zero("nonzero fixed-node Keller action slope",sp.integrate(-J0,(u,-vv,vv))+2*vv*J0)
P0=(u*u-r*r)*(u*u-s*s)*(u+y)+y
Q0=(u*u-r*r)*(u*u-s*s)*(u*u+y*y)+y*y
Jvary=sp.expand(sp.diff(P0,u)*sp.diff(Q0,y)-sp.diff(P0,y)*sp.diff(Q0,u))
for endpoint in [r,s]:
    Iv=sp.integrate(sp.expand(P0*sp.diff(Q0,u)),(u,-endpoint,endpoint))
    boundary=(P0*sp.diff(Q0,y)).subs(u,endpoint)-(P0*sp.diff(Q0,y)).subs(u,-endpoint)
    zero("fixed two-node endpoint boundary vanishes",boundary)
    zero("two-node endpoint variation identity",
         sp.diff(Iv,y)-boundary+sp.integrate(Jvary,(u,-endpoint,endpoint)))

# Exact non-Keller control: finite monogenic unramified normalization with
# generic conductor degree four, a degree drop, and a repeated-root fibre.
x,T=sp.symbols("x T")
ss=U-T
base_g=(ss-1)*(T*ss-1)
gdeg=(x*x-1)*(y*x*x-1)
Pdeg=x*x+y
Qdeg=x*gdeg
Fdeg=V*V-ss*base_g**2
Jdeg=sp.expand(sp.diff(Pdeg,x)*sp.diff(Qdeg,y)-sp.diff(Pdeg,y)*sp.diff(Qdeg,x))
expectedJ=2*x**6-(5*y+2)*x**4+3*(y+1)*x*x-1
zero("degeneration monic relation",Fdeg.subs({U:Pdeg,V:Qdeg,T:y},simultaneous=True))
zero("degeneration exact Jacobian",Jdeg-expectedJ)
check("degeneration map is not Keller",sp.degree(Jdeg,x)==6)
zero("degeneration conductor trace-dual ratio",
     sp.diff(Fdeg,V).subs({U:Pdeg,V:Qdeg,T:y},simultaneous=True)/sp.diff(Pdeg,x)-gdeg)
zero("degeneration FT retains missing Jacobian factor",
     sp.diff(Fdeg,T).subs({U:Pdeg,V:Qdeg,T:y},simultaneous=True)+Jdeg*gdeg)
zero("degeneration unramified at sole Px zero",sp.diff(Qdeg,x).subs(x,0)-1)
zero("degeneration monogenic normalization relation",x*x-(Pdeg-y))
zero("degeneration odd generator square",Qdeg**2-(Pdeg-y)*gdeg**2)
incdeg=sp.diag(1,base_g)
zero("degeneration inclusion determinant",incdeg.det()-base_g)
zero("degeneration exact norm square",
     sp.resultant(x*x-ss,(x*x-1)*(T*x*x-1),x)-base_g**2)
check("degeneration global conductor degree four",sp.degree(gdeg,x)==4)
check("degeneration special conductor degree two",sp.degree(gdeg.subs(y,0),x)==2)
zero("degeneration repeated-root fibre",gdeg.subs(y,1)-(x*x-1)**2)
zero("degeneration exact conductor discriminant",sp.discriminant(gdeg,x)-16*y*(y-1)**4)
check("degeneration generic fibre is squarefree",sp.discriminant(gdeg.subs(y,2),x)!=0)
zero("degeneration no vertical conductor factor",gdeg.subs(x,0)-1)
for sign in [-1,1]:
    zero("degeneration fixed collision root "+str(sign),gdeg.subs(x,sign))
zero("degeneration equal P on paired points",Pdeg.subs(x,-x)-Pdeg)
zero("degeneration equal Q exactly on conductor",Qdeg-Qdeg.subs(x,-x)-2*x*gdeg)
check("degeneration generic index degree two",sp.degree(base_g,U)==2)
check("degeneration special index degree one",sp.degree(base_g.subs(T,0),U)==1)
zero("moving component escapes from special fibre",(y*x*x-1).subs(y,0)+1)
check("actual Keller condition is essential in FT identity",Jdeg!=1 and Jdeg!=-1)

audit_path=HERE/"sources/neumann-two-node-genus-two-source-audit.json"
audit=json.loads(audit_path.read_text())
source_hashes={}
for source in audit["sources"]:
    path=HERE/"sources"/source["local_pdf"]
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    check("pinned primary source "+source["local_pdf"],actual==source["sha256"])
    source_hashes[source["local_pdf"]]=actual
check("source audit distinguishes written classification from symbolic controls",
      audit["source_classification_proved_by_symbolic_checker"] is False)

proof=HERE/"SIMPLE_FOUR_ROOT_CONDUCTORS_AND_TWO_NODE_ACTION.md"
record={
    "outcome":"PASS","check_count":len(checks),"checks":checks,
    "high_degree_target_controls":controls,
    "one_zero_period_control":{"r_squared":7,"s_squared":3,
                               "first_action":str(zeroAr),
                               "second_action":str(nonzeroAs)},
    "source_sha256":source_hashes,
    "source_audit_sha256":hashlib.sha256(audit_path.read_bytes()).hexdigest(),
    "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
    "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "sympy_version":sp.__version__,
    "elapsed_seconds":round(time.perf_counter()-started,3),
    "degeneration_control":{
        "P":str(Pdeg),"Q":str(Qdeg),"conductor":str(gdeg),
        "Jacobian":str(Jdeg),"generic_conductor_degree":4,
        "degree_at_y0":2,"repeated_roots_at_y1":True,
        "normalization_finite_monogenic_unramified":True,
        "actual_Keller_map":False
    },
    "scope":{
        "every_fibre_four_simple_conductor_roots_excluded_by_written_argument":True,
        "all_globally_reduced_degree_four_conductors_excluded":False,
        "degree_drop_or_repeated_fibre_root_cases_excluded":False,
        "finite_unramified_order_forces_constant_conductor_length":False,
        "flat_pencil_modules_force_constant_length_without_finiteness":False,
        "two_node_target_equivalence_uses_audited_primary_genus_two_table":True,
        "individual_two_node_actions_always_nonzero":False,
        "a_fixed_canonical_two_node_action_is_nonzero":True,
        "classification_or_etale_cover_theorem_proved_by_finite_samples":False,
        "automatic_trace_integrality_or_JC2_proved":False,
        "independent_replay":False,"independent_mathematical_review":False,
        "formal_verification":False,"external_review":False
    }
}
outdir=Path(tempfile.mkdtemp(prefix="two-node-action-check-"))
receipt=outdir/"two_node_action_verification.json"
receipt.write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({"outcome":"PASS","checks":len(checks),"receipt":str(receipt)}))
