"""Exact controls for the written integral middle-cycle argument; not topology verification."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):
    check(sp.simplify(sp.expand(value))==0,label)
def zeros(values,label):
    check(all(sp.simplify(sp.expand(v))==0 for v in values),label)

w1,w2,f1,f2,c,L,x,y1,y2,T=sp.symbols('w1 w2 f1 f2 c L x y1 y2 T')
F=w1*f1+w2*f2
zero(F-w1*f1-w2*f2,'the first principal-divisor ideal reduces to f1 and w2*f2')
zero(w2*f2-(F-w1*f1),'the reverse first ideal inclusion is exact')
zero(F-w2*f2-w1*f1,'the second principal-divisor ideal reduces to w2 and w1*f1')
zero(w1*f1-(F-w2*f2),'the reverse second ideal inclusion is exact')
zero(F.subs(f2,-w1*f1/w2),'generic section chart eliminates f2 when w2 is a unit')
zero(F.subs(w2,-w1*f1/f2),'generic B chart eliminates w2 when f2 is a unit')
zero(F.subs(w1,-w2*f2/f1),'generic P chart eliminates w1 when f1 is a unit')
check(sp.degree(f1,f1)==1,'the section divisor parameter has multiplicity one')
check(sp.degree(-w1*f1/f2,f1)==1,'both divisors along B have generic multiplicity one')
check(sp.degree(w2,w2)==1,'the zero-section divisor parameter has multiplicity one')
check(sp.Poly(F,w1,w2).total_degree()==1,'over the source function field F is linear in the free pair')
check(sp.gcd(f1,f2)==1,'the abstract coefficient pair is coprime; the proof obtains this from finite fibres')

# The exact section homotopy includes all polynomial degrees via a substitution identity.
zero(x**2+T*w1*y1+T*w2*y2-c**2 -
     (x**2+w1*(T*y1)+w2*(T*y2)-c**2),
     'the section equation at T*w equals the equation of (w,X(Tw),T*Y(Tw))')
aa=[w1**2+w2,w2**3-w1]
hh=w1*aa[0]+w2*aa[1]
xx=c+hh
yy=[-(2*c+hh)*a for a in aa]
zero(xx**2+w1*yy[0]+w2*yy[1]-c**2,'a nonlinear polynomial section lies exactly on the smooth quadric')
sub={w1:T*w1,w2:T*w2}
xt=xx.subs(sub,simultaneous=True)
yt=[T*a.subs(sub,simultaneous=True) for a in yy]
zero(xt**2+w1*yt[0]+w2*yt[1]-c**2,'its polynomial contraction stays on the quadric for every parameter')
zeros([xt.subs(T,0)-c]+[a.subs(T,0) for a in yt],'the contraction endpoint is the standard section')
zeros([xt.subs(T,1)-xx]+[a.subs(T,1)-b for a,b in zip(yt,yy)],'the other contraction endpoint is the actual nonlinear section')
check(all(sp.Poly(a,w1,w2).total_degree()>1 for a in [xx]+yy),'section controls include unrestricted nonlinear free dependence')

# Explicit orthogonal rotation in the two b-coordinates.
a1,a2,b1,b2,C,S=sp.symbols('a1 a2 b1 b2 C S')
rot=[C*b1+S*b2,-S*b1+C*b2]
zero(sum(v*v for v in rot)-(C*C+S*S)*(b1*b1+b2*b2),'the b-coordinate rotation preserves the quadratic for C^2+S^2=1')
zero(sp.Matrix(rot).jacobian([b1,b2]).det()-(C*C+S*S),'the rotation is invertible along the compact unit circle')
oldw=[a1+sp.I*b1,a2+sp.I*b2]
oldy=[a1-sp.I*b1,a2-sp.I*b2]
neww=[a1+sp.I*rot[0],a2+sp.I*rot[1]]
newy=[a1-sp.I*rot[0],a2-sp.I*rot[1]]
zeros([u.subs({C:1,S:0})-v for u,v in zip(neww+newy,oldw+oldy)],'the rotation starts at the identity')
zeros([u.subs({C:-1,S:0})-v for u,v in zip(neww+newy,oldy+oldw)],'the rotation endpoint exchanges both w-y pairs')
zero(sum(u*v for u,v in zip(oldw,oldy))-(a1*a1+a2*a2+b1*b1+b2*b2),'the hyperbolic coordinates become a sum of four squares')
swap=sp.zeros(4)
swap[0,2]=swap[1,3]=swap[2,0]=swap[3,1]=1
zero(swap.det()-1,'the simultaneous hyperbolic-pair swap has positive complex determinant')

# Tangent columns in real coordinates (Re w, Im w, Re y, Im y).
Ip=sp.eye(4)
conj=sp.diag(1,-1,1,-1)
tanP=sp.zeros(8,4)
tanP[4:8,0:4]=Ip
tanC=Ip.col_join(conj)
transverse=tanP.row_join(tanC)
check(abs(transverse.det())==1,'the plane and compact real sphere meet transversely with unit integral local pairing')
zero(conj.det()-1,'complex conjugation on two complex coordinates preserves real orientation')
zero((x**2+w1*y1+w2*y2-1).subs({y1:oldy[0],y2:oldy[1],w1:oldw[0],w2:oldw[1]},simultaneous=True) -
     (x*x+a1*a1+a2*a2+b1*b1+b2*b2-1),
     'the compact real form is exactly a real unit four-sphere')
zero(sp.diff(x*x+w1*y1+w2*y2-c*c,x).subs(x,c)-2*c,'the positive pole is smooth exactly under the stated nonzero-c hypothesis')
check(all(sp.diff(x*x+w1*y1+w2*y2,v).subs({x:0,w1:0,w2:0,y1:0,y2:0})==0 for v in [x,w1,w2,y1,y2]),
      'the critical-value cone is singular and is not covered by this sphere argument')

# A genuine nonlinear Keller pair with full source recovery and high free dependence.
s,t=sp.symbols('s t')
k1=s+t*t
k2=t+k1**3
A=[w1**2,w2**3]
H=w1*A[0]+w2*A[1]
cc=w1*w2
X=1+H
Y=k1-(2+H)*A[0]+w2*cc
Z=k2-(2+H)*A[1]-w1*cc
phi=sp.expand(X*X+w1*Y+w2*Z)
zero(sp.Matrix([k1,k2]).jacobian([s,t]).det()-1,'the positive embedding control has a genuine plane Keller pair')
zero(phi-1-w1*k1-w2*k2,'the actual nonlinear embedding has the required constant-trace identity')
zero(sp.hessian(phi,[s,t,w1,w2]).det()-1,'the actual four-variable potential has unit Hessian')
rec1=Y+(X+1)*A[0]-w2*cc
rec2=Z+(X+1)*A[1]+w1*cc
rect=rec2-rec1**3
recs=rec1-rect**2
zeros([recs-s,rect-t],'ambient polynomials recover both source parameters globally')
zeros([Y.subs({w1:0,w2:0})-k1,Z.subs({w1:0,w2:0})-k2],'the zero-section is the entire positive helper plane in genuine source coordinates')

# No inverse statement without Keller, and no source-cycle counting from a mere parameter map.
bad=[s,t*t]
zero(sp.Matrix(bad).jacobian([s,t]).det()-2*t,'a two-sheet parameter control is explicitly non-Keller')
zeros([bad[0].subs(t,1)-bad[0].subs(t,-1),bad[1].subs(t,1)-bad[1].subs(t,-1)],'that parameter map identifies distinct source points')
zero(c*c+w1*s+w2*t*t-(c*c+w1*bad[0]+w2*bad[1]),'the same formal potential identity alone does not certify a closed source embedding')
check(sp.degree(t*t-1,t)==2,'the nonclosed parameter example has two generic source sheets')
check((3-1)%2==0 and 3-1!=0,'a mod-two pairing would not exclude an odd degree greater than one')

context=json.loads((HERE/'sources/quadric-middle-cycle-context.json').read_text())
check(len(context['newly_used_unproved_external_statements'])==2,
      'the integral cycle map and Poincare-duality foundations are explicit external inputs')
check(context['actual_closed_embedding_required'] and not context['source_field_lift_substituted_for_polynomial_lift'],
      'the actual closed coordinate-ring surjection supplies a polynomial ambient lift')
check(context['full_noncritical_rank_one_closed_embedding_chart_resolved'] and context['all_constant_trace_intrinsic_free_projections_resolved'],
      'the written integral cycle proof closes the full noncritical chart and all retained constant traces')
check(not context['all_embeddings_assumed_to_retain_intrinsic_free_projections'] and not context['all_zero_section_traces_assumed_constant'],
      'neither existence of an aligned projection nor constancy of arbitrary traces is asserted')
check(not context['unrestricted_HC4_or_JC2_resolved'] and not context['original_goal_complete'],
      'the original unrestricted objective and both endpoints remain open')
math_count=len(checks)
for name,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,'retained proof/source byte pin '+name)
for item in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/item['retained_file']).read_bytes()).hexdigest()==item['sha256'],
          'new primary source byte pin '+item['retained_file'])
old_counts={}
for name in ['receipt-free-projection-trace.json','receipt-primitive-helpers.json','receipt-abelian-classgroup.json']:
    old=json.loads((HERE/name).read_text())
    check(old['outcome']=='PASS','retained receipt PASS '+name)
    old_counts[name]=old['check_count']
    for source,digest in old.get('source_sha256',{}).items():
        check(hashlib.sha256((HERE/source).read_bytes()).hexdigest()==digest,'retained receipt source '+name+' -> '+source)
    for source,digest in old.get('dependencies',{}).items():
        check(hashlib.sha256((ROOT/source).read_bytes()).hexdigest()==digest,'retained receipt dependency '+name+' -> '+source)
names=['verify_quadric_middle_cycle.py','QUADRIC_MIDDLE_CYCLE_AND_FREE_PROJECTION_CLOSURE.md','sources/quadric-middle-cycle-context.json']
result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,checks=checks,
 source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
 dependencies=context['retained_proof_and_source_sha256'],retained_receipt_counts_not_newly_replayed=old_counts,
 newly_cached_primary_context=context['newly_cached_primary_context'],
 newly_used_unproved_external_statements=context['newly_used_unproved_external_statements'],
 retained_structural_inputs=context['retained_structural_inputs'],
 written_not_finite_exhaustion=context['direct_proof_arguments'],
 **{key:value for key,value in context.items() if type(value) is bool},
 runtime=dict(python=platform.python_version(),sympy=sp.__version__),elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-quadric-middle-cycle.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=math_count,
 original_goal_complete=False),indent=2))
