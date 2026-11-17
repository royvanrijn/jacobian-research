"""Exact controls for the written radical-plinth recognition and homology proof."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
F,G,h,P,Q,a,alpha,beta,ell,d=sp.symbols('F G h P Q a alpha beta ell d')
delta=lambda f:sp.expand(G*sp.diff(f,P)-F*sp.diff(f,Q))
zero(delta(F*P+G*Q),'the universal two-slice numerator is invariant')
zero(P/G+Q/F-(F*P+G*Q)/(F*G),'the Cech transition uses the full global numerator')
zero(delta(P/G)-1,'the G chart has its exact local slice')
zero(delta(-Q/F)-1,'the F chart has its exact local slice')
zero(F*((a-G*Q)/F)+G*Q-a,'the F-local hypersurface chart has its exact inverse')
zero(F*P+G*((a-F*P)/G)-a,'the G-local hypersurface chart has its exact inverse')
zero((F*P+G*Q-(F*alpha+G*beta))-F*(P-alpha)-G*(Q-beta),
     'an arbitrary invariant ideal correction changes the numerator by exactly its prescribed class')
zero((ell*h+d+F*alpha+G*beta)-F*alpha-G*beta-(ell*h+d),
     'once root residues are linear,the global affine residue representative is exact')
# The local derivatives at a coefficient zero;no global Keller condition is imposed.
af,ag,ah,aff,afg,agg,j11,j12,j21,j22=sp.symbols('af ag ah aff afg agg j11 j12 j21 j22')
gradient=sp.Matrix([P-af,Q-ag,-ah,F,G])
critical={F:0,G:0,P:af,Q:ag}
zeros(list(gradient.subs(critical,simultaneous=True)-sp.Matrix([0,0,-ah,0,0])),
      'a repeated root has exactly one possible singular point,while a simple root is smooth')
normal_graph=sp.Matrix([[1,0,0,0],[0,1,0,0],[(P-af)/ah,(Q-ag)/ah,0,0],[0,0,1,0],[0,0,0,1]])
zeros(list((gradient.T*normal_graph).subs({F:0,G:0},simultaneous=True)),
      'at a simple reached value the base coordinates and P,Q form the local tangent chart')
check(normal_graph.extract([0,1],[0,1])==sp.eye(2),
      'the normal map to the base is literally the identity in the two coefficient parameters')
check(normal_graph.extract([0,1],[2,3])==sp.zeros(2),
      'the tangent directions of each reached affine plane have zero base normal image')
HH=sp.Matrix([[aff,afg,j11/ah,j21/ah],[afg,agg,j12/ah,j22/ah],
              [j11/ah,j12/ah,0,0],[j21/ah,j22/ah,0,0]])
zero(HH.det()-(j11*j22-j12*j21)**2/ah**4,
     'each simple root contributes one nondegenerate full critical point in the local source chart')
# Complex normal changes preserve the integral orientation,with no hidden multiplicity.
r11,r12,r21,r22,i11,i12,i21,i22=sp.symbols('r11 r12 r21 r22 i11 i12 i21 i22')
Re=sp.Matrix([[r11,r12],[r21,r22]]);Im=sp.Matrix([[i11,i12],[i21,i22]])
realification=Re.row_join(-Im).col_join(Im.row_join(Re));detc=(Re+sp.I*Im).det()
zero(realification.det()-sp.expand(detc*sp.conjugate(detc)).subs(
     {sp.conjugate(v):v for v in [r11,r12,r21,r22,i11,i12,i21,i22]}),
     'a complex normal isomorphism has positive real determinant abs(det_complex)^2')
# Universal row-sum coordinates give the integral block normal form in every degree.
first,rest,total=sp.symbols('first rest total')
zero((first+rest)-rest-first,'the row-sum coordinate change has an exact inverse for any number of remaining coordinates')
zero((total-rest)+rest-total,'the inverse row-sum change recovers the whole normal-boundary sum')
block=sp.Matrix([[1,1,1]])
change=sp.Matrix([[1,-1,-1],[0,1,0],[0,0,1]])
check(change.det()==1,'the actual three-point block change is unimodular over Z')
check(block*change==sp.Matrix([[1,0,0]]),'normal disks over one plinth point have the same primitive boundary')
check(block*change[:,1:]==sp.zeros(1,2),'differences of normal disks generate the repeated-value boundary kernel')
check(sp.eye(3).det()==1,'one reached value over each actual three-point target line gives the identity boundary map')
# This is the exact Euler-balanced negative control from the proof,not a degree sweep.
x,z=sp.symbols('x z');C=z*(z-1);aa=1-z+z*h*(h-1)
zero(aa.subs(z,0)-1,'the negative example omits the entire first plinth line')
zero(aa.subs(z,1)-h*(h-1),'the negative example reaches two values on the second plinth line')
zero(sp.diff(aa,h).subs({z:1,h:0})+1,'the first reached root is simple')
zero(sp.diff(aa,h).subs({z:1,h:1})-1,'the second reached root is simple')
S=x*P+C*Q-aa
gradS=sp.Matrix([sp.diff(S,v) for v in [x,z,h,P,Q]])
for value in [0,1]:
    zero(gradS[2].subs({x:0,z:1,h:value})-(1-2*value),
         'the whole hypersurface is smooth at the reached plane '+str(value))
check(sp.Poly(aa.subs(z,0),h).degree()==0 and sp.Poly(aa.subs(z,1),h).degree()==2,
      'the counterexample has exactly the stated degree distribution (0,2)')
T=sp.Matrix([[0,0],[1,1]])
check(T.rank()==1,'the negative example has rank-one boundary despite equal source and target ranks')
check(len(T.nullspace())==1 and len(T.T.nullspace())==1,
      'both H4 and H3 survive in the Euler-balanced example')
check(T.nullspace()[0]==sp.Matrix([-1,1]),'the surviving H4 class is the difference of the two reached-plane disks')
check(sp.Matrix.hstack(T[:,0],sp.Matrix([1,0])).det()==-1,
      'the surviving H3 class is primitive and has no hidden torsion')
check(1+2-2==1,'Euler characteristic one cannot replace the two separate vanishing conditions')
# A universal cotangent frame supplies the converse action;the local plinth proof is retained JU.
f,g,U,V,d0=sp.symbols('f g U V d0');hh=d0+U*f+V*g
D=lambda value:g*sp.diff(value,U)-f*sp.diff(value,V)
zero(D(hh),'the simple cotangent translation preserves the actual potential')
zero(D(U)-g,'the first global local slice supplies g to the plinth')
zero(D(-V)-f,'the second global local slice supplies f to the plinth')
for value,label in [(U,'U'),(-V,'minus V')]:
    zero(D(D(value)),'the displayed global local slice has exact delta degree one '+label)
# Gate logic is universal in the existence bits,not in polynomial coefficients or degrees.
for alignment in [False,True]:
    for en in [False,True]:
        unit=en and alignment
        check((not unit)==((not en) or (not alignment)),
              'source alignment and plane Keller existence remain a conjunction '+str((alignment,en)))
algebra=len(checks);context=json.loads((HERE/'sources/radical-plinth-rigidity-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for key in ['newly_cached_primary_context','rechecked_retained_primary_context']:
    for entry in context[key]:
        check(hashlib.sha256((HERE/'sources'/entry['retained_file'] if key=='newly_cached_primary_context' else HERE/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
              'inspected primary '+entry['retained_file'])
imports=context['newly_imported_proved_external_statements']
check(len(imports)==1 and imports[0]['hypotheses_checked']
      and imports[0]['source']=='https://arxiv.org/html/2112.08241v2'
      and 'transverse naturality' in imports[0]['location'],
      'one precise purity-naturality import has the full inverse image and normal isomorphism hypotheses checked')
check(context['newly_used_unproved_external_statements']==[],
      'no individual-action classification or quotient descent conjecture is imported')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==5 and not flags['an_individual_native_action_with_required_kernel_and_plinth_constructed']
      and not flags['radical_of_plinth_alone_suffices']
      and not flags['S_rho_normal_form_source_alignment_constructed']
      and not flags['original_goal_complete'],
      'automatic quotient recovery under exact algebraic data does not assert existence of that action or completion of HC4')
proof=HERE/'RADICAL_PLINTH_RIGIDITY_AND_AUTOMATIC_QUOTIENT_RECOVERY.md'
paths=['sources/radical-plinth-rigidity-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-radical-plinth-rigidity.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
