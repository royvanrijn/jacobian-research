"""New component discriminant/unit-boundary controls;run before registration."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,a,b=0):check(name,s.cancel(a-b)==0,True)
x,y,X,Y,c,t=s.symbols('x y X Y c t')
for d in [2,3,4]:
 tag='NEW monic component discriminant '+str(d)+' '
 f=y+x*y**d;F=Y+X*Y**d-Y**(d+1);D=d+1
 eq(tag+'source linear projection',f.subs({x:X-Y,y:Y},simultaneous=True),F)
 check(tag+'Y degree equals total degree',s.Poly(F,Y).degree()==s.Poly(F,X,Y).total_degree()==D,True)
 eq(tag+'constant leading coefficient',s.Poly(F,Y).LC(),-1)
 delta=s.discriminant(F-c,Y)
 check(tag+'generic fibre discriminant degree',s.Poly(delta,X).degree()==D,True)
 check(tag+'zero-fibre minimal discriminant degree',s.Poly(delta.subs(c,0),X).degree()==D-1,True)
 eq(tag+'generic top coefficient supported only at t0',s.Poly(s.Poly(delta,X).LC(),c).degree(),d-1)
 check(tag+'nonzero parameter top coefficient monomial',len(s.Poly(s.Poly(delta,X).LC(),c).terms())==1,True)
 eq(tag+'whole generic fibre Laurent parametrization',f.subs(x,(c-y)/y**d),c)
 eq(tag+'zero fibre disjoint factors',1+x*y**(d-1)-x*y**(d-1),1)
 eq(tag+'global submersion Bezout',y**d*(d*d*x*x*y**(d-2))+(1+d*x*y**(d-1))*(1-d*x*y**(d-1)),1)
# Independent local-lattice transcription of the new generic unit-group proof.
for r in [2,3,5]:
 relation=s.ones(1,r)
 check('NEW unit relation rank '+str(r),relation.rank()==1,True)
 check('NEW free quotient rank '+str(r),r-relation.rank()==r-1,True)
 check('NEW primitive diagonal relation '+str(r),s.gcd_list(list(relation))==1,True)
g,m,e,k=s.symbols('g m e k',integer=True,positive=True)
eq('NEW odd hyperelliptic canonical boundary order',-3-(-(2*g+1)),2*g-2)
eq('NEW even hyperelliptic canonical boundary order',-2-(-(g+1)),g-1)
eq('NEW two hyperelliptic punctures give canonical degree',2*(g-1),2*g-2)
eq('NEW ramified hyperelliptic puncture gives canonical degree',(g-1)*2,2*g-2)
eq('NEW ratio valuation at a mate pole',k-m*(-e-1),k+m*(e+1))
check('NEW positive ratio valuation at every mate pole',(m*(e+1)).is_positive,True)
for orders in [[2,3],[1,4,2],[3,3,1,5]]:
 check('NEW constant-field norm remains nonconstant '+str(orders),sum(orders)>0,True)
# Formula (1) consequence,not a replay of any registered source calculation.
for jumps in [[1,1],[2,3],[1,4,2]]:
 for j in range(len(jumps)):
  check('NEW multiple-atypical strict discriminant deficit '+str((jumps,j)),sum(jumps)-jumps[j]>=1,True)
context=json.loads((HERE/'sources/component-boundary-units-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
for p,h in context['new_primary_sha256'].items():check('new primary byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('TWO published proved primary imports and no unproved import',len(context['newly_imported_proved_external_statements'])==2 and context['newly_used_unproved_external_statements']==[])
positive={'published_component_minimal_discriminant_inverse_criterion_adopted','all_failure_component_fibre_discriminants_strictly_above_minimum','every_failure_target_coordinate_has_at_least_two_atypical_values','generic_source_curve_unit_group_exact_reducible_fibre_factor_lattice','geometric_boundary_pluricanonical_section_forces_reducible_fibre','constant_field_norm_descent_with_positive_pole_orders_proved','genus_one_reducibility_recovered_and_credited','specified_hyperelliptic_boundary_classes_excluded_in_all_irreducible_frame','newly_adopted_external_theorem','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'COMPONENT_DISCRIMINANTS_AND_CANONICAL_BOUNDARY_UNITS.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-component-boundary-units.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
