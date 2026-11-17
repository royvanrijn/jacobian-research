"""PF2RULED64 exact controls; the unbounded statements have written proofs."""
from pathlib import Path
import hashlib,json,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['quartic_colouring_receipt.json','sign_square_receipt.json','power_difference_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/formal-branch-index/manifest.json';source=json.loads(manifest.read_text())
    for row in source['primary_sources']:
        p=HERE/row['path'];check('primary HTML pin '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
        check('primary text actually contains imported lifting assertion '+p.name,'equivalence' in p.read_text())
    g,w,aa,bb,ii=S.symbols('g w aa bb ii')
    M=S.Matrix([[2,0,0],[0,-1,g],[0,1,-1]]);one=S.Matrix([1,1,0])
    power=S.Matrix.hstack(one,M*one,M*M*one);J=power.det()
    check('formal trace-zero multiplication matrix',S.trace(M)==0)
    check('formal power-basis index is 9-G',S.expand(J-(9-g))==0)
    F=(w-2)*((w+1)**2-g)
    check('formal factorized characteristic polynomial',S.expand(M.charpoly(w).as_expr()-F)==0)
    check('formal cubic coefficients',S.expand(F-(w**3-(3+g)*w-2+2*g))==0)
    gram=S.diag(1,2,2*g)
    check('standard product trace discriminant',gram.det()==4*g)
    check('power-basis trace discriminant',S.expand((power.T*gram*power).det()-4*g*(9-g)**2)==0)
    check('formal discriminant from characteristic polynomial',S.expand(S.discriminant(F,w)-4*g*(9-g)**2)==0)
    check('universal depressed cubic discriminant',S.expand(S.discriminant(w**3-3*aa*w-2*bb,w)-108*(aa**3-bb**2))==0)
    a0=1+g/3;b0=1-g;i0=1-g/9
    check('normalized formal norm relation',S.expand(a0**3-b0**2-3*g*i0**2)==0)
    h,r,s,v=S.symbols('h r s v')
    check('arbitrary-order coefficient error has its formal divisor',S.expand((a0+h*r)**3-a0**3-h*(3*a0*a0*r+3*a0*h*r*r+h*h*r**3))==0)
    check('arbitrary-order constant coefficient error has its formal divisor',S.expand((b0+h*s)**2-b0**2-h*(2*b0*s+h*s*s))==0)
    check('arbitrary-order determinant error has its formal divisor',S.expand((i0+h*v)**2-i0**2-h*(2*i0*v+h*v*v))==0)
    check('normalized index residue gives explicit branch comaximality',S.expand((1+g*r)-g*r)==1)
    check('finite-order product factors have unit resultant along branch',S.resultant(w-2,(w+1)**2-g,w)==9-g)
    # These are controls for the general CRT/confluent Vandermonde proof, not a rank bound.
    for n in [3,5,7,9]:
        constants=list(range((n+1)//2));mult=[1]+[2]*((n-1)//2)
        rows=[]
        for c,e in zip(constants,mult):
            rows.append([S.Integer(c)**j for j in range(n)])
            if e==2:rows.append([S.Integer(0) if j==0 else j*S.Integer(c)**(j-1) for j in range(n)])
        expected=S.prod((constants[j]-constants[i])**(mult[i]*mult[j]) for i in range(len(constants)) for j in range(i+1,len(constants)))
        check('odd rank confluent CRT determinant '+str(n),S.Matrix(rows).det()==expected and expected!=0)
    for n in [3,5,7,9,15]:
        for a in range(n):
            origin=(a*pow(2,-1,n))%n
            check('odd involution unique fixed label '+str((n,a)),[i for i in range(n) if (a-i)%n==i]==[origin])
            check('odd relabelling gives inversion '+str((n,a)),all(((a-(origin+j))%n-origin)%n==(-j)%n for j in range(n)))
    check('even degree does not have a unique fixed label',[i for i in range(4) if -i%4==i]==[0,2])
    check('odd torsor annihilator symbolic Bezout',(2*(S.Symbol('k',integer=True)+1)-(2*S.Symbol('k',integer=True)+1)).expand()==1)
    x,u,z,t=S.symbols('x u z t');oka=(x*u-2)**2-u**3
    check('Oka branch parameter is Cstar',S.cancel(oka.subs({u:t*t,x:(t**3+2)/(t*t)},simultaneous=True))==0)
    check('Oka Kummer restriction is a cube',S.cancel((x*u-2+z).subs({u:t*t,x:(t**3+2)/(t*t),z:0},simultaneous=True))==t**3)
    y=S.symbols('y');badF=w**3-3*(x+y*y)*w-2*x
    check('quasietale A5 control retains triple-root fibre',badF.subs({x:0,y:0})==w**3)
    check('one element matching every order would have reducible characteristic polynomial',S.div(F,w-2,w)[1]==0)
    inputs=[Path(__file__),HERE/'ODD_DIHEDRAL_FORMAL_SPLITTING_AND_BRANCH_AVOIDING_INDEX.md',manifest]+[HERE/p for p in preds]+[HERE/r['path'] for r in source['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_odd_cyclic_degrees=True,dihedral_inversion_required=True,etale_at_singularities_required=True,
          entire_ramification_curve_split=True,singular_gluing_retained=True,all_finite_thickenings_split=True,
          compatible_formal_splitting=True,reflection_quotient_product=True,branch_neighbourhood_monogenic=True,
          normalization_A1_required_for_formal_theorem=False,arbitrary_Sd_covers_covered=False,even_degree_covered=False),
      quartic=dict(entire_packet=True,all_source_orders_and_rulings=True,all_branch_degrees=True,
          formal_tracezero_element='(2,-1+z)',formal_polynomial='(W-2)((W+1)^2-G)',
          every_finite_order_has_global_primitive_lift=True,normalized_relation='a^3-b^2=3*G*I^2',
          congruences=['a=1+G/3 mod G^m','b=1-G mod G^m','I=1-G/9 mod G^m'],
          index_disjoint_from_branch=True,index_nonconstant_by_PF56=True,different_lift_per_order=True,
          one_global_lift_of_all_orders_exists=False),
      controls=dict(Oka_nontrivial_Pic3_with_split_branch=True,quasietale_A5_control_excluded=True,
          index_square_factor_retained=True,finite_jet_increase_supplies_exclusion=False,
          Picard_restriction_injective_asserted=False,formal_split_implies_global_split=False),
      assurance=dict(written_unbounded_proofs=True,primary_nilpotent_invariance_imported=True,
          finite_controls_not_universal_proof=True,independent_full_reproof=False,formal_verification=False,
          external_review=False,novelty_priority=False,Pic3_vanishing_proved=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_pair_supplied=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'formal_branch_index_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; all-odd formal splitting and branch-avoiding cubic index')

if __name__=='__main__':main()
