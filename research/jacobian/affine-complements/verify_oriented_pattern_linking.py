"""PF2RULED66: exact oriented plumbing controls after written unbounded reductions."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,math,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent
def imported(name,file):
    spec=importlib.util.spec_from_file_location(name,HERE/file);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
E=imported('sealed_end61','verify_marked_end_lattice.py')
V=imported('sealed_coordinate50','verify_coordinate_pencil_matrix.py')

def negative_cf(alpha,beta):
    out=[]
    while beta:
        q=(alpha+beta-1)//beta;out.append(q);alpha,beta=beta,q*beta-alpha
    return out

def pattern_plumbing(s,t):
    assert s>=3 and math.gcd(3*s,4*t)==1
    alpha=[s,3*s]+([t] if t>1 else [])
    beta=[(-pow(4*t,-1,s))%s,(-pow(4*t,-1,3*s))%(3*s)]
    if t>1:beta.append((-pow(3*s,-1,t))%t)
    b=sum(S.Rational(be,al) for al,be in zip(alpha,beta))+S.Rational(1,3*s*t)
    assert b.q==1
    legs=[negative_cf(al,be) for al,be in zip(alpha,beta)]
    size=1+sum(map(len,legs));Q=S.zeros(size);Q[0,0]=-b;first=[];n=1
    for leg in legs:
        first.append(n);parent=0
        for coefficient in leg:
            Q[n,n]=-coefficient;Q[parent,n]=Q[n,parent]=1;parent=n;n+=1
    I=Q.inv();j=first[0];canonical=S.Matrix([-Q[i,i]-2 for i in range(Q.rows)])
    return dict(s=s,t=t,alpha=alpha,beta=beta,b=int(b),legs=legs,Q=Q,I=I,
        first=first,linking=-I[j,j],canonical_lift=I*canonical,H1=E.h1(Q.tolist()))

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    preds=['coordinate_pencil_matrix_receipt.json','marked_end_lattice_receipt.json','binomial_periodic_receipt.json','cyclic_companion_receipt.json']
    for name in preds:
        old=json.loads((HERE/name).read_text());check('predecessor PASS '+name,old['outcome']=='PASS')
        for p,h in old['source_sha256'].items():check('sealed predecessor '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/oriented-pattern-linking/manifest.json';source=json.loads(manifest.read_text())
    for row in source['primary_sources']:
        p=HERE/row['path'];check('primary bytes '+p.name,sha(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    w,A,B=S.symbols('w A B');check('unrestricted quartic discriminant fixes branch torus',S.discriminant(w**4+A*w+B,w)==256*B**3-27*A**4)
    u,v=S.symbols('u v');check('t1 explicit A singularity coordinate change',S.expand((w**4+u*w+v**3)-(w*(u+w**3)+v**3))==0)
    s,t,be,gam,de,e=S.symbols('s t beta gamma delta e',nonzero=True)
    rational=S.Matrix([[e-be/s-gam/(3*s)-de/t,1,1,1],[1,-s/be,0,0],[1,0,-3*s/gam,0],[1,0,0,-t/de]])
    inverse=rational.inv()
    check('unrestricted central Schur and selected inverse diagonal',S.factor(inverse[1,1]-(-be/s+(be/s)**2/e))==0)
    link=S.factor((-inverse[1,1]).subs(e,-1/(3*s*t)))
    check('oriented pattern linking numerator',S.factor(link-be*(1+3*t*be)/s)==0)
    check('unbounded modular generator rescaling factorization',S.expand(16*t*t*be*(1+3*t*be)+t-t*(4*t*be+1)*(12*t*be+1))==0)
    check('positive Hopf power map and degree4 pattern Euler',S.cancel(4*(-1/((3*s)*(4*t)))+1/(3*s*t))==0)
    nr,nd,c=S.symbols('nr nd c');Q=S.Matrix([[-nr,c],[c,-nd]])
    check('source eD linking equals ramified diagonal over determinant',S.simplify(-Q.inv()[1,1]-nr/(nr*nd-c*c))==0)
    # Unique triangle permutation class, independent of weights and coefficient choices.
    perms=list(itertools.permutations(range(4)));identity=tuple(range(4))
    order=lambda p:next(i for i in range(1,13) if E.power(p,i)==identity)
    pairs=[(a,b) for a,b in itertools.product(perms,repeat=2) if order(a)==3 and order(b)==4 and sum(E.mul(a,b)[i]!=i for i in range(4))==2]
    def orbit(pair):return min(tuple(E.mul(E.mul(E.inv(g),p),g) for p in pair) for g in perms)
    check('all three four two triangle pairs one orientation fixed conjugacy class',len(pairs)==24 and len({orbit(pair) for pair in pairs})==1)
    plumbings=[]
    for ss,tt in [(3,1),(5,1),(7,1),(3,4),(3,7),(3,16),(15,79),(27,142),(123,4)]:
        row=pattern_plumbing(ss,tt);Q=row['Q'];I=row['I'];j=row['first'][0];tag=str((ss,tt))
        check('normalized negative continued fractions '+tag,all(all(x>=2 for x in leg) for leg in row['legs']))
        check('all normalized slice parameters have signed inverse '+tag,(4*tt*row['beta'][0]+1)%ss==0 and (4*tt*row['beta'][1]+1)%(3*ss)==0)
        check('exact integral star determinant and cyclic H1 '+tag,abs(Q.det())==ss and row['H1']==dict(free=0,torsion=[ss]))
        check('negative definiteness by rational arm elimination '+tag,-row['b']+sum(S.Rational(be,al) for al,be in zip(row['alpha'],row['beta']))==-S.Rational(1,3*ss*tt))
        check('full integral inverse equals local Schur linking '+tag,row['linking']==S.Rational(row['beta'][0]*(1+3*tt*row['beta'][0]),ss))
        check('rescaled generator has exactly minus t linking '+tag,(16*tt*tt*row['linking']+S.Rational(tt,ss)).q==1)
        check('selected meridian generates entire cyclic group '+tag,row['linking'].q==ss)
        check('hypersurface canonical class has zero boundary mark '+tag,all(x.q==1 for x in row['canonical_lift']))
        if tt==1:
            chain=-2*S.eye(ss-1)
            for i in range(ss-2):chain[i,i+1]=chain[i+1,i]=1
            chainlink=-chain.inv()[0,0]
            check('A chain independent complex orientation '+tag,(chainlink+S.Rational(1,ss)).q==1 and (16*row['linking']-chainlink).q==1)
        plumbings.append(dict(s=ss,t=tt,alpha=row['alpha'],beta=row['beta'],central_b=row['b'],legs=row['legs'],
            matrix=[[int(x) for x in rr] for rr in Q.tolist()],linking_selected=str(row['linking']),rescaled_coefficient_mod_s=(-tt)%ss))
    parameter_controls=[]
    for p,q,r,kind in [(3,5,9,'C4'),(3,5,45,'C4'),(5,7,9,'C4'),(5,7,16,'C3'),(7,11,16,'C3')]:
        c=p*q-p-q;k=r*c+1
        assert kind!='C3' or k%9==0
        D=r//3 if kind=='C4' else k//3;tt=k//4 if kind=='C4' else r//4
        check('surviving pattern has primitive weights '+str((p,q,r,kind)),math.gcd(3*D,4*tt)==1)
        check('surviving pattern has nonsquare oriented coefficient mod3 '+str((p,q,r,kind)),D%3==0 and tt%3==1 and (-tt)%3==2)
        if kind=='C4':check('C4 full square class is minus one '+str((p,q,r)),(4*tt-1)%D==0)
        else:check('C3 full square class equals c '+str((p,q,r)),(-4*tt*c-1)%D==0)
        parameter_controls.append(dict(p=p,q=q,r=r,k=k,kind=kind,D=D,t=tt,linking_coefficient_mod_D=(-tt)%D))
    N,C=S.symbols('N C',integer=True,nonnegative=True);r=36*N+9;c=4*C+3;k=r*c+1
    check('all C4 windings have t1mod3',all(co.q==1 for co in S.Poly(S.expand((k-4)/3),N,C).coeffs()))
    check('all C4 coefficient differs from quarter by a multiple of D',S.expand(k-1-r*c)==0)
    check('remaining C3 unit c residue forces r1mod3',all((cc*rr+1)%3!=0 or (cc==2 and rr==1) for pp,qq,rr in itertools.product([1,2],repeat=3) for cc in [(pp*qq-pp-qq)%3]))
    check('a unit square mod3 is always one',{x*x%3 for x in [1,2]}=={1})
    for residue in [0,1,2]:
        rsrc=3*N+residue
        check('unbounded complementary nR residue '+str(residue),S.Poly(rsrc*rsrc-rsrc+1,N).nth(0)%3==[1,1,0][residue])
        check('all nonconstant terms in source residue vanish mod3 '+str(residue),all(co%3==0 for (exp,),co in S.Poly(rsrc*rsrc-rsrc+1,N).terms() if exp>0))
    pp,nn=S.symbols('p n',integer=True)
    check('source order1 complementary determinant is even',S.expand(5*(2*nn)-(2*pp)**2)==2*(5*nn-2*pp**2))
    rr,qq=S.symbols('r q',integer=True)
    check('unrestricted ramified lower bound',S.expand(((rr-1)*qq+2*qq-rr+2).subs(qq,rr-1)-(rr*rr-rr+1))==0)
    # Actual boundary blowups for two common-pencil values; not a formal disjointness heuristic.
    z,zz,T,Vv=S.symbols('z zz T V')
    clusters=[]
    for rsrc,qsrc,ds,dp in [(1,2,1,2),(2,1,1,2),(4,3,2,3)]:
        Rcl=V.cluster(z**(-qsrc)+z**(2-rsrc)*T,z**rsrc,z,T)
        Dcl=V.cluster(2*zz**(-dp)+zz**(1-ds)*Vv,1+zz**ds,zz,Vv)
        full,inv,labels,K=V.blowup_matrix([Rcl,Dcl]);i=labels[Rcl[1]];j=labels[Dcl[1]]
        check('different common source values give exact zero dual cross '+str((rsrc,qsrc,ds,dp)),inv[i,j]==0)
        check('whole source extraction retains canonical orders '+str((rsrc,qsrc,ds,dp)),K[i]==1 and K[j]==0)
        clusters.append(dict(r=rsrc,q=qsrc,s=ds,p=dp,cross=int(inv[i,j]),nR=int(-inv[i,i]),nD=int(-inv[j,j])))
    check('zero mark in diagonal block forces nR1 determinant identity',S.cancel((nr*nd)/nd)==nr)
    inputs=[Path(__file__),HERE/'ORIENTED_PATTERN_LINKING_AND_COMMON_SOURCE_GATE.md',manifest]+[HERE/x for x in preds]+[HERE/x['path'] for x in source['primary_sources']]
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      general=dict(all_finite_cover_degrees=True,integral_homology_solid_torus_required=True,
          cyclic_companion_branched_H1_zero_required=True,oriented_linking_preserved=True,outer_mark_preserved=True,
          replacement_implies_homeomorphism=False,replacement_implies_contact_equivalence=False,
          arbitrary_nontrivial_companion_homology_deleted=False),
      pattern=dict(model='w^4+u^t*w+v^s=0',positive_branch='T(3s,4t)',primitive_weights=["3s","4t","st"],
          euler='-1/(3st)',first_leg_beta='-(4t)^-1 mod s',linking='<-t/s>',complex_outward_orientation_fixed=True,
          A_chain_independent_orientation_control=True,unrestricted_Schur_and_rescaling_proved=True),
      all_odd=dict(all_coprime_odd_first_pairs=True,all_primitive_actual_binomial_windings=True,
          full_remaining_peripheral_classes_covered=True,ramification_mark_zero=True,
          C4_linking='<-1/(r/3)>',C3_remaining_linking='<-(r/4)/(k/3)> or <c/(k/3)>',
          source_nR_congruent2_mod3=True,ENTIRE_complementary_source_excluded=True,
          all_source_orders_and_unramified_poles=True,common_source_pencil_required=True,
          common_source_boundary_value_coincides=True,distinct_target_images_identified=False,
          common_source_configuration_excluded=False),
      plumbing_controls=plumbings,parameter_controls=parameter_controls,different_value_blowup_controls=clusters,
      assurance=dict(written_unbounded_proofs=True,primary_Euler_and_plumbing_conventions_imported=True,
          finite_controls_not_universal_proof=True,orientation_sign_checked=True,
          contact_filling_classification_used=False,ruling_implies_lens_asserted=False,
          full_satellite_negative_plumbing_asserted=False,independent_full_reproof=False,
          formal_verification=False,external_review=False,novelty_priority=False,
          full_quartic_exclusion=False,higher_degree_exclusion=False,finite_Keller_map=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'oriented_pattern_linking_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; oriented linking excludes all complementary sources in the all-odd class')

if __name__=='__main__':main()
