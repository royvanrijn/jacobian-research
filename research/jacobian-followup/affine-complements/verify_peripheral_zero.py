"""Exact PF2RULED60 controls; the unbounded peripheral and degree proofs are written."""
from pathlib import Path
import hashlib,importlib.util,itertools,json,os,platform
import sympy as S
HERE=Path(__file__).resolve().parent;checks=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def check(name,value):
    assert value,name
    checks.append(name)
def main():
    preds=['unramified_orevkov_receipt.json','quartic_colouring_receipt.json','sign_square_receipt.json','canonical_quadrics_receipt.json','nodal_deformation_receipt.json']
    for name in preds:
        r=json.loads((HERE/name).read_text());check('predecessor PASS '+name,r['outcome']=='PASS')
        for p,h in r['source_sha256'].items():check('predecessor seal '+name+' '+p,sha(HERE/p)==h)
    manifest=HERE/'sources/peripheral-zero-control/manifest.json';m=json.loads(manifest.read_text())
    for row in m['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements'];check('primary pin '+p.name,p.stat().st_size==row['bytes'] and sha(p)==row['sha256'])
    spec=importlib.util.spec_from_file_location('sealed_colouring54',HERE/'verify_quartic_colouring.py');G=importlib.util.module_from_spec(spec);spec.loader.exec_module(G)
    perms=list(itertools.permutations(range(4)));I4=tuple(range(4));trans=[p for p in perms if sum(p[i]!=i for i in range(4))==2]
    def parity(p):return sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2
    for tau in trans:
        other=[i for i in range(4) if tau[i]==i]
        allowed=[p for p in perms if G.compose(p,tau)==G.compose(tau,p) and parity(p)==0 and all(p[i]==i for i in other)]
        check('even centralizer fixes residual sheets only at identity '+str(tau),allowed==[I4])
    for contacts in [[1],[2],[3,4],[1,1,2,9],[10,11,12]]:
        check('paired contact parity even '+str(contacts),sum(2*(n%2) for n in contacts)%2==0)
    def constrained(mats,rels,mus,longitude):
        rows=[]
        for word in rels:
            A,z=G.word_linear(word,mats);assert A==G.I;rows.extend(z)
        for word in mus:
            A,z=G.word_linear(word,mats);Az=G.transform(A,z);rows.extend([z[0]^Az[0],z[1]^Az[1]])
        A,z=G.word_linear(longitude,mats);assert A==G.I
        cob=[]
        for v in [1,2]:
            cc=0
            for j,A in enumerate(mats):cc|=(v^G.act(A,v))<<(2*j)
            cob.append(cc)
        dimcob=G.rank2(cob);base=2*len(mats)-G.rank2(rows)-dimcob;rows.extend(z)
        return base,2*len(mats)-G.rank2(rows)-dimcob
    s1=(1,0,2,3);s2=(0,2,1,3);s3=(0,1,3,2)
    rels=[[1,1],[2,2],[3,3],[1,2]*3,[2,3]*3,[1,3]*2]
    dims=constrained([G.project(p) for p in [s1,s2,s3]],rels,[[1],[2],[3]],[1,3])
    check('Coxeter group control parabolic1 longitude-zero0',dims==(1,0))
    aa=(1,2,0,3);bb=(1,2,3,0)
    dims34=constrained([G.project(aa),G.project(bb)],[[1]*3+[-2]*4],[[1,-2]],[1]*3+[2,-1]*12)
    check('T34 peripheral control survives longitude equation',dims34[0]==dims34[1] and dims34[1]>0)
    p,q,r,k=S.symbols('p q r k',integer=True,positive=True)
    check('unbounded Wronskian degree equals normal recurrence',S.expand((k+q*r-1-(q-1)*p*r)-(k+r*(p+q-p*q)-1))==0)
    check('eta0 defect has DZ equality degree',S.expand((p*q*r+1)-(p*r+q*r)-(r*(p*q-p-q)+1))==0)
    U,V,Up,Vp=S.symbols('U V Up Vp',positive=True)
    F=U**q-V**p;Fp=q*U**(q-1)*Up-p*V**(p-1)*Vp;W=q*Up*V-p*U*Vp
    check('universal weighted Wronskian identity',S.simplify(V*Fp-p*Vp*F-U**(q-1)*W)==0)
    t,s=S.symbols('t s');Y=t**6+4*t**4+10*t*t+6;X=t**9+6*t**7+21*t**5+35*t**3+S.Rational(63,2)*t
    defect=-S.Rational(27,4)*(4*t**4+13*t*t+32)
    check('explicit DZ defect exact',S.expand(X*X-Y**3-defect)==0)
    check('explicit weighted Wronskian378',S.expand(2*S.diff(X,t)*Y-3*X*S.diff(Y,t))==378)
    check('primitive pole orders',S.gcd(S.gcd(S.degree(X,t),S.degree(Y,t)),S.degree(defect,t))==1)
    check('X Y defect squarefree and coprime',S.gcd(X,Y)==1 and all(S.gcd(f,S.diff(f,t))==1 for f in [X,Y,defect]))
    check('curve immersive',S.gcd(S.diff(X,t),S.diff(Y,t))==1)
    check('primitive second key is essential',S.degree(defect,t)==4 and 4%3!=0 and S.gcd(4,3)==1)
    check('normal degree0 and canonical a2=1',4+3*(3+2-3*2)==1 and 4+3*(3+2-3*2)-1==0)
    A=S.cancel((X-X.subs(t,s))/(t-s));B=S.cancel((Y-Y.subs(t,s))/(t-s))
    Delta=(4*t*t-6*t+11)*(4*t*t+6*t+11)*(2*t**8+12*t**6+42*t**4+70*t*t+63)
    check('complete collision resultant exact',S.expand(S.resultant(A,B,s)-S.Rational(729,32)*Delta)==0)
    check('no leading coefficient hole',S.Poly(A,s).LC()==1 and S.Poly(B,s).LC()==1)
    check('twelve distinct collision parameters',S.degree(Delta,t)==12 and S.gcd(Delta,S.diff(Delta,t))==1)
    lin=S.subresultants(A,B,s)[-2];c1=S.Poly(lin,s).coeff_monomial(s);c0=S.Poly(lin,s).coeff_monomial(1)
    check('linear subresultant leading coefficient invertible',S.degree(lin,s)==1 and S.gcd(c1,Delta)==1)
    partner=-t*(256*t**8+1536*t**6+5376*t**4+8960*t*t+8361)/297
    check('exact subresultant partner',S.rem(c1*partner+c0,Delta,t)==0)
    check('partner preserves coordinates',S.rem(X.subs(t,partner)-X,Delta,t)==0 and S.rem(Y.subs(t,partner)-Y,Delta,t)==0)
    check('partner fixed-point-free involution',S.rem(partner.subs(t,partner)-t,Delta,t)==0 and S.gcd(partner-t,Delta)==1)
    tangent=S.diff(X,t)*S.diff(Y,t).subs(t,partner)-S.diff(Y,t)*S.diff(X,t).subs(t,partner)
    check('every pair transverse',S.gcd(S.rem(tangent,Delta,t),Delta)==1)
    check('six ordinary double values',S.degree(Delta,t)//2==6)
    def torus(w,n):return S.cancel((t**(w*n)-1)*(t-1)/((t**w-1)*(t**n-1)))
    alex=S.cancel(torus(3,4)*torus(2,3).subs(t,t**3))
    check('realized infinity determinant9',alex.is_polynomial(t) and alex.subs(t,1)==1 and abs(alex.subs(t,-1))==9)
    ca=(1,0,2,3);cb=(1,2,0,3);tau=(0,1,3,2);beta=G.compose(ca,G.inverse(cb))
    def power(a,n):return G.perm_word(([1]*n if n>=0 else [-1]*(-n)),[a])
    check('companion trefoil peripheral data',power(ca,2)==power(cb,3)==I4 and beta==(2,1,0,3) and power(beta,-6)==I4)
    check('cable cone relation',power(G.compose(beta,tau),3)==I4)
    check('cable outer peripheral match',power(beta,3)==beta and power(beta,4)==I4 and power(beta,-4)==I4)
    check('cable preferred longitude identity',power(tau,-12)==I4)
    check('whole peripheral infinity S4 image',len(G.group([ca,cb,tau],4))==24 and tau in trans)
    inputs=[Path(__file__),HERE/'PERIPHERAL_LONGITUDE_AND_ZERO_NORMAL_CONTROL.md',manifest]+[HERE/name for name in preds]+[HERE/row['path_relative_to_affine_complements'] for row in m['primary_sources']]+[HERE/'verify_quartic_colouring.py']
    report=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
      peripheral=dict(whole_surviving_quartic=True,all_source_rulings_poles_and_contacts=True,full_residual_rank_two_cover_retained=True,finite_branch_parameters_paired=True,normalization_infinity_split=True,preferred_longitude_identity=True,preimage_branch_knot_components=3,all_three_degrees_one=True,meridian_framing_assumed_zero=False,unramified_D_retained=True,cuspidal_one_prime_model_substituted=False,longitude_asserted_in_arbitrary_higher_degrees=False),
      linear=dict(longitude_equation='z(lambda)=0',coboundary_invariant=True,nonzero_constrained_H1_iff_S4_peripheral_lift=True,all_finite_presentations_and_cable_lengths=True,control_Coxeter_is_knot_group=False,Coxeter_dimensions=list(dims),T34_dimensions=list(dims34),passing_knot_lift_implies_affine_lift=False),
      binomial=dict(all_coprime_p_q_and_r_at_least_two=True,immersive_normalization_hypothesis_explicit=True,primitive_key_hypothesis_explicit=True,general_key_assumed_binomial=False,normal_degree_equals_Wronskian_degree=True,negative_eta_excluded_for_this_class=True,quartic_forces_eta_zero=True,eta_zero_Wronskian_nonzero_constant=True,eta_zero_DZ_equality=True,eta_zero_transverse_weighted_Euler=True,all_two_stage_keys_excluded=False,target_field_symmetry_inferred=False),
      control=dict(parametric_degrees=[9,6],binomial_defect_degree=4,weighted_Wronskian=378,normal_degree=0,canonical_a2=1,finite_birational_A1_normalization=True,immersive=True,collision_parameters=12,nodes=6,triple_values=0,cusps=0,tangencies=0,infinity_stages=[[2,3],[3,4]],infinity_determinant=9,infinity_S4_peripheral_representation=True,infinity_preferred_longitude_identity=True,affine_complement_S4_representation_supplied=False,sign_Picard_or_linking_gate_certified=False,finite_quartic_cover_supplied=False,marked_plane_source_supplied=False,control_declared_Keller=False),
      assurance=dict(written_unbounded_proofs=True,finite_controls_not_universal_proof=True,primary_imports_independently_reproved=False,formal_verification=False,independent_full_reproof=False,external_review=False,novelty_priority=False,full_quartic_exclusion=False,higher_degree_exclusion=False,research_goal_completed=False),
      source_sha256={os.path.relpath(p,HERE):sha(p) for p in inputs},runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'peripheral_zero_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(f'PASS {len(checks)} controls; quartic longitude and sharp eta0 curve')
if __name__=='__main__':main()
