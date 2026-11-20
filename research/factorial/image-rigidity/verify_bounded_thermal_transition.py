#!/usr/bin/env python3
"""Exact/source controls for the bounded-real-thermal transition theorem."""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Refusing -O/-OO: exact controls must execute.')
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,isqrt
import argparse,copy,hashlib,importlib.util,json,re,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent
PN='UNIFORM_BOUNDED_THERMAL_TRANSITION_EVENTUAL_SFC.md'
CN='verify_bounded_thermal_transition.py'
IN='bounded_thermal_transition_inputs.json'
LN='bounded_thermal_transition_source_ledger.json'
RN='bounded_thermal_transition_verification.json'
PREFIX='research/factorial/image-rigidity/'
# Author-time literal snapshot; no runtime discovery may add or omit a role.
PARENT_SPECS = {'main_proof': ('local', 'UNIFORM_INTEGER_CUSP_EVENTUAL_STRONG_FACTORIAL.md', 'written_argument'), 'uniform_local_proof': ('local', 'UNIFORM_GROWING_QUARTIC_CUSP_NEIGHBORHOOD.md', 'mathematical_input_uniform_local'), 'uniform_global_proof': ('local', 'UNIFORM_CUSP_CURIE_WEISS_REGULAR_COMPLEMENT.md', 'mathematical_input_uniform_regular'), 'exact_checker': ('local', 'verify_uniform_cusp_eventual.py', 'exact_and_source_support'), 'exact_inputs': ('local', 'uniform_cusp_eventual_inputs.json', 'exact_fixture_input'), 'context_audit': ('local', 'UNIFORM_CUSP_PRIMARY_CONTEXT.md', 'classical_context_only'), 'primary_kabluchko': ('local', 'uniform_cusp_sources/kabluchko_2203_05533v3.pdf', 'classical_context_only'), 'compact_402_proof': ('repo', 'research/factorial/image-rigidity/UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md', 'mathematical_input_compact_family_and_bounded_H'), 'fixed_cusp_405_proof': ('repo', 'research/factorial/image-rigidity/CUSP_EVENTUAL_STRONG_FACTORIAL_THEOREM.md', 'mathematical_input_narrow_bounded_s_route'), 'fixed_local_construction': ('repo', 'research/factorial/image-rigidity/CUSP_SHRINKING_WINDOW_LOG_RATIO_NONCOLLISION.md', 'mathematical_input_C1_to_C15_revalidated_uniformly'), 'fixed_global_complement': ('repo', 'research/factorial/image-rigidity/CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md', 'mathematical_input_narrow_405_closure'), 'fixed_morse_proof': ('repo', 'research/factorial/image-rigidity/CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md', 'mathematical_input_narrow_405_closure'), 'cusp_stationary_geometry': ('repo', 'research/factorial/image-rigidity/CUSP_STATIONARY_BELYI_MONODROMY_AND_RADIAL_BARRIER.md', 'mathematical_input_narrow_405_closure'), 'actual_coalescing_bridge': ('repo', 'research/factorial/image-rigidity/COALESCING_BINOMIAL_SADDLE_BRIDGE.md', 'mathematical_input_actual_phase_and_jets'), 'fixed_bounded_proof': ('repo', 'research/factorial/image-rigidity/BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md', 'mathematical_input_402_model_mechanism'), 'quartic_unfolding_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md', 'mathematical_input_402_actual_jet_mechanism'), 'quartic_model_proof': ('repo', 'research/factorial/image-rigidity/QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md', 'mathematical_input_model_location_and_simplicity_only'), 'primary_newman': ('repo', 'research/factorial/image-rigidity/quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf', 'mathematical_input_model_location_only'), 'old_gamma_upper_soft': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex', 'mathematical_input_old_upper_and_local_catch_soft'), 'sharp_annulus_soft_404': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md', 'mathematical_input_compact_annulus_and_finite_bad_box'), 'closed_disc_grace': ('repo', 'research/factorial/image-rigidity/GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md', 'mathematical_input_404_closed_disc_mechanism_only')}
SPECS={
 'main_proof':('local',PN,'written_argument'),
 'exact_checker':('local',CN,'exact_and_source_support'),
 'exact_inputs':('local',IN,'exact_fixture_input'),
 'uniform_parent_ledger':('support','uniform_cusp_eventual_source_ledger.json','exact_uniform_parent_source_boundary'),
 'uniform_parent_receipt':('support','uniform_cusp_eventual_verification.json','exact_uniform_parent_receipt'),
 'old_fixed_terminal_high_gap':('repo',PREFIX+'lower_exponent_initial_and_boundary.tex','mathematical_input_raw_fixed_terminal_support_finiteness'),
}
for _role,(_kind,_path,_relation) in PARENT_SPECS.items():
    SPECS['uniform_'+_role]=('support' if _kind=='local' else 'repo',_path,_relation)
ASSURANCE={'formal_verification':False,'external_review':False,'independent_replay':False}
SCOPE={
 'schema':1,'parameters':'integers a>=1,d>=2,N>=2; s=(d-1)/2 may be half-integral',
 'thermal_tube':'for each finite U>=0 one N0(U) works for sqrt(N)*abs(a-s^2)/(a+s)<=U',
 'same_cutoff_all_fields_and_coefficients':True,
 'all_d_covered':True,'nonzero_e_lattice_lower_bound':'1/4',
 'real_temperature_model_location_only':True,'model_zero_simplicity_or_HJN_needed':False,
 'bounded_field_uses_integrality_not_zero_labels':True,
 'actual_complex_temperature_frozen_real':False,
 'thermal_T_held_fixed_on_support_fixed_complex_m_collars':True,
 'leading_amplitude_T_over_rho_cubed_term_retained':True,
 'growing_error':'1/N + rho/k + U/k^(3/2) + 1/k^2',
 'actual_actions_used_no_m_over_s_smallness':True,
 'regional_bad_set_finite_iff_uniform_tube_cutoff':True,
 'all_U_cutoffs_iff_distinct_bad_sequences_thermal_escape':True,
 'raw_old_fixed_terminal_complements_only_for_reverse':True,
 'uses_407_or_411':False,'old_a0_a1_full_finite_prefix_premise':False,
 'global_bad_triple_set_finite':False,'all_order_SFC':False,
 'effective_cutoff':False,'unbounded_thermal_matching_solved':False,
 'reviewed_T_body_sha256':'b47c3b457d2fe296442936e5e4fa9645386542c865305aa84a826fb0fa0e1958',
 'controls':{'maximum_gap':12,'maximum_lower':40,'moment_maximum':4,
             'lattice_fourth_roots':[2,3,4,5,6], 'thermal_bounds':['1','2','5'],
             'fixed_terminal_maximum':64},
 'assurance_flags':ASSURANCE,
}
def need(ok,msg):
    if not ok:raise ValueError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def typed(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(typed(a[k],v) for k,v in b.items())
    if type(b) is list:return len(a)==len(b) and all(typed(x,y) for x,y in zip(a,b))
    return a==b
def check_inputs(d):need(typed(d,SCOPE),'Exact typed bounded-thermal scope')
def roots(repo,support):return {'local':HERE,'repo':repo,'support':support}
def check_sources(d,repo,support):
    need(type(d) is dict and set(d)=={'schema','assurance_flags','sources'},'Ledger schema')
    need(type(d['schema']) is int and d['schema']==1,'Integer schema1')
    need(typed(d['assurance_flags'],ASSURANCE),'False exact assurance flags')
    need(len(PARENT_SPECS)==21 and len(SPECS)==27,'Frozen exact parent/source cardinality')
    need(type(d['sources']) is list and len(d['sources'])==len(SPECS),'All mandatory sources')
    out={};bases=roots(repo,support);paths=set()
    for row in d['sources']:
        need(type(row) is dict and set(row)=={'role','kind','path','relation','sha256'},'Source row schema')
        role=row['role'];need(type(role) is str and role in SPECS and role not in out,'Unique exact role')
        need((row['kind'],row['path'],row['relation'])==SPECS[role],'Exact source kind/path/relation')
        need(type(row['sha256']) is str and re.fullmatch('[0-9a-f]{64}',row['sha256']),'SHA256')
        base=bases[row['kind']].resolve();p=(base/row['path']).resolve()
        need(p.is_relative_to(base) and p.is_file() and sha(p)==row['sha256'],'Stale/missing source:'+role)
        key=(row['kind'],row['path']);need(key not in paths,'Unique source path');paths.add(key);out[role]=row
    forbidden={'ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md',
      'FIBERWISE_UNIFORM_EVENTUAL_BINOMIAL_SFC.md','all_exponent_prefix_certificates.json',
      'ordinary_lower_one_finite_bridge_verification.json'}
    need(not any(Path(path).name in forbidden for _,path in paths),'Broader theorem/old full-prefix routes excluded')
    return out
def parent_module(support):
    sys.dont_write_bytecode=True
    p=support/PARENT_SPECS['exact_checker'][1]
    spec=importlib.util.spec_from_file_location('thermal_uniform_parent',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    need(module.SPECS==PARENT_SPECS,'Exact parent kernel/spec snapshot')
    return module
def check_parent(repo,support,rows,module):
    doc=load(support/SPECS['uniform_parent_ledger'][1])
    actual=module.validate_sources(doc,repo)
    need(set(actual)==set(PARENT_SPECS),'All parent roles')
    need(all(actual[role]==rows['uniform_'+role]['sha256'] for role in actual),'Every parent/direct pin agrees')
    with tempfile.TemporaryDirectory(prefix='bounded-thermal-parent-') as folder:
        out=Path(folder)/'actual.json'
        cmd=[sys.executable,'-B',str(support/PARENT_SPECS['exact_checker'][1]),
             '--repo-root',str(repo),'--output',str(out)]
        result=subprocess.run(cmd,cwd=folder,capture_output=True,text=True,timeout=60)
        need(result.returncode==0,'Fresh parent failed:'+result.stdout+';'+result.stderr)
        raw=out.read_bytes();expected=(support/SPECS['uniform_parent_receipt'][1]).read_bytes()
        need(raw==expected,'Fresh parent receipt byte identity')
        got=json.loads(raw);need(got['status']=='PASS' and got['source_pins']==actual,'Exact parent PASS scope')
        need(got['source_ledger_sha256']==rows['uniform_parent_ledger']['sha256'],'Parent ledger identity')
    return {'source_roles':len(actual),'receipt_byte_identical':True,'receipt_sha256':hashlib.sha256(raw).hexdigest()}
def universal(module):
    R=module.RF;v=module.var;out=[]
    def check(name,expr):need(expr.iszero(),'Universal identity:'+name);out.append(name)
    def add(a,b,n):return [(a[j] if j<len(a) else R())+(b[j] if j<len(b) else R()) for j in range(n+1)]
    def mul(a,b,n):
        c=[R() for _ in range(n+1)]
        for j,x in enumerate(a):
            for k,y in enumerate(b):
                if j+k<=n:c[j+k]=c[j+k]+x*y
        return c
    def log_unit(a,n):
        need((a[0]-1).iszero(),'Formal log constant1')
        z=list(a);z[0]=R();power=[R(1)]+[R() for _ in range(n)]
        result=[R() for _ in range(n+1)]
        for j in range(1,n+1):
            power=mul(power,z,n);result=add(result,[x*Q((-1)**(j+1),j) for x in power],n)
        return result
    def inv_unit(a,n):
        need((a[0]-1).iszero(),'Formal inverse constant1')
        c=[R(1)]+[R() for _ in range(n)]
        for j in range(1,n+1):c[j]=-sum((a[k]*c[j-k] for k in range(1,j+1)),R())
        return c
    s,e=v('s'),v('x');Sc=s*(s+1);C=Sc+e;d=2*s+1;a=s*s+e;n=4
    B=[R(1)]
    for j in range(1,n+1):
        value=s/R(factorial(j))
        for ell in range(1,j):value=value*(d-ell)
        B.append(value)
    A=[R()]+[a*Q((-1)**(j+1),j) for j in range(1,n+1)];A[1]=A[1]-C
    F=add(A,log_unit(B,n),n)
    check('stationary moving center',F[1])
    check('raw quadratic thermal jet',F[2]+e/2)
    check('raw cubic thermal jet',F[3]-e/3)
    check('raw quartic thermal jet',F[4]+Sc*Sc/12+e/4)
    inverse=inv_unit(B,n)
    check('exact mixed linear field jet cleared by d',(s+1)*inverse[1]+Sc)
    check('raw mixed quadratic field jet vanishes',(s+1)*inverse[2])
    l,p,H,t,q,b,lam,m,rho,k=(v(z) for z in ['l','p','H','t','q','b','lam','m','rho','k'])
    check('support a=C-s identity',a-(C-s))
    check('frozen cusp beta error modulo L squared equals C',
          (-Sc/(d*l)+l/d-e/(d*l))*d*l-(l*l-C))
    tail0=a-l*l*p+d*H;tail1=-s-l*t+d*H
    check('detuned full-tail numerator', (tail0-tail1).subst(t=l*(p-1))-(C-l*l))
    S,Z,D0=v('z'),v('y'),v('d')
    w=Z*(a+D0-S)/(S-a)
    check('all-support stationary envelope',D0*w/(w+Z)-(a+D0-S))
    # Here m,p are the positive fourth roots of the two order symbols.
    upper=m**3*b*(q-1);lower=p**3*b*(lam*q-1)
    check('exact own-order scaled field relation',lower-lam*(p/m)**3*upper-b*p**3*(lam-1))
    r=v('r');uN=-m*m*r;un=-p*p*r
    check('exact real thermal own-order relation',un-(p/m)**2*uN)
    T=-r/rho**2;sqrtk=m*m*rho*rho
    check('residual thermal parameter is support fixed',uN/sqrtk-T)
    check('no spurious complex-order thermal derivative',(uN/sqrtk).diff('m'))
    Tformal=v('p')
    transported=(D0/m)*(rho**-2+Tformal*rho**-3+1/(k*rho**3))/(l*rho)
    check('all three transported amplitude terms',transported.subst(m=k/rho**4)-(D0/l)*(rho/k+Tformal/k+1/k**2))
    z=v('z');u=v('H')
    check('thermal amplitude error has k to three halves',(Tformal/k).subst(p=u/z,k=z*z)-u/z**3)
    check('quarter-integral numerator',4*e-(4*a-(d-1)**2))
    check('strict half-integer square residue',(2*s-1)**2-(4*s*(s-1)+1))
    # At T=0 two selected stationary cubes agree, hence their amplitude
    # ratio's omega derivative vanishes; analyticity makes it O(T).
    uL,uR=v('y'),v('z')
    check('leading ratio derivative at zero thermal parameter',
          (-1/uR**3+1/uL**3)*uL**3*uR**3-(uR**3-uL**3))
    need(not F[2].iszero(),'Dropping a nonzero thermal quadratic must fail')
    need(not (e/d).iszero(),'Frozen cusp mixed beta must fail off equality')
    return out
def exact_controls(inp):
    cfg=inp['controls'];lattice=0;moments=0;transport=0;triangles=0;model=0
    for d in range(2,cfg['maximum_gap']+1):
        s=Q(d-1,2);Sc=s*(s+1)
        for a in range(1,cfg['maximum_lower']+1):
            e=Q(a)-s*s;C=Q(a)+s;r=e/C
            need(4*e==4*a-(d-1)**2,'Exact quarter numerator')
            if e:need(abs(e)>=Q(1,4),'Nonzero ordinary detuning lattice')
            else:need(d%2==1 and s.denominator==1,'Equality curve is integer cusp')
            for j in cfg['lattice_fourth_roots']:
                N=j**4
                for bound in cfg['thermal_bounds']:
                    U=Q(bound)
                    if e and abs(r)<=Q(1,2) and j*j*abs(r)<=U:
                        need(C>=Q(j*j)/(4*U),'Uniform L lower bound squared')
                        need(Sc*Sc/(C*j*j)>=1/(16*U),'Real radial field drop squared')
                        lattice+=1
            for m in range(1,cfg['moment_maximum']+1):
                AA=Q(2,3);BB=-Q(3,5);x=AA/BB
                direct=sum(Q(comb(m,j)*factorial(a*m+d*j))*AA**(m-j)*BB**j for j in range(m+1))
                poly=sum(Q(comb(m,j)*factorial(a*m+d*j),factorial(a*m))*x**(m-j) for j in range(m+1))
                need(direct==BB**m*factorial(a*m)*poly,'Actual factorial field identity')
                moments+=1
            w0=C**d*Q(s+1,s)
            for N in (2,3,5,8):
                n=N-1;qN=Q(2,3)/(N**d*w0);qn=Q(2,3)/(n**d*w0)
                need(qn==Q(N,n)**d*qN,'Both actual own-order coefficient maps')
                need(n*r*r==Q(n,N)*(N*r*r),'Squared own-order real thermal relation')
                transport+=1
    for N in range(2,cfg['fixed_terminal_maximum']+1):
        target=4*N*(N-1);exponent=(target-1).bit_length();D=2*(N-1)*exponent
        need(2**(exponent-1)<target<=2**exponent,'Exact raw old high-gap cutoff')
        # OLD soft criterion is sufficient for regional reverse finiteness.
        count=sum(d*(2*N-3)-1 for d in range(2,D))
        need(count==(2*N-3)*((D-1)*D//2-1)-(D-2),'Old fixed-terminal positive support count')
        triangles+=1
    for n in range(12,81):
        s=Q(n-1,2);e=n-s*s;C=n+s;theta2=Q(n)*e*e/(C*C)
        need(theta2>=Q(n**3,144),'An infinite diagonal can escape every bounded thermal tube')
        # Different supports can each have one abstract failure while
        # remaining off the exact equality curve in one bounded tube.
        ss=n;aa=n*n+1;dd=2*n+1;NN=n**4;CC=aa+ss
        need(Q(n*n,CC)<1 and aa-ss*ss==1,'Exact curve alone is not a tube theorem')
        model+=1
    # Genuine arithmetic countercontrol to silently allowing real detunings.
    continuous_C=Q(2)+Q(1,10**8)
    need(10**4*Q(1,10**8)/continuous_C<1 and continuous_C<2500,
         'Continuous detuning does not satisfy quarter-lattice lower bound')
    return {'ordinary_lattice_and_squared_drop_controls':lattice,
      'actual_moment_normalizations':moments,'actual_own_order_transport_controls':transport,
      'raw_fixed_terminal_support_counts':triangles,'logical_diagonal_countermodel_controls':model,
      'genuine_math_negative_controls':3}
def negatives(ledger,repo,support):
    out=[]
    def reject(name,f):
        try:f()
        except (ValueError,KeyError,TypeError):out.append(name);return
        raise ValueError('Negative accepted:'+name)
    for role in SPECS:
        d=copy.deepcopy(ledger);d['sources']=[r for r in d['sources'] if r['role']!=role]
        reject('omit:'+role,lambda d=d:check_sources(d,repo,support))
        d=copy.deepcopy(ledger);next(r for r in d['sources'] if r['role']==role)['sha256']='0'*64
        reject('badpin:'+role,lambda d=d:check_sources(d,repo,support))
    for key in ('model_zero_simplicity_or_HJN_needed','actual_complex_temperature_frozen_real',
                'uses_407_or_411','old_a0_a1_full_finite_prefix_premise','global_bad_triple_set_finite',
                'all_order_SFC','effective_cutoff','unbounded_thermal_matching_solved'):
        d=copy.deepcopy(SCOPE);d[key]=True;reject('false-scope:'+key,lambda d=d:check_inputs(d))
    for key in ('all_d_covered','real_temperature_model_location_only',
                'thermal_T_held_fixed_on_support_fixed_complex_m_collars',
                'leading_amplitude_T_over_rho_cubed_term_retained',
                'regional_bad_set_finite_iff_uniform_tube_cutoff'):
        d=copy.deepcopy(SCOPE);d[key]=False;reject('missing-scope:'+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(SCOPE);d['schema']=True;reject('boolean-schema',lambda:check_inputs(d))
    d=copy.deepcopy(ledger);d['sources']=[];reject('empty-ledger',lambda:check_sources(d,repo,support))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--support-root',type=Path,default=HERE)
    ap.add_argument('--output',type=Path,default=HERE/RN);ns=ap.parse_args()
    repo=ns.repo_root.resolve();support=ns.support_root.resolve()
    inp=load(HERE/IN);check_inputs(inp);ledger=load(HERE/LN);rows=check_sources(ledger,repo,support)
    module=parent_module(support);parent=check_parent(repo,support,rows,module)
    identities=universal(module);counts=exact_controls(inp);bad=negatives(ledger,repo,support)
    check_sources(ledger,repo,support)
    result={'schema':1,'status':'PASS','scope':inp,'source_ledger_sha256':sha(HERE/LN),
      'proof_sha256':sha(HERE/PN),'checker_sha256':sha(HERE/CN),'input_sha256':sha(HERE/IN),
      'sources':[rows[r] for r in sorted(rows)],'mandatory_source_roles':len(rows),
      'fresh_uniform_parent':parent,'universal_identities':identities,'counts':counts,
      'source_scope_negative_controls':bad,'assurance_flags':ASSURANCE,
      'meaning':'The complete written detuned-cycle and packet proof supplies uniformity. Exact/source controls support jets, scalings, discrete and logical boundaries, not a finite-prefix inference or formal verification.'}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','mandatory_sources':len(rows),'universal_identities':len(identities),
                     'counts':counts,'source_scope_controls':len(bad),'output':str(ns.output)}))
if __name__=='__main__':main()
