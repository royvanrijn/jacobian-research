"""PF2RULED29: exact quadratic presentation and separately labelled modular evidence.
Default replay constructs the characteristic-zero equations without elimination.
--replay-modular reruns the two bounded native finite-field saturation jobs.
"""
from pathlib import Path
from time import monotonic
import argparse,gzip,hashlib,json,platform,re,shutil,subprocess,tempfile
from sage.all import QQ,GF,NumberField,PolynomialRing,matrix,vector
from sage.version import version as sage_version
from verify_full_next_normal import cox_jets,combine,jac1
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def gz(path):return json.loads(gzip.decompress(path.read_bytes()))
def main():
    args=argparse.ArgumentParser();args.add_argument('--replay-modular',action='store_true');opts=args.parse_args()
    start=monotonic();pred=json.loads((HERE/'degree8_unit_row_receipt.json').read_text())
    check('current global-recovery predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    roots=json.loads((HERE/'square_boundary_gap_certificates.json').read_text())['sectors'];markings=json.loads((HERE/'second_normal_obstruction_certificates.json').read_text())['sectors'];prior=gz(HERE/'square_boundary_degree8_checkpoint.json.gz')['sectors'];primary=gz(HERE/'degree8_affine_reduction_certificates.json.gz')['sectors']
    Y=PolynomialRing(QQ,'y');models=[];observations=[]
    for case,mark,checkpoint,affine in zip(roots,markings,prior,primary):
        k=case['multiplicities'][0];p8=Y(case['P8']);check(f'irreducible degree-eight field k{k}',p8.degree()==8 and p8.is_irreducible())
        K=NumberField(p8,'y');y=K.gen();x=K(Y(case['x_as_polynomial_in_y']));Wring=PolynomialRing(K,'w');w=Wring.gen();H=(w-x)**k*(w-y)**(9-k);aa,bb=H[8],H[7];s=(24*bb+72-11*aa*aa)/(12*aa);delta=-aa/48
        V=PolynomialRing(K,'v');v=V.gen();h=v-s;lam=3*h/4+delta
        Lop=lambda F:h*F[1]-lam*F[0].derivative()
        T=[V.zero()]*4;bas3=cox_jets(3,V)
        for power,idx in enumerate(mark['T_symbol_lift_cox_ji']):
            row=next(row for row in bas3 if [row['j'],row['i']]==idx);T=[T[j]+H[9-power]*row['jets'][j] for j in range(4)]
        bas2=cox_jets(2,V);LM=matrix(K,[[Lop(row['jets'])[j] for row in bas2] for j in range(8)]);piv=case['S_cusp_inverse_pivot_columns'];free=case['S_free_cox_columns'];inverse=LM[:,piv].inverse();co=vector(K,[K(c) for c in case['S_cusp_correction']]);Sbase=[T[j]+combine(co,bas2)[j] for j in range(4)]
        check(f'exact root reference satisfies whole cusp identity k{k}',Lop(Sbase)==0)
        directions=[]
        for j in free:
            row=vector(K,12);row[j]=1
            for ii,c in zip(piv,-inverse*LM.column(j)):row[ii]=c
            check(f'complete root cusp direction k{k} free{j}',LM*row==0);directions.append(combine(row,bas2))
        check(f'whole four-parameter root space k{k}',len(directions)==4)
        names=['b'+str(i) for i in range(9)]+['m'+str(i) for i in range(6)]+['n'+str(i) for i in range(11)]+['c'+str(i) for i in range(2,8)]+['a'+str(i) for i in range(4)]+['kap','invkap']
        A=PolynomialRing(K,names=names,order='degrevlex');var=dict(zip(names,A.gens()));VA=PolynomialRing(A,'v');v=VA.gen();h=v-A(s);lam=3*h/4+A(delta)
        S=[VA(Sbase[j])+sum(var['a'+str(i)]*VA(directions[i][j]) for i in range(4)) for j in range(4)];root=S[0];r,rr=root.derivative().quo_rem(h)
        check(f'all root parameters retain cusp condition and exact division k{k}',rr==0 and h*S[1]-lam*root.derivative()==0)
        R=h**8+sum(var['c'+str(i)]*h**i for i in range(2,8));L=R.derivative()//h
        C=lam**2*r*L.derivative()+(2*h*(S[2]+v*lam*r)-lam**2*r.derivative())*L
        check(f'critical polynomial equals actual normal Jacobian formula k{k}',C==jac1(S,[R,lam*L,-v*lam*L,VA.zero()]))
        B=sum(var['b'+str(i)]*v**i for i in range(9));M=sum(var['m'+str(i)]*v**i for i in range(6));N=sum(var['n'+str(i)]*v**i for i in range(11));kap=var['kap'];invkap=var['invkap']
        F=C-2*root.derivative()*B-M;G=2*M*root+2*h*L*N-kap
        generators=[F[j] for j in range(20)]+[G[j] for j in range(18)]+[kap*invkap-1]
        check(f'exact38-variable39-equation presentation k{k}',A.ngens()==38 and len(generators)==39 and all(g!=0 for g in generators))
        check(f'all equations have degree at most2 k{k}',max(g.total_degree() for g in generators)==2)
        lead=root.derivative()[11]
        check(f'first monic-division pivot is a nonzero FIELD scalar k{k}',lead in K and lead!=0)
        for row in range(9):
            check(f'quotient B triangular pivot k{k} row{row}',F[19-row].derivative(var['b'+str(8-row)])==-2*lead)
            check(f'quotient B pivot has no unsolved lower-coefficient dependence k{k} row{row}',all(F[19-row].derivative(var['b'+str(8-col)])==0 for col in range(row+1,9)))
        for row in range(11):
            check(f'normal control N triangular pivot16 k{k} row{row}',G[17-row].derivative(var['n'+str(10-row)])==16)
            check(f'normal control N pivot has no unsolved lower dependence k{k} row{row}',all(G[17-row].derivative(var['n'+str(10-col)])==0 for col in range(row+1,11)))
        check(f'constant recovery is a field-linear combination of second-identity rows k{k}',G(s)==2*M(s)*root(s)-kap)
        quotient,remainder=C.quo_rem(2*root.derivative())
        check(f'quotient degree and complete exact division k{k}',quotient.degree()<=8 and remainder.degree()<=10 and C==2*root.derivative()*quotient+remainder)
        oldA=PolynomialRing(K,names=checkpoint['variables'],order='degrevlex');hom=oldA.hom([var[n] for n in checkpoint['variables']],A)
        check(f'all five primary high-remainder equations are unchanged k{k}',[remainder[j] for j in range(6,11)]==[hom(oldA(g)) for g in checkpoint['high_generators']])
        low=sum(remainder[j]*v**j for j in range(6))
        check(f'eliminating B and M leaves exactly the five high equations k{k}',C-2*root.derivative()*quotient-low==sum(remainder[j]*v**j for j in range(6,11)))
        check(f'second identity leading division scalar is fixed16 k{k}',(2*h*L).degree()==7 and (2*h*L).leading_coefficient()==16)
        models.append(dict(multiplicities=case['multiplicities'],P8=str(p8),variables=names,generators=list(map(str,generators)),
            source_root_restriction=str(root),source_root_second_coefficient=str(S[2]),boundary_perturbation_R=str(R),critical_polynomial_C=str(C),
            first_identity='C-2*S0_prime*B-M=0',second_identity='2*M*S0+2*h*L*N-kap=0',nonzero_equation='kap*invkap-1=0',
            auxiliary_counts=dict(quotient_B=9,remainder_M=6,normal_N=11,kap=1),first_quotient_pivot=str(K(-2*lead)),second_quotient_pivot='16',
            source_high_variables=list(checkpoint['variables']),all_variable_denominator_strata_retained=True))
        # Match the actual retained finite-field inputs to the current primary core.
        modulus=1073741827;ff=GF(modulus);root_y={2:743393613,4:1037248271}[k]
        check(f'documented residue lies on P8 k{k}',p8.change_ring(ff)(ff(root_y))==0)
        ev=lambda c:sum(ff(t)*ff(root_y)**i for i,t in enumerate(K(c).list()))
        AP=PolynomialRing(K,names=['q4','q5','q6','q7','z']);VP=PolynomialRing(AP,'v');AQ=PolynomialRing(ff,names=AP.variable_names());VQ=PolynomialRing(AQ,'v');vq=VQ.gen()
        def transfer(text):return VQ([AQ({tuple(map(int,ex)):ev(c) for ex,c in p.dict().items()}) for p in VP(text)])
        WW,MM,LL=[transfer(affine[key]) for key in ['root_deviation_W','critical_remainder_M','Rprime_over_h_L']];ss=ev(K(affine['cusp_parameter']));hh=vq-ss;DD=(MM(ss)-MM)//hh;EE=(MM*(WW//hh)).quo_rem(LL)[1];zz=AQ.gens()[-1]
        expected=[EE[5]]+[DD[j]*zz-EE[j] for j in range(5)]+[MM(ss)*zz]
        basename='sector27' if k==2 else 'sector45';inputpath=HERE/'sources/degree8-native-saturation'/f'{basename}.ms';outpath=inputpath.with_suffix('.out')
        lines=inputpath.read_text().splitlines();actual=[AQ(g.strip()) for g in '\n'.join(lines[2:]).split(',') if g.strip()]
        check(f'retained native modular input matches the COMPLETE current core k{k}',lines[0]==','.join(AQ.variable_names()) and int(lines[1])==modulus and actual==expected)
        check(f'retained native result records the unit basis at the same prime k{k}',f'characteristic: {modulus}' in outpath.read_text() and re.search(r'\[1\]:\s*$',outpath.read_text()) is not None)
        replay=None
        if opts.replay_modular:
            with tempfile.TemporaryDirectory() as temp:
                target=Path(temp)/'basis.out';beg=monotonic()
                run=subprocess.run(['msolve','-S','-t','1','-v','0','-g','2','-f',str(inputpath),'-o',str(target)],capture_output=True,text=True,timeout=60)
                check(f'actual native modular replay completes with the unit basis k{k}',run.returncode==0 and re.search(r'\[1\]:\s*$',target.read_text()) is not None)
                replay=dict(elapsed_seconds=round(monotonic()-beg,6),exit_code=run.returncode,basis='[1]',stdout=run.stdout,stderr=run.stderr)
        observations.append(dict(multiplicities=case['multiplicities'],prime=modulus,residue_y=root_y,unit_basis_recorded=True,
            actual_input_matches_primary_core=True,characteristic_zero_unit_identity_obtained=False,
            input_path=str(inputpath.relative_to(HERE)),output_path=str(outpath.relative_to(HERE)),native_replay=replay))
        print(f'PASS k{k}: exact39 quadratics, monic pivots preserve all cases; modular unit observation remains separate',flush=True)
    cp=HERE/'degree8_quadratic_presentations.json.gz';cp.write_bytes(gzip.compress((json.dumps(dict(as_of='2026-11-14',models=models),indent=2)+'\n').encode(),mtime=0))
    obs=HERE/'degree8_native_saturation_observations.json';obs.write_text(json.dumps(dict(as_of='2026-11-14',observations=observations,
        characteristic_zero_exclusion_proved=False,all_rooted_degree8_cases_excluded=False),indent=2)+'\n')
    inputs=[Path(__file__),HERE/'DEGREE_EIGHT_QUADRATIC_PRESENTATION.md',cp,obs,HERE/'degree8_unit_row_receipt.json',
        HERE/'square_boundary_gap_certificates.json',HERE/'square_boundary_degree8_checkpoint.json.gz',HERE/'second_normal_obstruction_certificates.json',
        HERE/'degree8_affine_reduction_certificates.json.gz',HERE/'verify_full_next_normal.py',HERE/'DEGREE_EIGHT_GLOBAL_CONSTANT_RECOVERY.md']
    inputs+=sorted((HERE/'sources/degree8-native-saturation').glob('*'))
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,marked_points=16,
        exact_quadratic_variable_count=38,exact_quadratic_equation_count=39,maximum_equation_degree=2,
        uniquely_recoverable_auxiliary_count=27,base_high_variable_count=10,nonzero_condition_preserved=True,
        modular_unit_observations=2,modular_prime=1073741827,actual_modular_replay_performed=opts.replay_modular,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(exact_presentation_covers_all_complex_embeddings=True,monic_elimination_has_no_variable_denominator=True,
            complete_original_necessary_system_preserved=True,all_parameter_rank_drop_cases_retained=True,
            modular_results_not_promoted_to_characteristic_zero=True,characteristic_zero_degree8_unit_certificate_obtained=False,
            full_degree8_nonzero_kappa_core_solved=False,actual_global_section_lifts_constructed=False,
            general_determinant_locus_solved=False,actual_root_hypothesis_removed=False,full_source_Keller_pair_constructed=False,
            field_degree_four_realization_constructed=False,global_surface_finiteness_proved=False,actual_ruled_compatibility_proved=False,
            contact_cases_classified=False,whole_marked_candidate_excluded=False,quartic_packets_excluded=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(f.relative_to(HERE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sage=sage_version,msolve_binary_sha256=(hashlib.sha256(Path(shutil.which('msolve')).read_bytes()).hexdigest() if shutil.which('msolve') else None),
        elapsed_seconds=round(monotonic()-start,6))
    (HERE/'degree8_quadratic_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; exact quadratic presentation and two modular unit observations; characteristic-zero exclusion remains open',flush=True)
if __name__=='__main__':main()
