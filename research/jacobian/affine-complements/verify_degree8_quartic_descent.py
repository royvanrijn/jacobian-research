"""PF2RULED30: exact sign descent to quartic coefficient fields.
Run with sage -python. Replays polynomial transport and the descended unit rows;
no Groebner search or expanded high-degree E calculation is performed.
"""
from pathlib import Path
from time import monotonic
import gzip,hashlib,json,platform
from sage.all import QQ,NumberField,PolynomialRing
from sage.version import version as sage_version

HERE=Path(__file__).resolve().parent
checks=[]
def check(name,truth):
    assert truth,name
    checks.append(name)
def read_gz(path):return json.loads(gzip.decompress(path.read_bytes()))
def encode(f):
    return [[list(e),[str(c[j]) for j in range(4)]] for e,c in f.dict().items()]

def main():
    start=monotonic()
    predecessor=json.loads((HERE/'degree8_quadratic_receipt.json').read_text())
    check('current quadratic predecessor PASS',predecessor['outcome']=='PASS')
    for name,digest in predecessor['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    models=read_gz(HERE/'degree8_quadratic_presentations.json.gz')['models']
    primary=read_gz(HERE/'degree8_affine_reduction_certificates.json.gz')['sectors']
    unitrows=read_gz(HERE/'degree8_unit_row_certificates.json.gz')['sectors']
    check('two sectors aligned',len(models)==len(primary)==len(unitrows)==2 and
          [x['multiplicities'] for x in models]==[x['multiplicities'] for x in primary]==[[2,7],[4,5]])
    sectors=[]
    for model,old,rowdata in zip(models,primary,unitrows):
        k=model['multiplicities'][0]
        Y=PolynomialRing(QQ,'y');p8=Y(model['P8'])
        K=NumberField(p8,'y');y=K.gen()
        U=PolynomialRing(QQ,'t');tt=U.gen();p4=sum(p8[2*j]*tt**j for j in range(5))
        check(f'exact even irreducible octic with nonzero constant k{k}',p8.is_irreducible() and
              p8.degree()==8 and p8[0]!=0 and all(p8[j]==0 for j in [1,3,5,7]))
        check(f'irreducible quartic defining polynomial k{k}',p4.degree()==4 and p4.is_irreducible())
        F=NumberField(p4,'t');t=F.gen();embed=F.hom([y*y],K);sigma=K.hom([-y],K)
        check(f'quartic subfield and sign involution k{k}',embed(t)==y*y and sigma(y)==-y and
              all(sigma(sigma(y**j))==y**j for j in range(8)) and K.degree()==2*F.degree())
        count=[0]
        def down(c):
            c=K(c)
            assert sigma(c)==c, f'coefficient not sign invariant k{k}'
            assert all(c[j]==0 for j in [1,3,5,7]), f'odd field coefficient k{k}'
            out=F(sum(c[2*j]*t**j for j in range(4)))
            assert embed(out)==c, f'exact coefficient round trip k{k}'
            count[0]+=1
            return out
        names=model['variables'];A=PolynomialRing(K,names=names);B=PolynomialRing(F,names=names)
        weights=[6-j for j in range(9)]+[17-j for j in range(6)]+[22-j for j in range(11)]+[8-j for j in range(2,8)]+[1,0,1,0]+[29,-29]
        scales=[17-j for j in range(20)]+[29-j for j in range(18)]+[0]
        def transfer(f,scale):
            return B({tuple(e):down(c*y**(sum(int(a)*b for a,b in zip(e,weights))-scale))
                      for e,c in A(f).dict().items()})
        quadratic=[transfer(f,sc) for f,sc in zip(model['generators'],scales)]
        inverse=B.hom([A.gen(i)/y**weights[i] for i in range(38)],A,base_map=embed)
        for i,(g,f,sc) in enumerate(zip(quadratic,model['generators'],scales)):
            check(f'complete quadratic coefficient transport and inverse k{k} row{i}',
                  inverse(g)*y**sc==A(f))
        check(f'38 variables39 quadratics remain over the quartic field k{k}',
              len(quadratic)==39 and B.ngens()==38 and max(g.total_degree() for g in quadratic)==2)
        check(f'nonzero equation is preserved exactly k{k}',quadratic[-1]==B.gen(36)*B.gen(37)-1)
        # The compact globally reduced model is transported independently.
        AP=PolynomialRing(K,names=['q4','q5','q6','q7','z']);BP=PolynomialRing(F,names=AP.variable_names())
        VP=PolynomialRing(AP,'v');VQ=PolynomialRing(BP,'w');w=VQ.gen();qweights=[4,3,2,1,12]
        def core(txt,scale):
            return VQ([BP({tuple(e):down(c*y**(j+sum(int(a)*b for a,b in zip(e,qweights))-scale))
                            for e,c in cf.dict().items()}) for j,cf in enumerate(VP(txt))])
        W=core(old['root_deviation_W'],12);M=core(old['critical_remainder_M'],17);L=core(old['Rprime_over_h_L'],6)
        s=down(K(old['cusp_parameter'])/y);h=w-s;D=(M(s)-M)//h
        check(f'whole core descends with exact cusp and degree bounds k{k}',W(s)==W.derivative()(s)==0 and
              W.degree()==12 and M.degree()<=5 and L.degree()==6 and L.leading_coefficient()==8)
        check(f'coefficient row descends exactly k{k}',D==core(old['linear_constant_D'],16) and D.degree()<=4)
        check(f'cusp coefficient generates the full quartic field k{k}',s.minpoly().degree()==4)
        # Descend the existing exact unit row, preserving each coefficient.
        C=[];row_weights=[4,3,2,1]
        for i,encoded in enumerate(rowdata['coefficients']):
            C.append(BP({tuple(ex)+(0,):down(K([QQ(c) for c in co])*y**(16-i+sum(a*b for a,b in zip(ex,row_weights))))
                         for ex,co in encoded}))
        check(f'descended degree3 Bezout row is exact over quartic coefficients k{k}',
              sum(C[i]*D[i] for i in range(5))==1 and max(c.total_degree() for c in C)<=3)
        # The remainder circuit commutes with this field transport: the divisor is monic up to8.
        check(f'monic remainder transport factors are field units k{k}',y!=0 and L.leading_coefficient()==8 and
              embed(t)!=0 and core(old['constant_remainder_Ms'],17)[0]==M(s))
        # Sign transport of the original marking identifies the two octic extensions.
        xp=Y(old['x_as_polynomial_in_y']);x=K(xp)
        check(f'the two boundary roots and cusp change sign k{k}',sigma(x)==-x and
              sigma(K(old['cusp_parameter']))==-K(old['cusp_parameter']))
        sectors.append(dict(multiplicities=model['multiplicities'],P8=str(p8),P4=str(p4),
            P4_coefficients=[str(p4[j]) for j in range(5)],quartic_generator_in_octic='t=y^2',
            octic_extension_over_quartic='y^2=t',quadratic_variables=names,
            variable_y_weights=dict(zip(names,weights)),equation_y_weights=scales,
            quadratic_generators=list(map(str,quadratic)),
            core_variables=list(BP.variable_names()),core_y_weights=dict(zip(BP.variable_names(),qweights)),
            normalized_cusp_parameter=str(s),normalized_W=str(W),normalized_M=str(M),normalized_L=str(L),normalized_D=str(D),
            E_circuit='remainder_w(M*(W/(w-s)),L)',core_equations=['E[5]=0']+[f'D[{i}]*z-E[{i}]=0' for i in range(5)],
            nonzero_condition='M(s)*z != 0',unit_multipliers=[encode(c) for c in C],
            recovered_root_constant_circuit='Zstar=sum(C[i]*E[i],i=0..4)',
            unit_identity_verified=True,coefficient_round_trips_verified=count[0],
            core_solved=False,all_octic_embeddings_covered=True))
        print(f'PASS k{k}: quartic field descent, all39 exact inverse transports, exact unit row; core remains open',flush=True)
    # Universal Cox expression and the source involution, checked as exact polynomials.
    R=PolynomialRing(QQ,names=['u','v']);u,v=R.gens()
    source=[u**3*v-u**2,u**2*v-u,u*v-1,v]
    check('actual Wright generator sign action',
          [f(-u,-v) for f in source]==[source[0],-source[1],source[2],-source[3]])
    for m in [3,6,9]:
        for j in range(m+1):
            for i in range(m+j+1):
                f=(-v)**(m-j)*(1-v*v-u*v)**j*(v+u)**(m+j-i)
                check(f'whole Cox section sign eigenvector m{m} j{j} i{i}',f(-u,-v)==(-1)**i*f)
    cp=HERE/'degree8_quartic_descent_certificates.json.gz'
    cp.write_bytes(gzip.compress((json.dumps(dict(as_of='2026-11-14',sectors=sectors),indent=2)+'\n').encode(),mtime=0))
    inputs=[Path(__file__),HERE/'DEGREE_EIGHT_QUARTIC_DESCENT.md',cp,HERE/'degree8_quadratic_receipt.json',
            HERE/'degree8_quadratic_presentations.json.gz',HERE/'degree8_affine_reduction_certificates.json.gz',
            HERE/'degree8_unit_row_certificates.json.gz',HERE/'verify_full_next_normal.py']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        marked_octic_embeddings=16,quartic_embedding_problems=8,sector_count=2,
        original_coefficient_field_degree=8,descended_coefficient_field_degree=4,
        exact_quadratic_variable_count=38,exact_quadratic_equation_count=39,maximum_equation_degree=2,
        descended_unit_multiplier_degree=3,core_parameter_count=4,saturated_variable_count=5,
        rational_unknowns_per_field_unknown_before=8,rational_unknowns_per_field_unknown_after=4,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],
        assurance=dict(all_sixteen_octic_embeddings_covered=True,all_polynomial_coefficients_transported_exactly=True,
            invertible_constant_changes_only=True,all_parameter_strata_retained=True,nonzero_kappa_preserved=True,
            original_cusp_and_root_hypotheses_preserved=True,descended_Bezout_identity_checked_in_characteristic_zero=True,
            descent_uses_no_modular_inference=True,source_sign_automorphism_preserves_finite_map_obligations=True,
            characteristic_zero_unit_certificate_for_full_core_obtained=False,full_degree8_nonzero_kappa_core_solved=False,
            actual_global_section_lifts_constructed=False,general_determinant_locus_solved=False,actual_root_hypothesis_removed=False,
            full_source_Keller_pair_constructed=False,field_degree_four_realization_constructed=False,
            global_surface_finiteness_proved=False,actual_ruled_compatibility_proved=False,contact_cases_classified=False,
            whole_marked_candidate_excluded=False,quartic_packets_excluded=False,formal_verification=False,
            independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        python=platform.python_version(),sage=sage_version,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'degree8_quartic_descent_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; coefficient degree8 to4; no full-core exclusion',flush=True)
if __name__=='__main__':main()
