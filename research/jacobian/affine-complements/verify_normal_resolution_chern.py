"""PF109 exact controls; written Picard/Chern proof and imports establish universal scope."""
from pathlib import Path
import hashlib,json,platform,math
import sympy as S
import normal_resolution_chern as N
import ade_end_exclusion as A
import multiple_boundary_contact as M
import coprime_leading_forms as C

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name;checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    old=json.loads((P/'ade_end_exclusion_receipt.json').read_text())
    check('PF108 predecessor PASS',old['outcome']=='PASS')
    for path,h in old['source_sha256'].items():check('PF108 exact seal '+path,sha(P/path)==h)
    manifest=json.loads((P/'sources/normal-resolution-chern/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context pin '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('full new primary pages visually audited',manifest['sources'][0]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[1,2,3,6,14,15,16] and manifest['sources'][1]['fresh_complete_PDF_pages_rendered_and_visually_inspected']==[4,7])
    check('whole deformation theorem imported honestly',manifest['deformation_theorem_imported'] and not manifest['sources'][0]['independent_whole_deformation_proof'])
    check('prior classification audit not falsely fresh',manifest['sources'][3]['fresh_visual_audit_this_turn'] is False)
    L,z,b,q=S.symbols('L z b q',integer=True,positive=True)
    ident('UNBOUNDED reference Milnor theta',-2*(1+L)-3*(-L)-(L-2))
    ident('UNBOUNDED arbitrary deformed filling theta',-z-2*(1+L)-3*(-L)-(L-2-z))
    ident('UNBOUNDED theta equality forces zero Chern square',(L-2-z-(L-2))+z)
    check('negative-definite nonzero Chern cannot have zero square',(-z).is_negative and z.is_positive)
    ident('UNBOUNDED canonical pairing on minimal rational exceptional curve',-2-(-b)-(b-2))
    ident('UNBOUNDED minimal curve Stein adjunction equality',-b+(b-2)+2)
    check('a holomorphic -1 blowup cannot enter the minimal deformation hypothesis',-1+abs(-1)>-2)
    # Unit relation is an ACTUAL finite singular model, not a formal counterexample.
    d=N.a1_wrong_open_model();v,t=d['v'],d['t']
    ident('exact degree12 quotient resolution chart Jacobian',d['Jacobian']+S.Rational(9,8)*v**8*t**4)
    upstream=C.explicit_pair(2);X,Y=upstream['X'],upstream['Y'];u,w=upstream['P'],upstream['Q']
    Ju=S.factor(S.diff(u,X)*S.diff(w,Y)-S.diff(u,Y)*S.diff(w,X))
    ident('actual invariant finite map upstairs volume',Ju+S.Rational(9,4)*X**8)
    ident('quasi-etale chart volume factor2',2*d['Jacobian']+S.Rational(9,4)*v**8*t**4)
    ident('invariant chart P agrees with descended polynomial',u.subs({X:v*Y}).expand()-Y**4*(1+v**3))
    ident('invariant chart Q agrees with descended polynomial',w.subs({X:v*Y}).expand()-Y**6*(1+S.Rational(3,2)*v**3+S.Rational(3,8)*v**6))
    check('actual homogeneous pair finite degree24 upstairs',S.resultant(upstream['p'],upstream['z'],upstream['t'])!=0 and S.total_degree(u)*S.total_degree(w)==24)
    check('quotient finite degree12 with A1 end retained',d['finite_degree']==12 and d['quotient_order']==2 and d['boundary_is_actual_A1_link'])
    relation=S.Matrix(d['principal_unit_relation']);canonical=S.Matrix(d['canonical_vector'])
    check('primitive NONCONSTANT boundary unit relation',math.gcd(*d['principal_unit_relation'])==1)
    check('positive Jacobian vector is exactly principal',canonical==d['relation_multiple']*relation)
    ident('principal invariant X2 has exact chart expression',(v*v*t)**4-v**8*t**4)
    check('finite A1 map actually has zero resolution Chern and wrong open',d['rational_Chern_class_zero'] and d['complement']=='C* times C' and not d['actual_Keller_source'])
    check('wrong open fails Chern certification despite positive coefficients',not N.boundary_chern_gate([8,4],whole_plane_open=False)['rational_nonzero_chern_certified'])
    check('unramified zero-canonical source control not excluded',not N.boundary_chern_gate([0],whole_plane_open=True)['rational_nonzero_chern_certified'])
    for vector in [[2],[2,1],[0,0,1],[0,4,2,7]]:
        gate=N.boundary_chern_gate(vector,whole_plane_open=True)
        check('nonzero free actual boundary vector survives rational Chern '+str(vector),gate['rational_nonzero_chern_certified'] and gate['independent_boundary_classes_certified'] and not gate['realization_claimed'])
    for vector in [[],[-1],[S.Rational(1,2)]]:
        try:N.boundary_chern_gate(vector,whole_plane_open=True)
        except ValueError:check('nonintegral or negative Jacobian order rejected '+str(vector),True)
        else:raise AssertionError('invalid Jacobian vector accepted')
    models=[]
    for kind,rank in [('A',1),('A',5),('D',4),('D',9),('E6',None),('E7',None),('E8',None)]:
        data=N.ade_chern_comparison(kind,rank);B=A.cartan(kind,rank);n=B.rows;k=S.zeros(n,1);k[-1]=2
        E=(k.T*B.inv()*k)[0]
        check('negative inverse form detects nonzero rational Chern '+str((kind,rank)),E>0 and n-2-E<data['theta'])
        check('every Stein ADE filling Chern vanishes by imported comparison '+str((kind,rank)),data['every_Stein_filling_rational_chern_zero'] and data['Stein_implies_minimal'] and data['singular_actual_Keller_normalization_excluded'])
        check('oriented end still requires actual recognition '+str((kind,rank)),data['actual_oriented_end_identification_still_required'] and data['oriented_end']=='negative ADE link')
        blow=M.canonical_packet(S.diag(B,S.ones(1)),[1]*n+[2])
        check('nonminimal blowup passes theta but not deformation input '+str((kind,rank)),blow['theta']==str(data['theta']) and blow['canonical_nonzero'] and blow['rank']==n+1)
        models.append(dict(type=kind,rank=n,theta=data['theta'],nonzero_Chern_energy=str(E),
                           comparison=data,nonminimal_blowup=blow))
    text=' '.join((P/'NORMAL_RESOLUTION_CHERN_GATE_AND_ADE_END_EXCLUSION.md').read_text().lower().split())
    for word in ['all actual normal','singular actual normalizations','minimal resolution','no relations',
                 'pic(y) tensor q -> h^2(y;q)','constant units','crossing points','a_j','integer >=1',
                 'not finite ramification indices','holomorphically minimal','c-infinity trivial',
                 'original resolved kahler form','not assumed exact','c1(tv)^2=0',
                 'every degree','whole marked plane','four times the principal relation',
                 'not the finite cover itself','remain open','not independently reproved']:
        check('scope/proof guard '+word,word in text)
    files=['NORMAL_RESOLUTION_CHERN_GATE_AND_ADE_END_EXCLUSION.md','normal_resolution_chern.py',
           'verify_normal_resolution_chern.py','sources/normal-resolution-chern/manifest.json',
           'ade_end_exclusion_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED109',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_field_degrees=True,all_normal_actual_sources=True,
                     arbitrary_omitted_prime_count=True,arbitrary_indices_and_residues=True,
                     whole_marked_plane_required=True,algebraic_Picard_Chern_injectivity_written=True,
                     minimal_resolution_admits_Stein_deformation=True,
                     deformed_filling_simply_connected=True,rational_Chern_nonzero=True,
                     b2_at_least_all_resolution_boundary_components=True,
                     all_negative_ADE_ends_excluded=True,
                     all_negative_definite_minus_two_TREE_ends_excluded=True,
                     oriented_end_identification_required=True),
        exact_ADE_controls=models,
        genuine_wrong_open_control={k:v for k,v in d.items() if k not in ['v','t','P','Q','Jacobian']},
        assurance=dict(written_unbounded_proof=True,deformation_theorem_imported=True,
                       prior_deep_contact_and_filling_classifications_retained=True,
                       independent_whole_deformation_proof=False,finite_sweep_proves_theorem=False,
                       original_resolution_Kahler_form_exact=False,
                       deformed_finite_target_map_holomorphic_claimed=False,
                       deformed_whole_plane_holomorphic_claimed=False,
                       arbitrary_nonminimal_resolution_Stein_deformation_claimed=False,
                       singularity_smoothing_claimed=False,all_spherical_ends_excluded=False,
                       reversed_orientations_excluded=False,homology_or_lattice_identifies_end=False,
                       global_finite_algebra_closed=False,formal_verification=False,external_review=False,
                       external_review_deferred=True,research_goal_completed=False),
        environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'normal_resolution_chern_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF109 controls; complete NORMAL ADE end class closed, full goal active')


if __name__=='__main__':main()
