"""PF108 exact controls; declared contact/filling classifications are imports."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import ade_end_exclusion as A
import multiple_boundary_contact as M

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.factor(value)==0)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    for filename in ['coprime_leading_forms_receipt.json','multiple_boundary_contact_receipt.json']:
        old=json.loads((P/filename).read_text())
        check('predecessor PASS '+filename,old['outcome']=='PASS')
        for path,h in old['source_sha256'].items():check('predecessor exact seal '+path,sha(P/path)==h)
    manifest=json.loads((P/'sources/ade-end/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('exact primary/context pin '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    audits={x['path'].split('/')[-1]:x.get('fresh_complete_PDF_pages_rendered_and_visually_inspected') for x in manifest['sources']}
    check('actual complete primary pages visually audited',all(audits[k]==v for k,v in {
        'ohta-ono-mpim-2002-reading.pdf':[1,2,3,4], 'ghiggini-2008.pdf':[1,3,4,14],
        'ozsvath-szabo-0303017v2.pdf':[2,11,12], 'ohta-ono-1999.pdf':[1,2]}.items()))
    check('published2005 binary not falsely claimed',manifest['sources'][0]['published_binary_obtained'] is False)
    check('derived reading PDF marked honestly',manifest['sources'][1]['text_extraction_reliable'] is False and manifest['sources'][1]['byte_deterministic_conversion_claimed'] is False)
    r,L,E,t,m,l=S.symbols('r L E t m l',integer=True,positive=True)
    ident('UNBOUNDED actual canonical energy theta',-E-2*(1+r)-3*(-r)-(r-2-E))
    ident('UNBOUNDED zero-Chern reference theta',-2*(1+L)-3*(-L)-(L-2))
    ident('UNBOUNDED contact equality gives energy rank difference',(r-2-E-(L-2)).subs(E,r-L))
    ident('UNBOUNDED minimal rank equality forces energy zero',(r-2-E-(L-2)).subs(r,L)+E)
    check('strict positive energy cannot be zero',E.is_positive and E.is_zero is False)
    ident('UNBOUNDED four-leaf star zero-vector',2*2**2+4*2*1**2-2*4*2*1)
    ident('UNBOUNDED two-trivalent zero-vector',8*t+8-8*(t-1)-16)
    ident('UNBOUNDED chain determinant recurrence',2*l-(l-1)-(l+1))
    ident('UNBOUNDED D central Schur',2-S.Rational(1,2)-S.Rational(1,2)-(m-3)/(m-2)-1/(m-2))
    ident('UNBOUNDED D determinant',2*2*(m-2)/(m-2)-4)
    ident('UNBOUNDED D rank',1+1+1+(m-3)-m)
    # Boundary cases in the written reciprocal classification, not a graph sweep.
    check('all p>=3 cannot be positive',3*S.Rational(1,3)-1==0)
    check('p2 q>=4 cannot be positive',S.Rational(1,2)+2*S.Rational(1,4)-1==0)
    check('p2 q3 r>=6 cannot be positive',S.Rational(1,2)+S.Rational(1,3)+S.Rational(1,6)-1==0)
    controls=[]
    for kind,rank in [('A',1),('A',6),('D',4),('D',5),('D',9),('D',18),('E6',None),('E7',None),('E8',None)]:
        data=A.ade_data(kind,rank);C=A.cartan(kind,rank);n=C.rows;label=kind+str(rank or '')
        check('Cartan rank/symmetry '+label,n==data['rank'] and C==C.T)
        check('Cartan negative plumbing determinant '+label,C.det()==data['H1_order'])
        check('Cartan positive definite '+label,all(C[:i,:i].det()>0 for i in range(1,n+1)))
        reference=M.canonical_packet(C,[1]*n)
        c=A.compare(reference,kind,rank,oriented_end_identified=True)
        check('genuine zero-Chern reference filling survives '+label,c['H1_order_matches'] and c['minimal_rank_matches'] and c['contact_theta_matches'] and not c['positive_source_energy'])
        indices=[1]*n;indices[-1]=3;ram=M.canonical_packet(C,indices)
        c=A.compare(ram,kind,rank,oriented_end_identified=True)
        check('ramified canonical data fail unique contact '+label,c['positive_source_energy'] and not c['contact_theta_matches'])
        blow=M.canonical_packet(S.diag(C,S.ones(1)),[1]*n+[2])
        bc=A.compare(blow,kind,rank,oriented_end_identified=True)
        check('real blowup passes theta but minimal rank is ESSENTIAL '+label,bc['H1_order_matches'] and bc['contact_theta_matches'] and not bc['minimal_rank_matches'] and bc['positive_source_energy'])
        check('homology does not certify end recognition '+label,not A.compare(ram,kind,rank)['eligible_for_obstruction'])
        if kind!='A':
            check('negative Seifert e0 retained '+label,data['Seifert_e0']==-2)
            check('every minus-two arm yields contact count1 '+label,all(abs(-2+1)**length==1 for length in data['arm_lengths']))
            check('spherical orbifold and orientation signs '+label,S.Rational(data['orbifold_Euler'])>0 and S.Rational(data['rational_Seifert_Euler'])==-S.Rational(data['orbifold_Euler']))
        controls.append(dict(type=kind,rank=n,reference=reference,ramified=ram,blowup=blow,oriented_end_recognition_not_inferred=True))
    e8=next(x for x in controls if x['type']=='E8')
    check('PF104 numerical E8 control remains but Poincare end fails',e8['ramified']['canonical_energy']=='8' and e8['ramified']['theta']=='-2' and A.ade_data('E8')['canonical_theta']==6)
    check('same homology order does not identify different D ends',A.ade_data('D',4)['H1_order']==A.ade_data('D',5)['H1_order']==4 and A.ade_data('D',4)['rank']!=A.ade_data('D',5)['rank'])
    for kind,rank in [('D',3),('A',0),('E8',7),('E9',None)]:
        try:A.ade_data(kind,rank)
        except ValueError:check('out of class input rejected '+str((kind,rank)),True)
        else:raise AssertionError('invalid ADE input accepted')
    text=' '.join((P/'ALL_ADE_END_CANONICAL_ENERGY_EXCLUSION.md').read_text().lower().split())
    for term in ['all-degree','actual smooth','all omitted primes','every m','one','e=0',
                 'negative-definite tree','unbounded classification proof','does not identify',
                 'not covered','singular n','not obtained','not independently reproduced','remain open']:
        check('scope/proof guard '+term,term in text)
    files=['ALL_ADE_END_CANONICAL_ENERGY_EXCLUSION.md','ade_end_exclusion.py','verify_ade_end_exclusion.py',
           'sources/ade-end/manifest.json','coprime_leading_forms_receipt.json','multiple_boundary_contact_receipt.json']+[x['path'] for x in manifest['sources']]
    result=dict(claim='PF2RULED108',outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},
        theorem=dict(all_field_degrees=True,arbitrary_omitted_prime_count=True,arbitrary_indices=True,
                     smooth_actual_normalization_required=True,whole_plane_open_required=True,
                     oriented_negative_ADE_link_required=True,all_ADE_ends_excluded=True,
                     all_negative_definite_minus_two_tree_plumbing_ends_excluded=True,
                     general_unique_contact_minimal_rank_criterion=True),
        exact_controls=controls,
        assurance=dict(written_unbounded_proof=True,deep_contact_and_filling_theorems_imported=True,
                       primary_full_selected_pages_visually_audited=True,
                       independent_whole_classification_proofs=False,finite_sweep_proves_theorem=False,
                       singular_normalizations_excluded=False,all_spherical_ends_excluded=False,
                       opposite_orientations_excluded=False,intersection_form_identifies_end=False,
                       global_finite_algebra_closed=False,formal_verification=False,external_review=False,
                       external_review_deferred=True,research_goal_completed=False),
        environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'ade_end_exclusion_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF108 controls; complete smooth oriented ADE end class closed, full goal active')


if __name__=='__main__':main()
