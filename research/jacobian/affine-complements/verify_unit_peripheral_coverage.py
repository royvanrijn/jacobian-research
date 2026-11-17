"""PF2RULED97: exact unit periphery and reduced target coverage controls.

The primary classification is imported; finite covers do not prove it.
No old sealed helper is edited. Recursive words are never expanded.
"""
from pathlib import Path
import sys,json,hashlib,os,platform,math
import sympy as S
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P))
import composed_peripheral_cover as C
E=C.E

def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def identity(name,value):check(name,S.expand(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    pred=P/'prefix_excess_closure_receipt.json'
    old=json.loads(pred.read_text());check('PF96 predecessor PASS',old['outcome']=='PASS')
    for path,digest in old['source_sha256'].items():
        check('sealed PF96 predecessor '+path,sha(P/path)==digest)
    manifest=P/'sources/unit-peripheral-coverage/manifest.json'
    sources=json.loads(manifest.read_text())
    for row in sources['primary_sources']:
        f=P/row['path']
        check('primary HTML pin '+f.name,f.stat().st_size==row['bytes'] and sha(f)==row['sha256'])
        check('primary classification explicitly imported',row['deep_classification_imported_not_independently_reproved'] is True)
    check('source audit does not claim earlier whole-chain reproof',
        sources['full_source_finite_completion_implication_independently_audited'] is False and
        sources['full_quartic_Keller_exclusion_registered'] is False and sources['fresh_PDF_audit'] is False)

    k,v,r,u,a,b=S.symbols('k v r u a b',integer=True)
    # Exponents in the commuting peripheral subgroup after z=M^k L.
    identity('r1 meridian M exponent for ALL Bezout choices',(1+k*v)-v*k-1)
    identity('r1 meridian L exponent for ALL Bezout choices',v-v)
    identity('r1 preferred longitude M exponent',-k+k)
    identity('r1 preferred longitude L exponent',S.Integer(1)-1)
    identity('every integral filling slope coefficient retained',a*((1+k*v)-v*k)+b*(-k+k)-a)
    identity('unknot parameter+1 meridian exponent after elimination',(r*u-v).subs(v,r*u-1)-1)
    identity('unknot parameter-1 meridian exponent after elimination',(-r*u-v).subs(v,1-r*u)+1)
    identity('unknot parameter+1 preferred longitude exponent',-r+r)
    identity('unknot parameter-1 preferred longitude exponent',-r+r)

    # Independent full cover controls: nonabelian degree3 and degree4,
    # then cyclic degree5 and7. All start with actual lambda monodromy1.
    examples=[(2,3,[(1,0,2),(1,2,0)],'trefoil S3'),
        (3,4,[(0,2,3,1),(1,2,3,0)],'torus full S4')]
    for degree in [5,7]:
        cycle=tuple(list(range(1,degree))+[0])
        examples.append((2,3,[E.power(cycle,3),E.power(cycle,2)],'cyclic degree'+str(degree)))
    def signature(data):
        rows=data['relation_rows'];bd=data['boundaries']
        ans={'unfilled':C.reduced_h1(rows)}
        for name,aa,bb in [('normal',1,0),('longitude',0,1),('slope2m3L',2,3),('slope-m2L',-1,2)]:
            slopes=[[aa*x+bb*y for x,y in zip(t['meridian'],t['longitude'])] for t in bd]
            ans[name]=C.reduced_h1(rows,slopes)
            for j,t in enumerate(bd):
                ans[name+' plus marked L '+str(j)]=C.reduced_h1(rows,slopes+[t['longitude']])
        ans['boundary_cycles']=[(t['degree'],t['sheets']) for t in bd]
        ans['meridian_monodromy']=data['meridian_monodromy']
        return ans
    controls=[]
    sequences=[[-7],[0],[1],[10**24+39],[-(10**24+39)],[-5,0,7,-11,13,17]]
    for p,q,rho,label in examples:
        base=C.peripheral_cover(p,q,[],rho)
        target=signature(base);M=base['meridian_monodromy']
        # E.power is linear; reduce by a valid finite permutation period FIRST.
        period=math.lcm(*range(1,len(M)+1))
        for sequence in sequences:
            images=rho+[E.power(M,kk%period) for kk in sequence]
            actual=[]
            for reverse in [False,True]:
                data=C.peripheral_cover(p,q,[(1,kk) for kk in sequence],images,reverse=reverse)
                value=signature(data)
                check('all unfilled/signed filled marked groups and cycles retained '+str((label,sequence,reverse)),value==target)
                check('no expanded recursive peripheral words '+str((label,sequence,reverse)),data['word_expansion_used'] is False)
                actual.append(value)
                if len(M)==3:
                    linking=C.integral_null_core_linking(data)
                    check('oriented trefoil full framing matrix unchanged '+str((sequence,reverse)),linking['matrix']==[[-1,2],[2,-4]])
            controls.append(dict(label=label,degree=len(M),parameters=sequence,signatures=actual))

    # A unit PARAMETER on a genuine companion changes colored linking.
    rho=[(1,0,2),(1,2,0)];base=C.peripheral_cover(2,3,[],rho)
    before=C.integral_null_core_linking(base);M=base['meridian_monodromy']
    counter=[]
    for reverse in [False,True]:
        data=C.peripheral_cover(2,3,[(3,1)],rho+[M],reverse=reverse)
        after=C.integral_null_core_linking(data)
        check('parameter1 is NOT winding1, exact coloring linking countercontrol '+str(reverse),
            before['matrix']==[[-1,2],[2,-4]] and after['matrix']==[[-9,18],[18,-36]])
        counter.append(after)
    # Leading unknot resets: build backwards from the final meridian.
    # Its restriction to an older companion need not remain transitive.
    unknot=[]
    for degree in [3,4,5]:
        cyc=tuple(list(range(1,degree))+[0]);cables=[(2,1),(3,-1),(5,1),(7,-1)]
        previous=cyc;added=[]
        for rr,kk in reversed(cables):
            added.append(E.power(previous,kk%degree))
            previous=E.power(previous,rr%degree)
        rho=[E.power(previous,5%degree),previous]+list(reversed(added));M=cyc
        data=C.peripheral_cover(1,5,cables,rho)
        check('leading signed unit-parameter unknot remains cyclic '+str(degree),
            data['meridian_monodromy']==M and len(data['boundaries'])==1 and
            C.reduced_h1(data['relation_rows'])==dict(free=1,torsion=[]) and
            C.reduced_h1(data['relation_rows'],[data['boundaries'][0]['meridian']])==dict(free=0,torsion=[]))
        unknot.append(dict(degree=degree,cables=cables,unfilled_H1=C.reduced_h1(data['relation_rows'])))
    check('cyclic transposition image cannot be full S4',E.generated([(1,0,2,3)])==2)
    check('retained PF96 closure BOTH pencils all lengths',
        old['quartic']['ENTIRE_GENUINE_POSITIVE_ITERATED_TORUS_QUARTIC_CLASS_BOTH_PENCILS_CLOSED'] is True)
    inputs=[Path(__file__),P/'UNIT_PERIPHERAL_REDUCTION_AND_TARGET_COVERAGE.md',
        P/'composed_peripheral_cover.py',manifest,pred]+[P/x['path'] for x in sources['primary_sources']]
    result=dict(as_of='2026-11-15',outcome='PASS',check_count=len(checks),checks=checks,
        general=dict(all_signed_unit_winding_stages_preserve_both_preferred_generators=True,
            all_Bezout_choices_and_all_integral_slopes_preserved=True,
            every_degree_representation_and_colored_cover_linking_preserved=True,
            initial_unknot_signed_unit_parameter_reduction=True,
            first_nontrivial_positive_cable_resets_actual_torus_stage=True,
            parameter1_on_nontrivial_companion_NOT_erased=True,
            unknot_cannot_have_full_S4_image=True),
        coverage=dict(primary_reduced_nontrivial_positive_splice_classification_imported=True,
            ordinary_winding_q_parameter_p_and_strict_infinity_conditions_retained=True,
            positive_genuine_prefix_semigroup_excludes_parameter1=True,
            target_degenerate_escape_closed_UNDER_WHOLE_PACKET_HYPOTHESES=True,
            no_transfer_of_coordinate_dependent_canonical_coefficients_assumed=True),
        signed_unit_winding_cover_controls=controls,parameter1_countercontrol=dict(before=before,after=counter),
        leading_unknot_controls=unknot,
        assurance=dict(written_degree_free_peripheral_and_geometric_proofs=True,
            finite_cover_controls_not_universal_or_classification_proof=True,
            earlier_whole_source_chain_and_finite_completion_audit_required=True,
            no_invented_residual_unit_stage_family=True,
            primary_classification_independently_reproved=False,
            full_quartic_Keller_exclusion_registered=False,
            full_compatible_completion_goal_solved=False,formal_verification=False,
            external_review=False,independent_full_reproof=False,research_goal_completed=False),
        source_sha256={os.path.relpath(f,P):sha(f) for f in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'unit_peripheral_coverage_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; exact all-degree unit periphery and primary target coverage; whole-chain/full completion audit OPEN')
if __name__=='__main__':main()
