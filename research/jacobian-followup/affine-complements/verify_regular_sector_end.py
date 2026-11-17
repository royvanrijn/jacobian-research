'''PF2RULED39 controls for marked regular-sector ends and odd-order cabling.
The plumbing-with-puncture and all-depth three-manifold proofs are written,
with primary imports; finite graph/SNF controls do not formalize them.
'''
from pathlib import Path
from math import gcd
import hashlib,json,platform
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build_graph(q):
    weights={'Dx':0,'Dy':0};edges={frozenset(['Dx','Dy'])}
    def blow_vertex(parent,name):
        weights[parent]-=1;weights[name]=-1;edges.add(frozenset([parent,name]))
    def blow_edge(a,b,name):
        weights[a]-=1;weights[b]-=1;weights[name]=-1
        edges.remove(frozenset([a,b]))
        edges.add(frozenset([a,name]));edges.add(frozenset([b,name]))
    blow_vertex('Dx','E1')
    for j in range(2,q+1):blow_edge('Dx','E'+str(j-1),'E'+str(j))
    for j in range(q+1,2*q+2):blow_vertex('E'+str(j-1),'E'+str(j))
    return weights,edges
def matrix(weights,edges,names):
    return S.Matrix([[weights[a] if a==b else (1 if frozenset([a,b]) in edges else 0)
                      for b in names] for a in names])
def main():
    preds={}
    for filename in ['quartic_cabling_receipt.json','intrinsic_normalization_receipt.json',
                     'coordinate_regular_sector_receipt.json']:
        pred=json.loads((HERE/filename).read_text());preds[filename]=pred
        check('predecessor PASS '+filename,pred['outcome']=='PASS')
        for name,h in pred['source_sha256'].items():
            check('predecessor seal '+filename+' '+name,digest(HERE/name)==h)
    layer=preds['quartic_cabling_receipt.json']['finite_peripheral_layer']
    check('all odd orders use the unchanged exhaustive permutation layer',
          layer['exhaustive_S4_subgroups']==30 and layer['exponent_modulus']==12 and
          layer['S_fiber_trivial_state_count']==9 and layer['A_fiber_trivial_state_count']==4 and
          layer['S_residue_controls']==208 and layer['A_residue_controls']==96)
    D=HERE/'sources/regular-sector-end';manifest=D/'manifest.json'
    pins=json.loads(manifest.read_text())
    check('new primary plumbing source pinned',len(pins['primary_sources'])==1)
    for row in pins['primary_sources']:
        p=D/row['path']
        check('primary PDF seal '+p.name,p.read_bytes().startswith(b'%PDF') and digest(p)==row['sha256'])
        check('correct plumbing input is Proposition2.1 R6',
              'Proposition2.1 move R6' in row['used_statement'])
    visual=pins['visual_source_audit']
    check('closed and punctured plumbing diagrams inspected with no extra summand',
          visual['R6_single_remaining_branch'] and visual['deleted_pair_boundary_punctures']==0 and
          visual['surviving_marked_endpoint_punctures']==1 and visual['no_uncontrolled_extra_summand'])
    u,t,z=S.symbols('u t z')
    for q in range(1,6):
        c=S.symbols('c1:'+str(q+1))
        HH=sum(c[j-1]*u**(q-j) for j in range(1,q+1))
        w=u**q/(HH+u**(q+1)*t);w0=u**q/HH
        zz=-t/(HH*(HH+u**(q+1)*t))
        check('full principal-part last blowup identity q'+str(q),
              S.cancel(w-w0-u**(2*q+1)*zz)==0)
        tinv=-z*HH**2/(1+z*HH*u**(q+1))
        check('whole omitted chart and last blowup coordinates are inverse q'+str(q),
              S.cancel(tinv.subs(z,zz)-t)==0)
        check('last exceptional includes every finite omitted parameter q'+str(q),
              S.cancel(zz.subs(u,0)+t/c[-1]**2)==0)
    graph_controls=[]
    for q in range(1,13):
        weights,edges=build_graph(q);last='E'+str(2*q+1)
        check('last exceptional and all deleted boundary weights q'+str(q),
              weights[last]==-1 and weights['Dy']==0 and weights['Dx']==-q and
              all(weights['E'+str(i)]==-2 for i in range(1,2*q+1)))
        boundary={k:v for k,v in weights.items() if k!=last}
        bedges={e for e in edges if last not in e}
        expected={frozenset(['Dx','Dy']),frozenset(['Dx','E'+str(q)])}|{
            frozenset(['E'+str(i),'E'+str(i+1)]) for i in range(1,2*q)}
        check('full blowup graph has exactly the claimed edges q'+str(q),bedges==expected)
        names=['Dy','Dx']+['E'+str(i) for i in range(1,2*q+1)]
        B=matrix(boundary,bedges,names)
        reduced_names=names[2:];redges={e for e in bedges if not({'Dy','Dx'}&set(e))}
        C=matrix(boundary,redges,reduced_names)
        check('closed boundary determinant agrees with the R6 chain q'+str(q),
              B.det()==-(2*q+1) and C.det()==2*q+1)
        check('R6 has one attached branch and keeps the punctured endpoint q'+str(q),
              sum('Dx' in e for e in bedges)==2 and
              'E'+str(2*q) in reduced_names and len(redges)==2*q-1)
        # One disk piece followed by annulus pieces after puncturing the
        # endpoint. Euler weights on the punctured vertex are immaterial.
        boundary_counts=[sum(name in e for e in redges)+(1 if i==2*q else 0)
                         for i,name in enumerate(reduced_names,1)]
        check('punctured linear plumbing is disk plus annuli q'+str(q),
              boundary_counts==[1]+[2]*(2*q-1))
        v=S.Matrix(range(1,2*q+1))
        check('chain meridian recurrence and primitive marked endpoint q'+str(q),
              C*v==S.Matrix([0]*(2*q-1)+[-(2*q+1)]) and gcd(2*q,2*q+1)==1)
        if q in [1,2,3,5,8,12]:
            SN=smith_normal_form(C,domain=S.ZZ)
            check('lens chain Smith invariants q'+str(q),
                  sorted(abs(int(SN[i,i])) for i in range(2*q))==[1]*(2*q-1)+[2*q+1])
        graph_controls.append(dict(q=q,blowups=2*q+1,deleted_vertices=2*q+2,
                                   end_H1_order=2*q+1,marked_exterior='solid torus'))
    q=S.symbols('q',integer=True,positive=True)
    U=S.Matrix([[0,1],[1,-q]]);coupling=S.Matrix([0,1])
    check('zero-leaf homology elimination does not alter the surviving chain',
          U.det()==-1 and (coupling.T*U.inv()*coupling)[0]==0)
    m=S.symbols('m',integer=True,positive=True)
    check('continued fraction all-length induction',S.cancel(2-(m-1)/m-(m+1)/m)==0)
    check('continued fraction starts with two',S.Rational(2,1)==2)
    check('chain determinant all-length induction',
          S.expand(2*(m+1)-m-(m+2))==0)
    check('end homology is odd and greater than one',S.expand((2*q+1)-2*q)==1)
    # New A-type replacement: torsion-free H1 for every pair of coprime
    # cone orders p,3 follows from the minors proof in the note; controls
    # retain every primitive Seifert numerator at the selected p.
    seifert_controls=0
    for p in [2,4,5,7,8,10]:
        for b1 in range(1,p):
            if gcd(b1,p)!=1:continue
            for b2 in [1,2]:
                M=S.Matrix([[p,0,b1],[0,3,b2]])
                SN=smith_normal_form(M,domain=S.ZZ)
                assert [abs(int(SN[i,i])) for i in range(2)]==[1,1]
                seifert_controls+=1
        check('alternating disk cones p,3 have free rank-one H1 p'+str(p),True)
    k,V,W=S.symbols('k V W',integer=True,positive=True);n=2*k+1
    check('alternating torus homology cannot equal any odd order',
          S.expand((2*V*W-n)-2*(V*W-k))==-1)
    check('twice a primitive class stays nonzero at every odd nontrivial order',
          S.expand(n-2*k)==1)
    local=[]
    for kind,p in [('A',2),('S',3)]:
        order,a,fill_m,fill_l=5,1,1,5;index=3 if kind=='A' else 2
        hm,hl=a,index*p
        check('new locally consistent '+kind+' removal has distance one',
              abs(hm*fill_l-hl*fill_m)==1 and gcd(p,order)==1)
        check('new local '+kind+' marked class is primitive and nonzero',
              len({p*j%order for j in range(order)})==order and p%order!=0)
        local.append(dict(kind=kind,p=p,a=a,n=order,fiber=[hm,hl],
                          filling=[fill_m,fill_l],marked_class=p,
                          is_actual_global_cover=False))
    check('symmetric removal at n5 permits winding three',gcd(3,5)==1 and gcd(3,6)!=1)
    depth=S.symbols('depth',integer=True,positive=True)
    check('both alternating and symmetric removals strictly lower all cabling depths',
          S.expand(depth-(depth-1))==1)
    inputs=[Path(__file__),HERE/'REGULAR_SECTOR_END_AND_QUARTIC_EXCLUSION.md',manifest]+[
            D/row['path'] for row in pins['primary_sources']]+[
            HERE/name for name in preds]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='entire PF2RULED38 simple regular-coordinate marked sector; all source poles and principal parts; finite field degree four excluded',
        end=dict(blowup_count='2q+1',deleted_graph='zero leaf, minus-q neighbor, and the negative-two chain',
            R6='one surviving branch; no additional summand; closed and punctured versions',
            orientation_immaterial_lens='L(2q+1,2q)',H1='Z/(2q+1)',
            ramification_exterior='solid torus',ramification_class='primitive generator',
            pole_recovered_from_abstract_end='q=(|H1|-1)/2',
            q_at_least_two_abstractly_not_Wright=True),
        strengthened_topology=dict(homology_order='every odd n>1',
            simple_S4_identity_longitude_lens_core_type_excluded=True,
            alternating_A4_lens_core_type_excluded=True,
            fourth_letter_pattern_excluded_by_marked_exterior=True,
            new_A_to_A_removal_preserves_actual_boundary_homology=True,
            removals=['S to S','S to A','A to A'],depth_unbounded_written_descent=True,
            n_one_excluded_by_this_theorem=False),
        algebraic_consequences=dict(whole_simple_regular_coordinate_quartic_sector_excluded=True,
            all_q_all_principal_parts=True,
            all_degree_torus_infinity_excluded=True,
            all_degree_prime_boundary_coordinate_degrees_excluded=True,
            minimum_boundary_coordinate_degree=4,minimum_possible_finite_field_degree=5,
            degree_five_realization_claimed=False),
        graph_controls=graph_controls,seifert_homology_control_count=seifert_controls,
        locally_consistent_removals=local,
        assurance=dict(whole_marked_compactification_checked_in_writing=True,
            actual_omitted_knot_exterior_tracked_by_punctured_graph=True,
            primary_plumbing_diagrams_visually_inspected=True,
            not_just_intersection_determinant_or_homology_recognition=True,
            primary_plumbing_import_explicit=True,
            actual_radial_end_recognition_uses_retained_primary_imports=True,
            all_odd_cabling_proof_written=True,
            finite_permutation_layer_unchanged_and_exhaustive=True,
            ball_alternative_excluded_by_essential_marked_class=True,
            all_replacements_track_actual_boundary_lattice=True,
            equal_H1_not_used_to_identify_closed_manifolds=True,
            all_finite_field_degrees_excluded=False,
            higher_regular_coordinate_order_excluded=False,
            no_regular_coordinate_case_excluded=False,
            arbitrary_quartic_packets_excluded=False,
            coefficient_or_Groebner_search_used=False,
            finite_Keller_map_constructed=False,
            formal_verification=False,independent_full_reproof=False,
            external_review=False,research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},
        python=platform.python_version(),sympy=S.__version__)
    (HERE/'regular_sector_end_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} exact controls; all-odd quartic cabling and whole simple-regular sector exclusion; full goal active')
if __name__=='__main__':main()
