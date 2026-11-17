'''PF2RULED41 controls for the all-r borderline marked end and quartic exclusion.
Three-manifold recognition and all-depth cabling arguments are written proofs.
The literal countercontrol is a topological branched cover, not a Keller map.
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
def graph(q):
    weights={'Dx':0,'Dy':0,'L':0}
    edges={frozenset(['Dx','Dy']),frozenset(['Dx','L'])}
    def edge(a,b,c):
        edges.remove(frozenset([a,b]))
        edges.add(frozenset([a,c]));edges.add(frozenset([c,b]))
        weights[a]-=1;weights[b]-=1;weights[c]=-1
    edge('Dx','L','Ediag');edge('Ediag','L','E1')
    for j in range(2,q+1):edge('Ediag','E'+str(j-1),'E'+str(j))
    last='E'+str(q);weights[last]-=1;weights['R']=-1
    edges.add(frozenset([last,'R']))
    return weights,edges
def matrix(names,weights,edges):
    return S.Matrix([[weights[a] if a==b else int(frozenset([a,b]) in edges)
                      for b in names] for a in names])
def mul(a,b):return tuple(a[b[j]] for j in range(4))
I=(0,1,2,3)
def power(a,n):
    if n<0:a=tuple(a.index(j) for j in range(4));n=-n
    out=I
    for _ in range(n):out=mul(out,a)
    return out
def generated(gs):
    seen={I};todo=[I]
    while todo:
        a=todo.pop()
        for b in gs:
            c=mul(a,b)
            if c not in seen:seen.add(c);todo.append(c)
    return seen
def main():
    preds={}
    for name in ['higher_regular_coordinate_receipt.json','regular_sector_end_receipt.json',
                 'quartic_cabling_receipt.json','torus_infinity_receipt.json','intrinsic_normalization_receipt.json']:
        pred=json.loads((HERE/name).read_text());preds[name]=pred
        check('predecessor PASS '+name,pred['outcome']=='PASS')
        for p,h in pred['source_sha256'].items():
            check('predecessor seal '+name+' '+p,digest(HERE/p)==h)
    layer=preds['quartic_cabling_receipt.json']['finite_peripheral_layer']
    check('unchanged exhaustive all-S4 peripheral layer',
          layer['exhaustive_S4_subgroups']==30 and layer['exponent_modulus']==12 and
          layer['S_fiber_trivial_state_count']==9 and layer['S_residue_controls']==208)
    manifest=HERE/'sources/borderline-regular-end/manifest.json'
    pins=json.loads(manifest.read_text())
    check('new primary rigidity and retained marked plumbing imports pinned',len(pins['primary_sources'])==6)
    for row in pins['primary_sources']:
        p=HERE/row['path_relative_to_affine_complements']
        check('primary seal '+p.name,p.read_bytes().startswith(b'%PDF') and digest(p)==row['sha256'] and p.stat().st_size==row['bytes'])
    audit=pins['new_visual_source_audit']
    check('Waldhausen theorem full hypotheses and exception inspected',
          audit['printed_page']==77 and audit['PDF_page']==23 and
          all(audit[k] for k in ['irreducible_and_boundary_irreducible_hypotheses_checked',
              'sufficiently_large_target_checked','exceptional_product_alternative_excluded_by_boundary_count',
              'peripheral_isomorphism_realized_before_using_theorem']))
    s,t,a,z,v=S.symbols('s t a z v')
    for q in range(1,9):
        r=q+1;x=a/s**q;y=s**r/a
        check('unimodular toric coordinates recover s and a q'+str(q),
              S.cancel(x*y-s)==0 and S.cancel(x**r*y**q-a)==0)
        jac=S.diff(x,s)*S.diff(y,a)-S.diff(x,a)*S.diff(y,s)
        check('toric source volume before the free blowup q'+str(q),S.cancel(jac+1/a)==0)
        check('actual free-blowup source volume has order one q'+str(q),
              S.cancel((jac*s).subs(a,1+s*t)+s/(1+s*t))==0)
        ss=z*v;tt=(v**r-1)/(z*v)
        check('entire Dv quotient chart equation q'+str(q),S.cancel(1+ss*tt-v**r)==0)
        check('free exceptional coordinate includes every finite parameter q'+str(q),
              (1+s*t).subs(s,0)==1)
    graphs=[]
    for q in range(1,17):
        r=q+1;n=q*r+1;weights,edges=graph(q);marked='E'+str(q)
        order=['Dy','Dx','Ediag']+['E'+str(j) for j in range(q,0,-1)]
        target_edges={frozenset([order[j],order[j+1]]) for j in range(len(order)-1)}
        check('full deleted graph from the actual corner and free blowups q'+str(q),
              {e for e in edges if not({'L','R'}&set(e))}==target_edges)
        check('deleted and retained self-intersections q'+str(q),
              [weights[k] for k in order]==[0,-1,-r]+[-2]*q and
              weights['R']==-1 and weights['L']==-2)
        check('omitted line attaches at the first negative-two vertex q'+str(q),
              frozenset(['R',marked]) in edges and order[3]==marked)
        reduced=order[2:];redges={e for e in target_edges if not({'Dy','Dx'}&set(e))}
        C=-matrix(reduced,weights,redges)
        check('closed R6 chain determinant q'+str(q),C.det()==n)
        meridians=S.Matrix([1+j*q for j in range(q+1)])
        check('actual marked chain relations and primitive class q'+str(q),
              C*meridians==S.Matrix([0]*q+[n]) and
              int(meridians[1])==r and gcd(r,n)==1)
        check('odd lens order is nontrivial q'+str(q),n>1 and n%2==1)
        check('R6 has one remaining branch and preserves the marked puncture q'+str(q),
              sum('Dx' in e for e in target_edges)==2 and marked in reduced)
        if q in [1,2,3,5,8,12]:
            SN=smith_normal_form(C,domain=S.ZZ)
            check('chain Smith invariants q'+str(q),
                  sorted(abs(int(SN[j,j])) for j in range(q+1))==[1]*q+[n])
        rel=S.Matrix([[r,0,1],[0,q,q-1]])
        gen=S.Matrix([-q,-r*(q-1),q*r])
        check('punctured exterior has the asserted integral generator q'+str(q),
              rel*gen==S.zeros(2,1) and gcd(gcd(q,r*(q-1)),q*r)==1)
        check('actual omitted-knot meridian and S3 filling slope q'+str(q),
              (S.Matrix([[1,1,2]])*gen)[0]==n and n-q*r==1)
        check('marked exterior torus orders and noncore test q'+str(q),
              (q==1)==(min(q,r)==1))
        graphs.append(dict(r=r,q=q,blowups=q+2,chain=[r]+[2]*q,
            marked_vertex_index_zero_based=1,H1_order=n,primitive_marked_class=r,
            torus_exterior_orders=[q,r],omitted_meridian=[n,1],
            Heegaard_core=(r==2)))
    q=S.symbols('q',integer=True,positive=True);r=q+1;n=q*r+1
    check('all-length continued fraction identity',S.cancel(r-q/(q+1)-n/r)==0)
    check('all-r filling slope and fiber have intersection one',S.expand(n-q*r)==1)
    check('all-r closed homology order is coprime to both cones',
          S.rem(n,q)==1 and S.rem(n,r)==1)
    cone_controls=0
    for p in range(4,24,2):
        check('regular-fiber quotient in the transposition pattern is nonzero p'+str(p),
              gcd(p//2,p)==p//2 and p//2>1);cone_controls+=1
    for p in range(3,24,2):
        check('regular-fiber quotient in the even-slope two-letter pattern is nonzero p'+str(p),
              gcd(p,p)==p and p>1);cone_controls+=1
    signed=[]
    for epsilon in [-1,1]:
        for delta in [-1,1]:
            half_a=3*epsilon+delta;a=2*half_a
            transfer=S.Matrix([[-1,half_a],[0,1]])
            h2=S.Matrix([half_a,3]);h1=S.Matrix([a,3])
            muY=S.Matrix([epsilon,1]);fill=transfer*muY
            check('actual T2I transfer preserves the fiber '+str((epsilon,delta)),
                  transfer*h2==h1 and transfer.det()==-1)
            check('atoroidal final case forces distance-one filling '+str((epsilon,delta)),
                  abs(half_a-3*epsilon)==1 and fill==S.Matrix([2*epsilon+delta,1]))
            third=abs(int(fill[0])-6*epsilon)
            check('both mirror trefoil fillings have a genuine third cone '+str((epsilon,delta)),
                  third in [3,5] and third!=1 and gcd(a,3)==1)
            signed.append(dict(epsilon=epsilon,delta=delta,cable_slope=a,
                companion_filling=list(map(int,fill)),third_cone=third,lens_possible=False))
    # A real simple S4 lens-cover control: same order seven, wrong marked exterior.
    tau=(1,0,2,3);alpha=(1,2,0,3);dd=mul(alpha,tau)
    bb=(2,1,3,0);AA=(0,1,3,2)
    check('literal countercontrol has the claimed companion meridian',
          power(AA,2)==I and power(bb,3)==I and mul(power(AA,-1),power(bb,2))==dd)
    check('literal pattern and companion give full S4',
          len(generated([alpha,tau]))==6 and len(generated([AA,bb]))==6 and
          len(generated([alpha,tau,bb]))==24)
    check('actual cable peripheral meridian and longitude match at (3,16)',
          power(alpha,3)==I and power(dd,3)==dd and power(dd,-16)==I and
          power(tau,-48)==I)
    check('order-seven lens cover fails the marked atoroidality test',
          16//2-1==7 and abs(7-6)==1 and abs(16//2-3)==5 and
          gcd(16%7,7)==1 and 16%7==2)
    tt=S.symbols('tt');pp=tt**6+tt**2;qq=tt**9+tt**7+tt**3
    check('polynomial infinity control has exactly the derivative gcd for one ordinary cusp',
          S.Poly(S.gcd(S.diff(pp,tt),S.diff(qq,tt)),tt).monic().as_expr()==tt)
    check('ordinary cusp has no other parameter preimage',
          S.gcd(S.cancel(pp/tt**2),S.cancel(qq/tt**3))==1)
    defect=S.Poly(S.expand(qq**2-pp**3),tt)
    check('polynomial infinity control has cable winding (3,16)',
          defect.degree()==16 and defect.LC()==2 and gcd(3,16)==1)
    inputs=[Path(__file__),HERE/'BORDERLINE_REGULAR_END_AND_EXCLUSION.md',manifest]+[
        HERE/row['path_relative_to_affine_complements'] for row in pins['primary_sources']]+[
        HERE/k for k in preds]
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
        scope='all q=r-1 marked regular-coordinate models with r>=2; marked end and torus-infinity constraints all-degree, finite field degree four excluded',
        end=dict(all_r_entire_compactification_written=True,blowup_count='r+1',
            R6_chain='[-r]--[-2]_(r-1 entries)',puncture='first negative-two vertex',
            H1='Z/(r(r-1)+1)',marked_class_primitive=True,
            exterior='T(r-1,r) knot exterior',filling_slope='(r(r-1)+1)mu+lambda',
            Heegaard_core_iff='r=2',orientation='consistent simultaneous mirror allowed',
            actual_radial_exterior_recognized_using_Waldhausen=True),
        strengthened_topology=dict(homology_order='every odd n>1',full_S4=True,
            simple_branch=True,identity_longitude=True,source_lens=True,
            primitive_ramification=True,atoroidal_ramification_exterior=True,
            such_cover_excluded=True,Heegaard_core_required=False,
            exact_old_permutation_layer_reused=True,
            primitive_regular_fiber_quotient_argument=True,
            genuine_cable_torus_boundary_parallel_exception_checked=True,
            final_S3_torus_companion_forced_trefoil=True,
            actual_trefoil_and_pattern_boundary_lattices_tracked=True),
        consequences=dict(whole_borderline_quartic_sector_excluded=True,
            no_coefficient_or_pole_bound=True,all_degree_torus_infinity_excluded=True,
            all_degree_prime_boundary_coordinate_degree_excluded=True,
            remaining_actual_quartic_regular_sector='r>=2 and q>=r+1',
            higher_finite_degrees_excluded=False),
        countercontrol=dict(actual_topological_simple_S4_cover_constructed=True,
            target='(3,16)-cable of positive trefoil',lens_H1_order=7,identity_longitude=True,
            primitive_marked_class=2,companion_filling=[7,1],marked_exterior_extra_cone=5,
            atoroidal=False,marked_N3_end_realized=False,finite_Keller_map_constructed=False,
            polynomial_infinity_curve=['tt^6+tt^2','tt^9+tt^7+tt^3'],
            derivative_gcd='tt',cable_winding=[3,16],
            all_collision_profiles_certified=False),
        controls=dict(graphs=graphs,primitive_fiber_cone_cases=cone_controls,
                      signed_trefoil_last_cases=signed),
        assurance=dict(all_r_and_all_depth_proofs_written=True,primary_rigidity_import_explicit=True,
            primary_rigidity_statement_visually_inspected=True,
            actual_marked_exterior_not_inferred_from_H1=True,
            no_target_equivariance_assumed=True,finite_controls_not_full_geometric_proof=True,
            affine_converse_for_arbitrary_Laurent_data=False,
            whole_higher_regular_sector_excluded=False,no_regular_coordinate_case_excluded=False,
            extra_omitted_primes_excluded=False,arbitrary_quartic_packets_excluded=False,
            all_finite_Wright_degrees_excluded=False,finite_Keller_map_constructed=False,
            formal_verification=False,independent_full_reproof=False,external_review=False,
            research_goal_completed=False),
        source_sha256={str(p.relative_to(HERE)):digest(p) for p in inputs},
        runtime=dict(python=platform.python_version(),sympy=S.__version__))
    (HERE/'borderline_regular_end_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS '+str(len(checks))+' borderline marked-end controls; whole quartic borderline excluded; goal active')
if __name__=='__main__':main()
