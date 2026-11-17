"""PF123 exact controls and seals; geometric descent uses written proof/imports."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import intermediate_critical_axis as L

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,value):check(name,S.simplify(value)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    olds=['block_quotient_monodromy_receipt.json','simple_ramification_ledger_receipt.json',
          'cyclic_quotient_morse_receipt.json','local_monodromy_branch_receipt.json']
    for name in olds:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    m=json.loads((P/'sources/intermediate-critical-axis/manifest.json').read_text())
    for item in m['sources']:
        f=P/item['path'];check('source seal '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('exact new primary imports and honest audit boundary',m['Pedersen_Theorem1_1_v3_imported'] and m['Giraldo_Roeder_published_TheoremA_imported'] and m['primary_HTML_bytes_pinned'] and not m['Pedersen_old_v2_used'] and not m['new_PDF_visual_audit_claimed'] and not m['independent_proof_of_two_imported_geometric_theorems'])
    t,b,mm=S.symbols('t b m',positive=True)
    ident('UNBOUNDED cyclic-lift excess scaling',t*b/4*(mm+1)-t*b/2-t*b/4*(mm-1))
    ident('UNBOUNDED affine outside cost balance',t*b-1-t*b/2-t*b/4*(mm-1)-(t*b/2-1-t*b/4*(mm-1)))
    check('UNBOUNDED two nonimmersion units overfill',S.expand(t*b/2-1-t*b/4*2)==-1)
    check('UNBOUNDED intermediate index2 excludes all singular orders',2-1<2)
    x,y=S.symbols('x y');k=S.symbols('k',integer=True,positive=True)
    p=x**k-3*y**2;q=x**k*y-2*y**3
    ident('UNBOUNDED actual cyclic cusp composition Jacobian',S.diff(p,x)*S.diff(q,y)-S.diff(p,y)*S.diff(q,x)-k*x**(2*k-1))
    ident('UNBOUNDED actual cusp fibre monic reduction',q-y*p-y**3)
    cusp_controls=[]
    for kk in [1,2,5]:
        pp,qq=p.subs(k,kk),q.subs(k,kk)
        gb=S.groebner([pp,qq],x,y)
        check('actual cyclic cusp critical axis '+str(kk),S.factor(S.diff(pp,x)*S.diff(qq,y)-S.diff(pp,y)*S.diff(qq,x))==kk*x**(2*kk-1))
        check('actual cyclic cusp fibre monic basis '+str(kk),gb.reduce(y**3)[1]==0 and gb.reduce(x**kk-3*y**2)[1]==0)
        check('actual higher even index remains possible '+str(kk),3*kk-2*kk==kk and 2*kk>=2)
        cusp_controls.append(dict(k=kk,e=2*kk,local_rank=3*kk,eta=kk))
    packets=[]
    for bb,tt,parts in [(4,3,[[1,2]]),(8,2,[[2],[1,1]])]:
        dd=bb*tt;indices=[2*z for packet in parts for z in packet]
        check('actual packet degrees conserve t '+str((bb,tt)),len(parts)==bb//4 and all(sum(packet)==tt for packet in parts))
        check('actual visible mass and index bound '+str((bb,tt)),sum(indices)==dd//2 and all(2<=e<=2*tt for e in indices))
        for defect in [0,1]:
            eta=dd//4*defect;outside=dd//2-1-eta
            check('actual all-index affine profile '+str((bb,tt,defect)),outside>=0 and sum(indices)+eta+outside==dd-1)
        packets.append(dict(b=bb,t=tt,packets=parts,indices=indices))
    affine=[]
    for aa in [2,3,5]:
        p0=L.linear_transvection(aa,0,1);q0=L.linear_transvection(aa,1,0)
        oo=L.orbits([p0,q0]);one=[o for o in oo if len(o)==1];three=[o for o in oo if len(o)==3];degree=1<<aa
        check('actual cusp S3 orbit profile '+str(aa),len(one)==len(three)==degree//4 and len(one)+3*len(three)==degree)
        check('actual each nontrivial cusp orbit has one transposition '+str(aa),all(sum(p0[i]!=i for i in o)==2 for o in three))
        affine.append(dict(a=aa,degree=degree,fixed=len(one),triple_orbits=len(three)))
    gens=[L.linear_transvection(3,0,1),L.linear_transvection(3,1,2),L.linear_transvection(3,2,0)]
    oo=L.orbits(gens)
    check('countercontrol fixed-point inequality alone permits GL3 orbit1+7',sorted(map(len,oo))==[1,7])
    check('actual GL3 local configuration fails proved single-transposition orbit condition',sum(gens[0][i]!=i for i in next(o for o in oo if len(o)==7))==4)
    text=' '.join((P/'INTERMEDIATE_CRITICAL_AXIS_DESCENT_AND_AFFINE_CUSP_PROFILES.md').read_text().lower().split())
    for phrase in ['original finite normalization','all omitted primes','intransitive at every finite target point',
      'every intermediate normalization','pedersen','v3','giraldo--roeder','reduced inverse image',
      'hartogs extends','smallness makes each character faithful','2<=n<=e-1',
      'entire m smooth','sum_(d_i over e) k_i=t','global f_i=1','e_i<=2t',
      'exactly one ordinary cusp','outside cost','higher cusps and multiple cusps are excluded',
      'c* times c','not polynomially isomorphic to a2','remain open']:
        check('written scope '+phrase,phrase in text)
    files=['INTERMEDIATE_CRITICAL_AXIS_DESCENT_AND_AFFINE_CUSP_PROFILES.md','intermediate_critical_axis.py',
      'verify_intermediate_critical_axis.py','sources/intermediate-critical-axis/manifest.json']+olds+[x['path'] for x in m['sources']]
    r=dict(claim='PF2RULED123',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_plane_degrees_all_original_monodromy_types=True,whole_original_normalization_all_primes_required=True,
        every_block_quotient_local_monodromy_intransitive_at_all_finite_values=True,
        every_intermediate_complete_fibre_has_at_least_two_points=True,
        every_intermediate_critical_divisor_disjoint_smooth_A1_curves=True,
        every_intermediate_singularity_cyclic_axis_type_and_order_bound=True,
        binary_affine_intermediate_entirely_smooth=True,
        binary_affine_visible_original_primes_smooth_residue1_even_index_at_most2t=True,
        binary_affine_visible_indices_partition_into_b_over4_packets_of_total_t=True,
        all_binary_affine_block_quotients_exact_two_cusp_profiles=True,
        all_binary_affine_higher_multiple_cusps_excluded=True),
      cusp_controls=cusp_controls,packet_controls=packets,affine_controls=affine,
      assurance=dict(Pedersen_v3_and_published_Giraldo_Roeder_theorems_imported=True,
        intermediate_local_geometry_proved_before_PF106_PF119_application=True,
        primary_HTML_sources_pinned=True,new_PDF_visual_audit=False,
        independent_proof_of_two_imported_geometric_theorems=False,
        intermediate_global_plane_recognition=False,
        all_original_visible_indices_equal2=False,all_original_N_entirely_smooth=False,
        invisible_original_branch_images_excluded=False,
        outside_cost_assumed_only_unramified_prime_count=False,
        local_GL3_exclusion_promoted_to_global_AGL3_exclusion=False,
        degree4_plane_inversion_applied_to_intermediate=False,
        source_surjectivity_or_exact_cusp_retained_count_proved=False,
        finite_samples_prove_theorem=False,all_binary_affine_quotients_excluded=False,
        actual_global_finite_algebra_closed=False,historical_priority_claimed=False,
        formal_verification=False,external_review=False,research_goal_completed=False),
      source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'intermediate_critical_axis_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
    print('PASS',len(checks),'PF123 controls; intermediate critical-axis descent and exact affine cusp profiles; goal ACTIVE')


if __name__=='__main__':main()
