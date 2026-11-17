"""PF119 small exact controls; unbounded results use the written proofs."""
from pathlib import Path
import hashlib,json,platform
import sympy as S
import simple_ramification_ledger as L

P=Path(__file__).resolve().parent


def main():
    checks=[]
    def check(name,value):
        assert value,name
        checks.append(name)
    def ident(name,z):check(name,S.simplify(z)==0)
    def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
    olds=['half_support_inertia_receipt.json','triple_boundary_residue_receipt.json',
          'cyclic_quotient_morse_receipt.json']
    for name in olds:
        r=json.loads((P/name).read_text());check('sealed predecessor PASS '+name,r['outcome']=='PASS')
        for f,h in r['source_sha256'].items():check('sealed predecessor '+name+' '+f,sha(P/f)==h)
    manifest=json.loads((P/'sources/simple-ramification-ledger/manifest.json').read_text())
    for item in manifest['sources']:
        f=P/item['path'];check('primary/source seal '+f.name,sha(f)==item['sha256'] and f.stat().st_size==item['bytes'])
    check('context hypotheses read, no classification imported',manifest['primary_HTML_hypotheses_read'] and not manifest['Kulikov_global_generic_cover_applied'] and not manifest['Orevkov_all_preimage_bijectivity_assumed'] and not manifest['fresh_PDF_audit_claimed'])
    x,y,t,z=S.symbols('x y t z');m=S.symbols('m',integer=True,positive=True)
    coefficient=S.integrate(1-z**m,(z,0,1))
    ident('UNBOUNDED critical image leading coefficient',coefficient-m/(m+1))
    check('UNBOUNDED coefficient cannot cancel',S.ask(S.Q.positive(m/(m+1))) is True)
    examples=[]
    def order(poly,var):return min(mon[0] for mon,c in S.Poly(S.expand(poly),var).terms())
    for mm in [2,3,7]:
        phi=y**mm*(1+2*y);V=1+x+3*y
        G=S.integrate(V*(x-phi),y);G=S.expand(G-G.subs(y,0))
        ident('actual unit-perturbed Jacobian '+str(mm),S.diff(G,y)-V*(x-phi))
        check('actual unit remains invertible '+str(mm),V.subs({x:0,y:0})==1)
        up=phi.subs(y,t);vp=S.expand(G.subs({x:up,y:t}))
        check('actual coprime critical-image orders '+str(mm),order(up,t)==mm and order(vp,t)==mm+1 and S.gcd(mm,mm+1)==1)
        unit=S.cancel(G.subs(x,0)/y**(mm+1))
        check('actual LOCAL fibre length, not global zero count '+str(mm),not S.denom(unit).has(y) and unit.subs(y,0)!=0)
        ident('actual image coefficient '+str(mm),S.Poly(vp,t).coeff_monomial(t**(mm+1))-S.Rational(mm,mm+1))
        examples.append(dict(multiplicity=mm,local_length=mm+1,excess=mm-1))
    r=S.symbols('r',integer=True,positive=True);k=S.symbols('k',integer=True,nonnegative=True)
    for M in [0,1]:ident('UNBOUNDED half-support exact balance '+str(M),2*r+r*M+(2*r-1-r*M)-(4*r-1))
    check('UNBOUNDED every defect>=2 impossible',(-1-r*k).is_negative)
    for rr in [1,2,8,10**20]:
        cases=L.half_support_cases(rr)
        check('two necessary cases '+str(rr),[(c['ordinary_cusps'],c['unramified_primes']) for c in cases]==[(0,2*rr-1),(1,rr-1)])
        check('two exact global ledgers '+str(rr),all(2*rr+c['total_excess']+c['unramified_primes']==c['degree']-1 for c in cases))
        check('no accidental binary affine realization/exclusion '+str(rr),all(not c['actual_plane_realization'] and not c['binary_affine_family_excluded'] for c in cases))
    double=L.source_fibre(32,[(8,[1,1])],[])
    cusp=L.source_fibre(32,[(8,[2])],[])
    check('double value consumes every sheet',double['source_points']==0 and len(double['ramified_local_lengths'])==16)
    check('cusp has eight etale normalization sheets before possible omissions',cusp['source_points']==8 and cusp['ramified_local_lengths']==[3]*8)
    check('other unramified image can remove fixed cusp sheets',L.source_fibre(32,[(8,[2])],[2,1])['source_points']==5)
    check('a cusp cannot share a double value',L.source_fibre(32,[(8,[2,1])],[])['source_points']<0)
    check('multiple-image multiplicity sum',L.source_fibre(10,[(1,[2]),(2,[1])],[1])['source_points']==2)
    # Actual affine S3 action: standard GL2(F2) on two bits, trivial on the rest.
    A=tuple(v ^ ((v&1)<<1) for v in range(32))
    B=tuple(v ^ ((v>>1)&1) for v in range(32))
    def mul(a,b):return tuple(a[b[i]] for i in range(32))
    I=tuple(range(32));group={I};todo=[I]
    while todo:
        g=todo.pop()
        for h in [A,B]:
            f=mul(g,h)
            if f not in group:group.add(f);todo.append(f)
    check('actual binary affine cusp image S3',len(group)==6 and mul(A,A)==I and mul(B,B)==I and mul(mul(A,B),A)==mul(mul(B,A),B))
    todo=set(range(32));orbits=[]
    while todo:
        a=min(todo);orbit={g[a] for g in group};todo-=orbit;orbits.append(len(orbit))
    check('actual binary affine cusp orbit profile',sorted(orbits)==[1]*8+[3]*8)
    text=' '.join((P/'SIMPLE_RAMIFICATION_MULTIPLICITY_LEDGER_AND_HALF_SUPPORT_CUSPS.md').read_text().lower().split())
    for phrase in ['all omitted primes','exact local degree','not an assertion','exact pointwise',
        'same parameters','no two ordinary cusps','arbitrary contact','fixed points can still be removed',
        'hypothesis fails','every irreducible component','l4','behaviour at infinity',
        'not exclude its entire binary affine family','no historical priority']:
        check('written scope '+phrase,phrase in text)
    files=['SIMPLE_RAMIFICATION_MULTIPLICITY_LEDGER_AND_HALF_SUPPORT_CUSPS.md','simple_ramification_ledger.py',
        'verify_simple_ramification_ledger.py','sources/simple-ramification-ledger/manifest.json']+olds+[i['path'] for i in manifest['sources']]
    result=dict(claim='PF2RULED119',outcome='PASS',check_count=len(checks),checks=checks,
      theorem=dict(all_simple_index2_plane_normalizations=True,normalization_smooth_and_all_residues1=True,
        local_length_equals_image_multiplicity_plus1=True,singular_image_pair_consecutive=True,
        local_symmetric_monodromy=True,exact_global_multiplicity_ledger=True,exact_pointwise_source_count=True,
        half_support_involutory_two_cases=True,half_support_at_most_one_ordinary_cusp=True,
        higher_and_multiple_cusps_excluded_in_half_support_case=True,half_support_cusp_local_group_S3=True),
      actual_unit_perturbed_controls=examples,half_support32=L.half_support_cases(8),
      cusp32_orbits=sorted(orbits),
      assurance=dict(PF105_local_calculation_credited=True,unbounded_proof_written=True,
        finite_samples_prove_theorem=False,global_symmetric_group_assumed=False,
        Kulikov_total_ramification_birationality_assumed=False,Orevkov_all_preimage_bijectivity_assumed=False,
        immersed_affine_curve_complement_assumed_cyclic=False,analytic_L4_gap_closed=False,
        binary_affine_family_excluded=False,all_degree6_maps_excluded=False,
        fixed_cusp_sheets_all_assumed_retained=False,source_surjectivity_assumed=False,
        actual_finite_Keller_realization_constructed=False,actual_global_finite_algebra_closed=False,
        historical_priority_claimed=False,formal_verification=False,external_review=False,research_goal_completed=False),
      source_sha256={f:sha(P/f) for f in dict.fromkeys(files)},environment=dict(python=platform.python_version(),sympy=S.__version__))
    (P/'simple_ramification_ledger_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',len(checks),'PF119 controls; all-simple multiplicity ledger and two half-support cusp cases; full goal ACTIVE')


if __name__=='__main__':main()
