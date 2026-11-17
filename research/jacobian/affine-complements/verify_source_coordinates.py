"""Exact controls for the two-chart all-height source-coordinate transfer.

Coverage, Cox geometry, all-height induction and the published degree
alternative are written/imported inputs, not formally certified here.
"""
from pathlib import Path
from time import monotonic
from math import gcd
import hashlib,json,platform
import sympy as S
HERE=Path(__file__).resolve().parent
start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(value)==0)
x,y,s,t,alpha,beta,a,b=S.symbols('x y s t alpha beta a b')
local={str(v):v for v in (x,y,s,t,alpha,beta,a,b)}
def parse(value):return S.sympify(value,locals=local)
def degree(value):return S.Poly(S.expand(value),x,y).total_degree()
def top(value):
    p=S.Poly(S.expand(value),x,y);d=p.total_degree()
    return S.expand(sum(c*x**i*y**j for (i,j),c in p.terms() if i+j==d))

data=json.loads((HERE/'source_coordinate_certificates.json').read_text())
check('producer source seal',data['source_sha256']==hashlib.sha256(
    (HERE/'compute_source_coordinates.py').read_bytes()).hexdigest())
check('exactly two contact-position charts',len(data['cases'])==2)
for row in data['cases']:
    name=row['name'];marked=name=='contact_at_marked'
    bb=list(map(parse,row['boundary_columns']['B']))
    ee=list(map(parse,row['boundary_columns']['E']))
    zero('contact-four determinant '+name,bb[0]*ee[1]-bb[1]*ee[0]-s**4)
    B=alpha*bb[0]+beta*bb[1];E=alpha*ee[0]+beta*ee[1]
    zero('stored B equation '+name,B-parse(row['cox_equations']['B']))
    zero('stored E equation '+name,E-parse(row['cox_equations']['E']))
    forward={local[k]:parse(v) for k,v in row['forward'].items()}
    inverse_y=parse(row['inverse_y'])
    zero('forward B=1 '+name,B.subs(forward,simultaneous=True)-1)
    zero('forward E=1 '+name,E.subs(forward,simultaneous=True)-1)
    zero('forward followed by inverse y '+name,inverse_y.subs(forward,simultaneous=True)-y)
    check('x is precisely the slice s coordinate '+name,forward[s]==x)
    R=parse(row['single_relation']);tv=parse(row['eliminated_t'])
    zero('eliminated relation is B equation '+name,
         R-(b*b if marked else b)*(B.subs(t,tv)-1))
    zero('eliminated t satisfies the E relation '+name,
         (E-1-(a*s+b*t)*(B-1)).subs(t,tv))
    for key in ('alpha','beta','t'):
        value=forward[local[key]].subs({x:s,y:inverse_y},simultaneous=True)
        expr=(value-local[key]).subs(t,tv)
        multiplier=parse(row['reverse_identity_multipliers'][key])
        zero('literal inverse certificate '+name+' '+key,expr-multiplier*R)
        check('inverse certificate only inverts b '+name+' '+key,
              not(S.denom(S.cancel(multiplier)).free_symbols-{b}))
    H=x**3*y
    tops=[H**4/b**3,H**3/b**2,-H**2/b,H]
    if not marked:tops=[-H**4/b**3,-H**3/b**2,H**2/b,H]
    raw=list(map(parse,row['ordered_section_basis']))
    for j,(section,expected) in enumerate(zip(raw,tops)):
        pulled=S.expand(section.subs(forward,simultaneous=True))
        zero('complete pulled section '+name+' '+str(j),pulled-parse(row['source_sections'][j]))
        zero('uniform top term '+name+' '+str(j),top(pulled)-expected)
        check('section degree '+name+' '+str(j),degree(pulled)==16-4*j)
        restriction=section.subs(beta,-alpha*s/t**2) if marked else section.subs(alpha,-beta*s*s/t)
        frame=raw[0].subs(beta,-alpha*s/t**2) if marked else raw[0].subs(alpha,-beta*s*s/t)
        sign=[1,-1,-1,-1][j] if marked else [1,1,-1,-1][j]
        zero('contact order and binary basis '+name+' '+str(j),restriction/frame-sign*(s/t)**j)
        check('no exceptional a-parameter denominator '+name+' '+str(j),
              not(S.denom(S.cancel(pulled)).free_symbols-{b}))

# Bounded controls for a separately written all-height factorization.
gens={(1,0,1,0),(0,1,1,0),(2,0,0,1),(1,1,0,1),(0,2,0,1)}
monomials={(0,0,0,0)};vanishing={0}
for m in range(1,13):
    monomials={tuple(v[i]+g[i] for i in range(4)) for v in monomials for g in gens}
    expected={(i,m+j-i,m-j,j) for j in range(m+1) for i in range(m+j+1)}
    check('Cox multiplication spans V_m control '+str(m),monomials==expected)
    check('section dimension control '+str(m),len(monomials)==(m+1)*(3*m+2)//2)
    vanishing={i+j for i in vanishing for j in range(4)}
    check('all binary monomial orders lift control '+str(m),vanishing==set(range(3*m+1)))
    check('contact-weight degree controls '+str(m),all(16*m-4*r>=4*m for r in vanishing))
check('Cox B,E degrees form a lattice basis',abs(S.det(S.Matrix([[1,1],[1,2]])))==1)

manifest_path=HERE/'sources/degree-transfer-manifest.json'
manifest=json.loads(manifest_path.read_text())
for entry in manifest['sources']:
    blob=(manifest_path.parent/entry['file']).read_bytes()
    check('primary pin bytes '+entry['file'],len(blob)==entry['bytes'])
    check('primary pin seal '+entry['file'],hashlib.sha256(blob).hexdigest()==entry['sha256'])
check('published theorem remains external',manifest['external_theorem'] is True and
      manifest['independent_whole_proof_replay'] is False)
check('conditional 2026 degree125 upgrade is not imported',
      manifest['conditional_2026_degree_125_exclusion_used'] is False)
tex=(HERE/'sources/2204.14178v1.tex').read_text()
check('pinned theorem retains the72,108 alternative',
      '$(\\deg(P),\\deg(Q))\\in\\{(72,108),(108,72)\\}$.' in tex)
check('all heights at most six have degree at most96',16*6==96<108)
for smaller in range(1,7):
    check('height-seven coprime leading-root case '+str(smaller),gcd(smaller,7)==1)
    check('smaller coordinate degree cannot be108 '+str(smaller),16*smaller<108)
check('nonvanishing seventh symbol gives112, not108 or125',16*7==112 and 108<112<125)
for multiplicity in (1,2,3):
    bound=max(16*6,16*7-4*7*multiplicity)
    check('vanishing seventh symbol degree bound '+str(multiplicity),bound==96)
open_pairs=[(m,8) for m in range(4,8) if 8%m!=0]
check('structural height-eight candidates before PF2RULED12',open_pairs==[(5,8),(6,8),(7,8)])
report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    source_charts=2,parameter_domain='Q[a,b,b^-1], b!=0',uniform_pole_to_degree_multiplier=16,
    section_degrees=[16,12,8,4],section_contact_orders=[0,1,2,3],
    written_results=['explicit A2 coordinates for every forced F1 contact-four boundary pair',
                     'all-height pole-to-degree bound and contact-sensitive improvement',
                     'no compatible quartic pole presentation of height<=7, using external2022 Theorem2.1'],
    excluded_pole_height_upper_bound=7,
    subsequent_all_height_pole_three_exclusion="PF2RULED8 (separate intrinsic proof and receipt)",
    frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[pair for pair in open_pairs if pair!=(5,8)],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
    structural_height_eight_candidates_before_PF2RULED12=open_pairs,
    assurance=dict(coverage='written two-chart classification',all_height_bounds='written section-factorization proof',
        plane_degree_theorem='external arXiv:2204.14178v1 Theorem2.1; not independently reproved',
        conditional_2026_degree125_upgrade_used=False,actual_ruled_compatibility_proved=False,
        unrestricted_higher_poles_excluded=False,quartic_packets_excluded=False,
        formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [
       Path(__file__),HERE/'SOURCE_COORDINATES.md',HERE/'POLE_SEVEN_EXCLUSION.md',
       HERE/'compute_source_coordinates.py',HERE/'source_coordinate_certificates.json',manifest_path]},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'source_coordinates_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; all compatible heights<=7 excluded with external degree input; height eight open')
