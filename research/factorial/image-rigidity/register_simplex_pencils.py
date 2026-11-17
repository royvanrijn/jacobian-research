"""Register the all-dimensional pencil progress with narrow, guarded edits."""
from pathlib import Path
import hashlib
import json
import os

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='factorial/image-rigidity/'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

# Historical continuation must not roll back later source-review scopes.
if any(r['id']=='FIRTWOANGLINEAR20261116' for r in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Later two-angular progress is registered; use register_two_angular.py.')


def update_array(path,key,changes,additions=()):
    before=path.read_bytes();text=before.decode();data=json.loads(text)
    known={r['id']:r for r in data[key]}
    encode=lambda row:'\n'.join('    '+line for line in json.dumps(row,indent=2).splitlines())
    for id,change in changes.items():
        old=known[id];row=change(dict(old))
        if row!=old:
            needle=encode(old)
            assert text.count(needle)==1,('record formatting changed',id)
            text=text.replace(needle,encode(row),1);known[id]=row
    new=[]
    for row in additions:
        if row['id'] in known:
            assert known[row['id']]==row,('existing pencil record differs',row['id'])
        else:new.append(row)
    if new:
        head,end=text.rsplit('\n  ]\n}',1);assert not end.strip()
        text=head+',\n'+',\n'.join(encode(r) for r in new)+'\n  ]\n}'+end
    if text.encode()!=before:
        assert path.read_bytes()==before,'Concurrent edit'
        temp=path.with_suffix(path.suffix+'.pencil-research/.scratch');temp.write_bytes(text.encode());os.replace(temp,path)


# A new continuation must validate actual previous artifacts, not its summary.
for name in ('verification.json','long_audit_verification.json'):
    receipt=json.loads((HERE/name).read_text())
    assert receipt['status']=='pass'
    for file,expected in receipt['files'].items():
        assert sha(HERE/file)==expected,('previous proof/control changed',file)

ledger_path=HERE/'source_ledger.json'
ledger=json.loads(ledger_path.read_text())
source='sources/deboor-0502036v1.html'
if source not in {r['file'] for r in ledger['cached_primary_sources']}:
    ledger['cached_primary_sources'].append(dict(file=source,url='https://arxiv.org/html/math/0502036v1',sha256=sha(HERE/source),
        review='Genocchi--Hermite and confluent divided-difference/contour formulas read; the current proof additionally derives the needed formula termwise on monomials, including repeated nodes.'))
for row in ledger['cached_primary_sources']:
    if row['file']=='sources/long-eit-vi-current.tex':
        row['review']='Theorem2.1 radial-limit proof checked in every dimension and adapted to fixed positive integer Gamma shapes. Other explicit higher-dimensional family claims and their Paper V nonvanishing inputs remain unaudited.'
for url in ['https://arxiv.org/abs/math/0502036','https://www.sciencedirect.com/science/article/pii/S0021869317300649']:
    if url not in ledger['additional_primary_references']:ledger['additional_primary_references'].append(url)
ledger.setdefault('published_input_pinpoints',{})['deBoor2005']='Genocchi--Hermite simplex formula and confluent contour divided differences; exact path weights are also independently derived on monomials in the pencil proof.'
ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')

proof=BASE+'SIMPLEX_PENCILS_AND_AFFINE_LIMITS.tex'
checker=BASE+'verify_simplex_pencils.py'
audit=dict(date='2026-11-16',previous_goal_turn='progress',
    previous_evidence='Both prior pass receipts were checked against every recorded current proof and checker hash before this continuation registration.',
    current_goal_turn='progress',canonical_proof=proof,
    mechanism=['finite sums of polynomial-weighted path integrals have no nonzero finite moment tail',
               'confluent Genocchi--Hermite with at least two distinct nodes produces path integrals and no point evaluations',
               'integer Dirichlet shapes and polynomial markers are repeated-node shifts',
               'homogeneous Gamma radial separation proves the full H(S,L) pencil class',
               'all-Q-algebra prime reduction retains actual coefficient nilpotence and ordinary torsion',
               'Lindemann plus a rational primitive pole-budget proves affine exponential nonvanishing',
               'the flat next-to-top affine layer permits arbitrary lower complex polynomial terms'],
    exact_producer='simplex_pencils.py: exact repeated-node chains, exponential coefficient groups, and membership identities for a supplied homogeneous radial pencil',
    boundaries=['general two-angular-coordinate simplex phases remain unresolved',
                'arbitrary inhomogeneous H(S,L) is outside the homogeneous theorem',
                'the radial sum S cannot be replaced by an arbitrary second linear form',
                'noninteger Dirichlet shapes are not covered by this chain theorem',
                'algebraic affine exponential periods may be algebraic and nonzero; transcendental nodes can give zero',
                'formal binary coefficients cannot replace actual composed coefficients at node collisions',
                'general SFC, full Paper V/VI family review and bibliography review remain open'],
    independent_replay=False,formal_verification=False,external_review=False,novelty_claim=False,whole_goal_complete=False)
(HERE/'simplex_pencil_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
common=dict(kind='theorem',state='proved',canonical_source=proof,checker=checker,proof_type='hybrid',
    independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_simplex_pencils.py'),
    software_lock=[proof,checker,BASE+'simplex_pencils.py',BASE+'moment_bridges.py',BASE+'simplex_pencil_verification.json',BASE+'simplex_pencil_examples.json',BASE+'simplex_pencil_audit.json',BASE+'source_ledger.json'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],consumers=[],invalidates_assumptions=[],
    replaced_by=[],priority='reference',programme_status='archived')
new=[
dict(common,id='FIRCHAINMOMENT20261116',title='Polynomial-weighted finite chain functionals have no nonzero finite moment tail along a nonconstant polynomial',
     scope='For Lambda(F)=sum_j integral_a_j^b_j r_j(z)F(z)dz with polynomial complex amplitudes and arbitrary complex endpoints, a nonzero polynomial H satisfying Lambda(H(p)*p^k)=0 for all k>=0 and nonconstant p forces Lambda(F(p))=0 for every entire F. In particular eventual moment zero forces all moments including the constant to vanish. Common inverse-tree paths, analytic integrable edge densities, connected complement and Cauchy jumps prove the assertion for the whole finite sum. No arithmetic assumption or common amplitude is imposed. A point evaluation is excluded and supplies a sharp negative control. Standard polynomial moment mechanism credited; no novelty claim.',
     dependencies=[]),
dict(common,id='FIRSIMPLEXLINE20261116',title='Every nonzero polynomial of one complex affine simplex coordinate has infinitely many nonzero moments for integer Dirichlet shapes',
     scope='Every dimension, every positive INTEGER Dirichlet shape vector, arbitrary complex p and complex affine linear form ell on the simplex. Eventual pure moment zero implies the actual phase p(ell) is zero. For any nonzero such phase and ANY polynomial marker w, eventual marked zero implies identically zero marked moments including exponent0. Confluent divided differences produce finitely many polynomial-weighted paths whenever ell is nonconstant; all-node coincidence is separately a point evaluation and constant phase case. Integer shape grouping and monomial shape shifts retain all original markers. The same polynomial-chain argument extends to uniform polytopes and polynomial densities of nonzero total mass. This is not a theorem for general phases in two independent angular coordinates or noninteger shapes.',
     dependencies=['FIRCHAINMOMENT20261116','external: de Boor math0502036v1 Genocchi--Hermite and confluent divided differences']),
dict(common,id='FIRPENCILGAMMA20261116',kind='corollary',title='Homogeneous factorial and integer-Gamma pencils H(sum Y,sum c_i Y_i) satisfy FC in every dimension and retain ordinary all-ring normalization',
     scope='All dimensions/degrees and complex node collisions. H is HOMOGENEOUS and S=sum Y is the ACTUAL radial sum. Every nonzero F=H(S,L), L=sum c_iY_i, has infinitely many nonzero moments for independent positive-integer-shape Gamma variables; alpha_i=1 is factorial FC. Exact homogeneous substitution gives the corresponding Gaussian even-block sector. Over ANY commutative Q-algebra, allowing H and L coefficients to vary, eventual pure zero, coefficientwise nilpotence of ACTUAL F, polynomial nilpotence and all-original-marker safety are equivalent. Scaled actual-period normalization equals Gamma_J(B)=union Ann(J^k), for the actual composed coefficient ideal J, including every collision/nilpotent stratum and non-Noetherian bases. Formal H coefficients give a wrong ideal on collisions. No arbitrary inhomogeneous H(S,L), arbitrary pair of linear forms, noninteger Gamma shape, first-moment bound, or unrestricted SFC conclusion.',
     dependencies=['FIRSIMPLEXLINE20261116','FIRGAMMAADAPTER20261116']),
dict(common,id='FIRAFFINESIMPLEX20261116',title='Affine exponential periods of every integer Dirichlet simplex are nonzero at algebraic coefficients, including collisions',
     scope='All positive integer Dirichlet shapes and algebraic affine node values. The period is a confluent exponential divided difference times the nonzero probability normalization. Lindemann--Weierstrass reduces a hypothetical zero to vanishing of every local residue coefficient of e^z/Q(z). Then a rational solution R of Rprime+R=1/Q would have pole budget N-s but require infinity decay N, impossible. Individual exponential group coefficients CAN vanish. Nonconstant triangle phase nodes(0,1,1) gives coordinate integral1/probability expectation2, so no general transcendence claim. Transcendental nodes(0,2pi*i) give zero and are rejected by the executable certificate. No arbitrary polynomial-phase simplex EIT assertion.',
     dependencies=['FIRSIMPLEXLINE20261116','external: Lindemann--Weierstrass theorem']),
dict(common,id='FIRFLATNEXTLINE20261116',kind='corollary',title='Every-dimensional flat factorial phases with affine next radial layer satisfy FC despite arbitrary lower complex terms',
     scope='Independent integer-shape Gamma variables, shape sum A, S=sum Y, D>=2, f=S^D+D*S^(D-2)*L+R with deg R<=D-2. For ALGEBRAIC coefficients of L and ARBITRARY complex lower R, E[f^m]/(A)_(Dm) tends to E_Dir exp(ell), which is nonzero; all sufficiently large moments survive. For arbitrary complex L, algebraic specialization proves FC; positive powers preserve the architecture, so infinitely many moments survive. A nonzero overall scalar may be removed. Degree1 f=S+c has direct limit exp(c). This does not assert a nonzero limit or eventual nonzero tail for every transcendental L. The all-layer radial proof is credited to Long--Le Blanc EIT VI Theorem2.1, with explicit integer Gamma adaptation; no Paper V geometric nonvanishing input is used here.',
     dependencies=['FIRAFFINESIMPLEX20261116','external: Long--Le Blanc EIT VI commit92db79e Theorem2.1 radial limit']),
dict(common,id='FIRPENCILEXEC20261116',proof_type='exact_symbolic',title='Exact confluent chain export and supplied homogeneous radial-pencil recognition produce finite identities',
     scope='Implemented simplex_pencils.py accepts exact algebraic nodes and positive integer shapes, groups collisions in an exact algebraic field, exports polynomial chain amplitudes with a separately typed point case, and groups affine exponential coefficients. A supplied-L recognizer transforms a homogeneous phase into radial/binary coordinates and verifies the reconstructed polynomial identity exactly. Accepted nonzero instances invoke the all-degree pencil theorem; no finite moment prefix is used as universal evidence. Outside-class answers make no conjecture decision. Unsupported noninteger shapes, inexact floats and nonalgebraic nodes for arithmetic exponential certificates are rejected. No first nonzero exponent or Strong Factorial window is asserted.',
     dependencies=['FIRPENCILGAMMA20261116','FIRAFFINESIMPLEX20261116'])]
new_ids=[r['id'] for r in new]
def open_scope(row):
    row['scope']='The full user goal remains open: general FC for arbitrary complex polynomial phases in dimension>=3, unrestricted SFC, complete primary bibliography and EIT VI/Paper V family review. Full rigidity and Long binary FC/EIT have internal written proof audits; unrestricted Image is false by SIC2C4. FIRCHAINMOMENT20261116/FIRSIMPLEXLINE20261116/FIRPENCILGAMMA20261116 now close EVERY homogeneous radial pencil H(S,L) in all dimensions and all positive integer Gamma shapes, including all original markers, node collisions and ordinary all-Q-algebra nilpotent normalization. FIRAFFINESIMPLEX20261116/FIRFLATNEXTLINE20261116 close the algebraic affine exponential and flat next-layer class with arbitrary lower terms. FIRPENCILEXEC20261116 implements exact supplied-pencils, with outside-class conclusions deliberately absent. General triangle phases depending on two independent angular coordinates remain unresolved, as do arbitrary noninteger shapes, general inhomogeneous H(S,L), and arbitrary two-linear-form substitutions. No finite-prefix, first-survival, SFC-window, nonzero-limit-at-every-transcendental-node, formal, external-review, bibliography-completeness or full-goal claim is made.'
    row['dependencies']=list(dict.fromkeys(row['dependencies']+new_ids));return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':open_scope},new)
known={r['id']:r for r in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][k] for k in ('state','title','scope')).encode()).hexdigest()[:16]
def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+new_ids))
    row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[proof,BASE+'simplex_pencil_audit.json']));return row
def work(row):
    row=refresh(row)
    row['next_step']='Seek a theorem for phases in TWO independent angular coordinates beyond the closed one-affine-coordinate simplex class, audit the EIT VI explicit families and exact Paper V inputs, finish Rocks2017 and the primary bibliography, and build original-unit period certificates for phases outside a supplied homogeneous radial pencil. Do not repeat closed H(S,L) degree/support searches or infer SFC from eventual nonvanishing.'
    row['prerequisites']='Read SIMPLEX_PENCILS_AND_AFFINE_LIMITS.tex, exact chain/recognition outputs and the completed Long/rigidity audits. Retain the typed point-evaluation branch at all-node coincidence, actual composed coefficient ideals on collision strata, noninteger-shape boundary, and algebraic-nonzero versus transcendental exponential claims. General two-angular-coordinate phases remain the mathematical gap; external review and formal verification are assurance fields, not a new permission blocker.'
    return row
def lesson(row):
    row=refresh(row)
    row['use']='For a homogeneous F=H(S,L), use the full all-dimensional integer-Gamma/simplex chain theorem and supplied-L identity recognizer. Repeated coefficients require confluent divided differences; at least two distinct nodes give polynomial paths, while all-node coincidence is a separately typed point evaluation. Polynomial marker shape shifts preserve the chain proof. Use the ACTUAL F coefficient ideal for all-ring nilpotence and Gamma_J normalization. Algebraic affine exponential nonvanishing supplies a flat next-layer class with arbitrary lower terms.'
    extra=' Do not replace the radial sum by an arbitrary linear form, drop homogeneity, extend integer-shape chains to arbitrary half-integer weights, infer transcendence from affine exponential nonvanishing, or treat an outside-supplied-pencil result as a general FC decision.'
    if extra not in row['avoid']:row['avoid']+=extra
    row['boundary']='All homogeneous radial pencils and one-affine-coordinate simplex phases with positive integer shapes are closed in every dimension. General FC3/two-angular-coordinate phases, unrestricted SFC, noninteger shapes, complete EIT V/VI family review and bibliography review remain open. Long binary FC/EIT and rigidity have internal written audits, with formal and external assurance false.'
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'simplex_pencils.py',checker]));return row
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
print(json.dumps({'registered':new_ids,'previous_goal_turn':'progress','goal_complete':False}))
