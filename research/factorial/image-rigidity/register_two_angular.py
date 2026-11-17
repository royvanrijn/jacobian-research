"""Register two-coordinate progress without rewriting unrelated ledger bytes."""
from pathlib import Path
import hashlib
import json
import os

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'factorial/image-rigidity/'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

# Historical continuation must not roll back later source/window scopes.
if any(row['id']=='FIRRECIPWINDOW20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Later window progress is registered; use register_moment_windows.py.')


def update_array(path, key, changes, additions=()):
    before = path.read_bytes()
    text = before.decode()
    known = {r['id']: r for r in json.loads(text)[key]}
    encode = lambda row: '\n'.join('    '+line for line in json.dumps(row, indent=2).splitlines())
    for id, change in changes.items():
        old = known[id]
        row = change(dict(old))
        if row != old:
            needle = encode(old)
            assert text.count(needle) == 1, ('record formatting changed', id)
            text = text.replace(needle, encode(row), 1)
            known[id] = row
    new = []
    for row in additions:
        if row['id'] in known:
            assert known[row['id']] == row, ('existing two-coordinate record differs', row['id'])
        else:
            new.append(row)
    if new:
        head, end = text.rsplit('\n  ]\n}', 1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(encode(r) for r in new)+'\n  ]\n}'+end
    if text.encode() != before:
        assert path.read_bytes() == before, 'Concurrent edit'
        temp = path.with_suffix(path.suffix+'.two-angular-research/.scratch')
        temp.write_bytes(text.encode())
        os.replace(temp, path)


# Validate actual retained sources before classifying previous progress.
for filename in ('verification.json', 'long_audit_verification.json',
                 'simplex_pencil_verification.json', 'two_angular_verification.json'):
    receipt = json.loads((HERE/filename).read_text())
    assert receipt['status'] == 'pass'
    for file, expected in receipt['files'].items():
        assert sha(HERE/file) == expected, ('proof/control source changed', file)

ledger_path = HERE/'source_ledger.json'
ledger = json.loads(ledger_path.read_text())
file = 'sources/long-eit-v-current.tex'
if file not in {r['file'] for r in ledger['cached_primary_sources']}:
    ledger['cached_primary_sources'].append(dict(file=file,
        url='https://raw.githubusercontent.com/long-mathematics/exponential-integral-theorem/92db79eedeb064694ff20e9201d73fc2bc0857f4/papers/05-polynomial-flags/polynomial-flag-exponential-periods.tex',
        sha256=sha(HERE/file),
        review='Linear criterion, actual polynomial Stokes construction, strict-extremum detector, normal form and ordinary-point lifting checked. The cubic rank-two case is independently proved here by a rational Riccati residue obstruction. General geometric/Galois independence, saddle/continued-density detectors and other examples are NOT fully audited.'))
for row in ledger['cached_primary_sources']:
    if row['file'] == 'sources/long-eit-vi-current.tex':
        row['review'] = 'Theorem2.1 radial limit checked in every dimension and adapted to integer-total-shape Gamma variables. Odd-family slicing constants and cubic phase checked; a standalone rank-two proof covers the needed a=3 linear input. General a geometric inputs, original elliptic continued-density detector, full block families and complete finite gcd range remain outside the completed audit. The h=0 good-reduction argument here is distinct from VI\'s h-variable mod-five matrix.'
pin = ledger['published_input_pinpoints']['Beukers2006']
if 'explicit rank-two cubic' not in pin:
    ledger['published_input_pinpoints']['Beukers2006'] += ' Also used for the explicit rank-two cubic system: actual boundary and interior E-functions have only finite pole0; functional injectivity forces every lifted relation to have zero interior coefficients.'
ledger_path.write_text(json.dumps(ledger, indent=2)+'\n')

proof = BASE+'TWO_ANGULAR_COORDINATES.tex'
checker = BASE+'verify_two_angular.py'
audit = dict(date='2026-11-16', previous_goal_turn='progress', current_goal_turn='progress',
    previous_evidence='Every recorded hash in all three previous pass receipts matches the retained proof/checker sources.',
    canonical_proof=proof,
    mechanisms=['free rank-two polynomial twisted quotient and actual Stokes primitives',
                'cyclic scalar equation and rational Riccati residue contradiction',
                'strict interior minimum separates actual periods from specified boundary moments',
                'ordinary-point Beukers linear lifting with redundant boundary coordinates',
                'h0 evaluation at a critical point with unit barycentric values handles every integer exponent triple',
                'Gamma grouping and radial limit give degree-four families with arbitrary lower complex terms',
                'finite exceptional sets, algebraic-shape specialization and power closure'],
    source_review='Required linear/extremum mechanism checked; cubic case proved directly. Full EIT V/VI review remains open.',
    boundaries=['failed valuation and zero remainder tests make no zero-period decision',
                'uniform arithmetic theorem requires a strictly negative valuation at a place above5',
                'three block shape sums must be positive integers; individual positive shapes may aggregate to them',
                'complex-scale specialization assumes algebraic individual shapes and excludes finite algebraic exceptions',
                'general FC>=3, unrestricted SFC, general coupled phases and full bibliography remain open'],
    coordinator_handoff='Inward coupled restrictions received as provisional. Binary feedback matches the already registered theorem; unrelated geometric/Mathieu exclusions have not been adopted without audit.',
    independent_replay=False, formal_verification=False, external_review=False,
    novelty_claim=False, whole_goal_complete=False)
(HERE/'two_angular_audit.json').write_text(json.dumps(audit, indent=2)+'\n')

common = dict(kind='theorem', state='proved', canonical_source=proof, checker=checker,
    proof_type='hybrid', independent_replay=False, formal_verification=False, external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_two_angular.py'),
    software_lock=[proof, checker, BASE+'two_angular.py', BASE+'moment_bridges.py',
                   BASE+'two_angular_verification.json', BASE+'two_angular_examples.json',
                   BASE+'two_angular_audit.json', BASE+'source_ledger.json'],
    supersedes=[], closes_problems=[], narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[], invalidates_assumptions=[], replaced_by=[], priority='reference', programme_status='archived')
new = [
dict(common, id='FIRTWOANGLINEAR20261116',
     title='The cubic triangle phase has an elementary irreducible rank-two actual linear comparison at every nonzero algebraic scale',
     scope='q=X^3-3X+Y^2, T={0<=X<=2, abs(Y)<=1-X/2}. The additive polynomial twisted quotient has free basis1,X and exact polynomial Stokes primitives in h=1/z. Connection A=[[-5/(6z),-2],[-2,-7/(6z)]] is irreducible: gauge v=zu gives vdoubleprime=(4-5/(36z^2))*v, and rational Riccati residues1/6 or5/6 at0 plus residue1 at other poles contradict zero infinity residue. A strict unique interior minimum separates the actual column from the specified boundary span; rank2 simplicity allows adjoining1. Closed interior/edge E-functions have only finite pole0. Beukers lifting gives numerical linear independence modulo algebraic specified-boundary values and1 at EVERY nonzero algebraic scale. Every algebraic polynomial amplitude with nonzero specialized remainder has a transcendental nonzero period. Zero remainder gives ONLY a boundary relation: positive amplitude1+2Y^2 at scale1 is exact and has a nonzero integral. The cubic case does not depend on unreviewed general Paper V geometric/Galois theorems. Mechanism credited to Long--Le Blanc; no novelty or external assurance claim.',
     dependencies=['FIREIT20261116', 'external: Beukers2006 Theorem3.2 and Lemma3.1', 'external: Long--Le Blanc EIT V commit92db79e linear/extremum mechanism']),
dict(common, id='FIRTWOANGSHAPES20261116',
     title='One arithmetic scale class makes every three-integer-shape cubic-triangle exponential period nonzero',
     scope='For the fixed cubic triangle, barycentric u=(2-X-2Y)/4,v=(2-X+2Y)/4,w=X/2, EVERY amplitude u^a*v^b*w^c with a,b,c>=0 has a transcendental nonzero integral against exp(xi*q) if algebraic xi has STRICTLY NEGATIVE valuation at some place above5. The exact remainder is 5-integral in h; its h0 evaluation at critical point(1,0) is the unit4^(-a-b)*2^(-c). Thus h reducing to0 excludes zero remainder for ALL triples, with no finite extrapolation or quotient-ring multiplication. xi=N*i/5 for any integer N prime to5 allows arbitrarily large oscillation. For each fixed triple the remainder-exceptional scale set is finite/algebraic and consists only of scales integral at every place above5; emptiness of every set is not proved. Failed valuation or zero remainder tests make no zero-period/FC decision.',
     dependencies=['FIRTWOANGLINEAR20261116']),
dict(common, id='FIRTWOANGGAMMA20261116', kind='corollary',
     title='Degree-four families in two angular coordinates have eventually nonzero Gamma/factorial moments for every integer shape vector',
     scope='Partition variables into three nonempty blocks with sums B1,B2,B3, S=sum Bj; set G=8*B3^3-6*B3*S^2+(B2-B1)^2*S. For D>=4 and algebraic xi with negative valuation at a place above5, EVERY F=S^D+D*xi*S^(D-4)*G+R, deg R<=D-2 with ARBITRARY COMPLEX lower coefficients, has eventually EVERY pure moment nonzero for independent Gamma variables whose three block shape sums are positive INTEGERS. All dimensions>=3 and positive integer individual shapes are covered, including positive noninteger shapes aggregating to integer blocks. At xi0 the normalized limit is1. The limit is Gamma(A)/(4*product Gamma(Aj)) times the weighted triangle period and is independent of R. xi=101*i/5,D4 supplies one fixed all-dimensional factorial family with arbitrary lower cross-variable couplings; three even Gaussian blocks follow. For fixed block shapes every algebraic scale outside the finite remainder-exception set has the same eventual conclusion. For ALGEBRAIC individual shapes, COMPLEX scales outside that set have infinitely many nonzero moments by specialization and power closure; no all-transcendental nonzero-limit claim. The angular Hessian determinant12X excludes a one-affine-coordinate phase. No exceptional-scale decision, arbitrary noninteger block sums or SFC window is asserted.',
     dependencies=['FIRTWOANGSHAPES20261116', 'FIRFLATNEXTLINE20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRTWOANGEXEC20261116', proof_type='exact_symbolic',
     title='Exact two-angular Stokes, valuation and structured-Gamma exporters retain actual boundaries and fail closed',
     scope='two_angular.py exports rational polynomial Stokes identities r=a(h)+b(h)*X+Dx(P)+Dy(Q), certifies a negative place above5 from the monic minimal polynomial of an exact algebraic scale, gives an all-integer-block-shape certificate without scanning amplitudes, and checks fixed shape amplitudes by exact specialized number-field normal forms. The family constructor checks partition, degree and integer shapes and constructs the actual F. A zero remainder or failed uniform test is not a zero-period conclusion; no first nonzero exponent or SFC window is emitted. Exact controls verify finite identities, actual boundary countercontrols, grouping and original Gaussian/Gamma constants. Universal quantifiers use the written rank-two/E-function/h0-unit argument, not a finite test.',
     dependencies=['FIRTWOANGLINEAR20261116', 'FIRTWOANGSHAPES20261116', 'FIRTWOANGGAMMA20261116'])]
new_ids = [r['id'] for r in new]


def open_scope(row):
    row['scope'] = 'Full goal open: arbitrary complex polynomial FC in dimension>=3, unrestricted SFC, complete bibliography and full EIT V/VI review. Rigidity and Long binary FC/EIT have internal written audits; unrestricted Image is false by SIC2C4. Homogeneous radial pencils H(S,L), one-affine-coordinate integer-Dirichlet phases, all original markers/collisions and ordinary all-Q-algebra normalization are closed; affine exponential nonvanishing supplies a flat next-layer class. FIRTWOANGLINEAR20261116/FIRTWOANGSHAPES20261116/FIRTWOANGGAMMA20261116 add a standalone cubic two-angular rank-two/Riccati comparison and h0 good-reduction theorem: degree-four factorial/Gamma families in every dimension and integer shape with arbitrary complex lower terms, at algebraic scales nonintegral at a place above5. Positive shapes aggregating to three integer blocks are included. Generic complex scales yield infinitely many moments under algebraic-individual-shape and finite-exception boundaries. Exact Stokes/arithmetic certificates retain actual boundaries. General coupled triangle phases, arbitrary noninteger block sums and unrestricted FC/SFC remain unresolved. No finite-prefix, first-survival, uniform SFC window, all-transcendental nonzero-limit, formal/external assurance, bibliography-completeness or full-goal claim.'
    row['dependencies'] = list(dict.fromkeys(row['dependencies']+new_ids))
    return row


update_array(ROOT/'MATH_STATUS.json', 'entries', {'OP-FACTORIAL-GENERAL20261116': open_scope}, new)
known = {r['id']: r for r in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger = lambda id: hashlib.sha256('\0'.join(known[id][k] for k in ('state', 'title', 'scope')).encode()).hexdigest()[:16]


def refresh(row):
    row['claims'] = list(dict.fromkeys(row['claims']+new_ids))
    row['reviewed_claims'] = {id: finger(id) for id in row['claims']}
    row['sources'] = list(dict.fromkeys(row['sources']+[proof, BASE+'two_angular_audit.json']))
    return row


def work(row):
    row = refresh(row)
    row['next_step'] = 'Extend the explicit rank-two comparison to broader coupled phases or prove all weighted cubic remainder exception sets empty; retain actual boundary terms. Complete untouched EIT V geometric/Galois/continued-density inputs, VI nonflat families, Rocks2017 and bibliography. General FC>=3, unrestricted SFC and arbitrary noninteger block sums remain the goal. Do not repeat closed radial-pencil searches or infer SFC from eventual nonvanishing.'
    row['prerequisites'] = 'Read TWO_ANGULAR_COORDINATES.tex and exact Stokes/valuation/Gamma outputs. The cubic has a standalone Riccati proof; no general unreviewed V theorem is required. Negative valuation is strict; zero remainder retains actual boundary contributions. Complex-scale specialization assumes algebraic individual shapes and excludes a finite algebraic set. Earlier pencil/collision and all-ring boundaries remain in force. The inward coordinator packet is provisional outside already verified binary feedback.'
    return row


def lesson(row):
    row = refresh(row)
    if 'two-angular cubic q=' not in row['use']:
        row['use'] += ' For the two-angular cubic q=X^3-3X+Y^2 on the specified triangle, use the rank-two/Riccati comparison, exact Stokes normal form and h0 unit evaluation at critical barycentric values(1/4,1/4,1/2). A negative valuation at a place above5 handles all integer shape triples and yields flat degree-four Gamma/factorial families with arbitrary lower complex terms.'
    extra = ' A zero specialized twisted remainder is only a boundary relation, not a zero period. Do not replace strictly negative valuation by an all-scale assertion, finite exception avoidance by a nonzero limit at transcendental scales, or integer block sums by arbitrary positive sums.'
    if extra not in row['avoid']:
        row['avoid'] += extra
    row['boundary'] = 'Closed classes include homogeneous radial pencils and a genuine two-angular cubic/flat-quartic family in every dimension and integer Gamma shape, including positive individual shapes aggregating to integer blocks. General FC>=3, arbitrary coupled phases, unrestricted SFC, arbitrary noninteger block sums, full EIT V/VI review and bibliography remain open. Formal, independent replay and external assurance remain false.'
    row['implementation'] = list(dict.fromkeys(row['implementation']+[BASE+'two_angular.py', checker]))
    return row


update_array(ROOT/'knowledge/work_items.json', 'actions', {'WORK-FACTORIAL-GENERAL': work})
update_array(ROOT/'knowledge/lessons.json', 'lessons', {'METHOD-FACTORIAL-RADIAL-PURITY': lesson})
print(json.dumps({'registered': new_ids, 'previous_goal_turn': 'progress', 'whole_goal_complete': False}))
