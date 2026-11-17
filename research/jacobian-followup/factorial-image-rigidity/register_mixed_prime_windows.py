"""Register mixed-prime lifting and degree-intersection progress narrowly."""
from pathlib import Path
import hashlib
import json
import os
import pypdf

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'jacobian-followup/factorial-image-rigidity/'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()


def update_array(path, key, changes, additions=()):
    before = path.read_bytes()
    text = before.decode()
    known = {row['id']: row for row in json.loads(text)[key]}
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
            assert known[row['id']] == row, ('existing mixed-prime record differs', row['id'])
        else:
            new.append(row)
    if new:
        head, end = text.rsplit('\n  ]\n}', 1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(encode(row) for row in new)+'\n  ]\n}'+end
    if text.encode() != before:
        assert path.read_bytes() == before, 'Concurrent edit'
        temporary = path.with_suffix(path.suffix+'.mixed-window-tmp')
        temporary.write_bytes(text.encode())
        os.replace(temporary, path)


if (HERE/'hypergeometric_window_audit.json').exists():
    raise SystemExit('Historical mixed-prime registration: use register_hypergeometric_windows.py to preserve later scopes and reviews.')

for filename in ('verification.json', 'long_audit_verification.json', 'simplex_pencil_verification.json',
                 'two_angular_verification.json', 'moment_window_verification.json',
                 'prime_power_window_verification.json', 'mixed_prime_window_verification.json'):
    receipt = json.loads((HERE/filename).read_text())
    assert receipt['status'] == 'pass'
    for file, expected in receipt['files'].items():
        assert sha(HERE/file) == expected, ('retained source changed', file)

ledger_path = HERE/'source_ledger.json'
ledger = json.loads(ledger_path.read_text())
url = 'https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf'
for file, review in [
('sources/filaseta-rocks-2016-published.pdf',
 'Author-hosted8-page journal-typeset Online First copy of the2016 Colloquium article, not claimed to be the final paginated file. Lemma2.3 single-edge statement/proof matched and rederived with coordinate j, including strict off-edge inequalities. Current note computes its full residual. Source typo saying n has leading digit mk is not used; n=k*p^r has digit k, while m*n has digit mk. Number-theoretic density proof remains mapped, not independently reconstructed.'),
('sources/filaseta-rocks-2016-published.txt',
 'Local layout text extraction with pypdf '+pypdf.__version__+'; PDF is authoritative. Used for the literal Lemma2.3 statement and proof, with all needed identities independently rederived.')]:
    if file not in {row['file'] for row in ledger['cached_primary_sources']}:
        ledger['cached_primary_sources'].append(dict(file=file, url=url, sha256=sha(HERE/file), review=review))
if url not in ledger['additional_primary_references']:
    ledger['additional_primary_references'].append(url)
ledger['published_input_pinpoints']['FilasetaRocks2016Mixed'] = 'Journal Online First Lemma2.3: N=k*p^r and mk<p give a single Newton edge. Current proof rederives the strict coefficient inequalities and computes its residual f_k,m(epsilon*T). Factorial prime-power insertion then handles d=m*p^a, Gamma integer shapes<=p^a. No all-index or priority extrapolation from the source.'
ledger_path.write_text(json.dumps(ledger, indent=2)+'\n')

proof = BASE+'MIXED_PRIME_RESIDUAL_WINDOWS.tex'
checker = BASE+'verify_mixed_prime_windows.py'
audit = dict(date='2026-11-16', previous_goal_turn='progress', current_goal_turn='progress',
    previous_evidence='Every recorded source hash in all six earlier pass receipts matches, including the4031-check all-index prime-power packet.',
    canonical_proof=proof,
    progress=['exact all-N,m,a prime-power insertion preserves every first Newton edge and residual factor profile',
              'small residual kernel proves an unbounded mixed-exponent/indexed SFC family for every a,r>=1',
              'square-divisor degree sets combine information across primes even when neither residual test suffices alone',
              'literal exponent-six terminal7/22 certificates require combining primes2/3',
              'terminal21/23 remain unresolved under just2/3 and are separately certified using prime5',
              'ordinary mixed-pair complete intersections retain length N(N-1) and sharp all-ring exponent2N-2',
              'positive rational Gamma shapes allowed for exact finite-input local certificates, with half-shape original Gaussian controls',
              'journal Online First primary source located; needed single-edge proof checked and full residual rederived'],
    boundaries=['indexed windows at N=k*p^r are not all starting indices',
                'a squarefree base polynomial alone cannot replace a squarefree base residual',
                'nontrivial degree intersections and repeated finite-field residuals are UNKNOWN, not SFC counterexamples',
                'no uniform useful-prime bound or all-index exponent-six conclusion follows from the finite probe',
                'universal insertion Gamma range alpha<=p^a is distinct from positive rational finite-input shapes',
                'full FC>=3, unrestricted SFC and remaining source audits remain open'],
    independent_replay=False, formal_verification=False, external_review=False,
    novelty_claim=False, whole_goal_complete=False)
(HERE/'mixed_prime_window_audit.json').write_text(json.dumps(audit, indent=2)+'\n')

common = dict(kind='theorem', state='proved', canonical_source=proof, checker=checker,
    proof_type='hybrid', independent_replay=False, formal_verification=False, external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_mixed_prime_windows.py'),
    software_lock=[proof, checker, BASE+'mixed_prime_windows.py', BASE+'prime_power_windows.py',
                   BASE+'moment_bridges.py', BASE+'mixed_prime_window_verification.json',
                   BASE+'mixed_prime_window_examples.json', BASE+'mixed_prime_d6_probe.json',
                   BASE+'mixed_prime_window_audit.json', BASE+'source_ledger.json'],
    supersedes=[], closes_problems=[], narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[], invalidates_assumptions=[], replaced_by=[], priority='reference', programme_status='archived')
new = [
dict(common, id='FIRNEWTONLIFT20261116',
     title='Prime-power insertion preserves the entire first Newton residual profile of factorial moment polynomials',
     scope='For EVERY N,m,a>=1 and prime p, d=m*p^a, b=m*(p^a-1)/(p-1), and integer1<=alpha<=p^a, the actual Gamma moment polynomial f_N,d,alpha(X/p^b) has exactly the Newton polygon of f_N,m(X) at p. At an edge of reduced slope h/e its residual equals R_base(gamma^e*T), gamma=(-1)^b mod p. The coefficient multiplier is ((m*p^a*j)!/(p^(b*j)*(m*j)!))*((alpha)_(d*j)/(d*j)!), a p-adic unit reducing to gamma^j. Thus ALL first residual irreducible degrees and multiplicities, including unresolved repeated branches, are preserved. This is all N,m,a, not a finite test. It does NOT infer actual squarefreeness of a lift merely from actual squarefreeness of the base when its residuals are repeated, and does not preserve uncomputed higher Newton data.',
     dependencies=['FIRPADICDEFORM20261116', 'external: Rocks2015 Lemma3.8 factorial-unit splitting, rederived']),
dict(common, id='FIRINDEXEDMIXEDSFC20261116',
     title='One small residual kernel certifies unbounded mixed-prime exponent and indexed two-moment families',
     scope='Fix positive m,k and a prime p>m*k. If f_k,m(T)=sum binom(k,j)*(m*j)!*T^j is SQUAREFREE MOD p, then for EVERY a,r>=1, d=m*p^a, N=k*p^r and integer1<=alpha<=p^a, no two Gamma(alpha) moments at indices N-1,N of any nonzero A+B*Y^d both vanish, for arbitrary COMPLEX A,B. The base polynomial f_N,m has a single Newton edge of reduced slope m*(p^r-1)/(p^r*(p-1)), and its full residual is f_k,m((-1)^(m*(p^r-1)/(p-1))*T) mod p. Insertion preserves its separability. For k1 the condition is automatic whenever p>m; examples d=2*3^a,N=3^r and d=2*7^a,N=2*7^r have both parameters unbounded. The m2,k2,p5 kernel is repeated and gives no decision. These are indexed SFC windows, NOT full SFC at all starting indices for those mixed exponents.',
     dependencies=['FIRNEWTONLIFT20261116', 'FIRPADICDEFORM20261116', 'external: Filaseta--Rocks2016 Online First Lemma2.3 single edge, rederived with residual']),
dict(common, id='FIRMULTIPRIMEWINDOW20261116',
     title='Intersected local residual multiplicity budgets exclude every rational square divisor and certify adjacent moment pairs',
     scope='For ANY rational polynomial P with nonzero constant/leading coefficients, each p-adic edge with reduced slope h/e and residual irreducible factor phi of multiplicity M permits only contributions e*deg(phi)*t, 0<=t<=floor(M/2), to the degree of a square divisor H^2|P over Q_p. The Minkowski sum over all branches is a necessary degree set D_p(P). If a finite-prime intersection is {0}, P is squarefree over Q and C. This can succeed when no individual prime has squarefree residuals: for f_7,6, D2={0,2},D3={0,1}; for f_22,6, D2={0,8},D3={0,1,2,3,4}. These certify the corresponding adjacent windows for EVERY complex coefficient pair. At terminal21/23 primes2/3 leave possible degree2 and are UNKNOWN; prime5 certifies those finite inputs. The criterion also applies to actual moment polynomials for positive RATIONAL Gamma shapes. Nonzero intersections are not counterexamples, and no uniform useful-prime bound or all-index exponent-six theorem is inferred.',
     dependencies=['FIRPADICDEFORM20261116']),
dict(common, id='FIRMIXEDPAIR20261116', kind='corollary',
     title='Certified mixed-prime adjacent pairs retain ordinary coefficient ideals and sharp all-ring nilpotence',
     scope='For EVERY adjacent moment pair at terminal N>=2 certified by FIRINDEXEDMIXEDSFC or FIRMULTIPRIMEWINDOW under its rational Gamma-shape hypotheses, the homogeneous coefficient ideal (mu_(N-1),mu_N) in Q[A,B] is a regular sequence of degrees N-1,N. The ORDINARY quotient has length N*(N-1), highest nonzero degree2N-3 of dimension1, and (A,B)^(2N-2) is contained in the ideal while the preceding power is not. Over EVERY Q-algebra the two zero moments therefore force the ACTUAL phase coefficient ideal to be nilpotent at exponent2N-2 sharply; every original marker vanishes at all later phase powers. Integral Gamma shapes give the identical ordinary even-block Gaussian profile. Two equations over coefficient rings imply nilpotence, not A=B0. This extends the pair-profile argument to the certified mixed/indexed and rational-shape inputs, without an all-index assertion.',
     dependencies=['FIRINDEXEDMIXEDSFC20261116', 'FIRMULTIPRIMEWINDOW20261116', 'FIRPRIMEPAIRI20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRMIXEDEXEC20261116', proof_type='exact_symbolic',
     title='Exact mixed-prime lifting, indexed-family and finite-input degree certificates retain unknown residual branches',
     scope='mixed_prime_windows.py verifies literal insertion identities, accepts an unbounded indexed family only after its small-kernel F_p gcd is1, and computes local square-divisor degree sets from the actual rational moment polynomial at supplied primes. Every residual factorization is reconstructed and each factor is separately checked by Frobenius/Rabin irreducibility remainders. Empty positive-degree intersection exports an adjacent-window certificate and ordinary pair profile; surviving degrees export UNRESOLVED with no counterexample claim. Positive rational Gamma shapes, including half-shapes, are allowed for finite-input certificates; they are not silently included in the universal integer insertion range. Tests retain exact terminal2..32 exponent-six profiles: only21/23 survive supplied primes2/3, and separate prime5 inputs close those particular windows. This finite table is not an all-index theorem. Original Gamma/Gaussian and ordinary Groebner checks support the written unbounded formulas.',
     dependencies=['FIRNEWTONLIFT20261116', 'FIRINDEXEDMIXEDSFC20261116', 'FIRMULTIPRIMEWINDOW20261116', 'FIRMIXEDPAIR20261116'])]
new_ids = [row['id'] for row in new]


def open_scope(row):
    addition = ' FIRNEWTONLIFT20261116/FIRINDEXEDMIXEDSFC20261116 now give exact prime-power insertion and unbounded MIXED-exponent indexed windows: d=m*p^a, terminal N=k*p^r, all a,r>=1, under p>mk and a squarefree finite kernel. They do not close ALL starts for d6. FIRMULTIPRIMEWINDOW20261116 combines residual square-divisor degrees across primes, certifying f_7,6 and f_22,6 where neither prime alone suffices. Literal terminal21/23 survive primes2/3 but are separately certified at prime5; no finite-table extrapolation. FIRMIXEDPAIR20261116 retains ordinary length N(N-1) and sharp all-ring nilpotence exponent2N-2. Rational Gamma finite-input certificates are separate from universal integer insertion shapes. General FC>=3, unrestricted and full exponent-six SFC, higher residual data and remaining bibliography/audits stay open.'
    if addition not in row['scope']:
        row['scope'] += addition
    row['dependencies'] = list(dict.fromkeys(row['dependencies']+new_ids))
    return row


update_array(ROOT/'MATH_STATUS.json', 'entries', {'OP-FACTORIAL-GENERAL20261116': open_scope}, new)
known = {row['id']: row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger = lambda id: hashlib.sha256('\0'.join(known[id][key] for key in ('state', 'title', 'scope')).encode()).hexdigest()[:16]


def refresh(row):
    row['claims'] = list(dict.fromkeys(row['claims']+new_ids))
    row['reviewed_claims'] = {id: finger(id) for id in row['claims']}
    row['sources'] = list(dict.fromkeys(row['sources']+[proof, BASE+'mixed_prime_window_audit.json']))
    return row


def work(row):
    row = refresh(row)
    row['next_step'] = 'For full exponent-six/all-starting-index SFC, seek higher Newton/residual information or a structural bound on common square-divisor degrees across primes. The exact prime-power insertion and small-kernel indexed families are closed; more successful finite inputs alone will not close the remaining all-index statement. Preserve unresolved first-order branches, including terminal21/23 under just primes2/3. General two-term/unrestricted SFC, higher-dimensional FC and remaining EIT/source review remain the full goal.'
    row['prerequisites'] = 'Read MIXED_PRIME_RESIDUAL_WINDOWS.tex. First residual data is preserved under p-power insertion, not arbitrary higher data or actual squarefreeness when residuals repeat. Kernel p>mk and separability are exact conditions; indexed terminal kp^r is not every start. Degree intersections combine one fixed global square-factor degree, and a nontrivial intersection is UNKNOWN. Distinguish rational Gamma finite-input certificates from universal integer shape<=p^a. Keep ordinary coefficient ideals and the sharp pair nilpotence exponent2N-2.'
    return row


def lesson(row):
    row = refresh(row)
    extra = ' For mixed exponents d=m*p^a, exact factorial unit scaling preserves the entire first Newton profile of the cofactor moment polynomial. A squarefree small kernel at p>mk certifies all a,r>=1 at terminal kp^r. When individual prime residuals repeat, intersect the necessary square-divisor degree sets from their residual multiplicities; this can exclude every global square divisor even when neither prime alone succeeds. Retain ordinary pair quotients and rational Gamma finite-input boundaries.'
    if extra not in row['use']:
        row['use'] += extra
    avoid = ' Do not infer all starting indices from a kp^r family or finite exponent-six table, actual squarefreeness preservation from a repeated first residual, or a counterexample from a surviving degree intersection. A third successful prime at a supplied index is not a uniform bound on future useful primes.'
    if avoid not in row['avoid']:
        row['avoid'] += avoid
    row['boundary'] = 'Prime-power all-index SFC remains closed. Mixed exponents now have unbounded indexed families and multiple-prime finite-input certificates, with ordinary pair profiles. Full exponent-six/all-starting-index and unrestricted SFC, arbitrary FC>=3, higher residual data, arbitrary Gamma shapes/markers and remaining source audits stay open; assurance/novelty flags remain false.'
    row['implementation'] = list(dict.fromkeys(row['implementation']+[BASE+'mixed_prime_windows.py', checker]))
    return row


update_array(ROOT/'knowledge/work_items.json', 'actions', {'WORK-FACTORIAL-GENERAL': work})
update_array(ROOT/'knowledge/lessons.json', 'lessons', {'METHOD-FACTORIAL-RADIAL-PURITY': lesson})
print(json.dumps({'registered': new_ids, 'previous_goal_turn': 'progress', 'whole_goal_complete': False}))
