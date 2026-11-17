"""Register the all-index prime-power two-term SFC progress narrowly."""
from pathlib import Path
import hashlib
import json
import os

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'factorial/image-rigidity/'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()

# Historical continuation must not roll back the later mixed-prime source/scope review.
if any(row['id']=='FIRNEWTONLIFT20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Later mixed-prime progress is registered; use register_mixed_prime_windows.py.')


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
            assert known[row['id']] == row, ('existing prime-power record differs', row['id'])
        else:
            new.append(row)
    if new:
        head, end = text.rsplit('\n  ]\n}', 1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(encode(row) for row in new)+'\n  ]\n}'+end
    if text.encode() != before:
        assert path.read_bytes() == before, 'Concurrent edit'
        temporary = path.with_suffix(path.suffix+'.prime-window-research/.scratch')
        temporary.write_bytes(text.encode())
        os.replace(temporary, path)


for filename in ('verification.json', 'long_audit_verification.json',
                 'simplex_pencil_verification.json', 'two_angular_verification.json',
                 'moment_window_verification.json', 'prime_power_window_verification.json'):
    receipt = json.loads((HERE/filename).read_text())
    assert receipt['status'] == 'pass'
    for file, expected in receipt['files'].items():
        assert sha(HERE/file) == expected, ('retained proof/control source changed', file)

ledger_path = HERE/'source_ledger.json'
ledger = json.loads(ledger_path.read_text())
ledger['published_input_pinpoints']['Rocks2015Newton'] = 'Author-preprint Lemma3.8 supplies the exact factorial unit congruence; Proposition3.22 gives the prime-power scaled Newton vertices. Both are rederived. Current proof computes every edge residual as a reversed truncated exponential and proves squarefreeness, extending the internally checked constant-plus-square case to every prime-power exponent. Full2017 text/priority and all unrelated two-monomial claims are not audited.'
for row in ledger['cached_primary_sources']:
    if row['file'] == 'sources/rocks-2015-author-preprint.pdf':
        addition = ' The prime-power factorial-unit and Newton-vertex inputs in Lemma3.8/Proposition3.22 are now also checked and rederived; the new note computes all edge residuals and independently proves the all-index prime-power two-term conclusion. This is not a claim that the unavailable final2017 article has the same scope.'
        if addition not in row['review']:
            row['review'] += addition
for url in ('https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf',):
    if url not in ledger['additional_primary_references']:
        ledger['additional_primary_references'].append(url)
ledger_path.write_text(json.dumps(ledger, indent=2)+'\n')

proof = BASE+'PRIME_POWER_TWO_TERM_SFC.tex'
checker = BASE+'verify_prime_power_windows.py'
audit = dict(date='2026-11-16', previous_goal_turn='progress', current_goal_turn='progress',
    previous_evidence='Every file hash in all five previous pass receipts matches the current retained sources, including the592-check finite-window audit.',
    canonical_proof=proof,
    progress=['weighted initial-product lemma: squarefree edge residuals exclude all square divisors over Q_p',
              'all-index falling-factorial unit-deformation theorem with geometric residue units',
              'exact digit vertices and every edge residual computed as a reversed truncated exponential of degree<p',
              'constant-plus-prime-power factorial SFC for every prime, exponent and positive starting index',
              'same two-moment theorem for integer Gamma shapes1..d and a specified rational low-degree marker class',
              'ordinary consecutive-pair ideal is a complete intersection of length n(n+1), with sharp coefficient nilpotence bound2n over every Q-algebra',
              'independent-linear adjunction produces all-dimensional SFC classes',
              'exact finite coefficient/hull/residual exporter and separate rational gcd/Bezout/countercontrols'],
    source_credit='Rocks2015 author-preprint factorial-unit and Newton-vertex inputs; all edge-residual and squarefree steps are written here. No priority assertion or reliance on the unavailable final2017 scope.',
    boundaries=['geometric-unit deformation is stated over Q_p, not arbitrary ramified coefficient fields',
                'test residual R(T), not its inseparable substitution R(Z^(p^i))',
                'prime-power exponent and integer Gamma shape<=d required by this theorem',
                'marker degree<d, nonzero mass and p-integral normalized rational divided-power coefficients required',
                'failed residual or outside marker test makes no opposite conjecture decision',
                'powers of the entire inhomogeneous phase do not inherit SFC windows automatically',
                'general non-prime-power two-term SFC, unrestricted SFC and general FC>=3 remain open'],
    independent_replay=False, formal_verification=False, external_review=False,
    novelty_claim=False, whole_goal_complete=False)
(HERE/'prime_power_window_audit.json').write_text(json.dumps(audit, indent=2)+'\n')

common = dict(kind='theorem', state='proved', canonical_source=proof, checker=checker,
    proof_type='hybrid', independent_replay=False, formal_verification=False, external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_prime_power_windows.py'),
    software_lock=[proof, checker, BASE+'prime_power_windows.py', BASE+'moment_bridges.py',
                   BASE+'prime_power_window_verification.json', BASE+'prime_power_window_examples.json',
                   BASE+'prime_power_window_audit.json', BASE+'source_ledger.json'],
    supersedes=[], closes_problems=[], narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[], invalidates_assumptions=[], replaced_by=[], priority='reference', programme_status='archived')
new = [
dict(common, id='FIRPADICDEFORM20261116',
     title='Every falling-factorial polynomial with geometric p-adic unit residues is squarefree in every degree',
     scope='For EVERY prime p and n>=1, if u_0,...,u_n are units in Z_p with u_k/u0 congruent to gamma^k mod p for some nonzero gamma in F_p, then P_n(X)=sum_(k0..n) n!/(n-k)!*u_k*X^k is squarefree over Q_p. For base-p digit a_i of n, its Newton edge has width a_i*p^i and height a_i*(p^i-1)/(p-1). Every lattice-point residual is H_(a_i)(epsilon_i*T)=sum a_i!/(a_i-j)!*(epsilon_i*T)^j, epsilon_i=(-1)^h_i*gamma^(p^i). This is a reversed truncated exponential of degree a_i<p and is squarefree. A written weighted-initial-product argument shows a repeated Q_p factor would give a repeated positive-degree residual factor, including wild slopes. The coefficient valuation group Z is retained; no blanket ramified-field extension is asserted. R(T) is tested before its inseparable R(Z^(p^i)) substitution. No finite degree scan proves the universal claim.',
     dependencies=['external: Rocks2015 author preprint Proposition3.22 Newton vertices, rederived here']),
dict(common, id='FIRSFCMONOPRIME20261116',
     title='Constant-plus-prime-power monomials satisfy Strong Factorial at every starting index, with Gamma-shape extensions',
     scope='For EVERY prime p, a>=1, d=p^a and arbitrary complex A,B with A+B*Y^d nonzero, no two consecutive positive factorial moments both vanish. More generally this exact two-moment conclusion holds under Gamma(alpha,1) for every INTEGER1<=alpha<=d. The parameter polynomials f_n,d,alpha=sum binom(n,k)*(alpha)_(d*k)*x^k become geometric-unit falling-factorial polynomials after x=X/p^b, b=(d-1)/(p-1). The Gamma multiplier product_(j1..alpha-1)(1+d*k/j) is a p-adic unit congruent1. Hence every parameter polynomial is squarefree, and n*f_(n-1)=n*f_n-x*f_nprime forbids adjacent common roots over C. Zero A or B is handled directly. This is all n,p,a, not fixed-d density or a bounded gcd range. Adding t independent nonzero linear factorial variables gives an all-dimensional window t+2 equal to the monomial count. d1 is already in the linear theorem. Non-prime-power d such as6, Gamma alpha>d, arbitrary powers of this inhomogeneous phase and unrestricted SFC remain outside this conclusion. Rocks inputs credited, no novelty or external assurance claim.',
     dependencies=['FIRPADICDEFORM20261116', 'FIRSFCAFFINE20261116', 'external: Rocks2015 Lemma3.8 factorial-unit congruence, rederived here']),
dict(common, id='FIRPRIMEMARKERS20261116',
     title='Specified low-degree original markers retain a two-moment window for constant-plus-prime-power phases',
     scope='For d=p^a and an exponential Y, let the RATIONAL marker be M(Y)=sum_(e0..d-1) a_e*Y^e/e!, with mass C=sum a_e nonzero and every a_e/C p-integral. For EVERY nonzero phase A+B*Y^d with arbitrary complex A,B and every positive start n, its two consecutive M-marked moments cannot both vanish. Normalized weighted monomial moments equal (d*k)! times a p-integral multiplier congruent1, so the parameter polynomials fall under FIRPADICDEFORM and the adjacent derivative identity. All monomial markers Y^e for e<d are covered after scalar normalization; these include Gamma shapes alpha=e+1<=d. Gaussian averaging transfers the specified markers, with odd coordinate markers separately identically zero. Failed marker criteria, zero mass, arbitrary higher marker degree and p-adic cancellation are outside this sufficient theorem, not zero-period decisions. No all-marker nonvanishing or general Gaussian-purity conclusion.',
     dependencies=['FIRPADICDEFORM20261116', 'FIRSFCMONOPRIME20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRPRIMEPAIRI20261116',
     title='Prime-power consecutive moment pairs have ordinary complete-intersection quotients and sharp all-ring coefficient nilpotence',
     scope='For d=p^a, integer1<=alpha<=d and every n>=1, let mu_k(A,B)=sum binom(k,j)*(alpha)_(d*j)*A^(k-j)*B^j. The ORDINARY pair ideal I_n=(mu_n,mu_(n+1)) in Q[A,B] is a regular sequence of degrees n,n+1, since the all-index two-moment theorem excludes every common projective zero. Its Hilbert series is (1+...+t^(n-1))*(1+...+t^n), length n(n+1), and highest nonzero degree2n-1 has dimension1. Thus (A,B)^(2n) is contained in I_n while (A,B)^(2n-1) is not. Over EVERY commutative Q-algebra, two consecutive zero moments force the ACTUAL coefficient ideal (A,B)^(2n)=0 and the phase power (A+B*Y^d)^(2n)=0, hence all original marked moments vanish for every exponent>=2n. The exponent is universally sharp in Q[A,B]/I_n because the phase-power coefficients generate exactly (A,B)^m. The Gaussian even-block substitution preserves actual coefficient ideal(A,B) and the same ordinary pair profile. This is not just a radical or reduced-point statement, and over coefficient rings the two equations need not force A=B0.',
     dependencies=['FIRSFCMONOPRIME20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRPRIMEWINDOWEXEC20261116', proof_type='exact_symbolic',
     title='Exact prime-power moment, Newton residual and all-index class exporters verify the actual polynomial identities',
     scope='prime_power_windows.py recognizes a prime-power exponent and Gamma integer shape1..d, exports the theorem all-index two-moment certificate, constructs any supplied finite moment polynomial, scales its parameter, computes a literal rational-coefficient Newton lower hull and normalized F_p residuals, and verifies equality with the proved digit formula. It also exports ordinary pair complete-intersection degrees, length and sharp nilpotence bound. A rational divided-power marker recognizer checks degree, normalized coefficient integrality and nonzero mass. Unsupported exponents/shapes are rejected; failed markers make no zero-period decision. Finite controls include generic unit deformations, edge lattice points, separate gcd/Bezout identities, original moments, repeated/wild residual countercontrols and adjunction. Ordinary Groebner controls verify lengths, all degree2n containments and a degree2n-1 survivor. Failed residual squarefreeness is only a failed sufficient test. Universal quantifiers use the written proof, not finite checks.',
     dependencies=['FIRPADICDEFORM20261116', 'FIRSFCMONOPRIME20261116', 'FIRPRIMEMARKERS20261116', 'FIRPRIMEPAIRI20261116'])]
new_ids = [row['id'] for row in new]


def open_scope(row):
    addition = ' FIRPADICDEFORM20261116/FIRSFCMONOPRIME20261116 now close ALL starting indices for constant-plus-prime-power monomials A+B*Y^(p^a): their parameter polynomials are squarefree by explicit Newton residuals, not irreducibility density. Gamma integer shapes1..p^a and the rational low-degree marker class FIRPRIMEMARKERS20261116 are included; independent linear adjunction gives all-dimensional SFC classes. The surviving constant-plus-monomial target has NON-PRIME-POWER exponent (starting d6); arbitrary Gamma shapes, arbitrary markers, powers of the inhomogeneous phase and general two-term/unrestricted SFC remain open. Exact exporters retain coefficient-field and residual-variable boundaries; formal/external assurance and novelty flags remain false.'
    if addition not in row['scope']:
        row['scope'] += addition
    ordinary = ' FIRPRIMEPAIRI20261116 retains the ordinary pair ideal: length n(n+1), highest nonzero degree2n-1, and sharp actual coefficient nilpotence bound2n over every Q-algebra. Two zero moments over rings imply this nilpotence, not coefficientwise zero.'
    if ordinary not in row['scope']:
        row['scope'] += ordinary
    row['dependencies'] = list(dict.fromkeys(row['dependencies']+new_ids))
    return row


update_array(ROOT/'MATH_STATUS.json', 'entries', {'OP-FACTORIAL-GENERAL20261116': open_scope}, new)
known = {row['id']: row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger = lambda id: hashlib.sha256('\0'.join(known[id][key] for key in ('state', 'title', 'scope')).encode()).hexdigest()[:16]


def refresh(row):
    row['claims'] = list(dict.fromkeys(row['claims']+new_ids))
    row['reviewed_claims'] = {id: finger(id) for id in row['claims']}
    row['sources'] = list(dict.fromkeys(row['sources']+[proof, BASE+'prime_power_window_audit.json']))
    return row


def work(row):
    row = refresh(row)
    row['next_step'] = 'Attack constant-plus-monomial SFC for NON-PRIME-POWER exponents d (first d6), or find a structural extension of the falling-factorial geometric-unit theorem that handles multiple prime factors. Do not repeat the now closed all-index prime-power case or substitute fixed-d irreducibility density for all indices. General two-term and unrestricted SFC, arbitrary FC>=3/coupled simplex phases, arbitrary Gamma shapes and the remaining EIT/bibliography audits remain the full goal.'
    row['prerequisites'] = 'Read PRIME_POWER_TWO_TERM_SFC.tex. Every Newton residual is a reversed truncated exponential of degree<p; the coefficient field is Q_p with integral valuation group. Differentiate the residual in T, not its inseparable substitution in Z. Gamma shape<=d, low marker degree and normalized rational coefficient integrality are exact boundaries. Outside a sufficient residual or marker criterion no opposite conclusion is allowed. Previous source-version, affine-window, common-prime, cubic, pencil and binary scopes remain in force.'
    return row


def lesson(row):
    row = refresh(row)
    extra = ' For prime-power two-term SFC use the all-index geometric-unit falling-factorial Newton theorem: edge residues are scaled reversals of truncated exponentials of digit degree<p, hence squarefree. Parameter polynomial squarefreeness plus the exact adjacent derivative identity excludes consecutive moment zeros. This covers Gamma shapes1..d, a stated rational low-degree marker class and independent-linear adjunction.'
    extra += ' Its ordinary two-moment ideal is a complete intersection of length n(n+1), forcing actual coefficient nilpotence exponent2n sharply over all Q-algebras; retain the nonreduced quotient instead of its radical.'
    if extra not in row['use']:
        row['use'] += extra
    avoid = ' Do not test R(Z^(p^i)) instead of the residual R(T), extend the Q_p deformation theorem to arbitrary ramified coefficient fields, replace prime-power d by all d, or infer SFC for powers of an entire inhomogeneous phase from SFC for that phase.'
    if avoid not in row['avoid']:
        row['avoid'] += avoid
    row['boundary'] = 'Constant-plus-prime-power monomial SFC is closed internally at EVERY index, including Gamma integer shapes1..d and the stated markers. The constant-plus-monomial survivor has non-prime-power exponent, starting d6. General two-term/unrestricted SFC, arbitrary FC>=3, arbitrary Gamma shapes/markers, full EIT and bibliography audits remain open. Independent/formal/external assurance and novelty flags remain false.'
    row['implementation'] = list(dict.fromkeys(row['implementation']+[BASE+'prime_power_windows.py', checker]))
    return row


update_array(ROOT/'knowledge/work_items.json', 'actions', {'WORK-FACTORIAL-GENERAL': work})
update_array(ROOT/'knowledge/lessons.json', 'lessons', {'METHOD-FACTORIAL-RADIAL-PURITY': lesson})
print(json.dumps({'registered': new_ids, 'previous_goal_turn': 'progress', 'whole_goal_complete': False}))
