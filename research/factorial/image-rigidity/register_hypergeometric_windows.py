"""Register uniform monomial windows with narrow, guarded ledger edits."""
from pathlib import Path
import hashlib
import json
import os

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'factorial/image-rigidity/'
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
            assert known[row['id']] == row, ('existing hypergeometric record differs', row['id'])
        else:
            new.append(row)
    if new:
        head, end = text.rsplit('\n  ]\n}', 1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(encode(row) for row in new)+'\n  ]\n}'+end
    if text.encode() != before:
        assert path.read_bytes() == before, 'Concurrent edit'
        temporary = path.with_suffix(path.suffix+'.hyper-window-research/.scratch')
        temporary.write_bytes(text.encode())
        os.replace(temporary, path)


if (HERE/'coupled_block_audit.json').exists():
    raise SystemExit('Historical monomial-window registration: use register_coupled_blocks.py to preserve later scopes and reviews.')

for filename in ('verification.json', 'long_audit_verification.json', 'simplex_pencil_verification.json',
                 'two_angular_verification.json', 'moment_window_verification.json',
                 'prime_power_window_verification.json', 'mixed_prime_window_verification.json',
                 'hypergeometric_window_verification.json'):
    receipt = json.loads((HERE/filename).read_text())
    assert receipt['status'] == 'pass'
    for file, expected in receipt['files'].items():
        assert sha(HERE/file) == expected, ('retained source changed', file)

proof = BASE+'HYPERGEOMETRIC_MONOMIAL_WINDOWS.tex'
checker = BASE+'verify_hypergeometric_windows.py'
audit = dict(date='2026-11-16', previous_goal_turn='progress', current_goal_turn='progress',
    previous_evidence='Every source hash in all seven earlier pass receipts matches, including the1154-check mixed-prime packet. Repeated residual branches and the remaining two-moment exponent-six gap are unchanged.',
    canonical_proof=proof,
    progress=['uniform all-starting-index D+1 window for every constant-plus-monomial Gamma phase and all positive real shapes',
              'every polynomial marker gives an exact R(t); zero R decides the whole zero sequence, otherwise window D+deg(R)+1',
              'Dirichlet extension for all positive real shapes: the extra singularity is checked and its indicial exponent N+beta exceeds degree N',
              'original Gaussian polynomial markers transfer by actual Reynolds averaging, including half-shapes and exact zero averages',
              'ordinary window ideal contains an exhibited binary regular subideal; all-ring nilpotence exponent at most2n+D-1',
              'actual expanded support proves full SFC for disjoint block-sum monomial products when C>=D+1, including every exponent with all block sizes>=2',
              'a+b*(Y1+Y2)^6 receives a seven-moment window within its eight actual monomials at every start'],
    boundaries=['D+1 windows do not establish two-monomial SFC when D>1; a+b*Y^6 still requires two instead of seven',
                'the equation controls root multiplicity up to D+deg(R), not simplicity',
                'constant phases use marker mass separately; the one-coordinate Dirichlet point law is excluded',
                'executable shapes are rational although the theorem allows all positive real shapes',
                'ordinary regular-subideal length is an upper bound for the full window quotient; no sharpness asserted',
                'no SFC for arbitrary powers of the whole phase or independent Gaussian linear adjunction is inferred',
                'full FC>=3, unrestricted SFC and remaining literature/source audits stay open'],
    independent_replay=False, formal_verification=False, external_review=False,
    novelty_claim=False, whole_goal_complete=False)
(HERE/'hypergeometric_window_audit.json').write_text(json.dumps(audit, indent=2)+'\n')

ledger_path = HERE/'source_ledger.json'
ledger = json.loads(ledger_path.read_text())
ledger.setdefault('internal_proof_extensions', {})['hypergeometric_monomial_windows'] = dict(
    canonical_source=proof,
    provenance='Derived directly from literal Gamma/Dirichlet moment coefficients, Euler identities and the ordinary binary complete-intersection Hilbert series. No newly cited external theorem or priority claim. Factorial independent-linear adjunction is the previously checked Rocks input.',
    singularity='Dirichlet x*=-D^D/prod(d_i^d_i); exceptional indicial exponent N+A_inactive+E-e+(s-1)/2>N for at least two coordinates.',
    boundary='Uniform D+e+1 window and expanded-support SFC; no two-moment claim for unrestricted mixed exponents.')
ledger_path.write_text(json.dumps(ledger, indent=2)+'\n')

common = dict(kind='theorem', state='proved', canonical_source=proof, checker=checker,
    proof_type='hybrid', independent_replay=False, formal_verification=False, external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_hypergeometric_windows.py'),
    software_lock=[proof, checker, BASE+'hypergeometric_windows.py', BASE+'moment_bridges.py',
                   BASE+'hypergeometric_window_verification.json', BASE+'hypergeometric_window_examples.json',
                   BASE+'hypergeometric_window_audit.json', BASE+'source_ledger.json'],
    supersedes=[], closes_problems=[], narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[], invalidates_assumptions=[], replaced_by=[], priority='reference', programme_status='archived')
new = [
dict(common, id='FIRGAMMAMONOWINDOW20261116',
     title='All positive Gamma shapes and all polynomial markers have exact monomial zero-sequence tests or uniform windows',
     scope='For independent Gamma(alpha_i,1) with EVERY positive REAL shape, d_i>=0, D=sum d_i>=1, Z=prod Y_i^d_i and F=a+bZ, arbitrary complex a and b!=0, a marker M=sum m_e Y^e has R(t)=sum m_e prod(alpha_i+d_i*t)_(e_i). R=0 iff ALL marked moments vanish including index0; otherwise EVERY D+deg(R)+1 consecutive positive indices contain a nonzero moment. In particular unmarked windows have length D+1, with no prime-power or integer-shape restriction. The coefficient differential equation has order D+deg(R)+1 and is ordinary at every nonzero parameter; falling Euler identities turn a zero window into an impossible root multiplicity. Constant phases use marker mass separately. Original Gaussian polynomial markers are first averaged by exact block Reynolds identities; half-shapes, odd and angular zero averages are retained. This is a degree-dependent window, NOT full two-term SFC for general D.',
     dependencies=['FIRGAMMAADAPTER20261116']),
dict(common, id='FIRDIRMONOWINDOW20261116',
     title='The uniform monomial and all-marker windows extend to every positive Dirichlet shape after checking the extra singularity',
     scope='For Dir(alpha) on at least TWO coordinates, EVERY positive REAL shape, Z=prod T_i^d_i, D>=1 and F=a+bZ with complex b!=0, homogenize every polynomial marker to degree E on sum T_i=1. R(t)=sum_|e|=E m_e prod(alpha_i+d_i*t)_(e_i) has degree e<=E or is zero. R=0 iff the whole marked moment sequence including0 is zero; otherwise every D+e+1 positive consecutive moments contain a nonzero value. The order D+e+1 Euler equation has one extra nonzero singularity x*=-D^D/prod d_i^d_i. A root of order at least D+e there would require order N+beta with beta=A_inactive+E-e+(s-1)/2>0, exceeding the polynomial degree N. Thus the singularity is actually handled, not treated as ordinary. The point law and constant phases are separate. Arbitrary coupled simplex phases and unrestricted SFC are not claimed.',
     dependencies=['FIRGAMMAMONOWINDOW20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRMONOWINDOWIDEAL20261116', kind='corollary',
     title='Uniform monomial windows retain ordinary coefficient ideals and an explicit all-ring nilpotence bound',
     scope='For either monomial law above with positive RATIONAL shapes, the unmarked ordinary window ideal I=(mu_n,...,mu_(n+D)) in Q[U,V] contains a regular sequence of degrees n,n+D. An actual coprime second generator G_t=sum_j=1..D t^(j-1)*U^(D-j)*mu_(n+j) exists for some integer0<=t<=n(D-1). Therefore (U,V)^(2n+D-1) is contained in I and the FULL window quotient has length at most n(n+D). The selected regular SUBIDEAL has exactly that CI length, not necessarily the full window ideal. Over EVERY Q-algebra, those zero moments force the ACTUAL phase coefficient ideal (a,b) nilpotent at exponent at most2n+D-1; phase powers and every original marker vanish thereafter. Dirichlet coordinate elimination and Gaussian even-block lift preserve the actual ideal. No radical replacement or sharpness of this bound is asserted.',
     dependencies=['FIRGAMMAMONOWINDOW20261116', 'FIRDIRMONOWINDOW20261116']),
dict(common, id='FIRBLOCKSUMSFC20261116', kind='corollary',
     title='Actual expanded block-sum support gives all-index Strong Factorial classes for arbitrary exponents',
     scope='For disjoint nonempty blocks B_i of factorial variables, L_i=sum_(j in B_i) Y_j, r_i=|B_i| and arbitrary positive integer exponents d_i, F=a+b*prod L_i^d_i has exactly C=1+prod binom(r_i+d_i-1,d_i) actual monomials when a,b!=0. If C>=1+sum d_i, the uniform Gamma window proves FULL SFC for this class at EVERY start. In particular every exponent is allowed when all r_i>=2; a+b*(Y1+Y2)^6 has eight monomials and a seven-moment window. Disjoint factorial linear adjunctions preserve the support/window inequality, with degenerate a0/b0 bases handled directly. Replacing block sums by half squared Gaussian norms gives the same base support-count Gaussian analogue, with every original marker controlled by the preceding theorem. No independent Gaussian-linear adjunction or arbitrary powers of the whole phase is claimed. The two-monomial mixed-exponent SFC gap remains.',
     dependencies=['FIRGAMMAMONOWINDOW20261116', 'FIRSFCAFFINE20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRHYPERWINDOWEXEC20261116', proof_type='exact_symbolic',
     title='Exact hypergeometric marker and ordinary-ideal exporter checks the simplex singular point and preserves support boundaries',
     scope='hypergeometric_windows.py accepts positive rational shapes and exact algebraic polynomial markers. It computes literal Gamma/Dirichlet Q,S,R, either an exact whole zero-sequence identity or the all-starting-index D+deg(R)+1 window, and the positive extra-singularity beta for Dirichlet. An unmarked exporter exhibits the finite-selection coprime regular subideal and returns bounds for the ordinary window ideal without asserting it is a complete intersection or sharp. A block-sum exporter certifies SFC only when actual expanded support covers the proved window; one-variable exponent6 returns insufficient budget, not a counterexample. Controls compare direct original Gamma/Dirichlet/Gaussian moments, independently expand ordinary differential coefficients and the indicial leading term, retain short initial zero prefixes and constant zero-mass exceptions, and verify ordinary Groebner inclusions. The written proof, not finite instances, establishes the universal statements.',
     dependencies=['FIRGAMMAMONOWINDOW20261116', 'FIRDIRMONOWINDOW20261116', 'FIRMONOWINDOWIDEAL20261116', 'FIRBLOCKSUMSFC20261116'])]
new_ids = [row['id'] for row in new]


def open_scope(row):
    addition = ' FIRGAMMAMONOWINDOW20261116/FIRDIRMONOWINDOW20261116 now give EVERY-start D+1 unmarked windows for constant-plus-monomial phases and EVERY positive REAL Gamma/Dirichlet shape, plus exact all-marker zero-sequence tests or D+deg(R)+1 windows. The simplex extra singularity is checked by indicial exponent N+beta>N. FIRMONOWINDOWIDEAL20261116 retains ordinary all-ring nilpotence bound2n+D-1 and an actual regular SUBIDEAL, without claiming the full window ideal is CI or sharp. FIRBLOCKSUMSFC20261116 proves full SFC for expanded disjoint block-sum products whose actual support covers D+1, including all exponents when every block has >=2 variables. For a+b*Y^6 the general window is SEVEN, not the required TWO; unrestricted and full two-term mixed-exponent SFC remain open. Executable shapes are rational, theorem shapes are all positive real.'
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
    row['sources'] = list(dict.fromkeys(row['sources']+[proof, BASE+'hypergeometric_window_audit.json']))
    return row


def work(row):
    row = refresh(row)
    row['next_step'] = 'Uniform all-start D+1 monomial windows, all positive shapes, all-marker classification, simplex singularity, ordinary ideal bounds and support-qualified block-sum SFC are closed by HYPERGEOMETRIC_MONOMIAL_WINDOWS.tex. To reach full two-term SFC, improve root multiplicity from <=D to1 for unmarked mixed exponents (first d6), using higher Newton data or a stronger structural property; the present ODE does not do that. Preserve repeated residual branches and indexed-family boundaries. General coupled simplex phases, FC>=3, unrestricted SFC and remaining EIT/source review stay in the full goal.'
    row['prerequisites'] = 'Read the hypergeometric proof and earlier mixed-prime proof. The Gamma ODE is ordinary away0; Dirichlet has an extra singular point whose exceptional indicial exponent exceeds polynomial degree. Markers use actual R(t), with Dirichlet homogenization and original Gaussian Reynolds averaging. Constant phases use marker mass, not the nonconstant criterion. D+1 is a weaker window than two-term SFC; expanded support must actually cover it. Ordinary regular-subideal length is an upper bound for the full window ideal, not an equality or sharpness statement.'
    return row


def lesson(row):
    row = refresh(row)
    extra = ' For constant-plus-monomial phases, derive the exact Euler differential equation from c_(k+1)/c_k=Q(k)/S(k). Every marker is encoded by R(k), so an ordinary-point multiplicity bound yields a uniform D+deg(R)+1 window for every positive Gamma shape. For Dirichlet, homogenize first and check the additional singularity: its indicial exponent N+A_inactive+E-deg(R)+(s-1)/2>N rules out a high-order polynomial root. This gives an all-shape result without extending p-adic separability beyond its range. Full SFC follows for expanded block sums only after counting actual support.'
    if extra not in row['use']:
        row['use'] += extra
    avoid = ' Do not equate a degree-dependent window with a support-size SFC window, assume an order-D+1 ODE implies simple roots, ignore the extra simplex singularity, or report a selected regular subideal Hilbert series as the full window quotient. A marker can have a short identically zero prefix while R is nonzero; constant zero-mass phases and original Gaussian zero averages are separate cases.'
    if avoid not in row['avoid']:
        row['avoid'] += avoid
    row['boundary'] = 'All positive real Gamma/Dirichlet shapes and arbitrary polynomial markers are now handled for NONCONSTANT constant-plus-monomial phases by exact R and uniform windows. Actual support-qualified block sums satisfy SFC at every start. The general monomial window D+1 does not close two-term SFC for mixed exponents, including d6; general coupled simplex phases, FC>=3, unrestricted SFC and remaining source audits stay open. All novelty/outside/formal assurance flags remain false.'
    row['implementation'] = list(dict.fromkeys(row['implementation']+[BASE+'hypergeometric_windows.py', checker]))
    return row


update_array(ROOT/'knowledge/work_items.json', 'actions', {'WORK-FACTORIAL-GENERAL': work})
update_array(ROOT/'knowledge/lessons.json', 'lessons', {'METHOD-FACTORIAL-RADIAL-PURITY': lesson})
print(json.dumps({'registered': new_ids, 'previous_goal_turn': 'progress', 'whole_goal_complete': False}))
