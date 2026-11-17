"""Register the all-shape block limit and coupled angular product family."""
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
            old = known[row['id']]
            if old != row:
                needle = encode(old)
                assert text.count(needle) == 1, ('record formatting changed', row['id'])
                text = text.replace(needle, encode(row), 1)
                known[row['id']] = row
        else:
            new.append(row)
    if new:
        head, end = text.rsplit('\n  ]\n}', 1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(encode(row) for row in new)+'\n  ]\n}'+end
    if text.encode() != before:
        assert path.read_bytes() == before, 'Concurrent edit'
        temporary = path.with_suffix(path.suffix+'.coupled-block-research/.scratch')
        temporary.write_bytes(text.encode())
        os.replace(temporary, path)


if any(row['id']=='FIRREALROOTPARITY20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Historical coupled-block registration: use register_real_place_windows.py to preserve the later SFC scopes and reviews.')

for filename in ('verification.json', 'long_audit_verification.json', 'simplex_pencil_verification.json',
                 'two_angular_verification.json', 'moment_window_verification.json',
                 'prime_power_window_verification.json', 'mixed_prime_window_verification.json',
                 'hypergeometric_window_verification.json', 'coupled_block_verification.json'):
    receipt = json.loads((HERE/filename).read_text())
    assert receipt['status'] == 'pass'
    for file, expected in receipt['files'].items():
        assert sha(HERE/file) == expected, ('retained source changed', file)

proof = BASE+'COUPLED_BLOCK_PERIODS.tex'
checker = BASE+'verify_coupled_block_periods.py'
audit = dict(date='2026-11-16', previous_goal_turn='progress', current_goal_turn='progress',
    previous_evidence='Every recorded source pin in all eight earlier pass receipts matches, including the577-check uniform Gamma/Dirichlet monomial-window packet. Its weaker multiplicity bound and unresolved two-term exponent-six SFC boundary are retained.',
    canonical_proof=proof,
    progress=['complete pinned EIT VI block-radial limit and algebraic block-separable product step audited',
              'block limit extended to all positive real Gamma shapes, including block totals below1, and every polynomial marker',
              'inactive variables retain fixed Gamma laws; compact boundary and unbounded tail estimates justify arbitrary strictly lower-degree active-inactive couplings',
              'markers use active degree and exact inactive Gamma expectations, rather than setting inactive marker variables to zero',
              'singleton blocks give eventual nonzero for every single nonconstant top monomial, allowing zero exponents and arbitrary complex lower couplings in additional variables; X^3+Q<=2 has limit exp(Q2(3,0,0)/27)',
              'coupled product-uniform angular period equals J_(r-1)(xi)/xi for every r, with actual E-function system and a rational-residue functional independence proof',
              'ordinary-point Beukers lifting makes the coupled unit period transcendental/nonzero at every nonzero algebraic scale',
              'every algebraic polynomial cube marker has an exact Laurent reduction to1,exp,J1,...,J_(r-1), giving a complete algebraic-scale single-period zero test',
              'new genuinely coupled block families have eventual nonzero factorial moments at algebraic scales and FC/infinitely many nonzero moments at every complex scale',
              'power closure and ordinary universal tail unit ideals retain all-Q-algebra no-tail and zero scalar-normalization conclusions',
              'original Gaussian markers, zero averages, zero leading-period controls and exact source normalizations retained'],
    boundaries=['coupled products on products of intervals are not arbitrary coupled phases on a triangle',
                'a zero leading asymptotic coefficient is not a whole marked zero sequence',
                'algebraic-scale eventual nonzero is distinct from complex-scale FC/infinitely many moments',
                'single-period arithmetic exporter requires algebraic marker coefficients and exact algebraic scales',
                'ordinary tail unit ideals give existence of finite Bezout identities, not a supplied universal cutoff or SFC window',
                'the block-limit and product-step audit does not complete the remaining EIT V geometric/elliptic or higher odd-degree inputs',
                'full FC>=3, unrestricted SFC and full two-term exponent-six SFC remain open'],
    independent_replay=False, formal_verification=False, external_review=False,
    novelty_claim=False, inactive_variable_extension_complete=True, general_FC_and_SFC_closed=False)
(HERE/'coupled_block_audit.json').write_text(json.dumps(audit, indent=2)+'\n')

ledger_path = HERE/'source_ledger.json'
ledger = json.loads(ledger_path.read_text())
for row in ledger['cached_primary_sources']:
    if row['file'] == 'sources/long-eit-vi-current.tex':
        row['review'] = 'Radial Theorem2.1 and the COMPLETE block-radial theorem proof, including compact coordinate-boundary and unbounded-tail estimates, now internally checked. Its block-separable angular substitution and product step are checked; usable nonzero inputs here are the independently checked affine and cubic classes. The block limit is extended to EVERY positive real Gamma shape and all polynomial markers, and the new coupled product-uniform angular family is proved independently by explicit E-functions. Odd-family general geometric inputs, original elliptic continued-density detector and complete finite gcd range remain outside the completed audit. No blanket audit of all source angular inputs is claimed.'
ledger['published_input_pinpoints']['LongEITVIBlock'] = 'Pinned EIT VI section Block-radial limits and nonflat families: entire Block-radial limit theorem and the algebraic product step of Several active angular blocks, with original coordinate-volume constants checked. Current proof adds positive real Gamma shapes, active-degree polynomial markers, inactive variables with arbitrary strictly lower-degree couplings and genuinely coupled product-uniform angular periods. Remaining odd-degree/elliptic geometric inputs are distinct.'
ledger['published_input_pinpoints']['BeukersCoupledProduct'] = 'Theorem3.2 and Lemma3.1 applied to the explicit actual E-function system1,e^z,J1,...,J_(r-1), whose sole finite system pole is0. Functional linear independence is proved by induction using the impossibility of a rational derivative1/z. Thus all nonzero algebraic points give independent values, enabling a complete exact algebraic-marker single-period test.'
ledger.setdefault('internal_proof_extensions', {})['coupled_block_products'] = dict(
    canonical_source=proof,
    inactive_variable_extension='Any finite tuple of fixed positive-shape Gamma variables may occur in strictly lower total-degree terms. Their next-layer couplings disappear at Z=0, while active-degree marker amplitudes retain exact inactive moments. The unique scaled maximum is at active d and inactive0, with explicit compact and tail estimates.',
    provenance='Explicit product-uniform period J_(r-1)(z)/z, exact Laurent reduction of every polynomial marker, and coupled block FC families derived here; credited block-limit and arithmetic lifting inputs retained. No novelty claim.',
    boundary='Product-of-interval coupled phase, algebraic-scale eventual nonzero, complex-scale FC/power closure. Zero leading periods and original Gaussian zero averages remain distinct.')
ledger_path.write_text(json.dumps(ledger, indent=2)+'\n')

common = dict(kind='theorem', state='proved', canonical_source=proof, checker=checker,
    proof_type='hybrid', independent_replay=False, formal_verification=False, external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_coupled_block_periods.py'),
    software_lock=[proof, checker, BASE+'coupled_block_periods.py', BASE+'moment_bridges.py',
                   BASE+'coupled_block_verification.json', BASE+'coupled_block_examples.json',
                   BASE+'coupled_block_audit.json', BASE+'source_ledger.json',
                   'jacobian/cross-programme-connections/check_inactive_variable_factorial_extension_20261116.py',
                   'jacobian/cross-programme-connections/check_inactive_variable_factorial_extension_20261116.json',
                   BASE+'sources/long-eit-vi-current.tex', BASE+'sources/beukers-math0405549.html'],
    supersedes=[], closes_problems=[], narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[], invalidates_assumptions=[], replaced_by=[], priority='reference', programme_status='archived')
new = [
dict(common, id='FIRBLOCKRADIAL20261116',
     title='The full block-radial limit allows inactive variables, every positive Gamma shape and every polynomial marker',
     scope='Partition a nonempty subset of independent positive-shape Gamma variables into active blocks with shape sums A_j, F0=prod S_j^d_j, d_j>=1, D=sum d_j, Cd=prod d_j^d_j. All remaining variables Z have unchanged independent fixed Gamma shapes beta_l>0. For F=cF0+G_(D-1)(X,Z)+R_(<=D-2), c!=0, arbitrary COMPLEX lower coefficients, and a polynomial marker H of ACTIVE degree e, the normalized moment E[H F^m]/(m^e*c^m*prod(A_j)_(d_j*m)) tends to E_(T,Z)[H_[e](dT,Z)*exp(G(dT,0)/(c*Cd))]. EVERY positive REAL shape is allowed, including shapes below1. The proof controls the scaled inactive boundary at v0 and unbounded tails without dividing by F0 there; Stirling normalization includes m^(-r/2-B) for total inactive shape B. A nonzero leading period gives eventual marked nonzero; a zero limit makes no whole-sequence decision. Setting inactive variables to zero is valid for G in the exponent, but marker amplitudes require their actual unchanged Gamma expectations. The complete pinned EIT VI all-active block-limit proof and product step retain their coordinate-volume constants. Remaining geometric/elliptic angular inputs are not thereby audited.',
     dependencies=['FIRGAMMAADAPTER20261116', 'external: Long--Le Blanc EIT VI pinned block-radial theorem, fully checked and extended']),
dict(common, id='FIRMONOMIALTOPFC20261116', kind='corollary',
     title='Every single nonconstant top monomial has eventual nonzero Gamma moments even with zero exponents and arbitrary lower couplings',
     scope='For independent Gamma(alpha_i,1) of EVERY positive REAL shape, any exponent tuple a_i>=0 with sum D>=1, and c!=0, every polynomial whose highest ordinary homogeneous part is the SINGLE monomial c*Y^a has ALL sufficiently large pure moments nonzero, with arbitrary COMPLEX strictly lower-degree coefficients in ALL variables. Let I={i:a_i>0}, Ca=prod_(i in I) a_i^a_i, and G=F_(D-1). Then E[F^m]/(c^m*prod_(i in I)(alpha_i)_(a_i*m)) tends to exp(G(a)/(c*Ca)), always nonzero, with zero coordinates in a kept literally. Variables absent from the leading monomial can occur freely in lower terms. For F=X^3+Q<=2(X,Y,Z), L(F^m)/(3m)! tends to exp(Q2(3,0,0)/27); for X^2Y+Q<=2(X,Y,Z,W), the denominator is (2m)!m! and exponent Q2(2,1,0,0)/4. For a marker of active degree e, multiply the unit limit by the actual inactive Gamma expectation of H_[e](a_I,Z). A zero marker amplitude is no zero-sequence decision: Z-1 has zero limit but E[(Z-1)*(X^3+XZ)]=1 for factorial shapes. Gaussian orthogonal block lifts retain all original averaged markers. No algebraic-scale restriction or short all-start SFC window is asserted.',
     dependencies=['FIRBLOCKRADIAL20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRCOUPLEDPERIOD20261116',
     title='Product-uniform coupled angular periods are transcendental and nonzero at every nonzero algebraic scale in every dimension',
     scope='For EVERY r>=1, independent uniforms U_j and nonzero algebraic xi, integral_[0,1]^r exp(xi*prod U_j)=J_(r-1)(xi)/xi is transcendental and nonzero, with J_k(z)=sum_n>=1 z^n/(n!*n^k), J0=e^z-1. At xi0 the period is1. The actual E-function list1,e^z,J1,...,J_(r-1) has a rational differential system with only finite pole0. It is linearly independent over C(z): induction would force the rational derivative1/z. Exponential common-denominator bounds establish the E-function hypotheses, then Beukers ordinary-point linear lifting gives independence of values at EVERY nonzero algebraic xi. For r>=2 the phase is genuinely coupled across angular blocks; no factorization into independent periods is substituted. No arbitrary coupled simplex phase theorem is claimed.',
     dependencies=['FIREIT20261116', 'external: Beukers2006 Theorem3.2 and Lemma3.1, checked linear form']),
dict(common, id='FIRCOUPLEDMARKER20261116',
     title='Every algebraic polynomial marker on the coupled uniform cube has an exact algebraic-scale period zero test',
     scope='For arbitrary P(U)=sum p_a U^a with algebraic coefficients, its coupled exponential period has Taylor weights W(k)=sum p_a/prod_i(k+a_i+1). Verified partial fractions A_(j,h)/(k+j)^h, h<=r, and signed Stirling numbers give exact Laurent coefficients in the function basis1,e^z,J1,...,J_(r-1). At every nonzero algebraic xi the period is zero IFF all specialized coefficient entries are zero; a nonzero nonconstant-basis entry certifies transcendence, otherwise the remaining constant is its exact algebraic value. At xi0 use the direct marker mass. This is a SINGLE-period test: P=1-2u+u(1-u)v has zero exponential period at xi1 but nonzero mass1/12. Exact W=0 is separately a cube marked-power zero identity. Original Gaussian averaging and leading Gamma marker limits are retained without upgrading a zero leading period to a whole Gamma zero sequence.',
     dependencies=['FIRCOUPLEDPERIOD20261116', 'FIRBLOCKRADIAL20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRCOUPLEDBLOCKFC20261116', kind='corollary',
     title='Coupled nonflat block families retain FC and algebraic-scale eventual nonzero under arbitrary lower inactive-variable couplings',
     scope='For r disjoint factorial pairs (X_j,Y_j), S_j=X_j+Y_j, d1>=2 and d_j>=1, D=sum d_j, adjoin any finite tuple Z of independent fixed positive-shape Gamma variables. Every F=prod S_j^d_j+d1*xi*S1^(d1-2)*prod_(j>=2)S_j^(d_j-1)*prod_j X_j+K(X,Y,Z)+R_(<=D-2) has ALL sufficiently large pure moments nonzero at EVERY algebraic xi, for arbitrary COMPLEX coefficients in K,R, where K is homogeneous of total degree D-1 and K(X,Y,0)=0. Its unit limit is J_(r-1)(xi)/xi, or1 at0, unaffected by inactive couplings. For every complex xi and fixed ALGEBRAIC inactive shapes, FC and infinitely many surviving moments follow by algebraic specialization and exact power closure: d->k*d and K->k*F0^(k-1)*K at the SAME xi. This includes arbitrary factorial inactive variables. Previously proved block-separable affine/cubic nonzero limits likewise survive whenever G(X,0) is unchanged. Each active arm may aggregate positive real shapes totaling1 for algebraic-scale eventual nonzero; fixed algebraic individual shapes also permit complex-scale specialization. Gaussian lifts retain original markers and any inactive half-integer-shape coordinates. No all-complex-scale eventual nonzero or short SFC window is claimed.',
     dependencies=['FIRBLOCKRADIAL20261116', 'FIRCOUPLEDPERIOD20261116']),
dict(common, id='FIRCOUPLEDBLOCKRINGS20261116', kind='corollary',
     title='Universal ordinary tail ideals remain unit after inactive-variable lower couplings over all rational algebras',
     scope='For each fixed coupled block architecture, finite inactive tuple with positive RATIONAL shapes, and start N>=1, the ORDINARY ideal of all pure moments m>=N in Q[xi,coefficients_of_K,coefficients_of_R] is the unit ideal, where K has degree D-1 and vanishes at inactive0. The identical conclusion holds for each normalized single-top-monomial family with nonnegative exponents, fixed positive RATIONAL shapes and arbitrary lower coefficients in ALL variables. A proper tail ideal would have a complex point contradicting infinitely many nonzero moments. Unit radical implies1 lies in the ordinary ideal itself, so a finite Bezout identity exists for each tail and persists over EVERY Q-algebra. Fixed unit top families have no eventual pure zero tail over any nonzero coefficient ring, and actual scalar normalization for all marked periods is ZERO. Unit overall phase scalars are allowed. No computed cutoff, sharp nilpotence bound or short SFC window follows. The retained restricted r2,d(2,1),K0,R0 example has mu2-(8xi+24)*mu1=432.',
     dependencies=['FIRCOUPLEDBLOCKFC20261116', 'FIRMONOMIALTOPFC20261116']),
dict(common, id='FIRCOUPLEDEXEC20261116', proof_type='exact_symbolic',
     title='Exact block and coupled-period exporter handles inactive variables, active-degree markers and zero monomial exponents',
     scope='coupled_block_periods.py verifies a literal top c*prod(block sums)^d on a NONEMPTY SUBSET of all actual variables, positive rational shapes, and arbitrary lower exact algebraic coefficients. Remaining variables are inactive. It extracts the marker ACTIVE degree, takes exact unchanged inactive Gamma moments of its leading part, sets inactive0 only in the next-layer exponential, and retains simplex constraints. monomial_top_certificate permits zero exponents with at least one positive. The coupled-family exporter accepts inactive variables/shapes and homogeneous next-layer couplings vanishing at inactive0, preserving exact power closure and leading-marker period decisions. Exact cube partial fractions and independent E-basis reduction decide every algebraic-marker period at algebraic scales, with origin separate. Controls compare original Gamma/Dirichlet/Gaussian moments, active-degree/inactive-marker amplitudes, X^3 and X^2Y limit exponents, true inactive couplings, powers, small shapes, outside-domain rejections and zero leading amplitudes with nonzero actual moments. The retained 50 draft transcription checks also pass. Numerical samples are diagnostics only; the universal all-shape limit uses the written boundary and tail proof.',
     dependencies=['FIRBLOCKRADIAL20261116', 'FIRMONOMIALTOPFC20261116', 'FIRCOUPLEDPERIOD20261116', 'FIRCOUPLEDMARKER20261116', 'FIRCOUPLEDBLOCKFC20261116', 'FIRCOUPLEDBLOCKRINGS20261116'])]
new_ids = [row['id'] for row in new]


def open_scope(row):
    addition = ' FIRBLOCKRADIAL20261116 now fully audits EIT VI block-radial limit/product algebra and extends the limit to EVERY positive real Gamma shape and every polynomial marker, with small-shape boundary estimates retained. FIRCOUPLEDPERIOD20261116 proves genuine product-uniform coupled angular nonvanishing for arbitrary r and ALL nonzero algebraic scales by an explicit independent E-function chain. FIRCOUPLEDMARKER20261116 decides every algebraic polynomial cube marker period exactly; a zero leading period is not a zero moment sequence. FIRCOUPLEDBLOCKFC20261116 adds nonflat coupled block families in2r factorial/4r Gaussian variables: eventual pure nonzero for algebraic scales, FC and infinitely many surviving moments for every complex scale, arbitrary complex lower terms and exact power closure. FIRCOUPLEDBLOCKRINGS20261116 makes each universal ordinary tail ideal unit over EVERY Q-algebra and the actual all-marker scalar normalization ideal zero, without a supplied effective cutoff. Arbitrary coupled triangle phases, remaining EIT geometric/elliptic inputs, general FC>=3, unrestricted SFC and full two-term exponent-six SFC remain open.'
    if addition not in row['scope']:
        row['scope'] += addition
    monomial = ' FIRMONOMIALTOPFC20261116 also gives eventual pure nonzero for every single nonconstant top monomial, including zero exponents and arbitrary lower couplings in inactive variables, every positive real Gamma shape and arbitrary complex lower layers; XYZ+arbitrary quadratic/affine/constant terms is an explicit FC3 class, with no algebraic-scale restriction.'
    row['scope'] = row['scope'].replace(' FIRMONOMIALTOPFC20261116 also gives eventual pure nonzero for every single top monomial involving all variables, every positive real Gamma shape and arbitrary complex lower layers; XYZ+arbitrary quadratic/affine/constant terms is an explicit FC3 class, with no algebraic-scale restriction.', monomial)
    if monomial not in row['scope']:
        row['scope'] += monomial
    row['scope'] = row['scope'].replace('every single top monomial involving all variables,', 'every single nonconstant top monomial, including zero exponents and arbitrary lower inactive-variable couplings,')
    extension = ' The inactive-variable extension is CLOSED in COUPLED_BLOCK_PERIODS.tex: active radii scale as m, inactive Gamma variables keep fixed laws, G is restricted at inactive0, and leading marker parts retain exact inactive moments. X^3+arbitrary Q<=2 in any number of extra variables and coupled block families with all strictly lower-degree inactive couplings are included. This closes the requested extension, not general FC or SFC.'
    if extension not in row['scope']:
        row['scope'] += extension
    row['dependencies'] = list(dict.fromkeys(row['dependencies']+new_ids))
    return row


update_array(ROOT/'MATH_STATUS.json', 'entries', {'OP-FACTORIAL-GENERAL20261116': open_scope}, new)
known = {row['id']: row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger = lambda id: hashlib.sha256('\0'.join(known[id][key] for key in ('state', 'title', 'scope')).encode()).hexdigest()[:16]


def refresh(row):
    row['claims'] = list(dict.fromkeys(row['claims']+new_ids))
    row['reviewed_claims'] = {id: finger(id) for id in row['claims']}
    row['sources'] = list(dict.fromkeys(row['sources']+[proof, BASE+'coupled_block_audit.json']))
    return row


def work(row):
    row = refresh(row)
    row['next_step'] = 'The full block-radial limit, inactive-variable extension with arbitrary lower couplings, positive-real-shape/active-degree-marker theorem and coupled product-uniform family are closed in COUPLED_BLOCK_PERIODS.tex. Seek nonvanishing for genuinely coupled phases beyond product uniforms, especially arbitrary triangle phases linked to FC3; preserve exact marker zero tests and nonzero-leading versus zero-leading distinctions. For full two-term SFC, the d6 simple-root/higher-Newton gap remains; the earlier D+1 ODE bound is not a solution. Finish the remaining EIT V geometric/elliptic inputs and bibliography. General FC>=3 and unrestricted SFC remain separate open programme obligations.'
    row['prerequisites'] = 'Read the coupled-block proof and current source ledger. The block limit allows inactive variables and all positive real shapes; use active marker degree and unchanged inactive Gamma expectations. Complex-scale specialization requires fixed algebraic shapes. The arithmetic cube certificates require algebraic scales and marker coefficients. Functional E-independence and an ordinary nonzero point are both needed for Beukers lifting. Products of uniforms are a genuine coupled class, not every triangle phase. Zero leading marker periods do not prove whole zero sequences. Complex-scale FC uses algebraic specialization plus exact power closure; all-ring ordinary tail-unit ideals give existence of finite Bezout witnesses, not an effective universal cutoff or SFC window.'
    return row


def lesson(row):
    row = refresh(row)
    extra = ' For block-monomial radial tops, control the compact coordinate boundary and the unbounded tail before dividing by the top. This extends the full limit to every positive real Gamma shape and all original averaged markers. Inactive Gamma variables keep their original laws and may occur freely in strictly lower total-degree terms. The scaled maximum lies at (active d,inactive0); its normalization retains m^(-r/2-B). Restrict G at inactive0, but take actual inactive expectations in active-degree marker amplitudes. The coupled product-uniform period has an explicit E-function chain1,e,J1,... with a short rational-residue independence proof; arithmetic lifting then gives every algebraic-scale unit value and a complete exact polynomial-marker period test. This closes arbitrarily many coupled angular blocks without a general geometric irreducibility theorem. Exact power closure gives complex-scale FC, and ordinary universal tail-unit ideals persist over all Q-algebras.'
    row['use'] = row['use'].replace(' For block-monomial radial tops, control the compact coordinate boundary and the unbounded tail before dividing by the top. This extends the full limit to every positive real Gamma shape and all original averaged markers. The coupled product-uniform period has an explicit E-function chain1,e,J1,... with a short rational-residue independence proof; arithmetic lifting then gives every algebraic-scale unit value and a complete exact polynomial-marker period test. This closes arbitrarily many coupled angular blocks without a general geometric irreducibility theorem. Exact power closure gives complex-scale FC, and ordinary universal tail-unit ideals persist over all Q-algebras.', extra)
    if extra not in row['use']:
        row['use'] += extra
    avoid = ' Do not infer a whole marked zero sequence from a zero leading period, all-complex-scale eventual nonzero from algebraic specialization, arbitrary coupled triangle nonvanishing from the product-uniform class, or a computable cutoff from a nonconstructive ordinary tail-unit ideal. Retain each block simplex constraint when simplifying the angular phase. Auditing the block-limit theorem does not audit its remaining geometric/elliptic angular inputs.'
    if avoid not in row['avoid']:
        row['avoid'] += avoid
    row['boundary'] = 'Complete block-limit/product algebra, inactive-variable lower-coupling extension and all-positive-shape/active-degree-marker limit are closed. Genuine coupled product-uniform block families now satisfy algebraic-scale eventual nonzero and complex-scale FC/power closure, with exact algebraic marker periods and ordinary all-ring tail-unit/zero-normalization results. Arbitrary coupled simplex phases, general FC>=3, unrestricted/two-term mixed-exponent SFC and remaining EIT geometric/elliptic/bibliography inputs stay open; outside/formal/novelty flags remain false.'
    row['implementation'] = list(dict.fromkeys(row['implementation']+[BASE+'coupled_block_periods.py', checker]))
    return row


update_array(ROOT/'knowledge/work_items.json', 'actions', {'WORK-FACTORIAL-GENERAL': work})
update_array(ROOT/'knowledge/lessons.json', 'lessons', {'METHOD-FACTORIAL-RADIAL-PURITY': lesson})
print(json.dumps({'registered': new_ids, 'previous_goal_turn': 'progress', 'inactive_variable_extension_complete': True, 'general_FC_and_SFC_closed': False}))
