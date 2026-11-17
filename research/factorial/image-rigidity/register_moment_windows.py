"""Register finite-window/source-audit progress with guarded narrow edits."""
from pathlib import Path
import hashlib
import json
import os
import pypdf

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'factorial/image-rigidity/'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()

# Historical continuation must not overwrite the later all-index prime-power scope.
if any(row['id']=='FIRPADICDEFORM20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Later prime-power progress is registered; use register_prime_power_windows.py.')


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
            assert known[row['id']] == row, ('existing window record differs', row['id'])
        else:
            new.append(row)
    if new:
        head, end = text.rsplit('\n  ]\n}', 1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(encode(row) for row in new)+'\n  ]\n}'+end
    if text.encode() != before:
        assert path.read_bytes() == before, 'Concurrent edit'
        temp = path.with_suffix(path.suffix+'.windows-research/.scratch')
        temp.write_bytes(text.encode())
        os.replace(temp, path)


for filename in ('verification.json', 'long_audit_verification.json',
                 'simplex_pencil_verification.json', 'two_angular_verification.json',
                 'moment_window_verification.json'):
    receipt = json.loads((HERE/filename).read_text())
    assert receipt['status'] == 'pass'
    for file, expected in receipt['files'].items():
        assert sha(HERE/file) == expected, ('retained proof/control changed', file)

ledger_path = HERE/'source_ledger.json'
ledger = json.loads(ledger_path.read_text())
sources = [
('sources/rocks-2015-author-preprint.pdf', 'https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf',
 'This is the September6,2015 author ARTICLE PREPRINT, not the dissertation. Core linear-power/reciprocal section, affine(corrected window), common-prime power, homogeneous lifting and disjoint-linear-adjunction proofs checked. Page9 affine display has an off-by-one typo; proof gives the correct m+1 window. Page10 Lemma3.9 uses algebraic closures, verified visually. Proposition3.11 retains n+m<p; other two-monomial/Newton polygon claims are not all audited. Final2017 publisher text not obtained; no assertion that the draft typo occurs in it.'),
('sources/filaseta-rocks-author-preprint.pdf', 'https://people.math.sc.edu/filaseta/papers/FilasetaRocksPaperSubmitted.pdf',
 'Primary author preprint of the publisher-confirmed2016 Colloquium article. Fixed-d lower density>=log2 irreducibility statement and proof architecture mapped; number-theoretic density inputs not independently reconstructed. The moment polynomial and adjacent derivative identity checked independently. Not a theorem for all n,d or general two-monomial SFC.')]
for file, url, review in sources:
    if file not in {row['file'] for row in ledger['cached_primary_sources']}:
        ledger['cached_primary_sources'].append(dict(file=file, url=url, sha256=sha(HERE/file), review=review))
    text_file = str(Path(file).with_suffix('.txt'))
    if text_file not in {row['file'] for row in ledger['cached_primary_sources']}:
        ledger['cached_primary_sources'].append(dict(file=text_file, url=url, sha256=sha(HERE/text_file),
            review='Local text extraction with pypdf '+pypdf.__version__+' in layout mode. PDF remains authoritative; ambiguous affine range and finite-field bars checked on page images.'))
for url in ('https://doi.org/10.1016/j.jalgebra.2017.01.039', 'https://doi.org/10.7936/K7BR8QB4',
            'https://doi.org/10.4064/cm6761-3-2016'):
    if url not in ledger['additional_primary_references']:
        ledger['additional_primary_references'].append(url)
ledger['published_input_pinpoints']['Rocks2017'] = 'Publisher metadata and abstract confirmed; proof audit uses the separately identified2015 author preprint Theorems3.3/3.7(corrected)/3.10/3.12 and Proposition3.4. Common-prime exponents and partial affine range retained. Full final journal text and full dissertation remain outside this audit.'
ledger['published_input_pinpoints']['FilasetaRocks2016'] = 'Publisher confirms Colloq.Math145(2016)307-314, DOIcm6761-3-2016. Author preprint Theorem1.1 gives lower density>=log2 in n for every FIXED d. Independent adjacent derivative identity converts squarefreeness to two-term window certificates, without all-index extrapolation.'
ledger_path.write_text(json.dumps(ledger, indent=2)+'\n')

proof = BASE+'SFC_WINDOWS_AND_ROCKS_AUDIT.tex'
checker = BASE+'verify_moment_windows.py'
audit = dict(date='2026-11-16', previous_goal_turn='progress', current_goal_turn='progress',
    previous_evidence='All four previous proof/checker receipts match their recorded source hashes, including the534-check cubic-family packet.',
    canonical_proof=proof,
    progress=['reciprocal r-section polynomial identity over every field gives all-index consecutive windows',
              'proper rational marked Gamma generating functions decide zero sequences and export certified windows',
              'integer sums at colliding nonzero nodes allow rational individual Gamma shapes including half-integers',
              'original simplex and Gaussian markers retained by homogenization and averaging',
              'Rocks linear-power/common-prime-power/disjoint-linear-adjunction proofs checked',
              'corrected source identity and the draft affine-window statement; visually verified algebraic closures',
              'mapped fixed-d irreducibility density and checked the exact adjacent-squarefree equivalence',
              'retained a failed cubic mod-five upgrade with12 modular exceptions rather than extrapolating125 rational gcds'],
    source_limits=['2015 author preprint is not the dissertation', 'full2017 publisher text was unavailable',
                   'full dissertation and all two-monomial claims not audited', 'Filaseta--Rocks density proof mapped, not independently reconstructed'],
    remaining=['general FC>=3', 'unrestricted SFC', 'all-index squarefreeness of f_n,d for general d>=3',
               'arbitrary noninteger active-node shape sums', 'full EIT V/VI and bibliography review'],
    independent_replay=False, formal_verification=False, external_review=False,
    novelty_claim=False, whole_goal_complete=False)
(HERE/'moment_window_audit.json').write_text(json.dumps(audit, indent=2)+'\n')

common = dict(kind='theorem', state='proved', canonical_source=proof, checker=checker,
    proof_type='hybrid', independent_replay=False, formal_verification=False, external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_moment_windows.py'),
    software_lock=[proof, checker, BASE+'moment_windows.py', BASE+'moment_bridges.py',
                   BASE+'moment_window_verification.json', BASE+'moment_window_examples.json',
                   BASE+'moment_window_audit.json', BASE+'two_angular_residue_probe.json', BASE+'source_ledger.json'],
    supersedes=[], closes_problems=[], narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[], invalidates_assumptions=[], replaced_by=[], priority='reference', programme_status='archived')
new = [
dict(common, id='FIRRECIPWINDOW20261116', kind='reproduction', proof_type='reproduction',
     title='Rocks reciprocal-polynomial sections have a nonzero coefficient in every degree-sized consecutive sampled window',
     scope='For EVERY field K, P in K[z] with P0=1 and degreeN>0, write 1/P=sum a_k z^k. For every n,r>=1, one of a_(rn),...,a_(r(n+N-1)) is nonzero. Factor P=product(1-c_i*z), define Q(t)=product(1-c_i^r*t), and use W(z)=Q(z^r)/P(z)=product sum_(j0..r-1)(c_i*z)^j. The r-section is numerator/Q with numerator constant1 and degree<N. An N-zero block propagates to an infinite zero tail, impossible for a nonzero proper rational function. No division by r or characteristic-zero assumption. This reproduces Rocks2015 author-preprint Proposition3.4/Corollary3.6; it concerns multiplicative inverse series, not Furter composition inverses.',
     dependencies=['external: Rocks2015 author preprint Proposition3.4 and Corollary3.6']),
dict(common, id='FIRMARKEDWINDOW20261116',
     title='Linear-power Gamma phases with integer active-node shape sums have exact consecutive-window and zero-marker certificates',
     scope='Independent positive-shape Gamma variables, L=sum c_iY_i nonzero, F=L^r with integer r>=1. If the shape sum at EACH distinct nonzero node c is an integer, unit moments have a nonzero value in every N=sum_(c_i!=0)alpha_i consecutive positive exponents. Integer individual shapes are included, as are half-integer or other positive individual shapes aggregating to integer nodes; zero nodes need no integrality. EVERY polynomial marker M has rational K_M(z)=E[M*exp(zL)]; its r-section is proper. The selected marked sequence is either identically zero INCLUDING exponent0, or every d consecutive positive moments contain a nonzero one, where d is the reduced denominator degree. A valid upper bound is sum_(c_i!=0)(alpha_i+deg_Yi M). All original simplex markers transfer after homogeneous lifting, and original Gaussian markers transfer by averaging. Actual zero phase is separate and can have a nonzero constant marked moment. For original factorial shapes, the unit bound is at most N(L^r), proving genuine SFC for powers of linear forms. General Gamma-shape SFC with the original monomial-count window, arbitrary noninteger active-node sums, optimal gap bounds, general polynomial phases and unrestricted SFC are NOT asserted.',
     dependencies=['FIRRECIPWINDOW20261116', 'FIRGAMMAADAPTER20261116']),
dict(common, id='FIRSFCAFFINE20261116',
     title='Correct affine windows and independent Gamma-variable adjunction preserve certified consecutive-moment bounds',
     scope='For a nonzero constant a and Gamma-linear L with integer active-node shape sums and denominator degreeN, normalized moments of a+L satisfy P(z)*series=exp(a*z), excluding every N+1 consecutive zero block. At a0 use the homogeneous N bound. For any nonzero phase F already having a window B, adding a new independent Gamma(shape alpha integer) term lambda*U gives the bound B+alpha by the exact degree-alpha difference identity of normalized moments. In the factorial case alpha1, a new nonzero variable term increases both the bound and actual monomial count by one, preserving SFC; iteration gives the disjoint-linear extension in Rocks Theorem3.12/Corollary3.13. The2015 preprint affine Theorem3.7 DISPLAY says m moments but its proof establishes m+1, the correct SFC bound; F=Y-1 refutes the shorter display. The page image confirms the typo, which is NOT attributed to the unavailable final2017 journal text. The linked2015 PDF is an article preprint, not the separate dissertation.',
     dependencies=['FIRRECIPWINDOW20261116', 'FIRMARKEDWINDOW20261116', 'external: Rocks2015 Theorem3.7 proof and Theorem3.12/Corollary3.13']),
dict(common, id='FIRSFCPRIME20261116', kind='reproduction', proof_type='reproduction',
     title='Rocks common-prime-power linear sums and all their positive powers satisfy Strong Factorial',
     scope='For a SINGLE fixed prime p, positive integers beta_i, G=sum_(i1..m)lambda_i*Y_i^(p^beta_i), every positive power G^r has a nonzero factorial moment in every m consecutive positive exponents. Normalize by n!, scale each coefficient by(-p)^(-b_i), b_i=(p^beta_i-1)/(p-1). The factorial ratio (p^beta_i*k)!/(p^(b_i*k)*k!) is integral and congruent to(-1)^(b_i*k); thus every normalized coefficient polynomial reduces to h_n. The sampled m-window has no nonzero common zero over ALGEBRAIC CLOSURE F_p by FIRRECIPWINDOW; homogeneous Nakayama/projective properness lifts it to characteristic0. The field-closure bars are visually verified in preprint Lemma3.9; rational finite-field point scans would be insufficient. The constant-added Proposition3.11 only covers starts n+m<p. This reproduces preprint Theorem3.10 with the publisher-confirmed2017 special-case scope; it does NOT cover mixed prime bases, all constant-added windows or unrestricted two-monomial SFC.',
     dependencies=['FIRRECIPWINDOW20261116', 'external: Rocks2015 author preprint Theorem3.10, Lemmas3.8-3.9']),
dict(common, id='FIRWINDOWEXEC20261116', proof_type='exact_symbolic',
     title='Exact rational section exporter retains original markers, node collisions and certified all-index moment windows',
     scope='moment_windows.py accepts exact algebraic nodes/marker coefficients and positive rational individual shapes whose sums at distinct nonzero nodes are integers. It constructs the actual marked Gamma rational generating function, exact polynomial r-section identity, cancels the final denominator and exports either an identically zero selected marked sequence or a certified window d. Root-of-unity aliases, repeated/zero nodes, half-shape aggregation and zero phase have explicit branches; no finite prefix determines the all-index claim. Reduced denominator degree is a bound, not an optimal zero-gap claim. Controls independently compare original Gamma/simplex/Gaussian moments, affine and adjunction differences, common-prime congruences and adjacent moment-polynomial derivative identities. The retained cubic probe has12 modular common-root cases among125 tested residues despite invertible fifth-power matrices, excluding a blanket mod5 upgrade; constant rational gcds in that finite box are not a universal theorem. Outside-class inputs fail closed.',
     dependencies=['FIRMARKEDWINDOW20261116', 'FIRSFCAFFINE20261116', 'FIRSFCPRIME20261116'])]
new_ids = [row['id'] for row in new]


def open_scope(row):
    addition = ' FIRRECIPWINDOW20261116/FIRMARKEDWINDOW20261116 now give all-index consecutive windows and exact zero-marker decisions for linear powers with integer active-node Gamma shape sums, including half-shape aggregation, all original markers and repeated-node/alias cancellations. FIRSFCAFFINE20261116/FIRSFCPRIME20261116 check the corrected affine and common-prime-power SFC cases and independent-variable construction from Rocks. The linked2015 source is an article preprint, not the dissertation; full final2017 text and other draft two-monomial claims remain unaudited. The general two-term constant-plus-monomial window is equivalent to all-index squarefreeness of f_n,d; fixed-d irreducibility density is not all-index closure. The failed cubic residue upgrade is retained rather than extrapolated. General FC>=3 and unrestricted SFC remain open.'
    if addition not in row['scope']:
        row['scope'] += addition
    row['dependencies'] = list(dict.fromkeys(row['dependencies']+new_ids))
    return row


update_array(ROOT/'MATH_STATUS.json', 'entries', {'OP-FACTORIAL-GENERAL20261116': open_scope}, new)
known = {row['id']: row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger = lambda id: hashlib.sha256('\0'.join(known[id][k] for k in ('state', 'title', 'scope')).encode()).hexdigest()[:16]


def refresh(row):
    row['claims'] = list(dict.fromkeys(row['claims']+new_ids))
    row['reviewed_claims'] = {id: finger(id) for id in row['claims']}
    row['sources'] = list(dict.fromkeys(row['sources']+[proof, BASE+'moment_window_audit.json']))
    return row


def work(row):
    row = refresh(row)
    row['next_step'] = 'For genuine SFC progress beyond the closed classes, study all-index squarefreeness of f_n,d=sum binom(n,j)*(d*j)!*x^j for general d>=3, equivalent to the constant-plus-monomial adjacent-window problem; distinguish fixed-d density from all indices. Continue arbitrary FC>=3/coupled simplex phases, complete EIT V/VI inputs and published/dissertation bibliography. Do not repeat closed linear-power or common-prime-power searches, mistake a marked zero sequence for a pure-null phase, or extend the failed blanket cubic mod5 argument.'
    row['prerequisites'] = 'Read SFC_WINDOWS_AND_ROCKS_AUDIT.tex, exact rational section outputs and original-marker identities. Keep integer ACTIVE-NODE sums, all-zero-phase branch, corrected affine m+1 window, common prime hypothesis, residue-field algebraic closure and partial n+m<p range. The2015 author article is separate from the thesis and the inaccessible final2017 text; later two-monomial draft claims are not automatically adopted. Earlier cubic, pencil and binary audit boundaries remain.'
    return row


def lesson(row):
    row = refresh(row)
    extra = ' For powers of a Gamma-linear phase with integer shape sums at each nonzero node, use the exact proper rational r-section: unit moments have a degree-sized consecutive window, and each polynomial marker is either identically zero or has a denominator-degree window. This retains collisions, root-of-unity aliases and rational individual shapes aggregating to integer nodes. Rocks common-prime sums and disjoint-linear adjunction supply genuine SFC classes; the affine draft proof needs m+1 moments.'
    if extra not in row['use']:
        row['use'] += extra
    avoid = ' Do not identify an article preprint with a dissertation or final publisher version, omit the COMMON prime in Rocks powers, replace algebraic-closure finite-field emptiness with rational point emptiness, infer all-index SFC from irreducibility density, or infer all-residue cubic nonvanishing from invertible fifth-power matrices.'
    if avoid not in row['avoid']:
        row['avoid'] += avoid
    row['boundary'] = 'Additional closed classes: all-index linear-power marked windows with integer active-node sums, corrected affine/disjoint-variable windows, and common-prime-power SFC. General FC>=3, unrestricted/two-term SFC, arbitrary noninteger active-node sums, full EIT V/VI and bibliography remain open. Full2017 text and other draft number-theoretic claims are not fully audited; all assurance flags remain false.'
    row['implementation'] = list(dict.fromkeys(row['implementation']+[BASE+'moment_windows.py', checker]))
    return row


update_array(ROOT/'knowledge/work_items.json', 'actions', {'WORK-FACTORIAL-GENERAL': work})
update_array(ROOT/'knowledge/lessons.json', 'lessons', {'METHOD-FACTORIAL-RADIAL-PURITY': lesson})
print(json.dumps({'registered': new_ids, 'previous_goal_turn': 'progress', 'whole_goal_complete': False}))
