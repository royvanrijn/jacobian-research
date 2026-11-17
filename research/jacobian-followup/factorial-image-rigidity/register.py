"""Historical initial registration; later scopes use register_long_audit.py.

Existing changed rows intentionally fail the equality guard. Never use this
initial snapshot to overwrite a later review or its source ledger.
"""
from pathlib import Path
import hashlib
import json
import os

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'jacobian-followup/factorial-image-rigidity/'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()


def append_records(path, key, rows):
    before = path.read_bytes()
    text = before.decode()
    data = json.loads(text)
    known = {r['id']: r for r in data[key]}
    new = []
    for row in rows:
        if row['id'] in known:
            if known[row['id']] != row:
                raise ValueError('Existing record differs; inspect before replacing: '+row['id'])
        else:
            new.append(row)
    if not new:
        return
    # All three established files end in the corresponding top-level array.
    tail = '\n  ]\n}'
    head, end = text.rsplit(tail, 1)
    if end.strip():
        raise ValueError('Unexpected registry trailer')
    payload = ',\n'.join('\n'.join('    '+line for line in json.dumps(row, indent=2).splitlines()) for row in new)
    after = (head+',\n'+payload+tail+end).encode()
    if path.read_bytes() != before:
        raise RuntimeError('Concurrent edit: '+str(path))
    temporary = path.with_suffix(path.suffix+'.factorial-tmp')
    temporary.write_bytes(after)
    os.replace(temporary, path)


sources = [
    ('long-fc2.tex', 'https://raw.githubusercontent.com/long-mathematics/factorial-conjecture-two-variables/8b1eb46ff530122ed25f588b78bf9fdc98cb9577/paper/factorial_conjecture_two_variables.tex',
     'Flat radial lemma and affine-next-layer subcase checked; full nonflat proof remains under review.'),
    ('long-fc2-readme.md', 'https://raw.githubusercontent.com/long-mathematics/factorial-conjecture-two-variables/8b1eb46ff530122ed25f588b78bf9fdc98cb9577/README.md',
     'Dependency and manuscript-status statements read.'),
    ('long-eit-pinned.tex', 'https://raw.githubusercontent.com/long-mathematics/exponential-integral-theorem/22d98c4f27e606b5b25fa673c33f9c56c53cc86b/papers/01-exponential-integral-theorem/exponential-integral-theorem.tex',
     'Invoked Corollary 1.4 matched exactly; selected technical sections inspected, full proof not verified.'),
    ('rigidity.typ', 'https://raw.githubusercontent.com/blueberryvertigo/polynomial-composition-rigidity/051c74aad9ef9bf8804c371964803aee099c5b63/rigidity.typ',
     'Sections 2-5 fully checked in the companion proof. Later scheme/polydegree consequences not independently audited.'),
    ('rigidity-readme.md', 'https://raw.githubusercontent.com/blueberryvertigo/polynomial-composition-rigidity/051c74aad9ef9bf8804c371964803aee099c5b63/README.md',
     'Claim and source location read.'),
]
manifest = {
    'date': '2026-11-16',
    'cached_primary_sources': [dict(file='sources/'+name, url=url, sha256=sha(HERE/'sources'/name), review=review)
                               for name, url, review in sources],
    'additional_primary_references': [
        'https://arxiv.org/abs/1008.3962v2',
        'https://arxiv.org/abs/1304.3956v2',
        'https://doi.org/10.1017/S0004972719000546',
        'https://doi.org/10.1112/jlms/jdu064',
        'https://arminstraub.com/downloads/pub/polydegree-conjecture.pdf',
        'https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf',
        'https://github.com/ipitchford/furter-r3-through-299',
        'https://sbseminar.wordpress.com/2026/07/20/the-new-counterexample-to-the-jacobian-conjecture/',
    ],
    'bibliography_complete': False,
    'external_review': False,
    'novelty_claim': False,
}
if not (HERE/'source_ledger.json').exists():
    (HERE/'source_ledger.json').write_text(json.dumps(manifest, indent=2)+'\n')

checker = BASE+'verify_bridges.py'
proof = BASE+'RIGIDITY_CRITICAL_VALUES.tex'
survey = BASE+'STATUS_AND_LITERATURE.tex'
common = dict(independent_replay=False, formal_verification=False, external_review=False,
              software_lock=[proof, survey, checker, BASE+'moment_bridges.py', BASE+'verification.json', BASE+'source_ledger.json'],
              supersedes=[], closes_problems=[], narrows_problems=[], consumers=[], invalidates_assumptions=[],
              replaced_by=[], priority='reference', programme_status='archived')


def claim(id, kind, title, scope, source, dependencies, proof_type='hybrid', state='proved'):
    return dict(common, id=id, kind=kind, state=state, title=title, scope=scope,
                canonical_source=source, dependencies=dependencies, checker=checker,
                proof_type=proof_type, artifact_hash='sha256:'+sha(HERE/'verify_bridges.py'))


rows = [
    claim('FIRRIGIDITY20261116', 'reproduction', 'The critical-value argument proves full Furter composition and inverse-coefficient rigidity',
          'All complex polynomial degrees and all positive inverse-window starts. For F=f_s composed through f_1 tangent to the identity and nonidentity, ord_0(F-z)<=1+sum #CV(f_i)<=1+sum(deg f_i-1). Perturbation creates k attracting fixed points; each full basin contains a critical point and distinct basins require distinct critical values. Complete internal written analytic verification of liqsweep September2026 Sections2-5, credited and pinned; no external refereeing/formal verification claim. Inverse truncation proves R(m)_n for every m,n>=1. A nonidentity polynomial a has no #CV(a) consecutive zero inverse coefficients above its linear coefficient. Finite symbolic tests check algebra/indices, not the analytic universal statement. General FC/SFC are not proved.',
          proof, ['external: liqsweep, polynomial-composition-rigidity, commit 051c74aad9ef9bf8804c371964803aee099c5b63'], 'reproduction'),
    claim('FIRSFCFAMILY20261116', 'corollary', 'Full rigidity transfers to the Edo--van den Essen factorial family and its simplex moments',
          'For f_mu=(product y_i)*sum mu_i*y_i in every dimension, L(f_mu^k)=(k+1)*(k!)^(r+1)*[z^(k+1)] inverse(z*product(1-mu_i*z)). Hence the nonzero family has no #CV(a_mu) consecutive zero factorial moments, in particular no N(f_mu) consecutive zeros. Zero parameters reduce by exact independent factorial factors. Homogeneity transfers the same assertion to p_mu on the uniform (r-1)-simplex, and block-radial substitution to Gaussian degree2(r+1). Exact bridge credited to Edo--van den Essen2014; rigidity argument credited to liqsweep2026. This is special-family SFC, not unrestricted SFC or triangle FC.',
          proof, ['FIRRIGIDITY20261116', 'external: Edo--van den Essen 1304.3956v2 Theorem2.25 and Lemma2.26']),
    claim('FIRGAMMAADAPTER20261116', 'theorem', 'Executable exact factorial and half-integer Gamma/Dirichlet adapters preserve original markers',
          'Exact polynomial moments, homogenization, block-Gaussian lifts, and arbitrary Gaussian monomial Reynolds factors are implemented. Serializer matches and is tested against the existing whole-source Gaussian producer QQ/simple-number-field parser. All-marker safety transfers through the previously written orthogonal-averaging proof. The serializer does not run a new Gaussian CAS campaign and does not decide purity on an unsafe phase. The finite exact controls include retained SIC2C4 formulas and characteristic-p negative controls. No universal FC/SFC or new Gaussian classification is asserted.',
          survey, ['GMCSOURCEHOMGEN20261116', 'SIC2C4'], 'exact_symbolic'),
    claim('FIRFLATAFFINE20261116', 'reproduction', 'The binary flat leading-part and affine next-layer factorial case needs no EIT input',
          'Every total degree D>=1, complex binary f with highest homogeneous part c*(x+y)^D for nonzero c. After normalization f(rt,r(1-t))=r^D+A1(t)*r^(D-1)+lower terms, assume deg A1<=1. An all-moment counterexample specializes algebraically retaining these rational equations. Uniform radial expansion gives L(f^n)/Gamma(Dn+2) tending to integral_0^1 exp(A1(t)/D)dt; the factorial-ratio bound is at most1 and supports dominated convergence. Hermite--Lindemann makes this affine exponential integral nonzero. Checked reproduction of Long pinned flat subcase, not a verification of Long full nonflat proof or general EIT. No new priority claim.',
          survey, ['external: Long FC2 commit 8b1eb46ff530122ed25f588b78bf9fdc98cb9577 flat case', 'external: Hermite--Lindemann theorem'], 'reproduction'),
]
open_row = dict(common, id='OP-FACTORIAL-GENERAL20261116', kind='open_problem', state='open',
                title='Resolve general Factorial and Strong Factorial conjectures and audit full binary FC claims',
                scope='The full thread programme remains open: prove or refute FC and SFC for arbitrary complex multivariate polynomials; complete current literature review and full Long FC2/EIT proof audit; exploit Gaussian/Gamma/Dirichlet certificates without assuming radial purity implies safety. Full rigidity is internally proved in FIRRIGIDITY20261116, special-family SFC in FIRSFCFAMILY20261116, and unrestricted Image is already disproved by SIC2C4. Homogeneous FC3 remains the general uniform-triangle pure-null question. Pointwise all-marker safety, finite zero prefixes and diagonal shadows do not settle these universal targets. No bibliography-completeness, independent-replay, formal or external-review assurance is claimed.',
                canonical_source=survey, dependencies=['FIRRIGIDITY20261116','FIRSFCFAMILY20261116','FIRGAMMAADAPTER20261116','FIRFLATAFFINE20261116','SIC2C4','FM1'],
                checker=None, proof_type='not-applicable', artifact_hash=None, priority='primary')
rows.append(open_row)
append_records(ROOT/'MATH_STATUS.json', 'entries', rows)

known = {r['id']: r for r in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
fingerprint = lambda id: hashlib.sha256('\0'.join(known[id][k] for k in ('state','title','scope')).encode()).hexdigest()[:16]
ids = [r['id'] for r in rows]+['SIC2C4','FM1']
action = dict(id='WORK-FACTORIAL-GENERAL', kind='research', area='geometry', programme_status='archived',
              title='Complete general FC/SFC and review Long full FC2/EIT after the rigidity closure',
              target=open_row['id'], claims=ids, sources=[survey,proof,BASE+'source_ledger.json'],
              next_step='Review the entire pinned Long nonflat continuation/monodromy argument and pinned EIT proof; build original-unit Gaussian period certificates for purity on unsafe phases; seek a structural theorem for the independent-block radial pure-null sector. Do not redo finite R3 windows or already disproved universal Image.',
              done_when='The full arbitrary-polynomial Factorial and Strong Factorial questions have proofs or counterexamples, Long claims and bibliography have completed scoped reviews, and transfers are backed by original-marker certificates. The existing rigidity proof and special-family SFC alone do not complete the objective.',
              prerequisites='Read the canonical survey, pinned source ledger, full critical-value proof, SIC2C4 and FM1. Retain unpublished versus established status, nonmultiplicativity of diagonal extraction, and the distinction between pure-nullity and eventual marked safety. No external messages or publishing are authorized.',
              compute='none', reviewed_claims={id:fingerprint(id) for id in ids})
append_records(ROOT/'knowledge/work_items.json', 'actions', [action])
lesson = dict(id='METHOD-FACTORIAL-RADIAL-PURITY', area='geometry', programme_status='archived',
              title='Use exact Gamma bridges while keeping purity, safety and moment windows distinct',
              when='Transferring Gaussian machinery to factorial or simplex moments and reporting current rigidity/Image status.',
              use='Homogenize with the degree retained, average every original Gaussian marker over independent orthogonal blocks, and keep the nonzero Gamma factor. Use the complete credited critical-value proof for full rigidity and its exact inverse-factorial identity for the special SFC family.',
              avoid='Treating an unsafe Gaussian phase as a purity refutation, a finite zero prefix as an infinite tail, an Image counterexample as a radial factorial polynomial, or a diagonal shadow of f^m as powers of one fixed polynomial. Do not reopen finite R3 computations when the universal analytic proof answers the same statement.',
              boundary='General FC/SFC and full Long FC2/EIT verification remain open here. Uniform triangles mean homogeneous FC3, not homogeneous binary FC. Internal analytic verification is separate from external refereeing and proof-assistant verification.',
              reopen_when='A valid original-unit period or universal radial-sector theorem closes the remaining pure-null implication, or a specific flaw is found in the credited analytic argument.',
              claims=ids, sources=[survey,proof], implementation=[BASE+'moment_bridges.py',checker],
              keywords=['factorial','Gamma','Dirichlet','rigidity','Image','pure-null','safety'], reviewed_claims={id:fingerprint(id) for id in ids})
append_records(ROOT/'knowledge/lessons.json', 'lessons', [lesson])
print(json.dumps({'registered':[r['id'] for r in rows], 'work':action['id'], 'lesson':lesson['id']}))
