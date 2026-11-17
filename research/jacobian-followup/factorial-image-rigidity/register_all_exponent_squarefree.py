"""Register the elementary cutoff and the transfer to every positive exponent."""
from pathlib import Path
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='jacobian-followup/factorial-image-rigidity/'
PROOF=BASE+'ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex'
CHECKER=BASE+'verify_all_exponent_squarefree.py'
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()

def update_array(path, key, changes, additions=()):
    before = path.read_bytes()
    text = before.decode()
    known = {row['id']:row for row in json.loads(text)[key]}
    encode = lambda row: '\n'.join('    '+line for line in json.dumps(row,indent=2).splitlines())
    for identifier, change in changes.items():
        old = known[identifier]
        row = change(dict(old))
        if row != old:
            needle = encode(old)
            assert text.count(needle)==1, ('record formatting changed',identifier)
            text = text.replace(needle,encode(row),1)
    extra = []
    for row in additions:
        if row['id'] in known:
            old=known[row['id']]
            if old!=row:
                needle=encode(old)
                assert text.count(needle)==1,('record formatting changed',row['id'])
                text=text.replace(needle,encode(row),1)
        else:
            extra.append(row)
    if extra:
        head,end = text.rsplit('\n  ]\n}',1)
        assert not end.strip()
        text = head+',\n'+',\n'.join(encode(row) for row in extra)+'\n  ]\n}'+end
    if text.encode()!=before:
        assert path.read_bytes()==before, 'Concurrent edit'
        temp=path.with_suffix(path.suffix+'.real-place-tmp')
        temp.write_bytes(text.encode())
        os.replace(temp,path)


receipt=json.loads((HERE/'all_exponent_squarefree_verification.json').read_text())
assert receipt['status']=='pass' and receipt['large_exponent_criterion_proved']
assert receipt['all_positive_exponents_share_one_exceptional_set'] and receipt['exponent_one_explicit']
assert receipt['MR_II_source_qualification_retained'] and not receipt['all_index_SFC_closed']
for file,expected in receipt['files'].items():assert sha(HERE/file)==expected,('stale all-exponent receipt',file)
common=dict(kind='theorem',state='proved',canonical_source=PROOF,checker=CHECKER,proof_type='hybrid',
    independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_all_exponent_squarefree.py'),
    software_lock=[PROOF,CHECKER,BASE+'all_exponent_squarefree.py',BASE+'all_exponent_squarefree_verification.json',
        BASE+'all_exponent_squarefree_audit.json',BASE+'all_exponent_squarefree_examples.json',BASE+'register_all_exponent_squarefree.py',
        BASE+'prime_power_windows.py',BASE+'PRIME_POWER_TWO_TERM_SFC.tex',BASE+'source_ledger.json',
        BASE+'SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex',BASE+'simultaneous_squarefree_verification.json',
        BASE+'simultaneous_squarefree_audit.json',BASE+'sources/rocks-2015-author-preprint.pdf',
        BASE+'sources/rocks-2015-author-preprint.txt'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],consumers=[],
    invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')
new=[dict(common,id='FIRLARGEM20261117',
    title='An elementary coefficient criterion makes every sufficiently large factorial exponent squarefree at a given index',
    scope='For ALL positive integers n,m, 2^m>=4*n^2 implies F_(n,m)(T)=sum binom(n,j)*(mj)!*T^j is squarefree over C and Q, ordinary Gamma shape1. Equality is allowed. The elementary factorial inequality u!/(v!)^(u/v)<=(u/v)^(u/2) follows from monotonicity of (k!)^(1/k)/sqrt(k), proved by AM-GM and (1+1/k)^(k+1)<=4. Normalize the constant and leading terms to1+z^n; intermediate positive coefficients b_j are bounded by binom(n,j)*q^min(j,n-j), q=2^(-m/2), using the concave logarithm chord. Symmetry bounds both weighted coefficient masses A=sum(j/n)*b_j and B=sum(1-j/n)*b_j by (1+q)^n-1<(exp(1/2)-1)<1 when q<=1/(2n). A repeated root of modulus<=1 contradicts B<1 via P-alpha*P_prime/n; one of modulus>=1 contradicts A<1 via alpha*P_prime/(n*alpha^n). No asymptotic analytic input, prime-power restriction or density extrapolation is used. Exact bit-length certificates avoid constructing2^m for an enormous supplied exponent. Other Gamma shapes, full all-index SFC and novelty/assurance upgrades are unclaimed.',
    dependencies=[]),
 dict(common,id='FIRLOGEXPCANDIDATES20261117',
    title='Every hypothetical bad exponent is logarithmically small and has at least two distinct prime divisors',
    scope='For EVERY positive integer index n, any nonsquarefree ordinary factorial F_(n,m) must satisfy STRICTLY m<2*log2(2n), equivalently2^m<4*n^2, and m has at least TWO DISTINCT prime divisors. The coefficient theorem excludes the complementary large-exponent range. All-index prime-power squarefreeness was ALREADY registered in FIRPADICDEFORM/FIRSFCMONOPRIME: its residual sigma=(-1)^(h+c*e) is exactly the retained geometric-unit deformation formula. Here m1 is explicit by r0,c0, every unit1 at prime2, so m1 is also excluded. This does not claim the prime-power theorem as a new contribution. The exporter retains the finite necessary candidate list below the integer ceiling cutoff, without declaring those candidates bad or mathematically unresolved; other registered index certificates may already exclude them. The first four indices are covered for all m by these two deterministic classes alone. Exponent6 remains possible in this necessary range for n>4, leaving full all-index SFC open. No other Gamma-shape or novelty/assurance claim.',
    dependencies=['FIRLARGEM20261117','FIRPADICDEFORM20261116','FIRSFCMONOPRIME20261116']),
 dict(common,id='FIRALLEXPLOG20261117',
    title='The published short-interval route gives one logarithmically sparse exceptional set for every positive factorial exponent',
    scope='Let D_all(X)=#{1<=n<=X:there exists ANY integer m>=1 with F_(n,m) not squarefree over Q}. Then D_all(X)<<X/(log X)^(1/50), with an ABSOLUTE implied constant, in ordinary Gamma shape1. There is NO upper cap on m. For every n>=2^144, 2*log2(2n)<=n^(1/16): put u=n^(1/16)>=512, g(u)=u-32log2(u)-2, g512=222 and g_prime>7/8 using log2>1/2. The large-exponent criterion then makes the all-exponent and capped bad sets agree at every such n. Thus D(X)<=D_all(X)<=D(X)+2^144, where D is the registered simultaneous m<=n^(1/16) count. Transfer FIRSIMULTLOG20261117, which uses PUBLISHED MR2016 independently of MR II. The explicit2^144 is ONLY an overlap threshold and finite additive count; it is not an extracted analytic starting threshold, numerical implied constant or pointwise audit of every smaller n. This logarithmic route alone does not imply convergent reciprocal exceptional sum. Full all-index SFC, finiteness of exceptions and other Gamma shapes remain unclaimed.',
    dependencies=['FIRLARGEM20261117','FIRSIMULTLOG20261117']),
 dict(common,id='FIRALLEXPPOWER20261117',
    title='The source-qualified power saving and harmonic sparsity extend simultaneously to all positive factorial exponents',
    scope='SOURCE-QUALIFIED all-exponent transfer: for D_all(X)=#{n<=X:exists ANY positive integer m with F_(n,m) not squarefree over Q}, D_all(X)<<X^(1-10^-18), with ABSOLUTE implied constant and ONE exceptional set for ALL m>=1. In ordinary shape1, the large-exponent theorem and exact overlap n>=2^144 imply D_all(X)<=D(X)+2^144, then FIRSIMULTPOWER20261117 supplies the rate. Consequently sum_(n>=1,exists m>=1 with F_(n,m) nonsquarefree)1/n<infinity; the finite initial index set does not affect convergence. The inherited power input is DENSE MR II PREPRINT Theorem9.2(ii), with the published sparse-sieve warning retained, exclusion of Section10/Lemma10.3, no Theorem1.7 advertised exponent, and a scoped application/dependency audit rather than full-preprint validation. The published logarithmic all-exponent route remains independent and weaker. Numerical implied constants/analytic starting thresholds and an empty or finite exceptional set are NOT established;2^144 is an explicit overlap threshold only. A convergent reciprocal sum is not all-index SFC, which stays open, especially exponent6. No other Gamma-shape, novelty or independent/formal/external assurance upgrade.',
    dependencies=['FIRLARGEM20261117','FIRSIMULTPOWER20261117'])]
ids=[row['id'] for row in new]
annulus_receipt=json.loads((HERE/'saddle_annulus_verification.json').read_text())
assert annulus_receipt['status']=='pass' and annulus_receipt['all_roots_in_annulus_above_cutoff']
assert not annulus_receipt['eventual_simplicity_proved'] and not annulus_receipt['finite_prefix_verified']
for file,expected in annulus_receipt['files'].items():assert sha(HERE/file)==expected,('stale saddle-annulus receipt',file)
ANNULUS_PROOF=BASE+'EFFECTIVE_SADDLE_ANNULUS.tex'
ANNULUS_CHECKER=BASE+'verify_saddle_annulus.py'
analytic_ids=['FIRANALYTICANNULUS20261117']
new.append(dict(common,id=analytic_ids[0],
    title='An effective annulus excludes both saddle-collision values from every large-index exponent-six root region',
    canonical_source=ANNULUS_PROOF,checker=ANNULUS_CHECKER,
    artifact_hash='sha256:'+sha(HERE/'verify_saddle_annulus.py'),
    software_lock=[ANNULUS_PROOF,ANNULUS_CHECKER,BASE+'saddle_annulus_audit.json',BASE+'saddle_annulus_verification.json',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json'],
    scope='For ordinary Gamma shape1, P_n(z)=sum binom(n,k)*(6k)!*z^(n-k) has every scaled root |z/n^6|<2800 for ALL n>=1; for EVERY n>=10^6 all roots lie in the fixed annulus1<|z/n^6|<2800. The outer exclusion is a coefficient bound; the inner exclusion is a complete analytic contour proof, including the vertical shift, local Gaussian error, compact complement and infinite tail, uniform on CLOSED |w|<=1. At its unique saddle |s0-6|<10^-3, a=1-5/s0 has Re a>1/7 and |a|<1/5. I_n differs from exp(n*Phi0)*sqrt(2pi/(n*a)) by less than exp(n*RePhi0)/sqrt(n), giving |I_n|>4*exp(n*RePhi0)/sqrt(n). Thus BOTH critical values0 and3125 of w=s^5(6-s) are excluded from the large-index root region. On1<=|w|<=2800, all six saddles have1/2<|s|<7, pairwise distance>10^-13, curvature modulus>1/4000, distinct saddle phase derivative differences>10^-18, and endpoint/saddle differences>1/40000. The exact discriminant6^6*w^4*(3125-w) supplies separation. These bounds are on SADDLE DATA, not polynomial-root separation. The10^6 cutoff concerns ZERO EXCLUSION ONLY, not simplicity or finite-prefix verification. Global contributing-contour selection, amplitude/differentiated-remainder bounds, multiple-contribution transitions and all-root simplicity remain OPEN. The proof also gives an exact two-contribution P/Pprime exclusion criterion but does NOT verify it globally. No MR/density input, numerical root scan, eventual-simplicity theorem, all-index SFC closure or novelty/assurance upgrade is claimed.',
    dependencies=[]))

simplicity_receipt=json.loads((HERE/'exponent_six_simplicity_verification.json').read_text())
assert simplicity_receipt['status']=='pass' and simplicity_receipt['finite_prefix_complete']
assert simplicity_receipt['all_index_exponent_six_two_term_SFC'] and not simplicity_receipt['unrestricted_SFC_closed']
for file,expected in simplicity_receipt['files'].items():assert sha(HERE/file)==expected,('stale exponent-six simplicity receipt',file)
SIMPLE_PROOF=BASE+'EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex'
SIMPLE_CHECKER=BASE+'verify_exponent_six_simplicity.py'
simple_ids=['FIREVENTUALSIMPLE620261117','FIRALLINDEXSFC620261117']
simple_common=dict(common,canonical_source=SIMPLE_PROOF,checker=SIMPLE_CHECKER,
    artifact_hash='sha256:'+sha(HERE/'verify_exponent_six_simplicity.py'),
    software_lock=[SIMPLE_PROOF,SIMPLE_CHECKER,BASE+'exponent_six_prefix_certificates.json',BASE+'exponent_six_simplicity_audit.json',BASE+'exponent_six_simplicity_verification.json',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json',ANNULUS_PROOF,BASE+'MIXED_PRIME_RESIDUAL_WINDOWS.tex',BASE+'QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex'])
new.extend([
 dict(simple_common,id=simple_ids[0],
    title='Every ordinary exponent-six reversed moment polynomial has simple roots above the effective cutoff10000',
    scope='For EVERY integer n>=10000, P_n(z)=sum binom(n,k)*(6k)!*z^(n-k) has only SIMPLE complex roots; every scaled root w=z/n^6 lies in60<|w|<240. Ordinary Gamma shape1 only. This is a uniform ALL-root analytic theorem, not a bulk zero law or density extrapolation. Entropy/coefficient bounds exclude |w|>=240; split I_n exactly into endpoint K_n on[0,3] and tail J_n on[3,infinity). On CLOSED |w|<=250, a complete contour shift through the unique saddle |s0-6|<1/25 and all local/vertical/compact/infinite tail estimates give J_n=U_n*(1+epsilon), |epsilon|<12/sqrt(n). Analyticity plus radius10 Cauchy estimates give |Jprime/J-n*(6-s0)/(6w)|<1 on |w|<=240. On |w|>=50, K_n=(w^n/n)*(1+D), with explicit D,Dprime=O(n^-5), both<1/1000. Inner comparison excludes |w|<=60. At every remaining cancellation, |Jprime/J-Kprime/K|>n/288-2>0. The EXACT interval split supplies the global two-contribution decomposition; no assumption about formal saddle selection or neglected multi-contribution regions is used. N0=10000 is an actual EFFECTIVE SIMPLICITY cutoff. No other exponent/Gamma shape, unrestricted FC/SFC, novelty or assurance upgrade is claimed.',
    dependencies=['FIRANALYTICANNULUS20261117']),
 dict(simple_common,id=simple_ids[1],
    title='Strong Factorial holds for every A+B Y^6 at every starting index in ordinary Gamma shape one',
    scope='For EVERY integer n>=1, F_(n,6)(t)=sum binom(n,k)*(6k)!*t^k and its reversal P_n have SIMPLE complex roots. Hence arbitrary complex A+B*Y^6 satisfies the two-term Strong Factorial assertion at EVERY starting index. The analytic theorem handles all n>=10000. The COMPLETE finite prefix consists of n1 linear and9998 exact modular Bezout certificates for2<=n<10000, reconstructed coefficientwise. For2<=n<30 use full-degree reduction at prime211>6n. For30<=n<10000 use preceding prime q and d=n-q; the sieve verifies q>6max(1,d), with maxd35. The registered shifted-Newton reduction forces any rational square divisor to reduce degree-preservingly into F_d modq. Each certificate has U*F_d+V*F_dprime=G with degG<=1, while the real-place theorem forces a positive rational square-divisor degree to be even and at least2, contradiction. Pprime_(n+1)=(n+1)P_n identifies the adjacent obstruction exactly. One-term coefficients are immediate. This is full closure of THIS ordinary single-variable exponent-six two-term family, NOT unrestricted Strong Factorial, general higher-dimensional Factorial, other Gamma shapes, or an implication from density alone. No numerical root scan or novelty claim; internal proof/finite certificates only, with independent/formal/external assurance false.',
    dependencies=[simple_ids[0],'FIRALLMLOCAL20261116','FIRREALROOTPARITY20261116'])])

uniform_receipt=json.loads((HERE/'uniform_two_term_sfc_verification.json').read_text())
assert uniform_receipt['status']=='pass' and uniform_receipt['all_positive_exponents_all_indices_squarefree']
assert uniform_receipt['complete_finite_prefix'] and not uniform_receipt['unrestricted_SFC_closed']
for file,expected in uniform_receipt['files'].items():assert sha(HERE/file)==expected,('stale uniform two-term receipt',file)
UNIFORM_PROOF=BASE+'ALL_UNIVARIATE_TWO_TERM_SFC.tex'
UNIFORM_CHECKER=BASE+'verify_uniform_two_term_sfc.py'
uniform_ids=['FIRUNIFORMTWOSIMPLE20261117','FIRALLUNIVTWOSFC20261117']
uniform_common=dict(common,canonical_source=UNIFORM_PROOF,checker=UNIFORM_CHECKER,
    artifact_hash='sha256:'+sha(HERE/'verify_uniform_two_term_sfc.py'),
    software_lock=[UNIFORM_PROOF,UNIFORM_CHECKER,BASE+'uniform_two_term_sfc.py',BASE+'all_exponent_prefix_certificates.json',BASE+'uniform_two_term_sfc_audit.json',BASE+'uniform_two_term_sfc_verification.json',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json',SIMPLE_PROOF,BASE+'exponent_six_simplicity_verification.json',PROOF,BASE+'all_exponent_squarefree_verification.json',BASE+'PRIME_POWER_TWO_TERM_SFC.tex',BASE+'QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex',BASE+'MIXED_PRIME_RESIDUAL_WINDOWS.tex'])
new.extend([
 dict(uniform_common,id=uniform_ids[0],
    title='Uniform effective simplicity for every factorial exponent6<=m<=sqrt(n) above index10000',
    scope='For ALL integer n>=10000 and6<=m<=sqrt(n), every root of P_(n,m)(z)=sum binom(n,k)*(mk)!*z^(n-k) is SIMPLE and satisfies1/2<|z/(n^m*(m/e)^m)|<21/10. Ordinary Gamma shape1 only. Uniform entropy/coefficient bounds give the exterior enclosure. Split the scaled integral EXACTLY at s=m/2. On the FULL disk |v|<=11/5, a unique saddle |s0-m|<1/25 and complete local/vertical/compact/infinite contour bounds give J=U*(1+epsilon), |epsilon|<22/sqrt(m*n). The analytic relative error and radius1/10 Cauchy bounds yield |Jprime/J-n*(m-s0)/(m*v)|<3 on |v|<=21/10. On |v|>=2/5 the endpoint term K=(A_m*v)^n/n*(1+D) has |D|,|Dprime|<1/1000. An inner comparison excludes |v|<=1/2. At EVERY remaining cancellation the log-derivative gap is>n/3-4>0. This supplies a GLOBAL exact two-piece decomposition and ALL-root separation from multiplicity, without an unproved saddle selection, bulk-zero inference or analytic-number-theory input. The tail estimate itself is uniform for all m>=6; m<=sqrt(n) is used only in coefficient and inner-comparison bounds. No other Gamma shape, arbitrary support, unrestricted FC/SFC, priority or assurance claim.',
    dependencies=['FIREVENTUALSIMPLE620261117']),
 dict(uniform_common,id=uniform_ids[1],
    title='Every ordinary A+B Y^m satisfies Strong Factorial at every index for every positive exponent',
    scope='For ALL positive integers n,m, F_(n,m)(t)=sum binom(n,k)*(mk)!*t^k is squarefree over Q and C. Thus every nonzero ordinary one-variable A+B*Y^m with arbitrary complex coefficients satisfies Strong Factorial at EVERY starting index, for EVERY positive integer m. The effective uniform large-index cutoff is10000: the coefficient criterion2^m>=4n^2 handles its range, and otherwise m<2log2(2n)<sqrt(n), so the uniform analytic theorem applies for m>=6; m1..5 use the retained prime-power theorem. Below10000, m>=30 is already excluded by the coefficient criterion, prime powers/m1 and the all-index m6 theorem are reused. ALL92935 remaining necessary pairs at m10,12,14,15,18,20,21,22,24,26,28 have exact reconstructed modular Bezout certificates with gcd degree<=1. There are92900 preceding-prime kernels, maxdegree35, with p>m*max(1,n-p), and35 full-degree reductions at p>m*n, largest index126. Degree preservation plus the proved real-place even-degree restriction exclude every rational square divisor. Every selected prime and every certificate is checked. In THIS ordinary constant-plus-one-variable-power sequence, E_m(X)=D_all(X)=0: the exceptional set is EMPTY, independent of both MR routes. Older density/source-qualified results remain retained evidence but no longer leave a possible exception here. This is NOT unrestricted Strong Factorial, general Factorial, both-nonconstant monomial supports, multivariate monomial supports, other Gamma shapes, or a reduction from arbitrary polynomial supports. No numerical root scan, priority or independent/formal/external assurance claim.',
    dependencies=[uniform_ids[0],'FIRLARGEM20261117','FIRSFCMONOPRIME20261116','FIRALLINDEXSFC620261117','FIRALLMLOCAL20261116','FIRREALROOTPARITY20261116'])])

six_boundary=(' CURRENT scoped update: FIRALLINDEXSFC620261117 now closes ordinary Gamma-shape1 A+B*Y^6 at EVERY starting index, by an effective analytic simplicity cutoff10000 plus a COMPLETE exact finite prefix. Earlier statements leaving THIS exponent-six family open are superseded by that result; this record retains its own original local/density/general-shape scope. Failure of an older sufficient exporter does not imply a mathematical exponent-six obstruction. Unrestricted FC/SFC and other unproved shapes/supports remain open.')
def refresh_six_boundary(row):
    text=row['scope'].lower()
    if row['id'] not in simple_ids+uniform_ids and ('open' in text or 'unresolved' in text or 'unproved' in text) and ('exponent6' in text or 'exponent-six' in text or 'especially exponent' in text or 'exponents starting6' in text or 'non-prime-power' in text):
        if six_boundary not in row['scope']:row['scope']+=six_boundary
    return row
pade_receipt=json.loads((HERE/'pade_monomial_separation_verification.json').read_text())
assert pade_receipt['status']=='pass' and pade_receipt['all_n_a1_through6_proved']
assert pade_receipt['unbounded_indexed_strip_family_proved'] and not pade_receipt['all_a_proved']
for file,expected in pade_receipt['files'].items():assert sha(HERE/file)==expected,('stale Pade separation receipt',file)
PADE_PROOF=BASE+'PADE_ADJACENT_MONOMIALS.tex'
PADE_CHECKER=BASE+'verify_pade_monomial_separation.py'
pade_ids=['FIRPADESEPARATION20261117','FIRADJACENTMONOSFC20261117']
pade_common=dict(common,canonical_source=PADE_PROOF,checker=PADE_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_pade_monomial_separation.py'),
    software_lock=[PADE_PROOF,PADE_CHECKER,BASE+'pade_monomial_separation.py',BASE+'pade_shift_positive_certificates.json',BASE+'pade_monomial_separation_audit.json',BASE+'pade_monomial_separation_verification.json',BASE+'pade_separation_diagnostics.json',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json',BASE+'sources/rocks-2015-author-preprint.pdf',BASE+'sources/rocks-2015-author-preprint.txt'])
new.extend([
 dict(pade_common,id=pade_ids[0],
    title='Exact Pade denominator strips and small auxiliary matrix certificates for adjacent nonconstant monomial supports',
    scope='For every integer M>=0,N>=1, the exponential Pade denominator Q_(M,N)(z)=sum(-1)^j*(M+N-j)!*N!/((M+N)!*(N-j)!*j!)*z^j has all roots M-N+2<=Re z<=M+N, STRICT at both ends for N>=2. An explicit complex symmetric tridiagonal matrix has this characteristic polynomial and real diagonal Hermitian part; the proof rederives the coefficient recurrence. For f=X^a(A+B*X), normalize B_(n,a)=sum binom(n,k)*(an+1)_k*t^k; its reciprocal is Q_(an,n), so adjacent moments reduce to COPRIMALITY of DIFFERENT-parameter Pade denominators, not a fixed-shape Appell simplicity test. The exact cross determinant equals a nonzero constant times z^((a+1)*(n-1)+1)*H_(a,n), degree a-1 for a>=2 and H1=2. A displayed rising/falling-factorial formula determines every coefficient. A positive Schwarz-tridiagonal reconstruction of H_(a,n)(z+(a-1)*n+2) excludes common zeros; exact nonnegative-coefficient numerator/denominator certificates in u=n-2 prove positivity for ALL n>=2 at any supplied a where they succeed. A second generic necessary obstruction is det(A+zB)=0 in dimension a-1: A upper bidiagonal with diag a*(n-1)+j+n and upper -(a*(n-1)+j+1); B lower bidiagonal with diag1 except last1/(a+1) and subdiag-1. It follows from exact generating-function recurrences, with a nonzero solution vector. The uniform spectral bound Re eigen(-B^-1*A)<(a-1)*n+2 is sufficient but UNPROVED. For a2..6 the pencil determinant and H agree by symbolic identities. No all-a, unrestricted SFC/FC, other Gamma shape, priority or assurance claim.',
    dependencies=[]),
 dict(pade_common,id=pade_ids[1],
    title='Strong Factorial holds at every index for X^a(A+B X) at a1 through6, with an unbounded indexed strip family',
    scope='For ordinary Gamma shape1 and arbitrary complex A,B giving a nonzero polynomial, X^a(A+B*X) satisfies Strong Factorial at EVERY starting index for each a=1,2,3,4,5,6. a0 was already covered by constant-plus-power closure. The retained Rocks2015 author draft Proposition3.14 already contains a0,1,2; this proof is independently algebraic and makes no priority claim for the a3..6 extension or final2017 scope. For a1 the Pade cross determinant auxiliary is constant2. For a2..6, exact symbolic Schwarz parameters reconstruct H_(a,n)(z+(a-1)*n+2); numerator/denominator polynomials in u=n-2 have nonnegative integer coefficients and positive constants, giving STRICT positivity for ALL n>=2, no finite-index extrapolation. A dissipative real tridiagonal matrix then places every auxiliary root strictly left of the Pade denominator strip. For ALL integers n>=2,a>=1 with a>=2*n-3, the strip bounds alone exclude simultaneous vanishing at orders n-1,n, giving an UNBOUNDED INDEXED family, not full all-starting-index closure for arbitrary fixed a. All coefficient-degenerate monomial cases are immediate. Generic finite shifted certificates also apply at supplied pairs without all-parameter extrapolation. Blanket real-rooted H and blanket unshifted Hurwitz stability are FALSE at (a,n)=(3,2),(5,2); they are not needed. General a>=7 beyond the indexed range, gaps d>1 in X^a(A+B*X^d), arbitrary monomial pairs, other Gamma shapes and unrestricted FC/SFC remain outside the proved scope. No numerical root scan, independent/formal/external assurance or novelty claim.',
    dependencies=[pade_ids[0]])])

parabola_receipt=json.loads((HERE/'pade_parabola_sfc_verification.json').read_text())
assert parabola_receipt['status']=='pass' and parabola_receipt['all_a_all_n_adjacent_monomial_SFC']
assert parabola_receipt['all_positive_Gamma_shapes'] and parabola_receipt['generalized_Pade_coprimality']
assert not any(parabola_receipt[k] for k in ('finite_scan_extrapolation','vertical_bound_required','unrestricted_SFC_closed','general_FC_closed','independent_replay','formal_verification','external_review','novelty_claim'))
for file,expected in parabola_receipt['files'].items():assert sha(HERE/file)==expected,('stale uniform Pade parabola receipt',file)
PARABOLA_PROOF=BASE+'UNIFORM_PADE_PARABOLA_SFC.tex'
PARABOLA_CHECKER=BASE+'verify_pade_parabola_sfc.py'
parabola_ids=['FIRPADEPARABOLA20261117','FIRUNIFORMADJACENTSFC20261117']
parabola_common=dict(common,canonical_source=PARABOLA_PROOF,checker=PARABOLA_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_pade_parabola_sfc.py'),
    software_lock=[PARABOLA_PROOF,PARABOLA_CHECKER,BASE+'pade_parabola_sfc.py',BASE+'pade_parabola_sfc_audit.json',BASE+'pade_parabola_sfc_verification.json',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json',BASE+'sources/saff-varga-1976-sharpness.pdf'])
new.extend([
 dict(parabola_common,id=parabola_ids[0],
    title='A shared parabola separates every adjacent-degree generalized Pade pair with a positive integer parameter shift',
    scope='For ALL N>=2, real l>-1 and integer a>=1, R_N(l+a,z) and R_(N-1)(l,z) are COPRIME over C, where R_N(t,z)=sum(-1)^k*binom(N,k)*(t+1)_k*z^(N-k). For EVERY real M>-1 all R_N(M,z) roots lie STRICTLY inside y^2<4*N*(x+N). A tridiagonal characteristic matrix has Re z=S-2*tau and (Im z)^2<=4*tau*(S-tau), with S=M+N>0,0<=tau<=N-1; the exact positive gap is(N-tau)*S+(N-tau)^2. A common root forces a zero of an explicit continuant of degree a-1. After normalization its recurrence has b_j*(1-b_(j-1)/c_j)>=N; a fully proved invariant Re mu_j<1 makes that continuant NONZERO in the same open parabola. The a1 case is direct parameter propagation. No unbounded vertical inequality is assumed. The elementary mechanism is compared with Saff-Varga1976 Theorem1.1 but proved self-contained; no whole-paper or novelty audit. No contour, numerical-root, density, parity, prime-gap, irreducibility or finite-prefix input. Arbitrary gaps and unrestricted FC/SFC remain outside scope; assurance flags false.',
    dependencies=[]),
 dict(parabola_common,id=parabola_ids[1],
    title='Every adjacent-monomial X^a(A+B X) satisfies Strong Factorial at every index and every positive Gamma shape',
    scope='For EVERY integer a>=0, positive real Gamma shape beta, complex A,B giving a nonzero polynomial, and terminal n>=1, the moments of f=X^a(A+B*X) at orders n-1,n CANNOT BOTH VANISH. Thus this family satisfies SFC at ALL starts, including ordinary factorials. For a>=1 use R_n(beta+a*n-1,z)=z^n*B_(n,a,beta)(-1/z), B=sum binom(n,k)*(beta+a*n)_k*t^k, with l=beta+a*(n-1)-1>-1,M=l+a. For a0 a fixed-shape recurrence propagates a common zero backward to R0=1. n1 has order-zero moment1; degenerate coefficients give monomials. Also product_j X_j^a_j*(A+B*X1) under independent positive Gamma shapes has the same assertion, because other variables contribute nonzero factors. This corollary requires monomials differing by one power in one variable. The earlier a>=7 gate is now CLOSED DIRECTLY for all a; the stronger vertical spectral target can remain unproved. No arbitrary gaps d>1, arbitrary monomial pairs, general FC or unrestricted SFC closure follows. Fully algebraic universal proof, not finite extrapolation; no priority or independent/formal/external assurance claim.',
    dependencies=[parabola_ids[0]])])
pade_current_boundary=(' CURRENT uniform adjacent-support update: FIRUNIFORMADJACENTSFC20261117 now proves X^a(A+B*X) SFC at EVERY start for EVERY integer a>=0 and EVERY positive Gamma shape. FIRPADEPARABOLA20261117 gives a self-contained shared-parabola separation, independent of contour/density/finite-prefix arguments. The earlier a>=7 support-family gate is closed; the stronger uniform VERTICAL inequality remains unproved and unnecessary. This record retains its earlier finite-a/strip scope. Arbitrary gaps d>1 and unrestricted FC/SFC remain open.')
for row in new:
    if row['id'] in pade_ids and pade_current_boundary not in row['scope']:row['scope']+=pade_current_boundary

dyadic_receipt=json.loads((HERE/'dyadic_binomial_sfc_verification.json').read_text())
assert dyadic_receipt['status']=='pass' and dyadic_receipt['dyadic_equal_lower_exponent_and_gap_all_indices']
assert dyadic_receipt['integer_Gamma_shapes_through_2q'] and dyadic_receipt['arbitrary_rational_unit_sequence_transfer']
assert dyadic_receipt['same_variable_scaling_at_consecutive_orders'] and dyadic_receipt['exact_resultant_valuation_proved']
assert not any(dyadic_receipt[k] for k in ('finite_parameter_extrapolation','other_lower_exponents_or_gaps_proved','all_positive_Gamma_shapes','general_FC_closed','unrestricted_SFC_closed','independent_replay','formal_verification','external_review','novelty_claim'))
for file,expected in dyadic_receipt['files'].items():assert sha(HERE/file)==expected,('stale dyadic binomial receipt',file)
DYADIC_PROOF=BASE+'DYADIC_BINOMIAL_SFC.tex'
DYADIC_CHECKER=BASE+'verify_dyadic_binomial_sfc.py'
dyadic_ids=['FIRBESSELUNITTRANSFER20261117','FIRDYADICBINSFC20261117','FIRLOWERONESIMPLE20261117','FIRCUBICPRIMECONF20261117']
dyadic_common=dict(common,canonical_source=DYADIC_PROOF,checker=DYADIC_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_dyadic_binomial_sfc.py'),
    software_lock=[DYADIC_PROOF,DYADIC_CHECKER,BASE+'dyadic_binomial_sfc.py',BASE+'dyadic_binomial_sfc_audit.json',BASE+'dyadic_binomial_sfc_verification.json',BASE+'dyadic_binomial_frontier_diagnostics.json',BASE+'prime_power_windows.py',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json',BASE+'sources/rocks-2015-author-preprint.pdf',BASE+'sources/rocks-2015-author-preprint.txt'])
new.extend([
 dict(dyadic_common,id=dyadic_ids[0],
    title='Bessel reductions transfer adjacent coprimality uniformly through arbitrary rational 2-adic unit sequences',
    scope='For EVERY rational 2-adic unit sequence u_L, h=0 or1, and n>=1, the two polynomials F_n(t)=sum binom(n,k)*(n+k+h)!*u_(n+k)*t^k and F_(n-1)(t) are COPRIME over C. Scale t=z/2 and divide by n!, identically at consecutive orders. Normalized coefficients belong to Z_(2), reduce to Y_n=sum (n+k)!/(2^k*k!*(n-k)!)*z^k for h0, or D_n=(n+1)*Y_n+z*Y_nprime for h1, and retain degree n because leading coefficients are odd. Both reductions satisfy C_n=z*C_(n-1)+C_(n-2), with initial pairs(1,1+z) or(1,z). Backward Euclidean propagation gives gcd1 for EVERY index, hence the rational resultant is a 2-adic unit and nonzero over C. Exact Bessel coefficient recurrences and integer units are proved directly. This does NOT transfer arbitrary characteristic-zero coprimality through every factorial substitution: full-degree coprime reduction is essential. h2 already has normalized reductions z^2 at both n2,n3. No finite-parameter extrapolation, contour/density/parabola input, external theorem import or priority/assurance claim.',
    dependencies=[]),
 dict(dyadic_common,id=dyadic_ids[1],
    title='Strong Factorial holds at every index for X^(2^r)(A+B X^(2^r)) and integer Gamma shapes through twice the gap',
    scope='For EVERY q=2^r,r>=0, EVERY integer Gamma shape1<=beta<=2q, EVERY complex A,B giving nonzero f=X^q(A+B*X^q), and EVERY terminal n>=1, the moments at n-1,n CANNOT BOTH VANISH. Ordinary factorial beta1 is included. For r>=1 the gap q>=2 and both terms are nonconstant, outside the adjacent-support family. Write beta-1=q*h+delta, h0 or1,0<=delta<q. The integers W_(L,delta)=(qL+delta)!/(2^((q-1)*L)*L!*delta!) are ODD, proved by dyadic doubling and the offset binomial valuation identity. The moment polynomial M_n(t)=sum binom(n,k)*(beta)_(q*(n+k))*t^k has normalized G_n=h!*W_(h,delta)*M_n(2^(-q)*z)/(n!*2^((q-1)*n)) with integral coefficients W_(n+k+h,delta)*(n+k+h)!/(2^k*k!*(n-k)!). Its full-degree reduction is Bessel Y_n or D_n, so Res(G_n,G_(n-1)) is odd. Consequently the exact valuation is v2 Res(M_n,M_(n-1))=(3q-2)*n*(n-1)+(n-1)*v2(n!)+n*v2((n-1)!), equivalently3q*n*(n-1)-(n-1)*s2(n)-n*s2(n-1), INDEPENDENT of beta in the proved range. Nonzero monomial cases are immediate. Monomial factors in OTHER independent variables preserve the assertion, with any positive Gamma shapes there. No other lower exponent/gap, all-positive-shape, general FC or unrestricted SFC closure. First-residue failures at gap2,a1,n3 and h2 are limitations of that certificate, not counterexamples. No novelty or independent/formal/external assurance upgrade.',
    dependencies=[dyadic_ids[0]])])
assert dyadic_receipt['conditional_cubic_prime_confinement_proved'] and dyadic_receipt['cubic_residual_transfer_proved']
assert dyadic_receipt['lower_one_multiple_root_equivalence_proved'] and dyadic_receipt['cubic_quarter_gap_square_degree_bound_proved']
assert not dyadic_receipt['cubic_all_index_SFC_closed']
new.extend([
 dict(dyadic_common,id=dyadic_ids[2],
    title='For every gap, a lower exponent of one turns the consecutive-moment obstruction into a multiple root',
    scope='For EVERY positive integer n,e, define A_(n,e)=sum binom(n,k)*(n+e*k)!*t^k/n! and T_(n,e)=sum binom(n,k)*(n)_(e*k)*t^k. Exact coefficient identities give A_(n,e)=T+(e/n)*theta*T and A_(n-1,e)=T-(1/n)*theta*T, theta=t*d/dt. Thus the two moment polynomials have a common root EXACTLY when T_(n,e) has a multiple root: T(0)=1 and the 2x2 coefficient matrix has determinant-(e+1)/n. Consequently ordinary X*(A+B*X^e) satisfies SFC at every start IFF every T_(n,e) is squarefree. The Gamma shape in T is n and varies with n; no fixed-shape Appell identity is assumed. This supplies a scalar structure for lower exponent ONE beyond the underdetermined order-e+1 parameter recurrence. It proves an equivalence, NOT the simplicity assertion or SFC for all gaps. More general lower exponents, unrestricted FC/SFC and assurance/priority upgrades remain unclaimed.',
    dependencies=[]),
 dict(dyadic_common,id=dyadic_ids[3],
    title='A prime above the index confines cubic common factors and square factors to explicit small residual kernels',
    scope='For ordinary f=X*(A+B*X^2), let C_n=sum binom(n,k)*(n+2k)!*t^k and J_n=gcd_Q(C_n,C_(n-1)). For EVERY n>=6 and prime n<q<=3*(n-1)/2, put d=q-n. Exact Newton polygons have disjoint nonzero slopes, including the even-n,d1 collinearity case. Hence all J_n roots are q-units and deg J_n<=floor((d-1)/2), sharpened by the real place to2*floor((d-1)/4). The horizontal reductions are H_d and H_(d+1), H_d=sum_(k<=floor((d-1)/2))(-1)^k*(d-1+k)!/(k!*(d-1-2k)!)*t^k. J_n reduces degree-preservingly to their gcd; residual gcd degree<=1 certifies coprimality. STRONGER via the lower-one equivalence: T_(n,2) has horizontal width floor(d/2) and primitive positive-rise edges, so any RATIONAL square divisor g^2 has deg g<=floor(d/4), and positive even degree gives deg g<=2*floor(d/8). In particular d<=7 certifies this cubic moment pair. Its exact residual is K_d=sum_(k<=floor(d/2))(-1)^k*d*(d+k-1)!/(k!*(d-2k)!)*t^k. Degree-preserving reduction sends g^2 into K_d modq. K_d=(H_d+2H_(d+1))/3 and theta*K_d=d*(H_(d+1)-H_d)/3, with q>=3d+3. Full residual square-factor degrees may be intersected with the even quarter-gap bound. This is a conditional all-parameter theorem, not an all-index prime-coverage theorem, an imported prime-existence assertion, or all-index SFC for the cubic. Finite exact windows only check bookkeeping. General FC/unrestricted SFC and all assurance upgrades remain open/unclaimed.',
    dependencies=[dyadic_ids[2]])])
frontier_boundary=(' CURRENT first-gap reduction: FIRLOWERONESIMPLE20261117 identifies X*(A+B*X^e) SFC exactly with squarefreeness of T_(n,e)=sum binom(n,k)*(n)_(e*k)*t^k. FIRCUBICPRIMECONF20261117 confines a rational square divisor of T_(n,2) to even degree<=2*floor((q-n)/8) for any prime n<q<=3*(n-1)/2, and transfers it degree-preservingly into the explicit residual K_(q-n). Uniform prime/kernel coverage remains unproved; this is not cubic or unrestricted all-index closure.')

dyadic_boundary=(' CURRENT nonadjacent dyadic family: FIRDYADICBINSFC20261117 proves X^q(A+B*X^q) SFC for EVERY q=2^r, EVERY start and EVERY integer Gamma shape1<=beta<=2q. FIRBESSELUNITTRANSFER20261117 gives the reusable arbitrary-unit transfer; a normalized resultant is odd, and an exact raw resultant valuation follows. This is a genuine all-index nonadjacent family for q>=2, not all gaps/lower exponents or unrestricted FC/SFC.')

tilted_receipt=json.loads((HERE/'tilted_cubic_frontier_verification.json').read_text())
assert tilted_receipt['status']=='pass' and tilted_receipt['uniform_shifted_polygon_and_residuals_proved']
assert tilted_receipt['all_prime_power_indices_irreducible'] and tilted_receipt['ten_uniform_shifted_index_families_proved']
assert tilted_receipt['exact_analytic_fold_data_proved']
assert not any(tilted_receipt[k] for k in ('finite_extrapolation_used','uniform_prime_coverage_proved','effective_simplicity_cutoff_proved','all_index_cubic_SFC_closed','general_FC_closed','unrestricted_SFC_closed','independent_replay','formal_verification','external_review','novelty_claim'))
for file,expected in tilted_receipt['files'].items():assert sha(HERE/file)==expected,('stale tilted cubic frontier receipt',file)
TILTED_PROOF=BASE+'TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex'
TILTED_CHECKER=BASE+'verify_tilted_cubic_frontier.py'
tilted_ids=['FIRTILTEDSHIFT20261117','FIRTILTEDPRIMEPOWER20261117','FIRCUBICFOLDS20261117']
tilted_common=dict(common,canonical_source=TILTED_PROOF,checker=TILTED_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_tilted_cubic_frontier.py'),
    software_lock=[TILTED_PROOF,TILTED_CHECKER,BASE+'tilted_cubic_frontier.py',BASE+'tilted_cubic_frontier_audit.json',BASE+'tilted_cubic_frontier_verification.json',BASE+'tilted_cubic_frontier_diagnostics.json',BASE+'prime_power_windows.py',BASE+'newton_slope_coverage.py',DYADIC_PROOF,BASE+'dyadic_binomial_sfc_verification.json',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json'])
new.extend([
 dict(tilted_common,id=tilted_ids[0],
    title='A uniform shifted Newton polygon has two distinct residual types for the tilted cubic moments',
    scope='For ALL r>=1,a>=1,i>=0 and prime p>3max(a,i), put P=p^r,n=aP+i,S=(P-1)/(p-1). The ordinary tilted T_n=sum binom(n,j)*(n)_(2j)*t^j has a horizontal edge of length i and an edge of width aP,slope2S/P with reduced denominator P. A levelwise Kummer/Legendre inequality proves the entire polygon including i0 negative-one factorial endpoints. The horizontal residual is T_i. The nonhorizontal residual is T_a=sum binom(a,v)*(a)_(2v)*Y^v if i0, but sum binom(a,v)*(a+1)_(2v)*Y^v for i>0. Wilson digit units give both formulas independently of r. A rational square divisor has degree u+vP,0<=2u<=i,0<=2v<=a, and even degree. No squarefreeness of those residuals is assumed, and general multipliers/offsets remain conditional. No extrapolation from finite r,n, other Gamma-shape theorem, priority or independent/formal/external assurance claim.',
    dependencies=[]),
 dict(tilted_common,id=tilted_ids[1],
    title='Every prime-power tilted cubic moment is irreducible, with ten additional uniform shifted index families',
    scope='For EVERY prime p and EVERY r>=1, T_(p^r,2)=sum binom(p^r,j)*(p^r)_(2j)*t^j is IRREDUCIBLE over Q_p and Q. For p>=5 the shifted theorem gives a primitive single edge. For p2 the exact single edge has widthP and height2P-1; for p3 it has heightP-1. Separate digit-sum and binomial valuation arguments put all interior coefficients above these chords, proving all r. In addition, for p>3max(a,i),r>=1, T_(a*p^r+i,2) is squarefree for (a,i)=(1,0),(1,1),(1,2),(1,3),(2,0),(2,1),(2,2),(3,0),(3,1),(4,0). Square occupancy plus even-degree parity removes all patterns except full rational squares in two cases; Bertrand supplies a valuation-one prime in leading(n)_(2n)=(3n-1)!/(n-1)!, excluding those. Thus ordinary X*(A+B*X^2) has no adjacent zero moments at EVERY terminal index in these unbounded families, for arbitrary complex coefficients. This is INDEXED coverage, not all-starting-index SFC for the entire cubic family. No general n, all multipliers/offsets, Gamma-shape, full FC/SFC or novelty/assurance upgrade.',
    dependencies=[tilted_ids[0],'FIRLOWERONESIMPLE20261117']),
 dict(tilted_common,id=tilted_ids[2],
    title='The exact scaled cubic phase has two ordinary folds with explicit nonzero derivatives',
    scope='For T_n=sum binom(n,j)*(n)_(2j)*t^j, T_n(w/n^2)=n^n/Gamma(n)*integral s^(-1)*exp(n*Phi)ds, Phi=log s+log(1+w*s^2)-s in local nonzero charts. The saddle cubic is F=w*s^3-3w*s^2+s-1, discriminant-4w*(1-9w+27w^2). Finite collisions are s=(3+-i*sqrt3)/2,w=1/6+-i*sqrt3/18, with |w|=1/(3sqrt3); w0 is a degree drop with unique finite saddle1, not a finite fold. At either collision Phi_s=Phi_ss=0, Phi_sss=-3*(s-1)^2/s^3 has modulus1/sqrt3 and Phi_sw=2*s/(1+w*s^2)^2 has modulus3sqrt3/2. The third saddle is3-2s. These are exact fold DATA only. No contour contributions, Airy estimates, global root enclosure, zero approach theorem, effective simplicity cutoff or all-index cubic SFC follows. Numerical root evidence is not used; priority and assurance flags false.',
    dependencies=[])])
cubic_receipt=json.loads((HERE/'cubic_adjoint_stability_verification.json').read_text())
assert cubic_receipt['status']=='pass' and cubic_receipt['all_index_cubic_SFC_closed']
assert cubic_receipt['all_positive_Gamma_shapes'] and not cubic_receipt['finite_prefix_required']
assert not any(cubic_receipt[k] for k in ('finite_scan_extrapolation','all_gaps_or_lower_exponents_proved','general_FC_closed','unrestricted_SFC_closed','independent_replay','formal_verification','external_review','novelty_claim'))
for file,expected in cubic_receipt['files'].items():assert sha(HERE/file)==expected,('stale cubic adjoint receipt',file)
CUBIC_PROOF=BASE+'CUBIC_ADJOINT_STABILITY_SFC.tex'
CUBIC_CHECKER=BASE+'verify_cubic_adjoint_stability.py'
cubic_ids=['FIRCUBICGAMMASIMPLE20261117','FIRCUBICALLSFC20261117']
cubic_common=dict(common,canonical_source=CUBIC_PROOF,checker=CUBIC_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_cubic_adjoint_stability.py'),
    software_lock=[CUBIC_PROOF,CUBIC_CHECKER,BASE+'cubic_adjoint_stability.py',BASE+'cubic_adjoint_stability_audit.json',BASE+'cubic_adjoint_stability_verification.json',BASE+'sources/wagner-multivariate-stable-polynomials-2011.pdf',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json'])
new.extend([
 dict(cubic_common,id=cubic_ids[0],
    title='An explicit adjoint and half-plane convolution prove tilted cubic simplicity uniformly in the index and Gamma parameter',
    scope='For EVERY n>=1 and EVERY real gamma>0 with gamma>=n-2, T_(n,gamma)=sum binom(n,k)*(gamma)_(2k)*t^k is squarefree over C. Reverse it as R(z)=z^(2n)*T(-z^-2), M=gamma+2n. Its monic third-order ODE has explicit adjoint solutions exp(-z)*z^(1-M)*B(z) and exp(z)*z^(1-M)*B(-z), where B(z)=E(z+U+V/2)^n for independent Gamma shapes gamma and n+1. A tridiagonal Gamma moment matrix bounds root real parts; the classical Grace-Walsh-Szego half-plane polarization theorem transfers this to B via Walsh additive convolution, giving Re roots<n-gamma-2 for n>=2. The exact Lagrange concomitants are constant exponential/power times rational functions, hence identically zero by their rational logarithmic derivatives. A repeated R root has multiplicity exactly2 and would force B(z)=B(-z)=0, contradicting strict half-plane separation. Noninteger powers are local gauges and introduce no rationality assumption. No asymptotic cutoff, finite-prefix extrapolation, parity, prime coverage or numerical root location is used. gamma below n-2, general gaps/supports, unrestricted FC/SFC, novelty and independent/formal/external assurance are unclaimed.',
    dependencies=['external: Grace-Walsh-Szego half-plane polarization; Wagner Theorem4.1, cached author PDF']),
 dict(cubic_common,id=cubic_ids[1],
    title='Strong Factorial holds for every cubic pair X(A+B X^2) at every starting index and every positive Gamma shape',
    scope='For EVERY beta>0, EVERY nonzero complex f=X*(A+B*X^2), and EVERY terminal n>=1, its Gamma-beta moments at n-1 and n cannot both vanish. Put gamma=n+beta-1 and T=T_(n,gamma). Exact coefficient identities give T_(n,gamma+1)=T+(2/gamma)*theta*T and T_(n-1,gamma)=T-(1/n)*theta*T; determinant-(gamma+2n)/(gamma*n) is nonzero. Since gamma>n-1, the uniform adjoint simplicity theorem excludes all adjacent moment failures. Nonzero monomial cases are immediate. Independent-variable monomial factors preserve the result for arbitrary positive Gamma shapes there. Ordinary beta1 is included. This closes the ENTIRE formerly indexed cubic frontier without a cutoff or finite prefix; earlier prime-power irreducibility and fold data remain valid independent results. It does not close arbitrary lower exponents/gaps, unrestricted SFC or general FC. No novelty, independent replay, formal verification or external review claim.',
    dependencies=[cubic_ids[0]])])
cubic_boundary=(' CURRENT whole-cubic closure: FIRCUBICGAMMASIMPLE20261117 proves T_(n,gamma) squarefree for ALL n>=1,gamma>0,gamma>=n-2, by explicit adjoint solutions and classical GWS half-plane convolution. FIRCUBICALLSFC20261117 therefore closes X*(A+B*X^2) SFC at EVERY start for EVERY positive Gamma shape. No arithmetic coverage, fold theorem, asymptotic cutoff or finite prefix is needed. Earlier cubic-open statements in this record describe only the preceding method stage and are superseded for this support. General gaps/lower exponents and unrestricted FC/SFC remain open; all assurance flags false.')
for row in new:
    if row['id'] in tilted_ids+dyadic_ids and cubic_boundary not in row['scope']:row['scope']+=cubic_boundary

gap_receipt=json.loads((HERE/'higher_gap_adjoint_rank_verification.json').read_text())
assert gap_receipt['status']=='pass' and gap_receipt['uniform_higher_gap_adjoint_proved']
assert gap_receipt['exact_repeated_root_rank_equivalence_proved'] and gap_receipt['full_Wronskian_reconstruction_proved']
assert gap_receipt['all_positive_shape_lower_one_bridge_proved'] and gap_receipt['gap_three_closed_sector_sufficiency_proved']
assert gap_receipt['Hurwitz_shortcut_refuted'] and gap_receipt['all_gap_initial_four_orders_all_positive_shapes_proved']
assert not any(gap_receipt[k] for k in ('gap_three_sector_zero_freeness_proved','higher_gap_all_index_SFC_closed','general_FC_closed','unrestricted_SFC_closed','finite_extrapolation','independent_replay','formal_verification','external_review','novelty_claim'))
for file,expected in gap_receipt['files'].items():assert sha(HERE/file)==expected,('stale higher-gap adjoint receipt',file)
GAP_PROOF=BASE+'HIGHER_GAP_ADJOINT_RANK.tex'
GAP_CHECKER=BASE+'verify_higher_gap_adjoint_rank.py'
gap_ids=['FIRGENERALGAPADJOINT20261117','FIRGENERALGAPRANK20261117','FIRALLGAPINITIAL20261117']
gap_common=dict(common,canonical_source=GAP_PROOF,checker=GAP_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_higher_gap_adjoint_rank.py'),
    software_lock=[GAP_PROOF,GAP_CHECKER,BASE+'higher_gap_adjoint_rank.py',BASE+'higher_gap_adjoint_rank_audit.json',BASE+'higher_gap_adjoint_rank_verification.json',BASE+'higher_gap_adjoint_rank_diagnostics.json',BASE+'cubic_adjoint_stability.py',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json'])
new.extend([
 dict(gap_common,id=gap_ids[0],
    title='Every monomial gap has explicit polynomial adjoint solutions and a uniform differential recurrence',
    scope='For EVERY d>=2,n>=1,gamma>0, put T=sum binom(n,k)*(gamma)_(d*k)*t^k, R=z^(dn)*T((-1)^(d+1)*z^-d), M=gamma+dn. Its monic order-d+1 Euler equation has d independent adjoint solutions exp(-omega*z)*z^(1-M)*B_(n,gamma,d)(omega*z), omega^d=1. The monic degree-n auxiliary has generating coefficients (1-x)^(-gamma)*[(1-(1-x)^d)/(d*x)]^(-n-1), equivalently exponential generating function (1-d*t)^(-(gamma+d-1)/d)*exp(z*(1-(1-d*t)^(1/d))). Formal residue/Lagrange substitution gives the exact all-n recurrence B_(n+1)=(z+dn+gamma+d-1)*B_n-z*B_nprime. A divergence-form adjoint calculation proves all gaps without finite extrapolation. The noninteger power gauge is local; no rationality or convergence away from the formal origin is assumed. This constructs the higher-gap structure, NOT a simplicity or SFC theorem for d>=3. No priority or independent/formal/external assurance upgrade.',
    dependencies=[]),
 dict(gap_common,id=gap_ids[1],
    title='Repeated roots are exactly a rank defect of the polynomial adjoint jets, with an explicit gap-three sector target',
    scope='For ALL parameters of FIRGENERALGAPADJOINT20261117 and all z0 nonzero, R has a repeated root IFF the d by d-1 matrix J_(j,k)=omega_j^k*((1-D)^k B)(omega_j*z0), k0..d-2, has rank<d-1. The d concomitants vanish by rational logarithmic derivatives, have full row rank on the complete original jet, and give the stated equivalence by restriction to derivatives Rsecond..R^(d). The square full matrix with k0..d-1 reconstructs R exactly, normalized by prod omega_j^n times their Vandermonde. Consecutive Gamma-beta moments of X*(A+B*X^d) correspond to a multiple root at gamma=n+beta-1, for EVERY beta>0. At gap3 the rank obstruction is a common three-point rotation orbit of a real degree-2n Wronskian W; zero-freeness on the CLOSED sector |arg z|<=pi/3 would suffice. Strict interlacing of any two cyclic sections W0,W1,W2 also suffices. Both UNIFORM conditions remain UNPROVED. Stronger Hurwitz stability of W is FALSE at n4,gamma4: exact seventh Hurwitz determinant -169684039754781133897728000000. The three sections at that index have gcd1, so it is not an SFC counterexample. Exact negative-root/interlacing controls at n=gamma4,6,8 are finite evidence only. Degenerate negative/complex Gamma parameters test actual repeated roots outside positive-shape scope. Higher-gap all-index SFC, general FC/unrestricted SFC and novelty/assurance upgrades remain unclaimed.',
    dependencies=[gap_ids[0]]),
 dict(gap_common,id=gap_ids[2],
    title='Every positive Gamma power moment is simple through degree four, uniformly in every monomial gap',
    scope='For EVERY d>=1,gamma>0 and n1..4, T_(n,gamma,d)=sum binom(n,k)*(gamma)_(d*k)*t^k is squarefree over C. For X=Y^d with Y Gamma-gamma, P_n(z)=E(z+X)^n has no real roots at even n and only simple real roots at odd n; conjugate multiplicity excludes repeats at n<=3. A repeated root at n4 would force the real monic quartic to be a square of a real quadratic and hence force zero third central moment. But a=(gamma)_d,b=(gamma+d)_d,c=(gamma+2d)_d satisfy b>a>0,c>=2b-a by increasing convexity of the rising factorial polynomial, so kappa3=a*b*c-3*a^2*b+2*a^3>=2*a*(b-a)^2>0. Thus every X*(A+B*X^d) has no consecutive zero moments through terminal order4, for ALL gaps and ALL positive Gamma shapes and arbitrary nonzero complex coefficients. This is an all-gap initial block, not all-index SFC. The proof is elementary and uniform in d,gamma; finite rising-factorial checks only replay it. General FC/unrestricted SFC, novelty and independent/formal/external assurance upgrades remain unclaimed.',
    dependencies=[gap_ids[1]])])
gap_boundary=(' CURRENT higher-gap structure: FIRGENERALGAPADJOINT20261117 constructs all d>=2 explicit adjoints via a degree-n polynomial recurrence; FIRGENERALGAPRANK20261117 gives an EXACT repeated-root rank equivalence and determinant reconstruction, with the all-positive-shape lower-one moment bridge. For gap3 a closed pi/3 zero-free sector or strict coprimality/interlacing of cyclic Wronskian sections would suffice; neither is proved uniformly. Full Hurwitz stability of that Wronskian is refuted exactly already at n=gamma4.  FIRALLGAPINITIAL20261117 proves all d>=1,gamma>0 at n1..4 by a strictly positive third central moment; only this initial block is closed uniformly. No higher-gap/unrestricted SFC or general FC closure follows.')

wr_receipt=json.loads((HERE/'gap_three_wronskian_verification.json').read_text())
assert wr_receipt['status']=='pass' and wr_receipt['uniform_Wronskian_shape_shift_and_recurrence_proved']
assert wr_receipt['bivariate_coefficient_positivity_proved'] and wr_receipt['no_common_real_section_root_proved']
assert wr_receipt['one_real_rooted_section_suffices_proved'] and wr_receipt['fifth_order_all_positive_Gamma_shapes_proved']
assert wr_receipt['separate_input_interlacing_induction_refuted']
assert not any(wr_receipt[k] for k in ('uniform_section_real_rootedness_proved','all_index_gap_three_SFC_closed','general_FC_closed','unrestricted_SFC_closed','finite_extrapolation','independent_replay','formal_verification','external_review','novelty_claim'))
for file,expected in wr_receipt['files'].items():assert sha(HERE/file)==expected,('stale gap-three Wronskian receipt',file)
WR_PROOF=BASE+'GAP_THREE_WRONSKIAN_RECURRENCE.tex'
WR_CHECKER=BASE+'verify_gap_three_wronskian.py'
wr_ids=['FIRGAPTHREEREC20261117','FIRGAPTHREEREAL20261117','FIRGAPTHREEFIFTH20261117']
wr_common=dict(common,canonical_source=WR_PROOF,checker=WR_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_gap_three_wronskian.py'),
    software_lock=[WR_PROOF,WR_CHECKER,BASE+'gap_three_wronskian.py',BASE+'gap_three_wronskian_audit.json',BASE+'gap_three_wronskian_verification.json',BASE+'gap_three_wronskian_diagnostics.json',GAP_PROOF,BASE+'higher_gap_adjoint_rank.py',BASE+'higher_gap_adjoint_rank_verification.json',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json'])
new.extend([
 dict(wr_common,id=wr_ids[0],
    title='A divided-difference identity gives a uniform gap-three Wronskian recurrence and coefficientwise shape positivity',
    scope='For ALL n>=0 and real gamma, the polynomial auxiliary obeys B_(n,gamma-1)=B_(n,gamma)-Bprime_(n,gamma). At alpha=exp(i*pi/3), W_n is the divided difference of consecutive B_(n+1),B_n at alpha*z and conjugate(alpha)*z. Its exact recurrence is W_(n+1,gamma)=a*(a-1+z)*W_(n,gamma)-a*z*Wprime_(n,gamma)+z^2*W_(n,gamma-1),a=3n+gamma+2. For F_n(z,b)=W_(n,n+b-1),a=4n+b+2, this is a triangular rational recurrence involving F_n(b+1) and F_n(b). Each coefficient of F_n as a polynomial in z,b is NONNEGATIVE, and every z coefficient is strictly positive at b>=0, by a uniform coefficient recurrence. Constant term is prod_(j0..n-1)(gamma+3j+1)*(gamma+3j+2); next-to-leading is n*(gamma+n+1). No radicals, sector theorem or root scan is needed for these identities. Positivity is not a real-rooted-section theorem. A proposed full-section interlacing induction from separate monic positive strictly log-concave inputs is refuted exactly: an explicit output section u^2+1772u+900000 has discriminant-460016. Inputs are not the actual common-Gamma pair; this does not refute its uniform section conjecture, the weaker one-section target, or SFC. No novelty/assurance upgrade.',
    dependencies=[gap_ids[0]]),
 dict(wr_common,id=wr_ids[1],
    title='One real-rooted Wronskian section suffices because positive Gamma moments exclude every common real section root',
    scope='For EVERY n>=1,gamma>0 at gap3, write W(z)=W0(z^3)+zW1(z^3)+z^2W2(z^3). The three real sections have NO common real root. Their constant-section value at zero is strictly positive. A common u!=0 is equivalent by the exact adjoint rank theorem to a repeated root of T_(n,gamma,3) at t=-1/u. At real t, even-n moments E(1+tY^3)^n are strictly positive and odd-n derivatives n*E[Y^3*(1+tY^3)^(n-1)] are strictly positive, for positive Gamma Y. Thus repeated real t is impossible. Consequently ANY NONZERO real constant-linear combination of the sections with only real roots certifies simplicity; in particular ONE nonzero real-rooted section is enough. The useful section/combination may depend on n,gamma. Pairwise strict interlacing of all sections and a closed-sector theorem are stronger sufficient targets, not required for this criterion. The monic common gcd is constant or has even degree and only nonreal roots. Uniform existence of a suitable real-rooted section/combination remains UNPROVED. No all-index higher-gap SFC, general FC/unrestricted SFC or novelty/assurance upgrade.',
    dependencies=[gap_ids[1],wr_ids[0]]),
 dict(wr_common,id=wr_ids[2],
    title='Every positive Gamma shape has simple fifth-order gap-three moments by an explicit positive discriminant',
    scope='For EVERY gamma>0, T_(5,gamma,3)=sum binom(5,k)*(gamma)_(3k)*t^k is squarefree over C. Its Wronskian section W2 is a quadratic whose A,B,C and discriminant are explicitly displayed in GAP_THREE_WRONSKIAN_RECURRENCE.tex. A,B,C have strictly positive coefficients in gamma; B^2-4AC has eleven strictly positive coefficients from1701*gamma^10 to66701465600, proving two distinct negative roots for ALL positive gamma. The one-real-section criterion applies. Together with FIRALLGAPINITIAL20261117 this closes terminal n1..5 at gap3 for ALL positive Gamma shapes: both X*(A+B*X^3) and A+B*X^3 have no consecutive zero moments through terminal5. Ordinary constant-plus-cubic was already closed at all starts; this is a stronger positive-shape INITIAL block, not a new ordinary all-index or generic higher-gap theorem. Finite parameter-edge/interlacing profiles are only controls, not extrapolation. General FC/unrestricted SFC, uniform higher-index simplicity and novelty/independent/formal/external assurance upgrades remain unclaimed.',
    dependencies=[wr_ids[1],'FIRALLGAPINITIAL20261117'])])
wr_boundary=(' CURRENT weaker gap3 gate: FIRGAPTHREEREC20261117 proves the exact two-parameter Wronskian recurrence and bivariate coefficient positivity. FIRGAPTHREEREAL20261117 excludes ALL common real section roots, so just ONE nonzero real-rooted section or real constant-linear combination certifies simplicity. A useful section may depend on n and gamma; proving all three strictly interlace is unnecessary. FIRGAPTHREEFIFTH20261117 proves terminal orders through5 at gap3 for ALL positive Gamma shapes by a symbolic positive discriminant. Uniform higher-index existence of such a section remains open. Separate-input interlacing/log-concavity does not by itself propagate full section interlacing under the recurrence; the exact counterexample is not an actual moment pair or SFC failure. General FC/unrestricted SFC remain open.')
for row in new:
    if row['id']==gap_ids[1] and wr_boundary not in row['scope']:row['scope']+=wr_boundary

pc_receipt=json.loads((HERE/'gap_three_polar_continuant_verification.json').read_text())
assert pc_receipt['status']=='pass' and pc_receipt['uniform_real_rooted_polar_companion_proved']
assert pc_receipt['companion_roots_simple_below_minus1_proved'] and pc_receipt['nonzero_pivot_backward_coprimality_criterion_proved']
assert pc_receipt['generic_inverse_transfer_to_real_rooted_individual_section_refuted'] and pc_receipt['zero_pivot_fail_closed']
assert pc_receipt['native_prime_squarefree_residual_shortcut_refuted']
assert not any(pc_receipt[k] for k in ('positive_weights_required_for_coprimality','uniform_pivot_nonvanishing_proved','uniform_section_existence_proved','all_index_gap_three_SFC_closed','general_FC_closed','unrestricted_SFC_closed','finite_extrapolation','independent_replay','formal_verification','external_review','novelty_claim'))
for file,expected in pc_receipt['files'].items():assert sha(HERE/file)==expected,('stale polar/continuant receipt',file)
PC_PROOF=BASE+'GAP_THREE_POLAR_AND_CONTINUANT.tex'
PC_CHECKER=BASE+'verify_gap_three_polar_continuant.py'
pc_ids=['FIRGAPTHREEPOLAR20261117','FIRGAPTHREEPIVOT20261117']
pc_common=dict(common,canonical_source=PC_PROOF,checker=PC_CHECKER,proof_type='exact_symbolic',
    artifact_hash='sha256:'+sha(HERE/'verify_gap_three_polar_continuant.py'),
    software_lock=[PC_PROOF,PC_CHECKER,BASE+'gap_three_polar_continuant.py',BASE+'gap_three_polar_continuant_audit.json',BASE+'gap_three_polar_continuant_verification.json',BASE+'gap_three_polar_continuant_certificates.json',BASE+'gap_three_polar_continuant_diagnostics.json',WR_PROOF,GAP_PROOF,BASE+'gap_three_wronskian.py',BASE+'higher_gap_adjoint_rank.py',BASE+'prime_power_windows.py',BASE+'sources/iskander-jain-2020-delta-laguerre-polya.pdf',BASE+'register_all_exponent_squarefree.py',BASE+'source_ledger.json'])
new.extend([
 dict(pc_common,id=pc_ids[0],
    title='A polar-derivative product gives a uniform simple real-rooted companion for the gap-three auxiliary',
    scope='For ALL n>=1 and real gamma>n-2, define p_n(x)=prod_(j0..n-1)(x+2+3j), B_n=sum b_j*z^j, J(y)=sum binom(n,k)*p_n(gamma-k)*y^k, and K(z)=sum n_falling(j)*b_j*z^j. BOTH J,K have n SIMPLE real roots strictly below-1. The exact identities are exp(-z)*B_n=sum p_n(gamma-k)*(-z)^k/k!, b_j=nabla^j p_n(gamma)/j!, J=prod_(j)(gamma+2+3j-yD)(1+y)^n, and J(y)=(1+y)^n*K(-y/(1+y)). Each polar constant exceeds n; an elementary rational-log-derivative proof gives root interlacing, multiplicity loss and strict separation below-1. This covers the needed gamma=n+beta-1,beta>0 at EVERY index. It is a real-rooted COMPANION, not a section or SFC theorem. A generic inverse-transform control at n8 with Kstar=prod_(j0..7)(1+(90+j)*z/108), all simple roots below-1, has three squarefree W sections of degrees5,5,4 with only3,3,0 real roots. Thus generic companion hyperbolicity does NOT imply a real-rooted individual section. The arithmetic-progression Gamma factors or a different valid transfer must be retained; possible real-rooted linear combinations are not refuted. Iskander-Jain Theorems1.1/1.2 were inspected for finite versus infinite zero-free hypotheses; the latter cannot be applied to the backwards Poisson sampling, whose forward ray contains roots. The proof here is self-contained. No uniform section/SFC/FC, novelty or assurance upgrade.',
    dependencies=[gap_ids[0]]),
 dict(pc_common,id=pc_ids[1],
    title='Nonzero backward three-step pivots give an exact coprimality certificate for the cyclic Wronskian sections',
    scope='For ANY real normalized polynomial triple q_N,q_(N-1),q_(N-2), constant coefficients1 and degrees<=floor(j/3), compute lambda_j=[u](q_j-q_(j-1)) and q_(j-3)=(q_j-q_(j-1))/(lambda_j*u) down through j3. If EVERY pivot is nonzero, the generated polynomials end at q0=q1=q2=1, and backward propagation excludes a common root over C. NO positivity or real-rootedness is required. Positive pivots additionally give a weighted nonoverlapping-length-three-interval path representation. At the actual gap3 family, normalize q_(2n-r)=W_r/W_r(0); coefficient positivity makes the constants nonzero in gamma=n+beta-1,beta>0. Uniform nonzero pivots have NOT been proved. Zero pivot is INCONCLUSIVE: a full-degree positive example with q6=1+u+u^2,q5=1+u,q4=1+2u has gcd1 but lambda6=0. Symbolic pivots remain formal nonzero conditions until their parameter ranges are certified; n2,3 have positive numerator/denominator coefficient certificates for ALL gamma>0 at those FIXED indices. Numeric ordinary n8,24 certificates are exact finite controls, with tamper rejection; n24,48 Sturm diagnostics are not extrapolated. Native p3 residual-squarefreeness is FALSE at n22 (1+Y^3 on edge19,27 to22,33), and p3 residual occupancy plus parity still permits degree2 at n63 (edge54,80 to60,89, denominator2, residual1+2Y^3). These are limits of specific sufficient methods, not SFC counterexamples or all-prime limitations. General FC/unrestricted SFC and all novelty/assurance upgrades remain unclaimed.',
    dependencies=[gap_ids[1],wr_ids[0]])])
pc_boundary=(' CURRENT structural gap3 additions: FIRGAPTHREEPOLAR20261117 proves uniform simple companion roots below-1 for ALL gamma>n-2, but a generic inverse-companion control disproves transfer to an individual real-rooted section without the actual arithmetic-progression structure. FIRGAPTHREEPIVOT20261117 supplies a backward coprimality certificate requiring NONZERO pivots only; positivity is optional and zero pivots are inconclusive. Uniform actual-family nonzero pivots or a valid real-section transfer remain open. Native p3 repeated residuals at22 and a surviving degree2 occupancy at63 rule out two stronger shortcuts; they are not SFC failures. General FC/unrestricted SFC remain open.')
for row in new:
    if row['id']==wr_ids[1] and pc_boundary not in row['scope']:row['scope']+=pc_boundary

tilted_boundary=(' CURRENT tilted-cubic indexed progress: FIRTILTEDSHIFT20261117 proves an all-r shifted polygon and both residual types. FIRTILTEDPRIMEPOWER20261117 proves T_(p^r,2) irreducible for ALL primes/r and squarefreeness in ten complete shifted families. FIRCUBICFOLDS20261117 supplies exact nondegenerate fold geometry, NOT a simplicity theorem. The specific prime set2..11 fails at n732 while13 completes that index; two tempting large-prime sufficient exclusions also fail in retained exact profiles. Uniform prime coverage, all-index cubic simplicity, general FC and unrestricted SFC remain open.')

uniform_boundary=(' CURRENT all-exponent family update: FIRALLUNIVTWOSFC20261117 now proves ordinary Gamma-shape1 A+B*Y^m Strong Factorial at EVERY starting index for EVERY positive integer m, with effective uniform cutoff10000 plus a complete finite prefix. For the sequence F_(n,m) the exceptional sets are EMPTY: E_m=D_all=0. Earlier possible-exception or open-composite-exponent boundaries for THIS family are superseded; older sufficient-exporter failures do not reopen it. The current record retains its original local/density/general-shape scope and source qualifications. General supports, both-nonconstant or genuinely multivariate monomials, other Gamma shapes and unrestricted FC/SFC remain outside this closure.')
def refresh_sequence_boundary(row):
    row=refresh_six_boundary(row)
    text=row['scope'].lower()
    if row['id'] not in uniform_ids+pade_ids+parabola_ids+dyadic_ids+tilted_ids+cubic_ids+gap_ids+wr_ids+pc_ids and any(word in text for word in ('f_(n,m)','f_(n,6)','two-term','factorial exponent','exponent-six','non-prime-power')):
        if uniform_boundary not in row['scope']:row['scope']+=uniform_boundary
    return row
new=[refresh_sequence_boundary(row) for row in new]
legacy_changes={row['id']:refresh_sequence_boundary for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']
    if row['id'].startswith('FIR') and row['canonical_source'].startswith(BASE) and row['id'] not in {item['id'] for item in new}}

def open_scope(row):
    extra=(' FIRLARGEM20261117 now proves the ELEMENTARY pointwise exclusion2^m>=4*n^2=>squarefree F_(n,m). FIRLOGEXPCANDIDATES20261117 combines this with the ALREADY proved prime-power theorem and explicit m1 unit case: any actual bad pair has m<2log2(2n) and at least two distinct prime divisors. The finite filter does not reopen indices excluded by other certificates. At n>=2^144, the all-positive-exponent bad set equals the prior capped m<=n^(1/16) set, giving D_all<=D+2^144. FIRALLEXPLOG20261117 therefore removes the exponent cap entirely with the independent PUBLISHED-2016 logarithmic rate; FIRALLEXPPOWER20261117 removes it with the source-qualified dense-MR-II power saving and convergent reciprocal exceptional sum. The sparse-sieve warning and limited preprint audit remain explicit. The2^144 threshold is overlap only, not an analytic threshold or a checked initial range. All-index SFC remains OPEN; the pointwise frontier is logarithmically small non-prime-power exponents such as6. Inactive-variable FC remains closed.')
    if extra not in row['scope']:row['scope']+=extra
    bridge=(' Density-to-all-index FAILURE PROPAGATION is an OPEN possible bridge: one genuine SFC failure would have to force enough DISTINCT further bad terminal indices, still in the F_(n,m) family counted by D_all, to contradict a retained bound. Positive upper density would contradict the published logarithmic route; a forced divergent reciprocal sum would contradict the source-qualified power/harmonic route. Infinitely many new failures, or many phases at one index, need not suffice. No such propagation theorem emerged from the comparison. Existing failure-preserving recomposition applies to pure/eventually pure unsafe seeds; short consecutive SFC zero blocks do not supply those hypotheses. The separate safety-reflection result also gives no such conversion. These established recomposition results do not currently bridge the two situations; no new theorem or all-index closure is claimed. Direct all-index arithmetic or simple-root proofs remain alternative routes.')
    if bridge not in row['scope']:row['scope']+=bridge
    analytic=(' PRIMARY focused exponent-six endpoint: prove P_(n,6)(z)=z^n F_(n,6)(1/z) has simple roots for EVERY n>=an EFFECTIVE N0, then verify the finite prefix exactly. Pprime_(n+1,m)=(n+1)P_(n,m) identifies the adjacent obstruction exactly. The canonical status note derives the integral scaling and saddle map, an elementary ALL-n root enclosure |z|<2800*n^6, excluding the nonzero saddle-collision value w3125, and distinct leading phase derivatives for distinct saddle/endpoint contributions. These auxiliary observations do NOT supply eventual simplicity or N0. Endpoint/transition regions, possible multiple comparable contributions, effective errors and ALL-root coverage remain open. Bulk zero distribution and density bounds do not exclude sparse double-root indices. Failure propagation remains a separate unproved alternative; no theorem state is upgraded.')
    if analytic not in row['scope']:row['scope']+=analytic
    extra=(' FIRANALYTICANNULUS20261117 now proves the effective large-index ROOT ANNULUS1<|w|<2800 for every n>=10^6, with complete contour/tail control on |w|<=1, and explicit uniform saddle-data separation on that annulus. Both saddle-collision values are excluded from every possible large-index root region. The10^6 threshold is zero-exclusion only; no simplicity cutoff, finite prefix or all-index SFC closure is claimed. The immediate missing gate is a global contributing-contour decomposition with effective amplitude/C1-error bounds at all cancellations, including possible multi-contribution transitions.')
    if extra not in row['scope']:row['scope']+=extra
    extra=(' COMPLETED exponent-six gate: FIREVENTUALSIMPLE620261117 proves every P_(n,6) simple for n>=the EFFECTIVE cutoff10000; FIRALLINDEXSFC620261117 supplies the COMPLETE exact finite prefix and closes ordinary A+B*Y^6 at ALL starts. The exact split at3 plus uniform analytic error/derivative estimates closes the former contour-selection/cancellation gate. The preceding annulus-only open roadmap is superseded. The overarching goal of resolving general FC and/or unrestricted SFC is still OPEN; closing this special family does not finish it.')
    if extra not in row['scope']:row['scope']+=extra
    extra=(' COMPLETED all ordinary constant-plus-power families: FIRALLUNIVTWOSFC20261117 proves F_(n,m) squarefree for EVERY positive n,m, hence A+B*Y^m SFC at ALL starts and ALL positive exponents. The uniform analytic cutoff is10000, combined with the large-exponent criterion; all92935 additional necessary finite pairs are exactly certified, reusing prime-power/m1 and m6 closure. The old D_all-counted exceptional set is now EMPTY, without MR inputs or failure propagation. The former open-composite-exponent and density-to-all-index gates for this counted family are closed directly. No reduction from unrestricted FC/SFC to constant-plus-one-variable-power supports is proved, so the OVERARCHING GOAL remains OPEN. The next scope must concern genuinely new polynomial supports or Gamma parameters, not reopening this completed family.')
    if extra not in row['scope']:row['scope']+=extra
    extra=(' COMPLETED new adjacent-support cases: FIRADJACENTMONOSFC20261117 proves ordinary X^a(A+B*X) SFC at ALL starts for a1..6, with a0..2 source-credit boundary retained and a3..6 proved by exact Pade/matrix separation. FIRPADESEPARATION20261117 supplies all-parameter denominator strips, exact degree-(a-1) auxiliaries, generic positive Schwarz certificates and an explicit necessary bidiagonal matrix pencil. The indexed condition a>=2*n-3 excludes every adjacent failure in that UNBOUNDED range. The residual a>=7,n>(a+3)/2 and the all-a spectral/shifted criterion remain open; the latter is stronger than the genuine common-root exclusion, not equivalent to SFC. No unrestricted FC/SFC resolution follows.')
    if extra not in row['scope']:row['scope']+=extra
    if pade_current_boundary not in row['scope']:row['scope']+=pade_current_boundary
    if dyadic_boundary not in row['scope']:row['scope']+=dyadic_boundary
    if frontier_boundary not in row['scope']:row['scope']+=frontier_boundary
    if tilted_boundary not in row['scope']:row['scope']+=tilted_boundary
    if cubic_boundary not in row['scope']:row['scope']+=cubic_boundary
    if gap_boundary not in row['scope']:row['scope']+=gap_boundary
    if wr_boundary not in row['scope']:row['scope']+=wr_boundary
    if pc_boundary not in row['scope']:row['scope']+=pc_boundary
    row['dependencies']=list(dict.fromkeys(row['dependencies']+ids+analytic_ids+simple_ids+uniform_ids+pade_ids+parabola_ids+dyadic_ids+tilted_ids+cubic_ids+gap_ids+wr_ids+pc_ids));return row

def capped_scope(row):
    row=refresh_sequence_boundary(row)
    target='FIRALLEXPLOG20261117' if row['id']=='FIRSIMULTLOG20261117' else 'FIRALLEXPPOWER20261117'
    extra=(' The exponent cap belongs to THIS analytic theorem. '+target+' now removes it using the independent elementary large-exponent bound and an exact2^144 overlap; no stronger analytic input or analytic starting threshold is claimed.')
    if extra not in row['scope']:row['scope']+=extra
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{**legacy_changes,'OP-FACTORIAL-GENERAL20261116':open_scope,
    'FIRSIMULTLOG20261117':capped_scope,'FIRSIMULTPOWER20261117':capped_scope},new)
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+ids+analytic_ids+simple_ids+uniform_ids+pade_ids+parabola_ids+dyadic_ids+tilted_ids+cubic_ids+gap_ids+wr_ids+pc_ids+['MVPUREFAILRECOMP20261116','MVRECOMPREFLECT20261116']));row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PC_PROOF,PC_CHECKER,BASE+'gap_three_polar_continuant_verification.json',BASE+'gap_three_polar_continuant_certificates.json',BASE+'gap_three_polar_continuant_diagnostics.json',BASE+'sources/iskander-jain-2020-delta-laguerre-polya.pdf',WR_PROOF,WR_CHECKER,BASE+'gap_three_wronskian_verification.json',BASE+'gap_three_wronskian_diagnostics.json',GAP_PROOF,GAP_CHECKER,BASE+'higher_gap_adjoint_rank_verification.json',BASE+'higher_gap_adjoint_rank_diagnostics.json',CUBIC_PROOF,CUBIC_CHECKER,BASE+'cubic_adjoint_stability_verification.json',BASE+'sources/wagner-multivariate-stable-polynomials-2011.pdf',PROOF,BASE+'all_exponent_squarefree_audit.json',BASE+'all_exponent_squarefree_verification.json',BASE+'STATUS_AND_LITERATURE.tex',ANNULUS_PROOF,ANNULUS_CHECKER,BASE+'saddle_annulus_verification.json',SIMPLE_PROOF,SIMPLE_CHECKER,BASE+'exponent_six_simplicity_verification.json',UNIFORM_PROOF,UNIFORM_CHECKER,BASE+'uniform_two_term_sfc_verification.json',TILTED_PROOF,TILTED_CHECKER,BASE+'tilted_cubic_frontier_verification.json',BASE+'tilted_cubic_frontier_diagnostics.json',DYADIC_PROOF,DYADIC_CHECKER,BASE+'dyadic_binomial_sfc_verification.json',BASE+'dyadic_binomial_frontier_diagnostics.json',PARABOLA_PROOF,PARABOLA_CHECKER,BASE+'pade_parabola_sfc_verification.json',PADE_PROOF,PADE_CHECKER,BASE+'pade_monomial_separation_verification.json',BASE+'pade_shift_positive_certificates.json','jacobian-followup/moment-vanishing-audit/POLYNOMIAL_RECOMPOSITION_REFLECTION_AND_UNIVERSAL_FAILURE.md']))
    return row

def lesson(row):
    row=refresh(row)
    extra=(' Use the deterministic coefficient bound2^m>=4n^2 before any arithmetic search. Normalize endpoint coefficients, majorize the middle coefficients by the symmetric binomial q^min envelope and bound BOTH weighted masses below1 to exclude repeated roots inside/outside the unit circle. The prime-power residual theorem is already registered; make m1 explicit via unit deformation at r0. Any remaining pair has2^m<4n^2 and at least two distinct prime divisors, so exponent candidates per n are logarithmically bounded. Exact bit lengths avoid building enormous powers. For n>=2^144 the remaining exponents fit below n^(1/16), transferring BOTH simultaneous analytic estimates to ALL positive exponents by D_all<=D+2^144. Harmonic convergence inherits the qualified power route.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not present the retained prime-power residual proof as a newly discovered contribution or a necessary-candidate filter as a counterexample/unknown-status test: other registered certificates can already exclude particular candidates. The explicit overlap2^144 does not extract any analytic constant/starting threshold or audit all smaller indices. Preserve MR II preprint/sparse-warning qualification for the all-exponent power and harmonic result; the published log route is independent. Removing the exponent cap asymptotically still does not close all-index SFC.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' If using the density results for all-index closure, formulate an actual failure-propagation theorem with enough DISTINCT bad terminal indices in the counted two-term family. Positive upper density would contradict the published log bound; forced divergent reciprocal mass would contradict the qualified power/harmonic bound. No such bridge is established by the comparison. The proved failure-preserving recomposition theorem starts from pure/eventually pure unsafe seeds, not a finite SFC zero block; safety reflection has different hypotheses and supplies no missing conversion either.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not apply all-order pure-failure recomposition to a short consecutive SFC obstruction, silently assume the block produces an all-order/eventual pure unsafe seed, or count recompositions outside F_(n,m) as new indices in D_all. Infinitely many failures can be geometrically sparse and remain compatible with a power-saving count and convergent reciprocal sum. This open propagation gate is not a proved theorem or a claim that no other all-index proof is possible.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' The PRIMARY exponent-six target is effective eventual simplicity: reverse to P_(n,m), use Pprime_(n+1,m)=(n+1)P_(n,m), and prove all roots simple for EVERY n>=effective N0 before checking the finite prefix. The scaled integral has saddle map w=s^5(6-s). The elementary root enclosure |w|<2800 excludes its nonzero collision value3125. Distinct saddle leading derivatives differ by(t-s)/(6w), and endpoint/saddle derivatives bys/(6w); use these only with justified contour contributions and effective differentiated errors. Seek disjoint one-zero Rouche contours covering all n roots, or an equivalent uniform P/Pprime exclusion, including endpoint and transition regions.')
    if extra not in row['use']:row['use']+=extra
    extra=(' A bulk zero law, most-root expansion, saddle collision, or formal difference of leading phase derivatives is not a proof of polynomial simplicity. That former no-cutoff boundary is now superseded for ordinary exponent six by N0=10000 and its complete finite prefix. The logarithmic phase is singular at(s,w)=(0,0); ordinary saddle reasoning does not apply there. Three or more comparable contributions require separate control. Root separation may shrink with n, but the certificate must cover every root and yield an effective cutoff. Eventual simplicity for exponent six would close that two-term family, not unrestricted SFC.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' FIRANALYTICANNULUS20261117 completes effective geometric localization: every exponent-six root lies in1<|w|<2800 for ALL n>=10^6. The full shifted-contour Gaussian proof, including all tails, excludes the inner disk; the elementary coefficient proof excludes the exterior. Both critical values0 and3125 are outside the remaining compact annulus. All formal saddles are uniformly separated and nondegenerate there, with explicit phase derivative gaps. The next gate is contour selection and effective amplitude/C1-remainder control across every cancellation on this annulus, including multiple comparable contributions. Use the exact U+V+E transversality criterion only after verifying its strict bound.')
    if extra not in row['use']:row['use']+=extra
    extra=(' The annulus cutoff10^6 is NOT an eventual-simplicity cutoff or a checked finite prefix. Saddle separation is not polynomial-root separation; an analytic decomposition into at most two contributions has not been proved. Excluding both saddle-collision values does not remove endpoint contributions at nonzero w or establish all-index SFC. The 34 finite controls verify constants and identities, while the contour note is the universal analytic evidence.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' COMPLETED ordinary exponent-six method: use the EXACT split I=K[0,3]+J[3,infinity). Entropy gives every root |w|<240 above n10000. A complete horizontal saddle contour estimates J uniformly on |w|<=250 with analytic relative error<12/sqrt(n); radius10 Cauchy estimates control its logarithmic derivative on |w|<=240. Watson bounds estimate K relative to w^n/n on |w|>=50. Inner comparison excludes |w|<=60, and at every cancellation the log-derivative gap is>n/288-2. This proves effective eventual simplicity with N0=10000. Every smaller index is then certified by small modular Bezout identities and parity, giving all-index A+B*Y^6. Retain the earlier annulus proof as a completed step, not the remaining frontier.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not keep exponent six marked open in ordinary Gamma shape1, confuse its actual simplicity cutoff10000 with the prior annulus-only10^6 threshold, or call the finite prefix a root sample: it covers every index below10000 exactly. This closes only the single-variable sixth-power two-term family. No valid reduction from arbitrary-polynomial FC/SFC to that family is supplied, and other shapes or exponents require their own proof. Do not turn an older exporter failing its sufficient local test into a mathematical unknown already settled by the global theorem.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['avoid']=row['avoid'].replace('No effective N0 or verified initial range has been obtained.','That former no-cutoff boundary is now superseded for ordinary exponent six by N0=10000 and its complete finite prefix.')
    extra=(' CURRENT completed uniform method: normalize z=n^m*A_m*v with A_m=(m/e)^m and split the integral EXACTLY at s=m/2. For n>=10000 and6<=m<=sqrt(n), every root has1/2<|v|<21/10. A uniform large-saddle proof on |v|<=11/5 gives relative error22/sqrt(m*n), Cauchy control and log-derivative error<3; the endpoint on |v|>=2/5 has D,Dprime<1/1000. Every cancellation has derivative gap>n/3-4. The elementary large-exponent criterion covers the remaining m, with the effective overlap m<2log2(2n)<sqrt(n). Complete finite certificates for92935 additional pairs plus reused prime-power/m1/m6 cases close ALL ordinary A+B*Y^m at ALL starts. The former density exceptional sets for F_(n,m) are EMPTY; MR source qualifications remain historical proof-route qualifications, not assumptions of the new direct theorem.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not keep any ordinary constant-plus-one-variable-power exponent open or run larger local coverage scans for that completed family. The uniform pointwise proof covers every n,m; the finite certificate grid is a complete prefix, not numerical asymptotic evidence. Do not infer unrestricted SFC, both-nonconstant monomial pairs, genuinely multivariate monomials or other Gamma shapes. For f=x^a(A+B*x^d), the Gamma shape is a*n+1 and varies with n, so the fixed-shape Appell derivative/squarefreeness equivalence used here does not transfer automatically. Check retained two-monomial literature before selecting that next gate.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'uniform_two_term_sfc.py',UNIFORM_CHECKER]))
    extra=(' CURRENT adjacent-support method: for X^a(A+B*X), use the exact reciprocal Pade denominator Q_(an,n). Its complex symmetric tridiagonal representation bounds ALL root real parts between M-N+2 and M+N, strictly for N>=2. Consecutive strips are disjoint whenever a>=2*n-3. Otherwise use the exact cross-determinant auxiliary H of degree a-1. For a2..6, positive Schwarz coefficients in n-2 reconstruct the shifted H at c=(a-1)*n+2; a dissipative tridiagonal matrix puts those roots strictly left of c for ALL indices. a1 has constant auxiliary2. This closes a1..6 at every start, with a0 already completed. The generic criterion can certify a supplied a symbolically or a supplied pair exactly. An explicit necessary pencil A+zB, dimension a-1 and affine in n, gives a concrete route toward the unbounded-a spectral target.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not leave ordinary X^a(A+B*X) at a1..6 open, infer all a from finitely many symbolic a certificates, or apply fixed-shape Appell simplicity to the varying shape an+1. Neither individual Pade root simplicity nor positive-real-part location alone proves adjacent coprimality; use the strip/auxiliary separation. H_(3,2) refutes blanket real-rootedness, and H_(5,2) refutes blanket unshifted Hurwitz stability. The SHIFTED spectral bound is sufficient and may be stronger than needed; if it fails, this does not refute SFC. Indexed a>=2*n-3 coverage is not all-starting-index coverage for an arbitrary fixed a. No all-powers or arbitrary recomposition closure is implied.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'pade_monomial_separation.py',PADE_CHECKER]))
    extra=(' CURRENT UNIFORM adjacent-support method: use a SHARED PARABOLA instead of the stronger vertical target. The R_N matrix gives x=S-2*tau,y^2<=4*tau*(S-tau); a positive exact gap puts every root inside y^2<4*N*(x+N). The necessary degree-(a-1) continuant is ZERO-FREE in that SAME parabola: normalized b,c satisfy b_j*(1-b_prev/c_j)>=N, and mu=z/(z+v), Re v>=N, preserves Re mu<1 by completing the square. This proves all integer shifts, all l>-1 and all indices, hence X^a(A+B*X) SFC for EVERY a>=0 and EVERY positive Gamma shape. No asymptotic, prefix, MR, parity or irreducibility input. Other-variable monomial factors contribute nonzero moments, extending to supports differing by one power in one variable.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not leave all-a adjacent supports open or keep the vertical bound as a required gate: shared-parabola separation closes the family independently. Finite replays check algebra; the written invariant and Cauchy inequalities are the universal proof. Do not infer arbitrary gaps d>1, unrestricted SFC or general FC. The parameter must satisfy l>-1 and its increment is a positive INTEGER. All positive Gamma shapes are covered; general complex shapes are not claimed. No novelty or formal/independent/external assurance upgrade.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' For gaps d>=2 retain U_N^(d)=sum(-1)^k binom(N,k)*(t+1)_(d*k)*z^(N-k), its exact Pascal/parameter-raising identities, and the sparse order-d+1 parameter recurrence in equations10-13 of the uniform proof. Two consecutive moments give only V_l=0 and z*V_M=(M+1)_d*V_(M+d). These two constraints no longer close a finite square pencil for d>1. Additional structure of the particular polynomial solution is required; importing the d1 zero-free argument without that input would leave a proof gap.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' CURRENT nonadjacent unit-transfer method: for q=2^r,1<=integer beta<=2q, X^q(A+B*X^q) satisfies SFC at ALL starts. Normalize moments by one common variable scaling2^(-q), and use W_(L,delta)=(qL+delta)!/(2^((q-1)*L)*L!*delta!) odd. The resulting degree-n coefficients reduce to Bessel Y_n or D_n=(n+1)*Y_n+z*Y_nprime; both retain odd leading coefficients and obey C_n=z*C_prev+C_before modulo2. Their consecutive gcd is1, so the rational resultant is a unit. This works for arbitrary rational unit sequences multiplying (L+h)!,h0 or1. It gives the exact raw resultant valuation3q*n*(n-1)-(n-1)*s2(n)-n*s2(n-1). Other independent-variable monomial factors contribute nonzero moments and are allowed.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not infer arbitrary factorial substitutions from characteristic-zero adjacent coprimality. Full-degree coprime reduction is what makes the dyadic transfer work. For h2, the natural normalized reductions at n2,n3 are both z^2; for gap2,a1,n3 the two Newton polygons share slope3/2 and residual1+Y. These are exact certificate limitations, not SFC counterexamples: the literal rational pairs are coprime. Neither arbitrary positive Gamma shapes nor every lower exponent at gap2 is proved by this unit argument. Do not reopen the completed dyadic diagonal family or extend it from finite q/n/shape replays.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' CURRENT first-gap scalar and arithmetic route: for every e,n>=1, T_(n,e)=sum binom(n,k)*(n)_(e*k)*t^k gives moment identities A_n=T+(e/n)*theta*T,A_prev=T-(1/n)*theta*T. Thus X*(A+B*X^e) SFC is EXACTLY all-index simplicity of these varying-shape polynomials, without a fixed-shape Appell assumption. For e2 and prime n<q<=3*(n-1)/2, T has primitive nonhorizontal edges, confining every rational square divisor g^2 to degree<=floor((q-n)/4), then even degree<=2floor((q-n)/8). It reduces degree-preservingly into a square divisor of explicit K_(q-n). The two raw moments have disjoint positive Newton slopes and reduce horizontally to H_d,H_(d+1); K=(H_d+2H_next)/3 and theta*K=d*(H_next-H_d)/3. Use full residual factor-degree sets and permit the prime to depend on the candidate degree.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not confuse the proved lower-one multiple-root equivalence with simplicity itself, import a fixed-shape derivative recurrence at Gamma shape n, or claim prime-window existence/coverage for every n from a conditional Newton theorem. The old two-boundary deficit is now supplied by a scalar critical-point characterization for lower exponent one, not for all lower exponents. A prime may have bad residual collisions yet another admissible prime may certify the index; demanding coprimality at every residual prime is unnecessary. Quarter-gap parity is for RATIONAL square divisors, not arbitrary complex linear factors. General cubic simplicity, all-index coverage and unrestricted FC/SFC remain unproved.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' CURRENT tilted-cubic prime-power mechanism: for n=a*p^r+i,p>3max(a,i), the all-level Kummer/Legendre estimate gives horizontal width i and reduced slope denominator p^r. The horizontal residual is T_i. IMPORTANTLY the sloped residual is T_a of shape a at i0, but the raw a-th cubic moment of shape a+1 at positive i. Wilson digit units prove both; copying a fixed-shape residual would be wrong. At prime-power INDICES, T_(p^r,2) is irreducible for EVERY prime. For2,3 separate digit bounds give primitive heights2P-1 and P-1. Parity plus a nonsquare leading coefficient proves all ten listed shifted index families without finite-r extrapolation. Use those completed indexed gates as inputs, not remaining unknowns.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not confuse prime-power INDICES of T_(n,2) with the older prime-power EXPONENT theorem. These indexed classes do not close the cubic at every start. The set2,3,5,7,11 leaves degree6 at n732, and prime13 removes it; do not infer that a convenient fixed small set covers all indices. The conditions p>2k+1 or p>3k with nmodp<2k do not themselves exclude degree k; exact counterprofiles at51 and118 are retained. Do not mistake the exact fold discriminant/derivatives for contour selection, Airy bounds, a root-location theorem or eventual simplicity. An analytic all-root argument must justify the contributing contours in potential fold neighborhoods or rigorously exclude them.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' The residual K_d kernels are NOT all real-rooted: K11 has degree5 and exactly3 positive real roots by exact Sturm count, while its gcd with its derivative is1. Alternating coefficients rule out negative roots. This kills a blanket real-rooted-kernel shortcut but gives no repeated-root or SFC counterexample; retain the exact diagnostic.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' CURRENT whole-cubic proof: build the Gamma-convolution auxiliary B and use the classical GWS half-plane polarization theorem with the rederived Gamma Jacobi bound. Construct both explicit adjoint solutions. Their Lagrange concomitants vanish by the rational logarithmic derivative argument; a double root would force B(z)=B(-z)=0 in disjoint half-planes. This proves all n,gamma>=n-2 with gamma>0 and closes X*(A+B*X^2) at every start/all positive Gamma shapes. The arithmetic/fold route is now independent evidence, not a remaining cubic obligation.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not assume convex averages preserve Hurwitz stability: use the actual Walsh additive convolution/polarization identity and verify the classical half-plane hypotheses. For noninteger gamma, the gauge is locally z^(1-M); only the stripped concomitant is rational. The simple-root theorem does not assert gamma below n-2 or arbitrary supports/gaps. Prior indexed certificates do not supply the new universal proof.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' HIGHER GAP: use B_(n+1)=(z+dn+gamma+d-1)*B_n-z*B_nprime, obtained by a formal residue/Lagrange change of variables. It constructs d independent exponential-polynomial adjoints at EVERY gap. A repeated root is EXACTLY rank<d-1 for the polynomial adjoint jet matrix; the full determinant reconstructs the original reversed moment. For d3, export the three cyclic sections of the normalized pair Wronskian. Uniform gcd1 or strict interlacing would suffice, as would zero-freeness in the CLOSED pi/3 sector. These are current open gates, not conclusions of a bounded scan.')
    if extra not in row['use']:row['use']+=extra
    extra=(' The cubic proof does NOT extend by setting each higher concomitant coefficient to zero: at gap d>=3 the surviving Rsecond..R^(d) terms can cancel. Use the full rank condition. The natural stronger Hurwitz assertion for the gap3 Wronskian is FALSE at n4,gamma4; its exact negative seventh Hurwitz determinant is retained. Finite strict interlacing at n4,6,8 is only evidence. A generalized Hurwitz total-nonnegativity theorem requires its own positivity proof and gives an OPEN sector, so boundary exclusion must also be addressed. Negative/complex Gamma failure controls are outside positive-shape SFC and are not counterexamples.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' UNIFORM INITIAL BLOCK: every Gamma power X=Y^d has positive third central moment. For r(x)=(x)_d, convexity gives c>=2b-a and strict increase b>a, hence kappa3>=2a(b-a)^2>0. Together with the real-root parity argument, this proves T_(n,gamma,d) simple for ALL gaps/shapes at n1..4. A real quartic moment square would force kappa3=0. This is reusable initial-range elimination, not an all-n argument.')
    if extra not in row['use']:row['use']+=extra
    extra=(' CURRENT gap3 mechanism: use the divided difference of consecutive auxiliaries and the exact Wronskian recurrence with Gamma shift gamma-1. Center gamma=n+b-1 to export F_(n+1)(b) from F_n(b+1),F_n(b) rationally; all z,b coefficients are nonnegative by a uniform induction. The real place rules out EVERY common real root of the three cyclic sections, so ONE nonzero real-rooted section or real constant-linear combination already certifies the index; allow the choice to depend on n/gamma. At terminal5 a quadratic section has a discriminant with eleven strictly positive gamma coefficients, proving every positive shape. Uniform higher-index real-rooted-section existence is the current gate.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not demand pairwise strict interlacing of all sections or a closed-sector theorem when one real-rooted section would suffice. Coefficientwise positivity alone does not prove real-rootedness. Separate monic positive strictly log-concave inputs, each with negative interlacing sections, can produce a nonreal section under the formal recurrence: the exact u^2+1772u+900000 counterexample has discriminant-460016. This refutes that full-interlacing induction shortcut, not the actual common-Gamma family, the weaker one-section target, or SFC. Retain Gamma-shift compatibility or another valid invariant. Finite edge/large-shape profiles are not a uniform theorem.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' UNIFORM polar companion: B_n has backward-difference coefficients of p_n(gamma)=prod(gamma+2+3j). The binomial companion J=prod(gamma+2+3j-yD)(1+y)^n and K=sum n_falling(j)*b_j*z^j have SIMPLE real roots below-1 for EVERY gamma>n-2, by a self-contained polar derivative argument and a real Mobius map. The full arithmetic progression remains relevant: generic inverse-companion hyperbolicity does not ensure even one real-rooted individual W section. An alternative exact certificate normalizes the sections and descends q_j-q_(j-1)=lambda_j*u*q_(j-3). Nonzero pivots alone prove gcd1; positivity gives a weighted length-three-interval interpretation but is not required. Uniform actual-family nonzero pivots or a valid companion-to-section transfer remain open.')
    if extra not in row['use']:row['use']+=extra
    extra=(' The finite-difference source has distinct FINITE interval and INFINITE zero-free-ray hypotheses. Iskander-Jain Theorem1.2 cannot be applied to p_n(gamma-s) sampled forward from0, since it has positive roots. Do not infer B or W hyperbolicity from the real companion. The n8 generic inverse control has section degrees5,5,4 but real root counts3,3,0; it is not an actual Gamma pair or SFC counterexample, and it says nothing against possible real-rooted linear combinations. Zero backward pivots are INCONCLUSIVE, not bad indices; keep symbolic pivots as conditions until parameter ranges are proved. Native p3 residuals already repeat at n22, and residual-plus-parity occupancy leaves degree2 at63; no all-prime limitation or actual SFC failure follows.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'gap_three_polar_continuant.py',PC_CHECKER]))
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'gap_three_wronskian.py',WR_CHECKER]))
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'higher_gap_adjoint_rank.py',GAP_CHECKER]))
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'cubic_adjoint_stability.py',CUBIC_CHECKER]))
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'tilted_cubic_frontier.py',TILTED_CHECKER]))
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'dyadic_binomial_sfc.py',DYADIC_CHECKER]))
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'pade_parabola_sfc.py',PARABOLA_CHECKER]))
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'all_exponent_squarefree.py',CHECKER]));return row

def work(row):
    row=refresh(row)
    row['next_step']='OVERARCHING GOAL ACTIVE: resolve general Factorial and/or unrestricted Strong Factorial rigorously. Completed support families, whole cubic, all-gap n1..4 and gap3 n1..5 remain closed in their recorded scopes. Gap3 identities, coefficient positivity and one-real-section criterion are proved. NEW uniform polar companions have simple roots below-1 for all gamma>n-2, but generic inverse transfer to W sections is refuted. NEXT exploit the actual arithmetic progression gamma+2+3j to prove a valid transfer yielding some real-rooted section/constant-linear combination, or prove actual-family nonzero backward three-step pivots. Do not require positive pivots when nonzero suffices, and treat zero pivots as inconclusive or use a justified look-ahead method. Neither uniform gate is proved; finite positive certificates are not an all-index theorem. Do not apply an infinite zero-free-ray theorem to backwards sampling that has positive roots. Native prime3 repeated residuals at22 and surviving degree2 occupancy at63 rule out those shortcuts only. Preserve routes toward arbitrary lower exponents, genuinely multivariate supports and unrestricted FC/SFC; no reduction from completed special families to the full conjectures is established. Goal stays active and assurance flags false.'
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
ledger_path=HERE/'source_ledger.json';before=ledger_path.read_bytes();ledger=json.loads(before)
ledger.setdefault('internal_proof_extensions',{})['all_exponents_from_large_exponent_cutoff']={
    'canonical_source':PROOF,'claims':ids,
    'provenance':'User-supplied factorial inequality, symmetric majorant, weighted-root exclusion, exact overlap and localization checked and integrated. The prime-power residual result is reconciled with the already registered unit-deformation proof, not claimed as new. No priority claim.',
    'evidence':{name:sha(HERE/name) for name in ('ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex','all_exponent_squarefree_audit.json','all_exponent_squarefree_verification.json')},
    'elementary_scope':'ALL n,m>=1,2^m>=4n^2 implies squarefree. Necessary remaining exponents are below2log2(2n) with at least two distinct prime divisors. m1 follows from the retained deformation at r0.',
    'analytic_transfer':'Capped and unrestricted bad sets agree at n>=2^144. Both absolute-constant rate bounds transfer to all positive exponents with finite additive count; harmonic convergence uses the qualified power route.',
    'source_boundary':'Published2016 log route independent of MR II. Dense MR II preprint9.2(ii) qualification, published sparse warning and Section10 exclusion retained for power/harmonic conclusion; no full-preprint validation.',
    'boundary':'Overlap threshold2^144 is not an analytic starting threshold or a pointwise initial-range audit. Infinite exceptions, all-index SFC and other Gamma shapes remain possible/open; all assurance flags false.'}
ledger['published_input_pinpoints']['PrimePowerResidualReconciliation20261117']='Rocks2015 author manuscript Lemma3.8/Proposition3.22 supplies the unit congruence/digit hull. The ALL-index residual squarefreeness theorem is ALREADY registered in PRIME_POWER_TWO_TERM_SFC.tex with sigma=(-1)^(h+c*e). Exponent one is the same deformation at r0,c0,units1. The new all-exponent note reuses this evidence, not a new literature or priority claim.'
ledger.setdefault('open_bridge_assessments',{})['short_block_failure_propagation_20261117']={
    'status':'OPEN; no propagation theorem established in this comparison',
    'canonical_source':BASE+'STATUS_AND_LITERATURE.tex',
    'context_claims':['FIRALLEXPLOG20261117','FIRALLEXPPOWER20261117','MVPUREFAILRECOMP20261116','MVRECOMPREFLECT20261116'],
    'missing_implication':'one genuine short consecutive SFC zero block => enough distinct bad indices in the counted F_(n,m) family to contradict an exceptional-set bound',
    'sufficient_targets':['positive upper natural density contradicts the published logarithmic bound',
        'divergent reciprocal mass contradicts the source-qualified power/harmonic bound'],
    'scope_boundary':'failure-preserving recomposition assumes a pure/eventually pure unsafe seed; safety reflection is also not a finite-block conversion. Infinitely many indices alone or new phases outside the counted family are insufficient.',
    'provenance':'User comparison assessment recorded; no new theorem, full literature nonexistence claim or all-index closure asserted'}
ledger.setdefault('open_analytic_targets',{})['effective_eventual_simplicity_exponent_six_20261117']={
    'status':'OPEN; no effective N0, eventual-simplicity theorem or completed finite prefix',
    'canonical_source':BASE+'STATUS_AND_LITERATURE.tex',
    'target':'P_(n,6) simple for every n>=effective N0, followed by exact finite verification; closes only the two-term exponent-six family',
    'exact_starting_identities':['P_(n,m)=z^n F_(n,m)(1/z), with nonzero constant term',
        'Pprime_(n+1,m)=(n+1)P_(n,m)',
        'P_(n,m)(n^m*w)=n^(mn+1)*integral_0^infty(w+s^m)^n*exp(-n*s)ds',
        'saddle map w=s^(m-1)*(m-s)'],
    'auxiliary_enclosure':'For ALL n>=1 every exponent-six root satisfies |z|<2800*n^6. Elementary coefficient bound a_k^(1/k)<1400*n^6 is proved in the canonical note. Nonzero saddle collision w3125 lies outside this disk; w0 is a singular logarithmic endpoint chart, not an ordinary saddle.',
    'remaining_gate':'Effective all-root coverage and error control in endpoint and distinct-saddle transition regions, including possible three-or-more contributions. Disjoint one-zero Rouche contours or an equivalent uniform P/Pprime exclusion would suffice. Formal leading derivative differences do not establish simplicity.',
    'scope_boundary':'Bulk zero distributions and density bounds permit a sparse double-root subsequence. No n-independent minimum root spacing is required. Failure propagation remains a separate open alternative, not a prerequisite for this analytic route.',
    'provenance':'User-selected analytic endpoint; exact reversal, derivative, scaling, saddle identities and elementary coefficient enclosure checked locally. No imported asymptotic theorem or numerical root campaign.'}
ledger.setdefault('internal_proof_extensions',{})['effective_saddle_annulus_exponent_six']={
    'canonical_source':ANNULUS_PROOF,'claims':analytic_ids,
    'evidence':{name:sha(HERE/name) for name in ('EFFECTIVE_SADDLE_ANNULUS.tex','verify_saddle_annulus.py','saddle_annulus_audit.json','saddle_annulus_verification.json')},
    'scope':'Every root satisfies |w|<2800 for all n>=1, and1<|w|<2800 for all n>=10^6. Full local Gaussian and all contour tails controlled uniformly, not a bulk zero law. Explicit separation/nondegeneracy concerns saddle data only.',
    'proof_inputs':'Elementary factorial/integral inequalities, Cauchy contour deformation, Rouche theorem and complex Gaussian integral; no imported asymptotic or density theorem.',
    'boundary':'The10^6 threshold is only for zero exclusion. Global contour contributions, uniform simplicity, an effective simplicity cutoff and finite-prefix verification remain open. No numerical root scan, priority or assurance upgrade.'}
ledger['open_analytic_targets']['effective_eventual_simplicity_exponent_six_20261117'].update({
    'completed_localization_claim':analytic_ids[0],
    'completed_localization':'ALL n>=10^6, every root is in1<|w|<2800, excluding both saddle critical values; all formal saddle data uniformly separated with explicit constants.',
    'remaining_gate':'Global justified contour decomposition on this compact annulus, effective amplitude and differentiated-remainder bounds at all root cancellations, and any region with three or more comparable contributions. Verify strict U+V+E transversality or disjoint one-zero contours with an effective simplicity cutoff, then exact finite prefix.',
    'cutoff_boundary':'10^6 is an inner zero-exclusion cutoff only, not eventual simplicity or a finite-prefix audit.'})
ledger.setdefault('internal_proof_extensions',{})['effective_eventual_simplicity_and_all_index_exponent_six']={
    'canonical_source':SIMPLE_PROOF,'claims':simple_ids,
    'evidence':{name:sha(HERE/name) for name in ('EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex','verify_exponent_six_simplicity.py','exponent_six_prefix_certificates.json','exponent_six_simplicity_audit.json','exponent_six_simplicity_verification.json')},
    'analytic_scope':'Every n>=10000; ALL roots in60<|w|<240; exact split at3, full Gaussian/contour tails, analytic relative error and Cauchy derivative bounds, endpoint Watson bounds, log-derivative gap>n/288-2.',
    'finite_prefix':'n1 linear, all9998 indices2..9999 certified by explicit reconstructed modular Bezout identities of degree<=1 plus registered parity/shifted-Newton lemmas. No numerical root scan.',
    'conclusion':'Ordinary Gamma shape1 A+B*Y^6 satisfies Strong Factorial at EVERY starting index. Actual effective simplicity cutoff10000, not merely an annulus cutoff.',
    'boundary':'Unrestricted SFC, general higher-dimensional FC, other exponents and Gamma shapes remain open. Internal proof only; no novelty, independent/formal/external assurance upgrade.'}
completed=ledger['open_analytic_targets'].pop('effective_eventual_simplicity_exponent_six_20261117')
completed.update({'status':'COMPLETED for ordinary exponent-six two-term family',
    'canonical_source':SIMPLE_PROOF,'claims':simple_ids,'effective_simplicity_cutoff':10000,
    'finite_prefix_complete':True,'remaining_gate':'No remaining obligation for this special-family theorem; unrestricted FC/SFC is the overarching open goal.'})
ledger.setdefault('completed_analytic_targets',{})['effective_eventual_simplicity_exponent_six_20261117']=completed
ledger['open_analytic_targets']['unrestricted_FC_SFC_after_exponent_six']={
    'status':'OPEN; overarching goal active',
    'next_gate':'Extend the exact endpoint/large-saddle decomposition to other composite exponents and seek a valid reduction or new argument for general polynomial supports.',
    'scope_boundary':'All-index ordinary exponent six is completed; even all constant-plus-one-monomial families would not alone close unrestricted SFC or general FC.'}
ledger['internal_proof_extensions']['effective_saddle_annulus_exponent_six']['boundary']='The10^6 threshold is annulus-only. The subsequent effective simplicity proof closes ordinary exponent six at ALL indices with actual cutoff10000 and complete exact prefix; other supports/shapes remain outside these claims.'
ledger.setdefault('internal_proof_extensions',{})['all_ordinary_univariate_constant_plus_power_SFC']={
    'canonical_source':UNIFORM_PROOF,'claims':uniform_ids,
    'evidence':{name:sha(HERE/name) for name in ('ALL_UNIVARIATE_TWO_TERM_SFC.tex','uniform_two_term_sfc.py','verify_uniform_two_term_sfc.py','all_exponent_prefix_certificates.json','uniform_two_term_sfc_audit.json','uniform_two_term_sfc_verification.json')},
    'analytic_scope':'n>=10000,6<=m<=sqrt(n): normalize by n^m*(m/e)^m, exact split at m/2, full-disk single-saddle error22/sqrt(m*n), Cauchy derivative bounds, endpoint Watson estimate and log-derivative gap>n/3-4. Complete all-root annulus1/2<|v|<21/10.',
    'all_exponent_bridge':'2^m>=4n^2 already gives squarefreeness; its complement has m<sqrt(n) for n>=10000. Prime powers and m1 retained.',
    'finite_prefix':'ALL92935 further necessary pairs below10000 at exponents10,12,14,15,18,20,21,22,24,26,28 have exact checked Bezout certificates.92900 kernels degree<=35 and35 full reductions index<=126. m6 and prime-power cases reused.',
    'conclusion':'Every ordinary A+B*Y^m satisfies SFC at all starts for every positive exponent. For this sequence the exceptional sets are EMPTY, E_m=D_all=0, independent of MR inputs.',
    'boundary':'Only ordinary one-variable constant-plus-power supports. Both-nonconstant pairs, multivariate monomials, other Gamma shapes and unrestricted FC/SFC remain separate. No priority/independent/formal/external assurance upgrade.'}
ledger.setdefault('completed_analytic_targets',{})['uniform_all_positive_univariate_power_exponents']={
    'status':'COMPLETED for ordinary A+B*Y^m at every exponent and starting index',
    'canonical_source':UNIFORM_PROOF,'claims':uniform_ids,'effective_uniform_simplicity_cutoff':10000,
    'finite_prefix_complete':True,'ordinary_exceptional_set_empty':True,
    'remaining_gate':'None for this special family; broader conjectures remain open.'}
ledger['open_analytic_targets'].pop('unrestricted_FC_SFC_after_exponent_six',None)
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families']={
    'status':'OPEN; overarching goal active',
    'next_gate':'Audit retained two-monomial results, then formulate adjacent obstructions for genuinely new supports, including varying-shape f=x^a(A+B*x^d), or another route to arbitrary-polynomial FC/SFC.',
    'scope_boundary':'All ordinary constant-plus-one-variable-power families are completed. This is not unrestricted SFC or general FC, and the fixed-shape Appell identity does not automatically extend to shape a*n+1.'}
ledger['open_bridge_assessments']['short_block_failure_propagation_20261117']['current_counted_family_status']='The subsequent uniform pointwise theorem now makes D_all=0 for F_(n,m). No propagation theorem was used. The old density-to-all-index gate for that counted family is closed directly; no conversion of a general SFC block into this family is established.'
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'retained_source_review':'Rocks2015 author manuscript Section3.4 Proposition3.14 covers a0,1,2 for X^a(A+B*X); Proposition3.15 and Remarks3.18-3.19 are partial for general a>=3. Not a full final2017 literature/priority audit.',
    'next_gate':'Prove and verify the Pade-denominator matrix strip and degree-(a-1) cross-determinant identity for adjacent B_(n,a), then seek all-index a3,4 closure using explicit Hurwitz auxiliary polynomials; no such upgrade registered yet.',
    'exact_reduction':'Q_(an,n)(z)=(-1)^n*(an)!/((a+1)n)!*z^n*B_(n,a)(-1/z), where B=sum binom(n,k)*(an+1)_k*t^k.',
    'failed_shortcuts':['H_(3,2)=z^2+14z+168 has no real roots, refuting blanket real-rootedness',
        'H_(5,2)=z^4+24z^3+360z^2+3960z+47520 fails the quartic Hurwitz determinant, refuting blanket unshifted Hurwitz stability'],
    'possible_uniform_target':'Every H_(a,n) root has real part<(a-1)*n+2; not proved.'})
ledger.setdefault('internal_proof_extensions',{})['Pade_adjacent_nonconstant_monomial_separation']={
    'canonical_source':PADE_PROOF,'claims':pade_ids,
    'evidence':{name:sha(HERE/name) for name in ('PADE_ADJACENT_MONOMIALS.tex','pade_monomial_separation.py','verify_pade_monomial_separation.py','pade_shift_positive_certificates.json','pade_monomial_separation_audit.json','pade_monomial_separation_verification.json','pade_separation_diagnostics.json')},
    'universal_scope':'Exact exponential Pade reduction, strict denominator strips for N>=2, cross determinant auxiliary degree a-1, positive Schwarz reconstruction criterion and an explicit necessary bidiagonal pencil, all proved algebraically.',
    'closed_supports':'Ordinary X^a(A+B*X), every starting index for a1..6. a0..2 already in retained source/previous results; a3..6 now internally proved. Also all indexed terminal pairs a>=2*n-3.',
    'verification':'Symbolic n-2 numerator/denominator positivity proves ALL n>=2, not a finite-index scan. Matrix recurrence and exact cross-determinant checks; finite diagnostics outside both registered ranges are only individual windows.',
    'source_boundary':'Rocks2015 author draft Section3.4 Proposition3.14 and Proposition3.15/Remarks3.18-3.19 inspected for scope. No full final2017 or priority audit; proof here rederives its own algebraic inputs.',
    'boundary':'General a beyond the indexed range, gaps d>1, arbitrary supports/Gamma shapes, general FC and unrestricted SFC remain open. Uniform spectral separation is sufficient, not equivalent to SFC. No independent/formal/external assurance or novelty upgrade.'}
ledger.setdefault('completed_analytic_targets',{})['adjacent_nonconstant_monomials_a1_through6']={
    'status':'COMPLETED for ordinary X^a(A+B*X), a1..6, every starting index',
    'canonical_source':PADE_PROOF,'claims':pade_ids,'all_a_claimed':False,
    'indexed_unbounded_addition':'a>=2*n-3 at terminal n>=2'}
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'status':'OPEN; overarching goal active; the a1..6 adjacent-support gate is now completed',
    'next_gate':'Prove a uniform separation/positivity argument for a>=7,n>(a+3)/2, using H_(a,n) or the explicit necessary pencil A+zB, rather than only extending finite scans. Continue toward arbitrary polynomial supports and general FC.',
    'proved_in_this_stage':'All-n a1..6 and every indexed a>=2*n-3. Exact symbolic certificates in n-2, no index extrapolation.',
    'uniform_matrix_target':'Re eigen(-B^-1*A)<(a-1)*n+2, with A upper bidiagonal diag a*(n-1)+j+n and upper -(a*(n-1)+j+1), B lower bidiagonal diag1 except last1/(a+1), subdiag-1; dimension a-1. Sufficient, unproved for all a.',
    'scope_boundary':'Neither the uniform shifted/spectral target nor generic all-a SFC is proved. Failed sufficient certificates do not refute SFC. Full FC/SFC and other supports remain open.'})
for cached in ledger.get('cached_primary_sources',[]):
    if cached['file']=='sources/rocks-2015-author-preprint.txt':
        review=' Section3.4 two-monomial definition, Proposition3.14 a0,1,2 and Proposition3.15/Remarks3.18-3.19 general-a partial scope inspected for the Pade route; the new proof rederives all algebraic inputs. No final2017 or priority claim.'
        if review not in cached['review']:cached['review']+=review
ledger.setdefault('internal_proof_extensions',{})['uniform_adjacent_monomials_shared_parabola']={
    'canonical_source':PARABOLA_PROOF,'claims':parabola_ids,
    'evidence':{name:sha(HERE/name) for name in ('UNIFORM_PADE_PARABOLA_SFC.tex','pade_parabola_sfc.py','verify_pade_parabola_sfc.py','pade_parabola_sfc_audit.json','pade_parabola_sfc_verification.json')},
    'conclusion':'ALL a>=0, all starts and all positive Gamma shapes for X^a(A+B*X); other monomial factors in independent variables permitted.',
    'universal_mechanism':'All R_N roots INSIDE the same parabola where every necessary shift continuant is ZERO-FREE. Complete algebraic proof is self-contained.',
    'verification_checks':parabola_receipt['total'],
    'next_gap_identities':'Equations10-13 give exact Pascal, raising and order-d+1 parameter recurrences; no nonadjacent SFC upgrade. Only two boundary constraints are available, so a direct finite-pencil extension is underdetermined.',
    'boundary':'Gaps d>1, arbitrary supports and unrestricted FC/SFC remain open. Independent/formal/external assurance false.'}
ledger.setdefault('completed_analytic_targets',{})['uniform_adjacent_monomials_all_shapes']={
    'status':'COMPLETED for every a>=0, every start and every positive Gamma shape',
    'canonical_source':PARABOLA_PROOF,'claims':parabola_ids,'all_a_proved':True,
    'vertical_bound_needed':False,'finite_prefix_needed':False}
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families']={
    'status':'OPEN; overarching goal active; all-a adjacent supports completed',
    'canonical_source':BASE+'STATUS_AND_LITERATURE.tex',
    'next_gate':'New gaps d>=2,a>=1 in X^a(A+B*X^d), or more general/multivariate supports. Derive recurrences and compare zero-free regions; check retained literature.',
    'completed_gate':'ALL a>=0 and positive Gamma shapes for X^a(A+B*X), by shared-parabola separation.',
    'retained_optional_target':'The stronger vertical spectral bound remains unproved, but is no longer an obligation for the completed family.',
    'scope_boundary':'No reduction of unrestricted FC/SFC to completed support families is established.'}
ledger['internal_proof_extensions']['Pade_adjacent_nonconstant_monomial_separation']['later_uniform_closure']=parabola_ids
parabola_source=dict(file='sources/saff-varga-1976-sharpness.pdf',
    url='https://www.math.kent.edu/~varga/pub/paper_98.pdf',
    sha256=sha(HERE/'sources/saff-varga-1976-sharpness.pdf'),
    review='Theorem1.1 displayed recurrence and parabola compared. Needed invariant fully reproved internally. No full-paper, dependency or priority audit.')
ledger['cached_primary_sources']=[r for r in ledger['cached_primary_sources'] if r['file']!=parabola_source['file']]+[parabola_source]
ledger.setdefault('internal_proof_extensions',{})['dyadic_binomial_unit_transfer']={
    'canonical_source':DYADIC_PROOF,'claims':dyadic_ids,
    'evidence':{name:sha(HERE/name) for name in ('DYADIC_BINOMIAL_SFC.tex','dyadic_binomial_sfc.py','verify_dyadic_binomial_sfc.py','dyadic_binomial_sfc_audit.json','dyadic_binomial_sfc_verification.json','dyadic_binomial_frontier_diagnostics.json')},
    'conclusion':'Every q=2^r,integer1<=beta<=2q gives X^q(A+B*X^q) SFC at ALL starts; arbitrary complex coefficients; independent-variable monomial factors allowed.',
    'universal_mechanism':'Degree-preserving mod2 Bessel/Fibonacci reductions, arbitrary rational unit sequences times (L+h)!,h0,1; normalized resultants are odd.',
    'exact_resultant_valuation':'3q*n*(n-1)-(n-1)*s2(n)-n*s2(n-1), independent of beta in the proved integer range',
    'verification_checks':dyadic_receipt['total'],
    'source_comparison':'Rocks2015 author-draft Lemma3.8 unshifted unit identity compared; all ingredients proved here; no full final-publication or priority audit.',
    'boundary':'Not all lower exponents/gaps/shapes. h2 and a1,d2,terminal3 first-residue failures are retained limitations, not counterexamples. General FC/unrestricted SFC and assurance upgrades unclaimed.'}
ledger.setdefault('completed_analytic_targets',{})['dyadic_equal_lower_exponent_and_gap']={
    'status':'COMPLETED for all q=2^r,all starts,integer shapes1..2q',
    'canonical_source':DYADIC_PROOF,'claims':dyadic_ids[:2],'finite_prefix_needed':False}
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'next_gate':'First remaining ordinary gap X*(A+B*X^2), then other nonadjacent pairs outside the completed dyadic diagonal. Seek a uniform higher-precision/unit or spectral invariant; first-mod2 residuals overlap at terminal3 and the order3 recurrence has only two known boundaries.',
    'additional_completed_gate':'X^q(A+B*X^q),q=2^r,all starts,integer Gamma shapes1..2q, by odd-resultant reduction.',
    'boundary_evidence':BASE+'dyadic_binomial_frontier_diagnostics.json',
    'scope_boundary':'No reduction of arbitrary FC/SFC to these completed classes; arbitrary shapes and remaining gaps are not inferred.'})
ledger['internal_proof_extensions']['dyadic_binomial_unit_transfer']['additional_frontier_theorems']={
    'claims':dyadic_ids[2:],
    'lower_one_equivalence':'Every gap e: A_n=T+(e/n)*theta*T and A_prev=T-(1/n)*theta*T. SFC for X*(A+B*X^e) iff each tilted Gamma-shape-n polynomial T_(n,e) is simple; simplicity itself unproved.',
    'cubic_prime_confinement':'Prime n<q<=3(n-1)/2,n>=6: all common-factor roots are q-units. A rational square divisor of T_(n,2) has even degree<=2floor((q-n)/8) and reduces degree-preservingly into K_(q-n).',
    'kernel_relation':'K_d=(H_d+2H_(d+1))/3,theta*K_d=d*(H_(d+1)-H_d)/3; residue collisions allowed, uniform coverage unproved.',
    'scope_boundary':'No all-index cubic SFC or prime-existence/coverage theorem claimed.'}
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'next_gate':'First remaining ordinary X*(A+B*X^2): prove all-index simplicity of T_(n,2)=sum binom(n,k)*(n)_(2k)*t^k, or uniform coverage of rational square-factor degrees by admissible prime/residual kernels. Use the quarter-gap even bound and full K_d factor-degree sets; allow the prime to vary with candidate degree.',
    'new_scalar_structure':'The lower-one exact two-moment/multiple-root equivalence works for every positive gap; no fixed-shape Appell inference. The old parameter boundary deficit still applies to general lower exponents.',
    'prime_window_boundary':'The confinement is conditional on prime n<q<=3(n-1)/2,n>=6. It proves no uniform existence/coverage assertion; finite window replays are bookkeeping.',
    'scope_boundary':'Dyadic diagonal and adjacent families completed; the cubic all-index simplicity/coverage and unrestricted FC/SFC remain open.'})
ledger.setdefault('internal_proof_extensions',{})['tilted_cubic_prime_power_polygons_and_folds']={
    'canonical_source':TILTED_PROOF,'claims':tilted_ids,
    'evidence':{name:sha(HERE/name) for name in ('TILTED_CUBIC_PRIME_POWER_AND_FOLD.tex','tilted_cubic_frontier.py','verify_tilted_cubic_frontier.py','tilted_cubic_frontier_audit.json','tilted_cubic_frontier_verification.json','tilted_cubic_frontier_diagnostics.json')},
    'universal_shift':'All a>=1,i>=0,r>=1,p>3max(a,i): exact horizontal/long edge; residuals T_i and either T_a at i0 or Gamma shape a+1 at i>0.',
    'completed_indexed_scope':'ALL prime-power terminal INDICES: irreducible over Qp. Ten complete shifted index families by parity and nonsquare leading coefficient. These are not all cubic starting indices.',
    'small_prime_detail':'p2 has primitive height2P-1, p3 heightP-1; separate all-r digit arguments.',
    'analytic_scope':'Exact saddle equation/discriminant, two finite cubic folds and their nonzero third/cross derivatives. No contour, Airy estimate, root proximity or effective cutoff proved.',
    'diagnostic_scope':'Fixed primes2..11 fail at732,13 completes that index; candidate-prime guesses fail at51/118. K11 has degree5 but exactly3 real/positive roots and is squarefree, disproving universal kernel real-rootedness. No all-primes limitation or SFC counterexample inferred.',
    'verification_checks':tilted_receipt['total'],
    'boundary':'General cubic simplicity/prime coverage and general FC/unrestricted SFC remain open; all priority/assurance upgrades unclaimed.'}
ledger.setdefault('completed_indexed_targets',{})['cubic_prime_power_indices_and_ten_shifts']={
    'status':'COMPLETED only for the exact unbounded index families stated',
    'canonical_source':TILTED_PROOF,'claims':tilted_ids[:2],
    'whole_cubic_all_starting_indices_closed':False}
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'completed_indexed_gate':'All prime-power indices and ten shifted families for the cubic; see FIRTILTEDPRIMEPOWER20261117. Not a whole-family SFC theorem.',
    'prime_method_boundaries':'Fixed2..11 fails at n732,13 works. p>2k+1 or p>3k with nmodp<2k is insufficient; retained exact counterprofiles.',
    'analytic_fold_data':'w+-=1/6+-i*sqrt3/18; Phi_sss and Phi_sw nonzero with explicit moduli. Exact geometry only; contour contributions and uniform all-root errors remain open.',
    'next_gate':'Uniform simplicity of T_(n,2) beyond the completed index families, via a valid arithmetic prime/residual coverage theorem or effective all-root contour/fold control plus a finite prefix. Preserve the full FC/SFC goal.'})
source_row={'file':'sources/wagner-multivariate-stable-polynomials-2011.pdf',
    'url':'https://www.math.uwaterloo.ca/~dgwagner/ceb_wagner.pdf',
    'sha256':sha(HERE/'sources/wagner-multivariate-stable-polynomials-2011.pdf'),
    'review':'Theorem4.1 statement and Section4 proof context inspected. Only classical half-plane polarization is imported, with full-degree/convex-region hypotheses checked for the explicit Walsh convolution. No full-paper or priority audit.'}
if source_row not in ledger['cached_primary_sources']:
    assert not any(row['file']==source_row['file'] for row in ledger['cached_primary_sources'])
    ledger['cached_primary_sources'].append(source_row)
ledger['published_input_pinpoints']['GWS_half_plane_for_cubic_adjoint_20261117']='Wagner Theorem4.1, PDF page11; open right half-plane, degree n, symmetric multiaffine polarization. All hypotheses explicit in CUBIC_ADJOINT_STABILITY_SFC.tex. No modern asymptotic theorem is used.'
ledger['internal_proof_extensions']['cubic_all_indices_by_adjoint_stability']={
    'canonical_source':CUBIC_PROOF,'claims':cubic_ids,
    'evidence':{name:sha(HERE/name) for name in cubic_receipt['files']},
    'simplicity_scope':'ALL n>=1,real gamma>0,gamma>=n-2',
    'SFC_scope':'ALL X*(A+B*X^2), every start and every positive Gamma shape; independent-variable monomial factors allowed',
    'mechanism':'Gamma Jacobi half-plane bound; explicit Walsh convolution auxiliary; classical GWS; exact third-order adjoint; both zero concomitants; disjoint B(z)/B(-z) roots.',
    'verification_checks':cubic_receipt['total'],'finite_prefix_required':False,
    'boundary':'General lower exponents/gaps, unrestricted SFC and general FC remain open. Internal proof only; all novelty/independent/formal/external assurance flags false.'}
ledger.setdefault('completed_analytic_targets',{})['whole_cubic_SFC_all_positive_shapes']={
    'status':'COMPLETED at all indices by algebraic adjoint/stability argument; no effective analytic cutoff needed',
    'canonical_source':CUBIC_PROOF,'claims':cubic_ids,'finite_prefix_needed':False,'all_positive_Gamma_shapes':True}
ledger['internal_proof_extensions']['tilted_cubic_prime_power_polygons_and_folds']['subsequent_whole_cubic_closure']=cubic_ids[1]+'; preceding open boundary is this earlier-stage method scope only.'
ledger['completed_indexed_targets']['cubic_prime_power_indices_and_ten_shifts']['later_whole_family_closure']=cubic_ids[1]
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'status':'OPEN; overarching FC/SFC goal active; whole cubic support now closed',
    'completed_cubic_gate':'Every X*(A+B*X^2), every start and positive Gamma shape; explicit adjoint/stability proof, no finite prefix.',
    'next_gate':'Genuinely new lower exponents/gaps: derive higher-order adjoint obstructions for X*(A+B*X^d),d>=3, or gap2 with other lower exponents beyond completed shapes; seek a valid reduction toward arbitrary supports.',
    'scope_boundary':'Arithmetic/fold coverage is not needed for the completed cubic theorem. It remains unproved as a separate stronger method target. No unrestricted FC/SFC conclusion.'})
ledger['internal_proof_extensions']['all_gap_explicit_adjoints_and_exact_rank']={
    'canonical_source':GAP_PROOF,'claims':gap_ids,
    'evidence':{name:sha(HERE/name) for name in gap_receipt['files']},
    'scope':'ALL d>=2,n>=1,gamma>0: polynomial adjoints, exact rank equivalence and full determinant reconstruction. ALL beta>0 lower-one consecutive moment bridge.',
    'mechanism':'Formal residue/Lagrange inversion gives a first-order auxiliary recurrence; divergence adjoint and local gauges; rational concomitants; companion cofactor identity.',
    'gap_three_conditional_target':'No common nonzero root of the three cyclic Wronskian sections. Strict real interlacing or CLOSED pi/3 sector exclusion would suffice; neither uniform statement established.',
    'failed_stronger_target':'Full Hurwitz stability of W refuted at n4,gamma4 by exact negative seventh determinant. This is a method limitation; the actual moment is squarefree.',
    'finite_evidence':'Exact Sturm/simple-negative-root/interlacing controls ONLY n=gamma4,6,8. Nonpositive/complex shape controls detect actual repeated roots outside positive Gamma scope.',
    'verification_checks':gap_receipt['total'],'all_gap_initial_four_orders':'ALL d>=1,gamma>0, n1..4 simple; strict third central moment via rising-factorial convexity',
    'boundary':'No all-index SFC at d>=3, general FC or unrestricted SFC closure. No imported literature theorem in this extension; all novelty/independent/formal/external assurance flags false.'}
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'completed_higher_gap_structure':'FIRGENERALGAPADJOINT20261117 and FIRGENERALGAPRANK20261117: all-gap adjoints, exact rank equivalence, determinant reconstruction and all-positive-shape lower-one bridge.',
    'next_gate':'At gap3, prove uniform coprimality/interlacing of cyclic Wronskian sections or CLOSED pi/3 sector zero-freeness at gamma=n+beta-1. The full Hurwitz shortcut is false already n4,gamma4. Continue toward arbitrary supports, not only finite certificates.',
    'scope_boundary':'Whole cubic remains closed. Higher-gap simplicity/rank exclusion and unrestricted FC/SFC remain open; finite interlacing controls are not an infinite family proof.'})
ledger.setdefault('external_followup_leads',{})['generalized_Hurwitz_sector_2015']={
    'url':'https://arxiv.org/abs/1506.07379',
    'authors':'Holtz, Khrushchev, Kushel',
    'status':'PRIMARY ABSTRACT LEAD ONLY; not an input to the internally proved adjoint/rank theorem',
    'reported_theorem':'Total nonnegativity of generalized Hurwitz H_M excludes roots in the OPEN sector |arg z|<pi/M.',
    'application_boundary':'No uniform H3 positivity for the current W has been verified. Closed-sector exclusion or a separate boundary argument would still be required. Not a full-paper or dependency audit.'}
ledger['internal_proof_extensions']['gap_three_Wronskian_recurrence_and_one_real_section']={
    'canonical_source':WR_PROOF,'claims':wr_ids,
    'evidence':{name:sha(HERE/name) for name in wr_receipt['files']},
    'all_parameter_structure':'B_(gamma-1)=B-Bprime; divided difference of consecutive B; exact fixed-shape W recurrence; centered rational recurrence and nonnegative z,b coefficients.',
    'weaker_uniform_target':'Common real section roots are impossible for positive Gamma. ONE nonzero real-rooted section or real constant-linear combination suffices; choice can vary with n/gamma. Existence at all indices is open.',
    'pointwise_conclusion':'ALL gamma>0, T_(5,gamma,3) simple by an explicit quadratic section and strictly positive discriminant coefficients. Reusing all-gap n<=4 gives gap3 initial n1..5 for all positive shapes.',
    'exact_method_boundary':'Separate monic positive strictly log-concave section-interlacing inputs do not guarantee full section interlacing after a formal recurrence step. Explicit section u^2+1772u+900000, discriminant-460016. Inputs are not the actual common-Gamma pair; no SFC counterexample.',
    'finite_controls':'Selected n6,10,12 near gamma=n-1 and at10n have exact negative-root/interlacing checks; bounded evidence only.',
    'verification_checks':wr_receipt['total'],
    'boundary':'No all-index gap3 simplicity, arbitrary support, general FC or unrestricted SFC closure. No imported external theorem; all novelty/independent/formal/external assurance flags false.'}
ledger['internal_proof_extensions']['all_gap_explicit_adjoints_and_exact_rank']['subsequent_gap_three_refinement']=wr_ids
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'completed_gap_three_recurrence_gate':'Exact Wronskian/Gamma-shift recurrences and bivariate coefficient positivity, plus the one-real-section criterion and all-positive-shape initial n<=5.',
    'next_gate':'For ALL n and positive beta, prove existence of one nonzero real-rooted cyclic section or real constant-linear combination of the gap3 Wronskian at gamma=n+beta-1. The choice can vary. All-section interlacing/closed-sector assertions are stronger than necessary. Exploit actual Gamma-shift compatibility under the new recurrence; separate-input interlacing/log-concavity is insufficient for full section preservation.',
    'scope_boundary':'Whole cubic support remains closed; gap3 is closed only through terminal5 plus individual certificates. Uniform higher-index section existence and general FC/unrestricted SFC remain open.'})
pc_source={'file':'sources/iskander-jain-2020-delta-laguerre-polya.pdf',
    'url':'https://arxiv.org/pdf/2002.09906',
    'sha256':sha(HERE/'sources/iskander-jain-2020-delta-laguerre-polya.pdf'),
    'review':'arXiv v1 Theorems1.1/1.2 and finite-interval versus infinite-ray zero-free hypotheses inspected. The real binomial companion is proved independently by polar derivatives. Infinite-ray theorem cannot be applied to backward sampling with positive roots. No full-paper/dependency/priority audit or imported simplicity conclusion.'}
if pc_source not in ledger['cached_primary_sources']:
    assert not any(row['file']==pc_source['file'] for row in ledger['cached_primary_sources'])
    ledger['cached_primary_sources'].append(pc_source)
ledger['internal_proof_extensions']['gap_three_polar_companion_and_backward_certificate']={
    'canonical_source':PC_PROOF,'claims':pc_ids,
    'evidence':{name:sha(HERE/name) for name in pc_receipt['files']},
    'uniform_companion_scope':'ALL n>=1,gamma>n-2: simple real roots below-1 for J and K, by exact finite-difference/polar/Mobius identities.',
    'missing_transfer':'Generic inverse-transform counterprofile at n8 has no individually real-rooted W section despite simple companion roots below-1. The actual arithmetic-progression Gamma structure or another invariant is still needed. Real-rooted combinations are not refuted.',
    'backward_certificate_scope':'Any normalized degree-bounded triple. Nonzero pivots alone prove gcd1 by descent to1. Positive pivots have a weighted interval representation. Zero pivots are inconclusive; formal symbolic pivots retain nonzero conditions.',
    'finite_controls':'Symbolic all-positive-gamma weights at FIXED n2,3, numeric ordinary n8,24 certificates with tamper rejection, and ordinary n24,48 Sturm diagnostics. No index extrapolation.',
    'native_prime_boundary':'At n22,p3 edge residual1+Y^3 repeats. At n63,p3 edge denominator2,residual1+2Y^3 permits square degree2 despite parity. No SFC failure or all-prime limitation.',
    'source_audit':'Iskander-Jain finite/infinite zero-free hypotheses distinguished. No external theorem imported for the new proof.',
    'verification_checks':pc_receipt['total'],
    'boundary':'Uniform actual pivot nonvanishing/section existence, all-index gap3 SFC, general FC and unrestricted SFC remain open; all novelty/independent/formal/external assurance flags false.'}
ledger['published_input_pinpoints']['Iskander_Jain_finite_infinite_zero_free_boundary']='Theorems1.1/1.2, cached arXiv v1 page3. Finite interval supports the binomial companion; backwards Poisson sampling has positive roots on its forward ray and violates the infinite-ray hypothesis. Needed polar theorem independently proved.'
ledger['open_analytic_targets']['general_supports_after_all_constant_power_families'].update({
    'completed_polar_companion_gate':'Uniform simple companion roots below-1 at gamma>n-2, and exact backward nonzero-pivot coprimality certificates. No actual-family transfer/pivot nonvanishing theorem.',
    'next_gate':'Exploit the arithmetic-progression polar factors to obtain a valid real-section/real-combination transfer, or prove actual-family backward pivots nonzero. Positivity is stronger than necessary; zero pivots are inconclusive. Keep one-section target available rather than requiring all-section interlacing.',
    'scope_boundary':'All prior completed support families remain closed. Uniform gap3 transfer or pivot theorem and general FC/unrestricted SFC remain open; generic inverse controls and native-p3 failures are method limitations, not conjecture counterexamples.'})
after=(json.dumps(ledger,indent=2)+'\n').encode()
if after!=before:
    assert ledger_path.read_bytes()==before,'Concurrent source edit'
    temp=ledger_path.with_suffix('.all-exponent-tmp');temp.write_bytes(after);os.replace(temp,ledger_path)
print(json.dumps({'registered':ids+analytic_ids+simple_ids+uniform_ids+pade_ids+parabola_ids+dyadic_ids+tilted_ids+cubic_ids+gap_ids+wr_ids+pc_ids,'checks':receipt['total']+annulus_receipt['total']+simplicity_receipt['total']+uniform_receipt['total']+pade_receipt['total']+parabola_receipt['total']+dyadic_receipt['total']+tilted_receipt['total']+cubic_receipt['total']+gap_receipt['total']+wr_receipt['total']+pc_receipt['total'],'uniform_polar_companion_proved':True,'uniform_nonzero_pivots_proved':False,'gap_three_one_real_section_criterion_proved':True,'gap_three_all_positive_shapes_fifth_order_proved':True,'higher_gap_adjoint_rank_proved':True,'higher_gap_uniform_simplicity_proved':False,'whole_cubic_all_starts_all_positive_shapes':True,'tilted_cubic_all_prime_power_indices_irreducible':True,'dyadic_equal_exponent_gap_all_starts':True,'uniform_all_a_all_positive_shapes':True,'adjacent_nonconstant_a1_through6_all_starts':True,'ordinary_all_positive_power_exponents_all_starts':True,'ordinary_exponent_six_all_starts':True,'effective_simplicity_cutoff':10000,'exponent_range':'ALL positive integers',
    'overlap':'2^144, not an analytic threshold','prime_power_result_reused':True,'all_index_SFC':'open'}))
