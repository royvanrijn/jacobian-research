"""Register the direct quantitative theorem for every fixed positive exponent."""
from pathlib import Path
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='factorial/image-rigidity/'
PROOF=BASE+'QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex'
CHECKER=BASE+'verify_quantitative_squarefree.py'
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
        temp=path.with_suffix(path.suffix+'.real-place-research/.scratch')
        temp.write_bytes(text.encode())
        os.replace(temp,path)


if any(row['id']=='FIRSIMULTPOWER20261117' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Use register_simultaneous_squarefree.py; retain the growing-exponent result and preprint/published route distinction.')

receipt=json.loads((HERE/'quantitative_squarefree_verification.json').read_text())
assert receipt['status']=='pass' and receipt['all_fixed_positive_integer_exponents']
assert not receipt['old_density_one_theorem_used'] and not receipt['bootstrap_used'] and not receipt['parity_used']
assert not receipt['all_index_SFC_closed'] and not receipt['numerical_constant_or_starting_threshold_extracted']
for file,expected in receipt['files'].items():assert sha(HERE/file)==expected,('stale quantitative receipt',file)
common=dict(kind='theorem',state='proved',canonical_source=PROOF,checker=CHECKER,proof_type='hybrid',
    independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_quantitative_squarefree.py'),
    software_lock=[PROOF,CHECKER,BASE+'quantitative_squarefree.py',BASE+'quantitative_squarefree_verification.json',
        BASE+'quantitative_squarefree_audit.json',BASE+'quantitative_squarefree_examples.json',BASE+'register_quantitative_squarefree.py',
        BASE+'prime_power_windows.py',BASE+'source_ledger.json',
        BASE+'sources/filaseta-rocks-2016-published.pdf',BASE+'sources/filaseta-rocks-2016-published.txt',
        BASE+'sources/matomaki-radziwill-2016-annals.pdf',BASE+'sources/matomaki-radziwill-2016-annals.txt'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],consumers=[],
    invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')
new=[dict(common,id='FIRALLMLOCAL20261116',
    title='Exact shifted Newton valuations and preceding-prime square-divisor bounds hold for every positive factorial exponent',
    scope='For EVERY positive integer m, n=a*p+i with a>=1,i>=0 and prime p>m*max(a,i) has m*n<p^2. For j=b*p+c,0<=c<p, the literal coefficient valuation is m*b+floor(m*c/p)+indicator(c>i), by one-level Legendre and Kummer. Comparing to m*(j-i)/p gives a horizontal edge of width i (absent at i0), then an edge of width a*p,slope m/p and reduced denominator p. The FULL reduction is F_(n,m) modp=F_(i,m) modp. Every rational g with g^2|F_(n,m) has degree u+v*p,2u<=i,2v<=a; degree<p forces wholly horizontal occupancy. Primitive Gauss normalization gives integer quotient and g(0)=+-1, hence leading p-adic unit and degree preservation on reduction. For the preceding prime q_n and d=n-q_n, q_n>m*max(1,d) therefore gives degree(g)<=floor(d/2), since the sloped width q_n is unavailable to a doubled divisor. This condition follows from d<sqrt(n) and n>(m+1)^2. A squarefree horizontal kernel then excludes every nonconstant rational square divisor when the global bound<p. No parity or exponent-six-only input, no claim that every residual is squarefree, no all-index SFC or other Gamma-shape conclusion.',
    dependencies=['FIRPADICDEFORM20261116',
        'external: standard Legendre/Kummer valuations, Gauss lemma and Dumas Newton multiplication, as in Filaseta-Rocks2016 Section2']),
 dict(common,id='FIRQUANTDENSITY20261116',
    title='Every fixed positive factorial exponent has at most X over log-power-one-fiftieth nonsquarefree moment indices',
    scope='For EVERY FIXED positive integer m, E_m(X)=#{1<=n<=X:F_(n,m) not squarefree over Q} satisfies E_m(X)<<_m X/(log X)^(1/50), unconditionally in ordinary factorial shape1. This DIRECT proof uses Filaseta-Rocks Theorem1.1 irreducible-offset lower density>=log2 for the original sequence, NOT the prior exponent-six density-one theorem or bootstrapping. Take H=floor(X^(1/8)),Y=sqrt(2mX). For irreducible offsets i<=H, R_i=Res(F_(i,m),F_(i,m)_prime) is nonzero, and Hadamard gives abs(R_i)<=i^i*(2^i*(m*i)!)^(2i-1), log abs R_i=O_m(i^2 log(i+1)). Counting EVERY progression n=i modp for p>Y dividing R_i bounds the entire bad-reduction union by ((X/Y+1)/logY)*sum log abs R_i=O_m(X^(7/8)); no uniform all-offset good-reduction or C-smoothness theorem is assumed. The smooth indicator has long mean tending to1-log2 for fixed m. Uniform Matomaki-Radziwill Theorem1 with growing H and delta1/4 gives at most O(X^(799/800)*(log X)^(1/3)+X/(log X)^(1/50)) short-interval exceptions. Pigeonhole against more than2H/3 irreducible offsets supplies p>Y at a good kernel. Filaseta-Rocks equation3.1 bounds the bad preceding-gap indices by O(X^(2/3)), including gaps crossing a dyadic endpoint. Outside these sets, the all-m local lemma forces degree(g)<sqrt(n)/2<p, degree-preserving horizontal reduction contradicts the squarefree residual. Boundary O(H) and every power-saving term are absorbed by the logarithmic bound; dyadic summation split at sqrt(X) proves the stated global rate. The exponent1/50 is explicit, but numerical implied constant and starting threshold are NOT extracted; constants depend on FIXED m, with no growing-m uniformity. Simultaneously for all complex nonzero A+B*Y^m, only O_m(X/(log X)^(1/50)) terminal indices<=X can fail the adjacent nonzero-moment window. Infinite exceptional sets, full all-index SFC and uniform all-index residual coverage remain possible/open. No irreducibility-density-one, general Gamma-shape, novelty or independent/formal/external assurance upgrade; deep published analytic inputs are credited imports, not consequences of finite tests.',
    dependencies=['FIRALLMLOCAL20261116',
        'external: Filaseta-Rocks2016 Theorem1.1 and equation3.1 page6, pinned journal Online First author copy',
        'external: Matomaki-Radziwill2016 Annals Theorem1 pp1015-1016, uniform h,delta,f with explicit logarithmic exponent1/50',
        'external: classical PNT, retained DLMF27.2.3 reference; disjoint large-prime multiples and partial summation'])]
ids=[row['id'] for row in new]

def old_density_scope(row):
    addition=(' FIRQUANTDENSITY20261116 now strengthens this exponent-six qualitative conclusion to E_m(X)<<_m X/(log X)^(1/50) for EVERY fixed positive integer m, using direct irreducible-offset input and a counted bad-prime union. The present fixed-H proof is retained as an independent qualitative route; the quantitative proof does not depend on this theorem. Full all-index SFC remains open.')
    if addition not in row['scope']:row['scope']+=addition
    return row

def open_scope(row):
    addition=(' FIRALLMLOCAL20261116 proves the exact r1 shifted valuation, reduction, occupancy and preceding-prime square-divisor bound for EVERY positive integer factorial exponent, with no parity input. FIRQUANTDENSITY20261116 now proves the direct unconditional quantitative bound E_m(X)<<_m X/(log X)^(1/50) for every FIXED positive integer m. Count bad resultant primes in their shift progressions to allow H=X^(1/8): the bad-reduction union costs O_m(X^(7/8)). Original Filaseta-Rocks irreducible-offset density plus uniform Matomaki-Radziwill at delta1/4 gives the logarithmic rate; equation3.1 supplies O(X^(2/3)) prime-gap exceptions. This does NOT bootstrap the earlier density-one theorem. No numerical implied constant or starting threshold is extracted, and uniformity for growing m is unclaimed. The possible quantitative sparse exceptional set and full all-index SFC remain the focused open gap; no all-C smoothness or universal-offset theorem is assumed. Inactive-variable FC stays closed.')
    if addition not in row['scope']:row['scope']+=addition
    row['dependencies']=list(dict.fromkeys(row['dependencies']+ids))
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':open_scope,
    'FIRDENSITYONE20261116':old_density_scope},new)
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+ids))
    row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PROOF,BASE+'quantitative_squarefree_audit.json',BASE+'quantitative_squarefree_verification.json']))
    return row

def lesson(row):
    row=refresh(row)
    extra=(' To allow polynomially growing offsets, count bad reduction primes instead of passing all of them: use integer R_i=resultant(F_i,F_i_prime), Hadamard abs(R_i)<=i^i*(2^i*(mi)!)^(2i-1), and the progression count X/p+1 for each p>Y dividing R_i. With H=X^(1/8),Y=sqrt(2mX), the union costs O_m(X^(7/8)). The general-m exact shifted valuation and preceding-prime bound use no parity. Original irreducible-offset density and uniform MR at delta1/4 then give the direct E_m(X)<<_m X/(log X)^(1/50) rate for every fixed positive integer m, without bootstrapping earlier density one.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Growing H(X) is now justified ONLY with the bad-prime progression union controlled, rather than an uncontrolled replacement in the older fixed-H proof. Keep m FIXED: no variable-m uniformity follows from fixed-m irreducible density or constants. An explicit logarithmic exponent does not mean a numerical implied constant or starting threshold has been extracted. For the rate, use quantitative gap equation3.1, not merely density-one equation3.3. Finite algebra checks do not numerically establish the asymptotic estimate.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'quantitative_squarefree.py',CHECKER]))
    return row

def work(row):
    row=refresh(row)
    row['next_step']='Current focused B gate: remove the possible sparse exceptional indices to prove full ALL-index constant-plus-monomial SFC, especially non-prime-power exponents such as6. The asymptotic result is now quantitative for EVERY FIXED positive integer m: E_m(X)<<_m X/(log X)^(1/50), directly from irreducible offsets, counted bad-resultant prime progressions, uniform MR and quantitative gap equation3.1. Do not reopen density-one coverage, restrict the new theorem to exponent6, or treat growing offsets as uncontrolled after the O_m(X^(7/8)) bad-prime union has been proved. The explicit exponent is closed; numerical implied constants and starting thresholds are unextracted, and growing-m uniformity is unclaimed. The complete bounded-multiplier families and explicit fixed-multiplier cutoffs remain closed. A uniform all-index residual-degree mechanism or another argument must still handle the possible infinite exceptional set. All-C positivity is equivalent to the characteristic-zero all-index problem in the exponent-six sequence, with all-C small-prime support stronger. Inactive-variable FC remains completed; general FC/SFC, nonunit Gamma-shape rate extensions and other programme obligations remain separate open work.'
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
ledger_path=HERE/'source_ledger.json';before=ledger_path.read_bytes();ledger=json.loads(before)
for row in ledger['cached_primary_sources']:
    if row['file'] in ('sources/filaseta-rocks-2016-published.pdf','sources/filaseta-rocks-2016-published.txt'):
        extra=' Quantitative extension also checks equation3.1 page6: sum of consecutive gaps>sqrt(left endpoint), left endpoint<=Z, is O(Z^(2/3)); the containment argument gives O(X^(2/3)) bad preceding-prime indices. Theorem1.1 applies to every FIXED positive integer m. Deep analytic inputs imported, not internally reproved.'
        if extra not in row['review']:row['review']+=extra
    if row['file'] in ('sources/matomaki-radziwill-2016-annals.pdf','sources/matomaki-radziwill-2016-annals.txt'):
        extra=' Quantitative application uses H=floor(X^(1/8)),delta1/4 and Y=sqrt(2mX): uniform h and f permit both dependencies, giving first exception power799/800 and the published log exponent1/50.'
        if extra not in row['review']:row['review']+=extra
ledger['published_input_pinpoints']['FilasetaRocksQuantitativeGaps']='Equation3.1 page6 in the pinned2016 journal Online First source: sum_(p_j<=Z, gap>sqrt(p_j)) gap << Z^(2/3), credited there to K.Matomaki2007. Counting contained bad preceding-prime indices, including an interval crossing the upper endpoint, yields O(X^(2/3)). Published statement checked; deeper proof imported.'
ledger['published_input_pinpoints']['MRQuantitativeAllExponents']='Theorem1 pp1015-1016, h,delta,f uniform, fixed delta1/4,H=floor(X^(1/8)): exception O(X^(799/800)*(log X)^(1/3)+X/(log X)^(1/50)). The original irreducible-offset lower density for each fixed m and a Hadamard/progression union O_m(X^(7/8)) justify growing H without bootstrapping density one.'
ledger.setdefault('internal_proof_extensions',{})['quantitative_all_fixed_exponents']={
    'canonical_source':PROOF,'claims':ids,
    'provenance':'User-supplied quantitative general-m refinement checked against pinned primary theorems. Exact local valuation and no-parity q-bound are rederived; Hadamard and bad-prime progression counting justify polynomially growing offsets. No novelty claim.',
    'evidence':{name:sha(HERE/name) for name in ('QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex','quantitative_squarefree_audit.json','quantitative_squarefree_verification.json')},
    'scope':'E_m(X)<<_m X/(log X)^(1/50) for every FIXED positive integer m, ordinary Gamma shape1; common adjacent nonzero-moment windows outside the same quantitative sparse set.',
    'independence':'Direct from original irreducible-offset density, NOT the previous density-one theorem or bootstrapping; no parity or all-offset C-smoothness assumption.',
    'boundary':'Numerical implied constant and starting threshold unextracted; growing-m uniformity, other Gamma shapes, irreducibility density one and full all-index SFC unclaimed. Infinite exceptional sets remain possible. All assurance flags false.'}
after=(json.dumps(ledger,indent=2)+'\n').encode()
if after!=before:
    assert ledger_path.read_bytes()==before,'Concurrent source edit'
    temp=ledger_path.with_suffix('.quantitative-research/.scratch');temp.write_bytes(after);os.replace(temp,ledger_path)
print(json.dumps({'registered':ids,'checks':receipt['total'],'rate':receipt['quantitative_rate'],
                  'old_density_one_used':False,'all_index_SFC':'open'}))
