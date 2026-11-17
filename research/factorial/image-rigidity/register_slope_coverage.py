"""Register finite slope coverage and proved infinite shifted-index windows."""
from pathlib import Path
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='factorial/image-rigidity/'
PROOF=BASE+'SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex'
CHECKER=BASE+'verify_slope_coverage.py'
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


if any(row['id']=='FIRPRIMEGAP20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Historical slope/large-degree registration: use register_prime_gap_windows.py to preserve the prime-gap and small-kernel scopes.')

receipt=json.loads((HERE/'slope_coverage_verification.json').read_text())
assert receipt['status']=='pass'
for file,expected in receipt['files'].items():
    assert sha(HERE/file)==expected,('stale slope verification',file)
base=json.loads((ROOT/'MATH_STATUS.json').read_text())
common=dict(kind='theorem',state='proved',canonical_source=PROOF,checker=CHECKER,
    proof_type='hybrid',independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_slope_coverage.py'),
    software_lock=[PROOF,CHECKER,BASE+'newton_slope_coverage.py',BASE+'prime_power_windows.py',
        BASE+'mixed_prime_windows.py',BASE+'MIXED_PRIME_RESIDUAL_WINDOWS.tex',
        BASE+'slope_coverage_verification.json',BASE+'slope_coverage_probe.json',
        BASE+'shifted_window_examples.json',BASE+'register_slope_coverage.py',BASE+'source_ledger.json',
        BASE+'sources/nagura-1952-prime-interval.pdf',BASE+'sources/nagura-1952-prime-interval.txt'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[],invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')
new=[dict(common,id='FIRSLOPECAP20261116',
    title='Slope-only square-divisor coverage has exact finite and candidate-dependent prime restrictions',
    scope='For the ordinary factorial polynomial P_(n,d)=sum binom(n,j)*(d*j)!*T^j, d>=1 and n>=2, a Newton edge of length L and reduced slope denominator e permits square-divisor contributions e*t, 0<=t<=floor(L/(2e)). The Minkowski sum over edges is the exact retained slope-only degree budget. All coefficient prime divisors are <=d*n, so every prime above d*n has one horizontal edge and admits every degree0..floor(n/2). Thus exhaust primes<=d*n and intersect with positive EVEN degrees for a complete slope-and-parity diagnostic at a supplied index: a survivor is a limitation of this method that no larger prime can fix, not an SFC counterexample. This cap does NOT apply to residual-multiplicity tests: at n7,d6,p43>42 the horizontal residual certifies squarefreeness despite unrestricted slopes. For a candidate degree k, only primes p<=2*d*k OR prime divisors of n(n-1)...(n-2*k+1) can exclude it by slopes. Otherwise the first2k coefficients are units and a horizontal edge of length>=2k admits k. For d6 quadratics this gives primes<=24 plus divisors of n(n-1)(n-2)(n-3). This is a proved necessary useful-prime restriction, not uniform existence of an excluding prime. Shape1 only for these bounds; no all-index exponent-six closure.',
    dependencies=['FIRPADICDEFORM20261116','FIRREALROOTPARITY20261116']),
 dict(common,id='FIRSHIFTEDINDEX20261116',
    title='Exact shifted-index Newton kernels prove an infinite exponent-six adjacent-window family',
    scope='For d,u>=1, i>=0, prime p>d*max(u,i), r>=1, e=p^r, c=(e-1)/(p-1) and N=u*e+i, the ordinary factorial polynomial has a horizontal edge of length i and edge (i,0)->(N,d*u*c) of reduced denominator e. Normalized residuals are P_(i,d)(T) modp and P_(u,d)((-1)^(d*c)*T) modp, by the integer factorial/binomial quotient and factorial digit units. Square-divisor degrees are h+e*t, h<=floor(i/2),t<=floor(u/2), before further residual/real/global restrictions. For exponent6, EVERY 1<=a<=3, 0<=i<=3, prime p>6*max(a,i), and r>=1 gives SQUAREFREE P_(a*p^r+i,6), hence an adjacent nonzero pair for every nonzero complex A+B*Y^6. This includes the user blocks a2,i0..2,p>12 and a3,i0..1,p>18, with a2,i2 closed by the ordinary factorial leading coefficient and constant1 full-square obstruction. The fixed discriminant gcd of P2,6 and P3,6 is2^10*3^4*5^2, independently closing a2,i3 and a3,i2 for p>18. The uniform degree<3n/10 theorem closes a3,i3 and additionally a4..6,i0..1,p>6a. Both residual kernels and parity are independent of r; symbolic global degree inequalities are checked for every r>=1. The former raw local pattern at n=2*13^r+2 is now excluded globally. A small pattern h=t1 still survives this one-prime/global test at N=4*137^r+3, so that test returns UNKNOWN, not a counterexample. All-ring/Gamma generalizations of the factorial leading or prime-band arguments are not asserted. No all-index or novelty claim.',
    dependencies=['FIRINDEXEDMIXEDSFC20261116','FIRPADICDEFORM20261116','FIRREALROOTPARITY20261116','FIRLARGESQUARE20261116']),
 dict(common,id='FIRSLOPECOVEREXEC20261116',proof_type='exact_symbolic',
    title='Exact slope-only diagnostics and shifted-family exporter keep finite coverage separate from all-index SFC',
    scope='newton_slope_coverage.py exports exact Legendre slope hulls, FULL slope cap p<=d*n, candidate-prime shortlists and shifted residual kernels at N=u*p^r+i. Raw local/even patterns remain visible when full-square and uniform large-degree arithmetic exclude them. The all-r large-degree rule checks (10t-3u)*p^r+10h-3i>=0 for EVERY r>=1 by a nonnegative coefficient and its r1 minimum, under the n>=9 shape1 Nagura theorem. large_degree_certificate returns a literal prime 3n<q<18n/5 and its exact two-edge hull. Independent coefficient/residual, small-kernel discriminant/resultant, factorial-square, Nagura prime, rational gcd and shifted-family controls validate transcription. Every finite n2..64 has complete slope coverage; no numerical range extension is used for a uniform claim. Uniform families a1..3,i0..3,p>6max(a,i), and a4..6,i0..1,p>6a are backed by the written proofs. Source pins retain the imported 1952 prime-interval theorem. The all-index exponent-six theorem and lower-degree candidate coverage remain open; assurance and novelty flags stay false.',
    dependencies=['FIRSLOPECAP20261116','FIRSHIFTEDINDEX20261116','FIRMULTIPRIMEWINDOW20261116','FIRLARGESQUARE20261116']),
 dict(common,id='FIRLARGESQUARE20261116',
    title='A uniform prime choice bounds every factorial square-divisor degree below three tenths and excludes full squares',
    scope='For EVERY n>=9 in the ordinary factorial exponent-six polynomial P_n, any rational H with H^2 dividing P_n satisfies deg(H)<3n/10. Nagura1952 supplies 3n<q<18n/5. Since q>n and 6n<2q, coefficient valuations are0 up to L=floor((q-1)/6) and1 after L, giving a horizontal edge then a height-one edge of width n-L and denominator n-L. Doubled occupancy cannot use that second edge, so deg(H)<=floor(L/2)<3n/10. Thus every candidate at least3n/10 has a UNIFORM excluding prime, leaving lower-degree coverage open. Independently, for every ordinary factorial P_(n,d), n>=2,d>=1, a square divisor of half the total degree is impossible: constant1 normalizes a full scalar square to a rational square, whereas (d*n)! is nonsquare by Bertrand. This closes the user-supplied 2p^r+2 block and, together with shifted degree patterns, the additional 2p^r+3 and 3p^r+2,+3 blocks. Shape1 is essential to the current proof; no replacement for general Gamma shapes is claimed. The Nagura theorem is a credited, pinned published input, not independently reproved here. No all-index SFC conclusion or novelty claim.',
    dependencies=['FIRPADICDEFORM20261116','external: Nagura1952 Theorem pp180-181, prime x<q<6x/5 for x>=25; Bertrand postulate recalled there'])]
ids=[row['id'] for row in new]

def open_scope(row):
    addition=(' FIRSLOPECAP20261116 now distinguishes COMPLETE finite slope-only coverage from residual tests: only primes<=6n can restrict slopes, and an even candidate k can be excluded only by primes<=12k or divisors of n(n-1)...(n-2k+1). FIRSLOPECOVEREXEC20261116 exhausts this finite cap at n2..64, with no survivors in that bounded range. FIRSHIFTEDINDEX20261116 proves an INFINITE shifted family: n=p^r+i, every prime p>18, r>=1 and offsets0..3, via one horizontal edge and one edge of denominator p^r plus parity. General N=u*p^r+i has exact residual kernels P_(i,6),P_(u,6) under p>6*max(u,i), giving a uniform candidate-degree exclusion and fixed-kernel all-r refinements. These are proved mechanisms, not all-index prime coverage or nearly-closed SFC. The all-index exponent-six coverage theorem and general SFC remain open; inactive-variable FC stays closed.')
    if addition not in row['scope']:row['scope']+=addition
    new_step=(' FIRLARGESQUARE20261116 now gives a UNIFORM prime choice excluding every repeated rational factor degree>=3n/10 for every n>=9, using the pinned Nagura interval. The user full-square obstruction closes a2,i2; the exact F2/F3 discriminant gcd closes a2,i3 and a3,i2, and the large-degree theorem closes a3,i3. All a1..3,i0..3 under p>6max(a,i), and a4..6,i0..1 under p>6a, are now uniform indexed families. The remaining coverage target is small even degrees below3n/10, not the former 2P+3 boundary.')
    if new_step not in row['scope']:row['scope']+=new_step
    row['dependencies']=list(dict.fromkeys(row['dependencies']+ids))
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':open_scope},new)
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+ids))
    row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PROOF,BASE+'slope_coverage_verification.json']))
    return row

def lesson(row):
    row=refresh(row)
    addition=(' Keep slope-only and residual-multiplicity budgets separate. For shape1, full slopes need only primes<=d*n, whereas horizontal residuals can help above that cap. For a candidate k, primes above2*d*k can restrict slopes only when dividing the first2k descending factors from n. Use exact shifted indices N=u*p^r+i with p>d*max(u,i): an integer factorial/binomial quotient transfers the unshifted single-edge inequality; the two residual kernels are P_i,d and P_u,d, independent of r up to a unit scaling. This proves all-r offset0..3 exponent-six windows near p^r for every p>18; doubled occupancy and real parity close them without a finite-index extrapolation.')
    if addition not in row['use']:row['use']+=addition
    addition=(' A survivor after primes<=d*n limits ONLY the slope-and-length method; do not cap residual tests there. An empty finite table is not all-index coverage. The shifted family is infinite but not every starting index, and a useful-prime necessary shortlist does not prove an excluding prime exists for every candidate.')
    if addition not in row['avoid']:row['avoid']+=addition
    addition=(' Ordinary factorial constant1 and a nonsquare factorial leading coefficient exclude a full square. For remaining large degrees, use Nagura at3n and the exact0/1 valuation hull: every square divisor has degree<3n/10 for n>=9. Keep raw local survivors visible when these GLOBAL rules remove them. For kernels F2 and F3 the discriminant gcd2^10*3^4*5^2 excludes simultaneous repeated residual factors for p>18, uniformly across all powers.')
    if addition not in row['use']:row['use']+=addition
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'newton_slope_coverage.py',CHECKER]))
    return row

def work(row):
    row=refresh(row)
    addition=(' Current focused B mechanism: for each even candidate k use only primes<=12k or divisors of n(n-1)...(n-2k+1) for slope exclusions. A shifted N=u*p^r+i with p>6*max(u,i) excludes floor(i/2)<k<p^r, and its fixed residual kernels handle further degrees. All n=p^r+i, p>18, offsets0..3, r>=1 are now closed indexed windows. Seek arithmetic cases covering the remaining (n,k), or a FULL-prime-cap slope survivor before refining its residuals; tested n2..64 alone give no uniform coverage.')
    if addition not in row['next_step']:row['next_step']+=addition
    addition=(' The former 2P+3 and 3P+2,+3 boundaries are closed. Focus on remaining even k<3n/10 using the candidate-prime shortlist and higher residual information. The raw h=t1 pattern at4*137^r+3 survives the current single-prime/global filters; it is UNKNOWN for that test, not a global obstruction or SFC counterexample.')
    if addition not in row['next_step']:row['next_step']+=addition
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
ledger_path=HERE/'source_ledger.json'
ledger=json.loads(ledger_path.read_text())
ledger.setdefault('internal_proof_extensions',{})['slope_coverage_and_shifted_windows']={
    'canonical_source':PROOF,'claims':ids,
    'provenance':'User-proposed finite slope cap and degree-specific coverage target; the exact candidate-prime restriction and shifted-index two-kernel formulas are proved here, retaining the credited Dumas/Filaseta--Rocks unshifted single-edge input. No novelty claim.',
    'unshifted_primary_input':'Filaseta--Rocks2016 Lemma2.3, already rederived in the mixed-prime proof; shifted transfer is independently derived by the integer factorial/binomial quotient.',
    'boundary':'Infinite shifted-index families and finite full-cap diagnostics do not give every-index exponent-six SFC. The cap is for slopes only, not residual polynomials; shape1 is explicit.'}
url='https://www.jstage.jst.go.jp/article/pjab1945/28/4/28_4_177/_pdf'
for name in ('sources/nagura-1952-prime-interval.pdf','sources/nagura-1952-prime-interval.txt'):
    row={'file':name,'url':url,'sha256':sha(HERE/name),
         'review':'Published1952 source. Theorem pp180-181 gives x<q<6x/5 for x>=25; current application uses x=3n>=27. Imported as a classical prime-interval result; proof not independently reconstructed. Text is an extraction of the pinned authoritative PDF.'}
    existing=next((x for x in ledger['cached_primary_sources'] if x['file']==name),None)
    if existing is None:ledger['cached_primary_sources'].append(row)
    else:assert existing==row
ledger['additional_primary_references']=list(dict.fromkeys(ledger['additional_primary_references']+[url]))
ledger['published_input_pinpoints']['NaguraLargeSquare']='1952 Theorem, pp180-181, n5 row and threshold25: x<q<6x/5. Applied at3n for every n>=9 to a0/1 factorial coefficient valuation hull. Imported, not internally reproved.'
ledger['internal_proof_extensions']['slope_coverage_and_shifted_windows']['leading_and_large_degree_extension']='User full-square argument and shifted blocks retained. Independent F2/F3 discriminant gcd removes their common repeated-residual pattern. Nagura gives a uniform exclusion for all square-divisor degrees>=3n/10, completing all a1..3,offset0..3 and a4..6,offset0..1 families. Small-degree all-index coverage remains open.'
ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
print(json.dumps({'registered':ids,'checks':receipt['total'],'infinite_shifted_family':'a1..3,i0..3,p>6max(a,i); a4..6,i0..1,p>6a; every r>=1','all_index_exponent_six':'open'}))
