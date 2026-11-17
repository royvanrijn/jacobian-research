"""Register effective cutoffs, the complete bounded-multiplier family and residual degrees."""
from pathlib import Path
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='factorial/image-rigidity/'
PROOF=BASE+'EFFECTIVE_CUTOFFS_AND_RESIDUAL_DEGREES.tex'
CHECKER=BASE+'verify_effective_residual_windows.py'
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


if any(row['id']=='FIRDENSITYONE20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Use register_density_one.py; ordinary natural density one is now proved and must remain in the current work scope.')

receipt=json.loads((HERE/'effective_residual_verification.json').read_text())
assert receipt['status']=='pass' and receipt['effective_fixed_multiplier_cutoff_computed']
assert receipt['complete_bounded_multiplier_family']['multipliers']==[1,24000]
assert not receipt['complete_bounded_multiplier_family']['finite_exceptions_remaining']
assert not receipt['all_index_SFC_closed'] and not receipt['uniform_nearby_prime_coverage_proved']
for file,expected in receipt['files'].items():assert sha(HERE/file)==expected,('stale effective/residual receipt',file)
common=dict(kind='theorem',state='proved',canonical_source=PROOF,checker=CHECKER,
    proof_type='hybrid',independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_effective_residual_windows.py'),
    software_lock=[PROOF,CHECKER,BASE+'effective_residual_windows.py',BASE+'effective_residual_verification.json',
        BASE+'effective_cutoff_examples.json',BASE+'residual_degree_examples.json',BASE+'register_effective_residual_windows.py',
        BASE+'initial_prime_gaps_500000.json',BASE+'multiplier_endpoint_certificate.json',
        BASE+'prime_power_windows.py',BASE+'prime_gap_windows.py',BASE+'mixed_prime_windows.py',BASE+'newton_slope_coverage.py',
        BASE+'PRIME_GAPS_AND_SMALL_KERNELS.tex',BASE+'small_kernel_certificates.json',BASE+'source_ledger.json',
        BASE+'sources/dusart-1002.0442.pdf',BASE+'sources/dusart-1002.0442.txt'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[],invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')
new=[dict(common,id='FIREFFECTIVECUTOFF20261116',
    title='An explicit Dusart cutoff makes every fixed-multiplier offset-through-twelve factorial family effective',
    scope='For ordinary factorial F_n and EVERY integer a>=15, put P0(a)=ceil(max(400000,exp(sqrt(a/50)))/(a-2)). For EVERY prime p>6a, r>=1 and offset i0..12, P=p^r>=P0(a) implies F_(aP+i) is squarefree. Dusart Proposition6.8 at x=(a-2)P+i>=400000 gives x<q<=x+x/(25log^2x)<x+2P=n, since log^2x>=a/50 and i<2P. The preceding-prime gap bounds any square-divisor degree<P, eliminating the shifted sloped contribution; the retained horizontal kernels give degree<=1, and real parity finishes. Multipliers1..14 were already uniform under p>6max(a,i). Exact cutoff values15:30770,100:4082,1000:401 are reproduced; p>6000 automatically covers ALL inputs at a1000. Each fixed multiplier has a finite explicitly bounded exception domain p^r<P0, not a counterexample set. Exported ceiling decisions use exact integer/rational logarithm bounds, not rounded transcendental numerics. Published prime-interval estimate is imported, not independently reconstructed. All-index SFC, all-offset kernels, general Gamma shapes and assurance upgrades are NOT asserted.',
    dependencies=['FIRPRIMEGAP20261116','FIRSMALLKERNEL20261116','FIRSHIFTEDINDEX20261116',
                  'FIRREALROOTPARITY20261116','FIRFOURTEENTH20261116','FIRCOUPLEDPERIOD20261116',
                  'external: Dusart2010 arXiv1002.0442 Proposition6.8, x>=396738, unconditional imported prime interval']),
 dict(common,id='FIRRESIDUALDEGREES20261116',
    title='Full nearby-prime residual square-part degree intersections give a weaker all-index coverage target',
    scope='For THIS ordinary factorial moment sequence, C_d=gcd(Res(F_d,F_d_prime),coefficients Sres1(F_d,F_d_prime))>0 IFF F_d is squarefree over Q: any repeated rational irreducible factor has even degree>=2, forcing resultant and first subresultant to vanish. This is an equivalence, NOT a proof of positivity for every d; all-d prime support<=6d is an even STRONGER modular target. For n=q+d and prime q>6max(1,d), equivalently q>6n/7 with q>6, a rational square divisor lies on the horizontal edge with degree-preserving reduction dividing K_(d,q)=product phi^floor(e_phi/2), where F_d modq=unit*product phi^e_phi. Its degree belongs to the FULL bounded sum of irreducible-factor degrees of K. Intersect these sets across any explicit finite nearby-prime list, global even parity and valid prime-gap/n14 bounds; only zero implies squarefreeness. A surviving degree makes no counterexample claim. The weaker scalar shortcut q not dividing C_d also suffices, offsets0/1 automatic. Uniform coverage for n>=N would prove eventual squarefreeness; ALL-index completion also needs the finite initial range. Synthetic irreducible-cubic square-part controls illustrate parity without being passed off as actual kernels. Fixed p113,offset13 gives a genuine ALL-r F_(113^r+13) certificate, without all-prime C13 smoothness or all-offset claims. Uniform residual-degree intersection coverage remains OPEN; finite prime inputs never imply it.',
    dependencies=['FIRSMALLKERNEL20261116','FIRPRIMEGAP20261116','FIRFOURTEENTH20261116',
                  'FIRPADICDEFORM20261116','FIRREALROOTPARITY20261116','FIRSHIFTEDINDEX20261116'])]
new.append(dict(common,id='FIRBOUNDEDMULTIPLIER20261116',
    title='Every admissible prime power is covered for factorial multipliers through twenty-four thousand and offsets through twelve',
    scope='For ordinary factorial F_n(T)=sum binom(n,j)*(6j)!*T^j, EVERY integer a1..24000, offset i0..12, prime p>6max(a,i) and integer r>=1 gives SQUAREFREE F_(a*p^r+i), with NO excluded initial prime powers. The retained uniform theorem handles a<=14. For a15..24000, P=p^r>6a and x=(a-2)P+i>6a(a-2). H(a)=50*log(6a(a-2))^2-a is strictly concave on real a>=3: H_second<=100*(2-log(6a(a-2)))*(a^-2+(a-2)^-2)<0. Exact positive endpoint bounds H15>1785,H24000>1225/72 imply log^2x>a/50 throughout the interval. Right endpoint uses e<2719/1000, e^(1/12)>13/12 and the retained integer rational comparison (2719/1000)^22*12/13<3455712000. If x>=400000, the credited Dusart interval gives n-2P<q<n. If x<400000, n<461539; a complete Eratosthenes prime list through500000, replayed by a separate Euler least-prime-factor sieve, has41538 primes, maximum consecutive gap114 at492113/492227, and maximum n-q_n=113, including the last-prime tail. Since P>90, both branches give degree(G)<P for any rational G^2|F_n. Shifted occupancy eliminates the sloped contribution; the registered horizontal kernels through12 give degree<=1 and parity finishes. The finite certificate concerns prime gaps, not a scan or squarefreeness assertion for every index<=500000. The bounded multiplier interval is complete, while unbounded multipliers/offsets and all-index SFC remain OPEN. Ordinary shape1 only; no novelty or independent/formal/external assurance upgrade.',
    dependencies=['FIREFFECTIVECUTOFF20261116','FIRGAPKERNELFAMILY20261116','FIRPRIMEGAP20261116',
                  'FIRSHIFTEDINDEX20261116','FIRSMALLKERNEL20261116','FIRREALROOTPARITY20261116']))
ids=[row['id'] for row in new]

def old_family(row):
    row['scope']=row['scope'].replace('This eventual statement uses the imported PNT and supplies NO effective cutoff;',
        'The retained PNT proof was qualitative; FIREFFECTIVECUTOFF20261116 now supplies the explicit Dusart cutoff for every a>=15;')
    row['dependencies']=list(dict.fromkeys(row['dependencies']+['FIREFFECTIVECUTOFF20261116']))
    return row

def old_gap(row):
    row['scope']=row['scope'].replace('No effective eventual cutoff or all-index SFC conclusion is asserted;',
        'This prime-gap lemma alone asserts no eventual cutoff or all-index SFC conclusion; its effective fixed-multiplier application is FIREFFECTIVECUTOFF20261116;')
    return row

def open_scope(row):
    row['scope']=row['scope'].replace('Every fixed multiplier and offset<=12 is also eventually squarefree qualitatively by PNT, with NO numerical cutoff claimed.',
        'The former qualitative fixed-multiplier effectiveness gap is now closed by FIREFFECTIVECUTOFF20261116.')
    addition=(' FIREFFECTIVECUTOFF20261116 proves the explicit unconditional P0(a)=ceil(max(400000,exp(sqrt(a/50)))/(a-2)) cutoff for every a>=15,offset0..12,p>6a, retaining only a finite bounded set of smaller prime powers. All-r a1000 is automatic from p>6000. FIRRESIDUALDEGREES20261116 retains the full square-part divisor-degree sets for nearby primes and intersects them with global even bounds. For this sequence C_d>0 is already equivalent to characteristic-zero squarefreeness; demanding all C_d have only primes<=6d is stronger still and must not be treated as a weaker independent route. The substantive all-index target is uniform coverage by these permissive residual-degree intersections, plus the initial range in any eventual argument; neither coverage nor all-index SFC is proved. Inactive-variable FC remains completed.')
    if addition not in row['scope']:row['scope']+=addition
    bounded=(' FIRBOUNDEDMULTIPLIER20261116 now discharges EVERY initial cutoff exception for a1..24000,offset0..12,p>6max(a,i),r>=1. Exact endpoint inequalities plus strict concavity give the whole logarithmic multiplier range; one retained finite prime-gap certificate through500000, replayed by two sieve algorithms, handles all small x. This is a COMPLETE bounded-multiplier family, not merely eventual or sampled coverage. The explicit P0 cutoff remains available beyond24000; a uniform mechanism for unbounded multipliers or offsets remains the substantive all-index gate.')
    if bounded not in row['scope']:row['scope']+=bounded
    row['dependencies']=list(dict.fromkeys(row['dependencies']+ids))
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':open_scope,
    'FIRGAPKERNELFAMILY20261116':old_family,'FIRPRIMEGAP20261116':old_gap},new)
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+ids))
    row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PROOF,BASE+'effective_residual_verification.json',BASE+'residual_degree_examples.json',
        BASE+'initial_prime_gaps_500000.json',BASE+'multiplier_endpoint_certificate.json']))
    return row

def lesson(row):
    row=refresh(row)
    row['use']=row['use'].replace('PNT gives only a qualitative eventual fixed-multiplier statement here.',
        'The retained PNT proof is qualitative; the Dusart cutoff now makes every fixed multiplier effective through offset12.')
    extra=(' Use the explicit Dusart cutoff P0(a) for a>=15 and preserve exact integer/rational log witnesses when exporting its ceiling; list smaller prime powers only as a finite exception domain of this sufficient test. On nearby-prime horizontal residuals retain full irreducible factor multiplicities and ALL square-part divisor degrees, intersected across explicit primes and global even bounds. This tolerates modular repeated kernels that the scalar subresultant shortcut cannot settle.')
    if extra not in row['use']:row['use']+=extra
    row['avoid']=row['avoid'].replace('or turn a PNT existence cutoff into a computed bound.',
        'or turn a PNT existence cutoff into a computed bound without the explicit prime-interval theorem.')
    extra=(' For this sequence unbounded C_d positivity is equivalent to all-index characteristic-zero squarefreeness; all-d smoothness is stronger. Do not present either as a weaker intermediate objective, infer uniform residual coverage from samples, or omit the initial finite range when an argument proves only eventual coverage. Synthetic residual fixtures are not moment-polynomial examples.')
    if extra not in row['avoid']:row['avoid']+=extra
    extra=(' For a1..24000 and offsets0..12 there are NO initial exceptions under p>6max(a,i): use strict concavity of H(a)=50*log(6a(a-2))^2-a and exact endpoint witnesses for large x, plus the single independently cross-algorithm-checked finite prime-gap table for x<400000. Keep finite gap evidence separate from the infinite squarefreeness deduction; no factorial-polynomial sweep is needed. Beyond this multiplier range retain the explicit cutoff and full residual degree intersections.')
    if extra not in row['use']:row['use']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'effective_residual_windows.py',CHECKER]))
    return row

def work(row):
    row=refresh(row)
    row['next_step']='Current focused B gate: a uniform coverage mechanism for UNBOUNDED multipliers or offsets, using full residual square-part divisor-degree intersections across nearby primes and the proven global even bounds. The complete family a1..24000,offset0..12,p>6max(a,i),EVERY r>=1 is CLOSED with no initial exceptions: exact concavity endpoints plus one two-algorithm prime-gap certificate discharge the whole initial range. Do not reopen or enumerate those finite exceptions. For larger fixed multipliers the explicit Dusart P0(a) cutoff still gives effective eventual coverage and a finite domain below the criterion. Kernels may have exceptional large primes or large derivative gcd; retain full divisor-degree sets. Unbounded C_d positivity is already equivalent to all-index characteristic-zero squarefreeness here, and all-C small-prime support is stronger. A future eventual all-index coverage proof must also handle its initial range; the closed bounded-multiplier family does not imply arbitrary-index coverage. Inactive-variable FC and completed shifted families stay closed. General FC/SFC and other programme obligations remain separate open work.'
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
ledger_path=HERE/'source_ledger.json';before=ledger_path.read_bytes();ledger=json.loads(before)
url='https://arxiv.org/pdf/1002.0442'
ledger['additional_primary_references']=list(dict.fromkeys(ledger['additional_primary_references']+[url]))
for name in ('sources/dusart-1002.0442.pdf','sources/dusart-1002.0442.txt'):
    record={'file':name,'url':url,'sha256':sha(HERE/name),
            'review':'Primary2010 arXiv v1, Proposition6.8 page8: for x>=396738 there is a prime x<q<=x*(1+1/(25log(x)^2)), unconditional. Statement read and application inequalities checked. Published numerical proof imported, not independently reconstructed. PDF authoritative; text is a local extraction.'}
    old=next((r for r in ledger['cached_primary_sources'] if r['file']==name),None)
    if old:old.update(record)
    else:ledger['cached_primary_sources'].append(record)
ledger['published_input_pinpoints']['DusartEffectiveFixedMultiplier']='Proposition6.8, primary arXiv1002.0442 v1 page8, threshold396738. At x=(a-2)P+i>=400000 and log^2x>=a/50, interval width<=2x/a<2P since i<2P. This makes the retained fixed-multiplier theorem effective; no new numerical prime-interval proof claimed.'
ledger.setdefault('internal_proof_extensions',{})['effective_cutoffs_and_residual_degrees']={
    'canonical_source':PROOF,'claims':ids,
    'provenance':'User-supplied explicit Dusart threshold, stronger-target distinction and full residual square-part degree refinement checked and integrated. Exact rational log brackets evaluate ceilings; explicit finite nearby-prime sets preserve factor-degree information. No novelty claim.',
    'external_inputs':'Pinned Dusart2010 Proposition6.8; retained Nagura/prime-gap and integer small-kernel proofs. Nonzero algebraic exponential transcendence, already credited in the coupled proof, ensures exact comparison termination.',
    'boundary':'The COMPLETE all-r family is now closed for a1..24000,offset0..12 under p>6max(a,i), with no initial exceptions. Effectiveness remains closed through offset12 for every larger fixed multiplier. Full all-index SFC and uniform nearby-prime residual-degree coverage remain open. C positivity equivalence is proved, not unbounded positivity; all-C smoothness is stronger. Fixedprime113/offset13 is an individual all-r certificate, not uniform C13 smoothness. Shape1 only; all assurance flags false.'}
if 'prime_gap_and_small_kernels' in ledger['internal_proof_extensions']:
    ledger['internal_proof_extensions']['prime_gap_and_small_kernels']['boundary']='All-index SFC and arbitrary-offset smoothness remain open. The retained uniform kernel prefix is0..12. Its former qualitative eventual statement is made effective by FIREFFECTIVECUTOFF20261116; the old proof and receipt remain evidence of the earlier stage. Shape1 only.'
ledger['internal_proof_extensions']['bounded_multiplier_initial_ranges']={
    'canonical_source':PROOF,'claims':['FIRBOUNDEDMULTIPLIER20261116'],
    'provenance':'User-supplied complete multiplier-through24000 corollary checked and integrated. Strict logarithmic concavity uses two exact rational endpoint comparisons; the small-x case uses a retained complete prime list and independent Eratosthenes/Euler agreement. No priority claim.',
    'evidence':{name:sha(HERE/name) for name in ('initial_prime_gaps_500000.json','multiplier_endpoint_certificate.json','effective_residual_verification.json')},
    'scope':'EVERY a1..24000,i0..12,p>6max(a,i),r>=1; no initial exceptions. Finite table proves prime gaps only; the written proof transfers the gap bound to the infinite squarefree family.',
    'boundary':'Unbounded multipliers, offsets and all-index SFC remain open. Two sieve algorithms are internal controls, not an independent reviewer/formal assurance upgrade.'}
after=(json.dumps(ledger,indent=2)+'\n').encode()
if after!=before:
    assert ledger_path.read_bytes()==before,'Concurrent source edit'
    temp=ledger_path.with_suffix('.effective-research/.scratch');temp.write_bytes(after);os.replace(temp,ledger_path)
print(json.dumps({'registered':ids,'checks':receipt['total'],'effective_cutoff':'exact P0(a), every fixed a>=15',
                  'complete_multiplier_range':[1,24000],'initial_exceptions_in_range':False,
                  'uniform_nearby_prime_coverage':'open','all_index_SFC':'open'}))
