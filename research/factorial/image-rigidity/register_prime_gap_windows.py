"""Register the exposed prime-gap edge and retained residual-kernel criteria."""
from pathlib import Path
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='factorial/image-rigidity/'
PROOF=BASE+'PRIME_GAPS_AND_SMALL_KERNELS.tex'
CHECKER=BASE+'verify_prime_gap_windows.py'
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


if any(row['id']=='FIREFFECTIVECUTOFF20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Use register_effective_residual_windows.py; the effective cutoff supersedes the earlier qualitative scope.')

receipt=json.loads((HERE/'prime_gap_verification.json').read_text())
assert receipt['status']=='pass' and receipt['all_kernel_remainders_one']
for file,expected in receipt['files'].items():assert sha(HERE/file)==expected,('stale prime-gap receipt',file)
common=dict(kind='theorem',state='proved',canonical_source=PROOF,checker=CHECKER,
    proof_type='hybrid',independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_prime_gap_windows.py'),
    software_lock=[PROOF,CHECKER,BASE+'prime_gap_windows.py',BASE+'prime_gap_verification.json',
        BASE+'small_kernel_certificates.json',BASE+'prime_gap_examples.json',BASE+'register_prime_gap_windows.py',
        BASE+'prime_power_windows.py',BASE+'mixed_prime_windows.py',BASE+'newton_slope_coverage.py',
        BASE+'SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex',BASE+'source_ledger.json',
        BASE+'sources/nagura-1952-prime-interval.pdf',BASE+'sources/nagura-1952-prime-interval.txt'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[],invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')
new=[dict(common,id='FIRPRIMEGAP20261116',
    title='An exposed primitive Newton edge bounds factorial square-divisor degrees by half the preceding prime gap',
    scope='For ordinary F_n(T)=sum binom(n,j)*(6j)!*T^j, any prime q>6 with n/2<q<=n and q^2>6n gives an EXPOSED edge of width q,height6,reduced denominator q. Literal coefficient valuations are floor(6j/q)+indicator(n-q<j<q). Adjusting by6j/q leaves unique paired minima at h,h+q, since multiplication by6 is injective modulo q on0..n-q. A doubled factor cannot occupy that edge, so any H in Q_q[T] with H^2|F_n has degree<=floor((n-q)/2). With q_n the largest prime<=n, Nagura ensures the hypotheses for EVERY n>=30, giving deg(H)<=floor((n-q_n)/2) and the sharper EVEN bound2*floor((n-q_n)/4) for rational H. The million-index example q_n999983 leaves only degrees2,4,6,8, without constructing a million-degree polynomial. The credited classical PNT implies n-q_n=o(n). No effective eventual cutoff or all-index SFC conclusion is asserted; shape1 only.',
    dependencies=['FIRPADICDEFORM20261116','FIRREALROOTPARITY20261116','external: Nagura1952 prime-interval theorem; PNT as DLMF27.2.3, classical imported inputs']),
 dict(common,id='FIRFOURTEENTH20261116',
    title='A uniform one-fourteenth factorial square-divisor bound extends multiplier blocks through fourteen',
    scope='For EVERY n>=30 and rational H with H^2|F_n in ordinary factorial shape1, deg(H)<n/14. Choose5n/6<q<n by Nagura. If q>6n/7, the primitive q-edge leaves width n-q<n/7. Otherwise L=floor(q/6) gives exact vertices(0,0),(L,0),(L+q,6),(n,7); the width-q,height6 edge and the final height-one edge are both primitive and unavailable to doubled occupancy. Hence degree<=L/2<q/12<=n/14. Combined with shifted h+p^r*t and parity, this proves squarefreeness for EVERY r>=1 at n=a*p^r+i, a1..14,offset0..3,p>6max(a,i); also a1..27,offset0..1,p>6a, and a28,offset0,p>168. Small indices below30 are already covered by the retained shifted/full-square families. This is a uniform prime-selection argument, not a finite scan. Full all-index SFC and general Gamma shapes remain outside the claim.',
    dependencies=['FIRPRIMEGAP20261116','FIRSHIFTEDINDEX20261116','FIRREALROOTPARITY20261116']),
 dict(common,id='FIRSMALLKERNEL20261116',
    title='Integer first-subresultant certificates exclude residual square-divisor degree above one through offset twelve',
    scope='For any d>=2 define R_d=Res(F_d,F_d_prime), full INTEGER first subresultant Sres1=A_d*T+B_d, and C_d=gcd(R_d,A_d,B_d)>=0. If C_d>0 has only prime divisors<=6d, then EVERY prime q>6d has deg gcd(F_d modq,F_d_prime modq)<=1. Integer determinant Bezout combinations specialize without division; a gcd of degree>=2 would force all of R,A,B to vanish. The exact packet proves this smoothness for EVERY d2..12, with all11 remainders1. C2=2^5*3^2*5 and C3=2^30*3^17*5^6*7^3*11^2*13*17 are reproduced. Every R,A,B,C, integer-gcd Bezout witness and small-prime factorization is retained; independent Sylvester and maximal-minor determinants verify the normalization. Residual discriminants MAY vanish: e.g d2,q13 and d3,q107 still have only a linear gcd. Degree-preserving horizontal reduction and real parity then prove EVERY F_(p^r+i) squarefree for i0..12,p>6max(1,i),r>=1; near-prime n=q+d is the r1 case. No arbitrary-offset or shape-generalization claim.',
    dependencies=['FIRPADICDEFORM20261116','FIRSHIFTEDINDEX20261116','FIRREALROOTPARITY20261116']),
 dict(common,id='FIRGAPKERNELFAMILY20261116',
    title='Prime-gap and retained small-kernel bounds combine into a fourteen-by-thirteen uniform shifted family',
    scope='Combining deg(H)<n/14 with the exact horizontal-kernel degree bound<=1 proves F_(a*p^r+i) SQUAREFREE for EVERY a1..14, i0..12, prime p>6max(a,i), and r>=1, hence every nonzero complex A+B*Y^6 has a nonzero adjacent moment at those terminal indices. Horizontal degree h<=1 and parity force any positive surviving h+P*t to be >=P+1; n/14<=P+12/14<P+1 excludes it. The few possible n<30 cases are already covered by a<=2,i<=3 retained families. The wider offsets0/1 multiplier blocks a<=27,plus a28i0,remain valid. For EVERY FIXED a>=1 and fixed i0..12 the same prime-conditioned family is also squarefree for all sufficiently large P=p^r: preceding prime gap o(n)=o(P) excludes any positive sloped contribution, then kernel degree<=1 and parity finish. This eventual statement uses the imported PNT and supplies NO effective cutoff; an unresolved all-r test is not promoted to all-r success. Exporters retain these separate statuses, e.g a28i1 and a15i12. All-index SFC and all-offset residual behavior remain open; ordinary shape1 and assurance limits retained.',
    dependencies=['FIRPRIMEGAP20261116','FIRFOURTEENTH20261116','FIRSMALLKERNEL20261116'])]
ids=[row['id'] for row in new]

def open_scope(row):
    addition=(' FIRPRIMEGAP20261116 now confines EVERY square-divisor degree, for n>=30, to half the preceding-prime gap and its even rounding; PNT gives an o(n) bound. FIRFOURTEENTH20261116 strengthens the explicit uniform range to degree<n/14 using both exact Nagura polygon branches. FIRSMALLKERNEL20261116 supplies replayable integer resultant/first-subresultant certificates for offsets2..12, universal over primes>6d, allowing repeated residual discriminants but only gcd degree<=1. FIRGAPKERNELFAMILY20261116 combines these into ALL-r uniform families a1..14,offset0..12,p>6max(a,i), plus the wider offset0/1 multiplier blocks. Every fixed multiplier and offset<=12 is also eventually squarefree qualitatively by PNT, with NO numerical cutoff claimed. These close the former wider-multiplier test examples, while all-index SFC and uniform residual behavior beyond the retained prefix remain open. Inactive-variable FC remains completed.')
    if addition not in row['scope']:row['scope']+=addition
    row['dependencies']=list(dict.fromkeys(row['dependencies']+ids))
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':open_scope},new)
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+ids))
    row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PROOF,BASE+'small_kernel_certificates.json',BASE+'prime_gap_verification.json']))
    return row

def lesson(row):
    row=refresh(row)
    extra=(' For ordinary factorial square divisors use the exposed width-q,height6 edge for n/2<q<=n,q^2>6n; adjusted heights prove uniqueness of the support pair. The preceding-prime gap and even rounding replace coarse linear degree bounds; both Nagura branches give the effective degree<n/14 bound. On a horizontal residual use the FULL integer resultant and first subresultant, not just the discriminant: smooth gcd(R,A,B) proves modular gcd degree<=1 for EVERY prime>6d. Retain integer Bezout and determinant witnesses; offsets2..12 are certified uniformly in the prime. This combines with the shifted polygon into a1..14,offset0..12 all-r families. PNT gives only a qualitative eventual fixed-multiplier statement here.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not discard integer subresultant content by monic or primitive rescaling, infer arbitrary-offset smoothness from the prefix through12, or turn a PNT existence cutoff into a computed bound. A vanished residual discriminant need not defeat the linear-gcd certificate. Former small one-prime survivors may be closed by the stronger prime-gap stage and must not be reported as current mathematical obstructions.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'prime_gap_windows.py',CHECKER]))
    return row

def work(row):
    row=refresh(row)
    row['next_step']='Current focused B gate: all-index exponent-six squarefreeness after the completed prime-gap and small-kernel theorems. For n>=30, remaining rational repeated factors have even degree<=2*floor((n-q_n)/4), with an explicit universal degree<n/14 bound. Use degree-preserving horizontal residuals and the integer resultant/first-subresultant criterion to understand kernels beyond the retained offsets0..12, or uniformly exclude the remaining small even candidates. All-r a1..14,offset0..12 and the wider offset0/1 blocks are closed. Every fixed multiplier with offset<=12 is eventually closed qualitatively by PNT, with no effective cutoff supplied. Do not reopen inactive-variable FC or the former2P+3 boundary. General FC/SFC and the remaining source/geometry programme obligations remain separate open work.'
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
ledger_path=HERE/'source_ledger.json'
ledger=json.loads(ledger_path.read_text())
url='https://dlmf.nist.gov/27.2.E3'
ledger['additional_primary_references']=list(dict.fromkeys(ledger['additional_primary_references']+[url]))
ledger['published_input_pinpoints']['PNTPrimeGap']='Classical PNT pi(x)~x/log x, DLMF27.2.3, checked online. It implies pi(n)-pi((1-epsilon)n)>0 eventually for each fixed epsilon, hence preceding-prime gap o(n). Credited asymptotic import; no numerical cutoff claimed.'
ledger.setdefault('internal_proof_extensions',{})['prime_gap_and_small_kernels']={
    'canonical_source':PROOF,'claims':ids,
    'provenance':'User-supplied adjusted-height prime-gap argument, two-branch n/14 bound and integer first-subresultant criterion independently checked. The exact eleven certificates reproduce C2/C3 and retain determinant and gcd witnesses. Combining them gives a1..14,offset0..12. No novelty claim.',
    'external_inputs':'Pinned Nagura1952 threshold25; classical PNT from the official DLMF27.2.3 reference. PNT is an imported asymptotic theorem, not an implemented effective bound.',
    'boundary':'All-index SFC and arbitrary-offset smoothness remain open. The retained prefix is0..12; every fixed multiplier there has eventual squarefreeness with no computed cutoff. Shape1 only.'}
ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
print(json.dumps({'registered':ids,'checks':receipt['total'],'uniform_rectangle':'a1..14,offset0..12','effective_eventual_cutoff':None,'all_index_SFC':'open'}))
