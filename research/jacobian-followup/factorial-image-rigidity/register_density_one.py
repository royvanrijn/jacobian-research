"""Register ordinary natural density one with scoped imported analytic evidence."""
from pathlib import Path
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='jacobian-followup/factorial-image-rigidity/'
PROOF=BASE+'DENSITY_ONE_SQUAREFREENESS.tex'
CHECKER=BASE+'verify_density_one.py'
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


if any(row['id']=='FIRQUANTDENSITY20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Use register_quantitative_squarefree.py; the quantitative all-fixed-exponent scope must remain current.')

receipt=json.loads((HERE/'density_one_verification.json').read_text())
assert receipt['status']=='pass' and receipt['ordinary_natural_density_one_proved']
assert receipt['fixed_H_then_X_then_H_order_checked'] and not receipt['all_index_SFC_closed']
for file,expected in receipt['files'].items():assert sha(HERE/file)==expected,('stale density-one receipt',file)
claim=dict(id='FIRDENSITYONE20261116',kind='theorem',state='proved',
    title='Exponent-six factorial moment polynomials are squarefree on a set of ordinary natural density one',
    scope='For F_n(T)=sum binom(n,j)*(6j)!*T^j in ordinary factorial shape1, the number of n<=X for which F_n is NOT squarefree over Q is o(X): ORDINARY NATURAL density one, not merely positive lower density or logarithmic density. Filaseta-Rocks Theorem1.1 supplies irreducible, hence squarefree, offsets I with lower density>=log2>2/3, and equation3.3 gives n-q_n<sqrt(n) at density one. Fix a sufficiently large H first; I_H is finite with size>2H/3 and finitely many nonzero discriminants, so all reductions are squarefree past a finite T_H. For Y=(2X)^(2/3), f_X(u)=indicator(Pplus(u)<=Y) is bounded real completely multiplicative. Disjoint large-prime multiples and PNT derive its long mean1-log(3/2)<3/5. Matomaki-Radziwill Theorem1 is explicitly UNIFORM in f,h,delta, so X-dependent f_X is permitted. With fixed delta1/100 and sufficiently large fixed H, at most eta(H)X+o_H(X) backward H-blocks have smooth count>=2H/3, eta(H)->0. Therefore some i in I_H gives a prime p>(2X)^(2/3) dividing n-i. For fixed H and large X, p>T_H,p>6max(a,i),n=ap+i. The global prime-gap square-divisor bound gives k<sqrt(n)/2<p. The shifted polygon forces a primitive integral g with g^2|F_n to be wholly horizontal, preserving degree on reduction; Gauss gives g(0)=+-1. Lucas/factorials give F_n modp=F_i modp, which is squarefree, so k=0. First let X tend to infinity with H fixed, then H tend to infinity; dyadic summation proves o(X). Simultaneously for every nonzero complex A+B*Y^6, a COMMON density-one set of terminal indices has a nonzero moment in the adjacent pair n-1,n. No fixed multiplier/offset cutoff or all-C smoothness hypothesis is used. This DOES NOT prove irreducibility density one, squarefreeness at every index, uniform all-index residual coverage, arbitrary Gamma shapes or full SFC. A possible density-zero exceptional set remains. Deep published inputs are credited imports, not independently reconstructed or proved by finite controls; no novelty or assurance upgrade.',
    canonical_source=PROOF,checker=CHECKER,proof_type='hybrid',
    independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_density_one.py'),
    dependencies=['FIRPRIMEGAP20261116','FIRSHIFTEDINDEX20261116','FIRREALROOTPARITY20261116',
        'external: Filaseta-Rocks2016 Theorem1.1 and Section3 equation3.3, pinned journal Online First author copy',
        'external: Matomaki-Radziwill2016 Annals Theorem1 pp1015-1016, explicit uniformity in h,delta,f',
        'external: classical PNT, retained DLMF27.2.3 reference; elementary partial summation'],
    software_lock=[PROOF,CHECKER,BASE+'density_one_windows.py',BASE+'density_one_verification.json',
        BASE+'density_one_audit.json',BASE+'density_one_local_examples.json',BASE+'register_density_one.py',
        BASE+'prime_power_windows.py',BASE+'prime_gap_windows.py',BASE+'source_ledger.json',
        BASE+'PRIME_GAPS_AND_SMALL_KERNELS.tex',BASE+'SLOPE_COVERAGE_AND_SHIFTED_WINDOWS.tex',
        BASE+'sources/matomaki-radziwill-2016-annals.pdf',BASE+'sources/matomaki-radziwill-2016-annals.txt',
        BASE+'sources/filaseta-rocks-2016-published.pdf',BASE+'sources/filaseta-rocks-2016-published.txt'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],
    consumers=[],invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')

def open_scope(row):
    addition=(' FIRDENSITYONE20261116 now proves ORDINARY NATURAL density-one squarefreeness for exponent6: #{n<=X:F_n not squarefree over Q}=o(X), hence a common density-one set of adjacent moment windows for every nonzero complex A+B*Y^6. The proof combines credited Filaseta-Rocks dense irreducible offsets and density-one small preceding gaps with Matomaki-Radziwill short averages, whose explicit uniformity permits an X-dependent smooth indicator. Its long mean is derived by disjoint large-prime multiples and PNT. Fix H before sending X to infinity to pass finitely many bad discriminant primes, then send H to infinity and sum dyadic blocks. This gives asymptotic coverage without a fixed multiplier or offset cutoff, not squarefreeness at every offset/index. No all-C smoothness theorem is needed. The remaining exponent-six gap is the possible density-zero exceptional set and uniform all-index coverage; full SFC remains OPEN. The theorem is not a numerical extrapolation, irreducibility-density-one claim or independent/formal/external assurance upgrade.')
    if addition not in row['scope']:row['scope']+=addition
    row['dependencies']=list(dict.fromkeys(row['dependencies']+['FIRDENSITYONE20261116']))
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':open_scope},[claim])
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+['FIRDENSITYONE20261116']))
    row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PROOF,BASE+'density_one_audit.json',BASE+'density_one_verification.json']))
    return row

def lesson(row):
    row=refresh(row)
    extra=(' Ordinary natural density-one exponent-six squarefreeness is now proved by analytic amplification: fix a finite offset set of lower density>2/3 from Filaseta-Rocks, apply uniform Matomaki-Radziwill to the X-dependent smooth indicator whose long mean is<3/5, and choose a large prime factor at a good offset. Pass the finite discriminant threshold while H remains fixed. The prime-gap bound forces degree<p, hence degree-preserving horizontal reduction contradicts a squarefree residual. Let X tend to infinity before H, then sum dyadic blocks. This does not require uniformly smooth C_d at all offsets.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Do not replace fixed-H then X then H limits by an uncontrolled growing H(X), treat the X-dependent multiplicative function as covered without explicit uniformity, or infer density from local finite controls. Ordinary natural density one can leave infinitely many density-zero exceptional indices; it does not close all-index SFC or prove irreducibility density one. Finite bad discriminant primes must be passed before enlarging the offset set.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'density_one_windows.py',CHECKER]))
    return row

def work(row):
    row=refresh(row)
    row['next_step']='Current focused B gate: remove the POSSIBLE DENSITY-ZERO exceptional set in exponent-six squarefreeness by a uniform all-index residual-degree mechanism or another argument. Ordinary NATURAL density-one squarefreeness is CLOSED by FIRDENSITYONE20261116 using fixed finite dense squarefree offsets, uniform Matomaki-Radziwill short averages, finite bad-discriminant avoidance and the prime-gap/local horizontal contradiction. No fixed multiplier or offset cutoff remains in that asymptotic theorem; it is not all-index SFC. Preserve the order fixed H, X to infinity, then H to infinity. The all-r bounded family a1..24000,offset0..12,p>6max(a,i) also remains closed with no initial exceptions, and explicit P0(a) remains available for larger fixed multipliers. Do not pursue all-C positivity as a weaker gate: it is equivalent to the characteristic-zero problem here, and all-C small-prime support is stronger. Inactive-variable FC stays completed. General FC/SFC, other Gamma shapes and other programme obligations remain separate open work.'
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
ledger_path=HERE/'source_ledger.json';before=ledger_path.read_bytes();ledger=json.loads(before)
url='https://annals.math.princeton.edu/wp-content/uploads/annals-v183-n3-p06-p.pdf'
ledger['additional_primary_references']=list(dict.fromkeys(ledger['additional_primary_references']+[url,'https://annals.math.princeton.edu/2016/183-3/p06']))
for name in ('sources/matomaki-radziwill-2016-annals.pdf','sources/matomaki-radziwill-2016-annals.txt'):
    record={'file':name,'url':url,'sha256':sha(HERE/name),
            'review':'Published Annals183(2016)1015-1056. Theorem1 pp1015-1016 and immediate uniformity statement checked: real bounded multiplicative f; h,delta,f may vary uniformly; exceptional estimate retains the (log h)^(1/3) factor. Applied to X-dependent smooth indicators with fixed H first. Deep proof imported, not independently reconstructed; PDF authoritative, text a local pypdf extraction.'}
    old=next((r for r in ledger['cached_primary_sources'] if r['file']==name),None)
    if old:old.update(record)
    else:ledger['cached_primary_sources'].append(record)
for row in ledger['cached_primary_sources']:
    if row['file'] in ('sources/filaseta-rocks-2016-published.pdf','sources/filaseta-rocks-2016-published.txt'):
        extra=' Density-one extension additionally checks exact polynomial normalization, Theorem1.1 lower density>=log2 and Section3 definition of B/equation3.3 giving density one for n-q_n<sqrt(n). These published analytic inputs are imported; their deeper prime-gap proof is not independently reconstructed.'
        if extra not in row['review']:row['review']+=extra
ledger['published_input_pinpoints']['MatomakiRadziwillDensityOne']='Theorem1 pp1015-1016, publisher PDF DOI10.4007/annals.2016.183.3.6, immediate uniformity in h,delta,f. Fixed delta1/100, fixed sufficiently large H; eta(H)=C*(log H)^(1/3)/(delta^2*H^(delta/25))->0. X-dependent smooth indicator is allowed. Imported unconditional theorem, not an internal analytic proof.'
ledger['published_input_pinpoints']['FilasetaRocksDensityOne']='Exact sequence is f_(n,6). Theorem1.1 gives lower density of irreducible offsets>=log2>2/3. Section3 B definition and equation3.3 give natural density one for n-q_n<sqrt(n). Credited published inputs, statements checked in pinned author-hosted journal Online First PDF.'
ledger.setdefault('internal_proof_extensions',{})['ordinary_natural_density_one']={
    'canonical_source':PROOF,'claims':['FIRDENSITYONE20261116'],
    'provenance':'User-supplied density-one argument checked against primary analytic inputs and registered local lemmas. The smooth mean is independently derived from disjoint large-prime multiples and PNT; Gauss normalization, degree preservation and limits are explicit. No novelty claim.',
    'evidence':{name:sha(HERE/name) for name in ('DENSITY_ONE_SQUAREFREENESS.tex','density_one_audit.json','density_one_verification.json')},
    'scope':'#{n<=X:F_n not squarefree over Q}=o(X), ordinary natural density one, exponent6 shape1; common density-one adjacent windows for all complex A,B.',
    'boundary':'A density-zero exceptional set may remain. No all-index SFC, irreducibility density one, all-C smoothness or uniform arbitrary-offset family claim. Deep external analytic results are imported, not reproved by finite controls. All assurance flags false.'}
after=(json.dumps(ledger,indent=2)+'\n').encode()
if after!=before:
    assert ledger_path.read_bytes()==before,'Concurrent source edit'
    temp=ledger_path.with_suffix('.density-tmp');temp.write_bytes(after);os.replace(temp,ledger_path)
print(json.dumps({'registered':'FIRDENSITYONE20261116','checks':receipt['total'],
                  'ordinary_natural_density_one':True,'exceptional_count':'o(X)','all_index_SFC':'open'}))
