"""Register simultaneous published-log and dense-preprint-power routes separately."""
from pathlib import Path
import hashlib,json,os
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE='jacobian-followup/factorial-image-rigidity/'
PROOF=BASE+'SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex'
CHECKER=BASE+'verify_simultaneous_squarefree.py'
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


if any(row['id']=='FIRALLEXPPOWER20261117' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Use register_all_exponent_squarefree.py; the exponent cap has been removed by the elementary large-exponent theorem.')

receipt=json.loads((HERE/'simultaneous_squarefree_verification.json').read_text())
assert receipt['status']=='pass' and receipt['published_2016_fallback_proved']
assert receipt['preprint_dense_power_application_checked'] and receipt['MR_II_preprint_source_qualification_retained']
assert not receipt['complete_preprint_validation_claimed'] and not receipt['sparse_Section10_used']
assert not receipt['analytic_exception_multiplied_by_M'] and not receipt['all_index_SFC_closed']
for file,expected in receipt['files'].items():assert sha(HERE/file)==expected,('stale simultaneous receipt',file)
common=dict(kind='theorem',state='proved',canonical_source=PROOF,checker=CHECKER,proof_type='hybrid',
    independent_replay=False,formal_verification=False,external_review=False,
    artifact_hash='sha256:'+sha(HERE/'verify_simultaneous_squarefree.py'),
    software_lock=[PROOF,CHECKER,BASE+'simultaneous_squarefree.py',BASE+'simultaneous_squarefree_verification.json',
        BASE+'simultaneous_squarefree_audit.json',BASE+'simultaneous_squarefree_parameters.json',
        BASE+'simultaneous_squarefree_examples.json',BASE+'register_simultaneous_squarefree.py',
        BASE+'quantitative_squarefree.py',BASE+'prime_power_windows.py',BASE+'source_ledger.json',
        BASE+'sources/mr-ii-2007.04290v1.pdf',BASE+'sources/mr-ii-2007.04290v1.txt',
        BASE+'sources/matomaki-2022-almost-primes-accepted-v3.pdf',BASE+'sources/matomaki-2022-almost-primes-accepted-v3.txt',
        BASE+'sources/matomaki-publications-20261117.html',BASE+'sources/filaseta-rocks-2016-published.pdf',
        BASE+'sources/filaseta-rocks-2016-published.txt',BASE+'sources/matomaki-radziwill-2016-annals.pdf',
        BASE+'sources/matomaki-radziwill-2016-annals.txt'],
    supersedes=[],closes_problems=[],narrows_problems=['OP-FACTORIAL-GENERAL20261116'],consumers=[],
    invalidates_assumptions=[],replaced_by=[],priority='reference',programme_status='archived')
new=[dict(common,id='FIRCOMMONOFFSETS20261117',
    title='One dense offset set is irreducible for every exponent up to the growing dyadic maximum',
    scope='For integer H,M>=1, every offset i<=H with i>(M+1)^2, i-q_i<sqrt(i) and Pplus(i)>sqrt(MH) has F_(i,m) irreducible over Q for EVERY1<=m<=M. At p=Pplus(i), p^2>MH>=mi and p|i make every factor degree divisible by p. At q_i the general shifted polygon gives a factor of degree>=q_i, with complement<=i-q_i<sqrt(i)<p; divisibility forces complement zero. For M=floor((2X)^(1/16)),H=floor(X^(1/5)), these common offsets J have size(log(32/21)+o(1))*H by disjoint large-prime multiples/PNT, deleting O(M^2)+O(H^(2/3)). This does NOT apply a fixed-m irreducibility-density theorem with growing m. The common smooth indicator cutoffY=sqrt(2MX) has long mean1-log(32/17)+o(1), leaving gap Delta=log(1024/357)-1>1/20, exactly certified with rational exponential bounds. The Hadamard resultant estimate has the ABSOLUTE form log abs R_(i,m)<=5*m*i^2*log(2mi), hence total<=5*M^2*H^3*log(2MH). If a common analytic argument supplies at leastH/50 usable offsets at an index with a good preceding gap, choosing ONE failing exponent and counting all its bad-prime incidences costs O(X^(159/160)), after division by H. Analytic exceptions are never multiplied by M. No parity, all-offset irreducibility, empirical-density or full all-index SFC claim.',
    dependencies=['FIRALLMLOCAL20261116',
        'external: classical PNT and elementary partial summation; standard Hadamard determinant inequality',
        'external: Filaseta-Rocks2016 equation3.1, absolute O(H^(2/3)) bad preceding-gap count']),
 dict(common,id='FIRSIMULTLOG20261117',
    title='Published 2016 short averages give one logarithmically sparse exceptional set for all exponents up to n to the one-sixteenth',
    scope='Define D(X)=#{1<=n<=X:there exists integer1<=m<=n^(1/16) with F_(n,m) not squarefree over Q}, ordinary Gamma shape1. Then D(X)<<X/(log X)^(1/50) with an ABSOLUTE implied constant. This simultaneous theorem uses the common offset set FIRCOMMONOFFSETS20261117, H=floor(X^(1/5)),M=floor((2X)^(1/16)),Y=sqrt(2MX), and the PUBLISHED MR2016 Theorem1 at delta1/100. The mean margin yields H/50 usable offsets outside one common analytic set, whose first polynomial exception saving is1/12500, plus X/(log X)^(1/50). Bad incidences costX^(159/160); common gap exceptions costX^(2/3), endpointsH. Dyadic summation gives the bound. The mathematical dependency route is INDEPENDENT of MR II, its Section10 sparse sieve, Theorem1.7 and the simultaneous power claim. Numerical absolute constant/starting threshold are unextracted. This logarithmic bound alone does NOT imply a convergent reciprocal exceptional sum. Full all-index SFC, arbitrary Gamma shapes and larger growing-exponent ranges remain unclaimed.',
    dependencies=['FIRCOMMONOFFSETS20261117',
        'external: published Matomaki-Radziwill2016 Annals Theorem1, uniform in h,delta,f']),
 dict(common,id='FIRSIMULTPOWER20261117',
    title='A dense-preprint analytic input gives a simultaneous power saving and harmonic sparsity for growing factorial exponents',
    scope='SOURCE-QUALIFIED deduction importing MR II arXiv2007.04290v1 Theorem9.2(ii), DENSE branch: D(X)=#{n<=X:exists integer1<=m<=n^(1/16) with F_(n,m) not squarefree over Q}<<X^(1-10^-18), with an ABSOLUTE implied constant and ONE common exceptional set. Ordinary Gamma shape1. Explicit fixed source parameters nu1=10^-5,nu2=3/20,eta=1/150,rho=1/10,delta=1/2000,J_MR=1 and prime intervals P1=X^(10^-8),Q1=P2=X^nu1,Q2=P3=X^sqrt(nu1*nu2),Q3=X^nu2 satisfy Section9 conditions eventually. Nonnegative f permits twist0. Theorem9.2(ii) applied to the smooth indicator and1 has error1/2000+10^5/(log X)^(1/120) and common exceptions O(X^(1-33/10^10)*(log X)^2+X^(1-5*10^-18)). Separate ordinary linear sieve Lemma4.3/equation78 gives complement loss<13/625+o(1); positivity counts it ONCE, total109/5000+o(1). The exact margin leaves at leastH/50 usable offsets. FIRCOMMONOFFSETS20261117 then gives bad incidences O(X^(159/160)), not a union of M analytic errors; common gap exceptions O(X^(2/3)) and endpoints O(H) are smaller. Every term is O(X^(1-10^-18)), and dyadic summation proves the bound. Partial summation yields sum_(n in the common exceptional set)1/n<infinity. Each fixed-m nonsquarefree set also has convergent reciprocal sum, since n>=m^16 enters the simultaneous range. Crucial SOURCE BOUNDARY: MR II remains listed as a PREPRINT; Matomaki2022 published warning flags an incorrect sparse-sieve calculation in its first arXiv version. The selected dense proof uses Proposition8.3/earlier inputs and the separate ordinary linear sieve, AVOIDING Section10/Lemma10.3. This is an application/dependency audit, NOT a reconstruction or validation of the entire preprint or a claim of obtaining a corrected/published MR II. Savings are derived directly from9.2, not Theorem1.7 advertised exponent. Numerical constants/starting thresholds, finiteness of exceptions, full all-index SFC, arbitrary Gamma shapes, larger exponent ranges and novelty/independent/formal/external assurance are unclaimed.',
    dependencies=['FIRCOMMONOFFSETS20261117',
        'external: MR II arXiv2007.04290v1 Theorem9.2(ii), Section9 pp42-51, dense PREPRINT input with scoped application/dependency audit',
        'external: ordinary linear upper sieve Lemma4.3 and equation78 in MR II; separate from Section10 sparse sieve'])]
ids=[row['id'] for row in new]

def open_scope(row):
    extra=(' FIRCOMMONOFFSETS20261117 now constructs one dense offset set irreducible for ALL m<=M=floor((2X)^(1/16)), avoiding misuse of fixed-m density with a growing exponent. Many usable offsets and incidence division by H give X^(159/160), with no factor M on analytic errors. FIRSIMULTLOG20261117 independently uses PUBLISHED MR2016 to prove simultaneous D(X)<<X/(log X)^(1/50), absolute constant. FIRSIMULTPOWER20261117 additionally imports the DENSE MR II PREPRINT Theorem9.2(ii), with directly audited explicit parameters, to obtain D(X)<<X^(1-10^-18) and a convergent reciprocal sum for the common exceptional set and each fixed-m exceptional set. The published sparse-sieve warning is retained; the selected route avoids Section10/Lemma10.3, and neither a full-preprint reconstruction nor Theorem1.7 advertised exponent is used. Numerical constants/thresholds and an empty or finite exceptional set remain unclaimed. The current all-index SFC gap is the possible infinite harmonically summable exceptional set, not asymptotic simultaneous coverage.')
    if extra not in row['scope']:row['scope']+=extra
    row['dependencies']=list(dict.fromkeys(row['dependencies']+ids));return row

def old_quant_scope(row):
    extra=(' The earlier no-growing-m boundary belongs to THIS fixed-m proof. FIRSIMULTLOG20261117 now separately proves a PUBLISHED-input simultaneous bound for every m<=n^(1/16); FIRSIMULTPOWER20261117 strengthens it using the explicitly qualified dense MR II preprint input. The common-offset construction and incidence division justify the new growing range; it is not obtained by merely making m vary in the present theorem.')
    if extra not in row['scope']:row['scope']+=extra
    return row
update_array(ROOT/'MATH_STATUS.json','entries',{'OP-FACTORIAL-GENERAL20261116':open_scope,
    'FIRQUANTDENSITY20261116':old_quant_scope},new)
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda id:hashlib.sha256('\0'.join(known[id][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh(row):
    row['claims']=list(dict.fromkeys(row['claims']+ids));row['reviewed_claims']={id:finger(id) for id in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PROOF,BASE+'simultaneous_squarefree_audit.json',BASE+'simultaneous_squarefree_verification.json']))
    return row

def lesson(row):
    row=refresh(row)
    extra=(' For growing m use a COMMON offset construction, not fixed-m irreducible density: i>(M+1)^2, a good preceding gap and Pplus(i)>sqrt(MH) force irreducibility for every m<=M. With M=(2X)^(1/16),H=X^(1/5), the common density/mean margin is log(1024/357)-1>1/20. Preserve at leastH/50 usable offsets and divide bad-incidence counts by H to get X^(159/160), rather than union-bounding analytic errors over m. Published MR2016 alone gives a simultaneous logarithmic bound. The separate power/harmonic result explicitly imports dense MR II9.2(ii) and ordinary linear sieve4.3/78; source parameters and tiny savings are exact rationals.')
    if extra not in row['use']:row['use']+=extra
    extra=(' Keep the MR II PREPRINT status and Matomaki2022 sparse-sieve warning explicit. Exclude Section10/Lemma10.3, do not use Theorem1.7 advertised saving, and do not turn the scoped application/dependency audit into full-preprint validation or published-input assurance. Positivity incurs the sieve complement only once. The published logarithmic fallback alone does not imply harmonic convergence; the power route does. A convergent reciprocal sum does not imply finitely many exceptional indices or all-index SFC. The older fixed-m proof still supplies no growing-m uniformity by itself; the new1/16 range needs the common-offset mechanism.')
    if extra not in row['avoid']:row['avoid']+=extra
    row['implementation']=list(dict.fromkeys(row['implementation']+[BASE+'simultaneous_squarefree.py',CHECKER]));return row

def work(row):
    row=refresh(row)
    row['next_step']='Current focused B gate: remove the possible infinite harmonically summable exceptional set to prove full ALL-index constant-plus-monomial SFC, especially non-prime-power exponents such as6. Simultaneous growing-exponent coverage m<=n^(1/16) is CLOSED: published MR2016 independently gives D(X)<<X/(log X)^(1/50), while the explicitly source-qualified dense MR II preprint9.2(ii) plus separate ordinary linear sieve gives D(X)<<X^(1-10^-18) and convergent reciprocal exceptional sums. Preserve the preprint/sparse-warning boundary: selected proof avoids Section10/Lemma10.3; no full-preprint validation or Theorem1.7 advertised saving. The mechanism is one common irreducible offset set, many usable offsets and incidence division by H, with no M factor on analytic errors. Numerical constants/thresholds, finiteness of exceptions and larger exponent ranges remain unclaimed. Do not reopen simultaneous asymptotic coverage, completed bounded-multiplier families or inactive-variable FC. General FC/SFC, other Gamma shapes and other programme obligations remain separate open work.'
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
ledger_path=HERE/'source_ledger.json';before=ledger_path.read_bytes();ledger=json.loads(before)
records=[('sources/mr-ii-2007.04290v1.pdf','https://arxiv.org/pdf/2007.04290v1',
 'MR II PREPRINT v1,8July2020,89pages. Section9 conditions pp42-43, dense Theorem9.2(ii) p44, proof pp45-51, ordinary linear sieve4.3 p21 and equation78 p62 checked for this application. Known sparse-sieve warning retained; Section10/Lemma10.3 and Theorem1.7 advertised exponent NOT used. No full-preprint reconstruction or corrected-version claim.'),
 ('sources/mr-ii-2007.04290v1.txt','https://arxiv.org/pdf/2007.04290v1',
 'Local pypdf layout extraction of pinned v1. Rotated marginal/arXiv text may be omitted; PDF authoritative. Used with direct PDF plain-text extraction for scoped dense-proof reference audit.'),
 ('sources/matomaki-2022-almost-primes-accepted-v3.pdf','https://arxiv.org/pdf/2012.11565v3',
 'Accepted AUTHOR version v3,18January2022, of the published JLMS106(2022)1061-1097 article, DOI10.1112/jlms.12592. Section6 footnote before Lemma6.1 and Remark6.2 retain the sparse-sieve warning; matched to publisher HTML online. This is not described as publisher typesetting.'),
 ('sources/matomaki-2022-almost-primes-accepted-v3.txt','https://arxiv.org/pdf/2012.11565v3',
 'Local layout extraction of the pinned accepted author version; PDF authoritative, rotated margin text may be omitted. Warning and missing-coprimality/off-diagonal explanation checked.'),
 ('sources/matomaki-publications-20261117.html','https://users.utu.fi/ksmato/publications/',
 'Dated17November2026 author-publication-page snapshot: Multiplicative functions in short intervals II is listed as Pre-print. Other similarly named higher-uniformity papers are distinct. This status check does not claim a journal MR II version was obtained.')]
for name,url,review in records:
    record={'file':name,'url':url,'sha256':sha(HERE/name),'review':review}
    old=next((row for row in ledger['cached_primary_sources'] if row['file']==name),None)
    if old:old.update(record)
    else:ledger['cached_primary_sources'].append(record)
ledger['additional_primary_references']=list(dict.fromkeys(ledger['additional_primary_references']+
    [url for _,url,_ in records]+['https://londmathsoc.onlinelibrary.wiley.com/doi/full/10.1112/jlms.12592']))
ledger['published_input_pinpoints']['MRII_Dense_Simultaneous']='SOURCE QUALIFICATION: MR II is a preprint. Dense Theorem9.2(ii) p44 and Section9 conditions pp42-43 with nu1=10^-5,nu2=3/20,eta1/150,rho1/10,delta1/2000,J1 give direct savings33/10^10 and5*10^-18. Ordinary linear upper sieve Lemma4.3/equation78 is separate. Theorem1.7 advertised exponent and sparse Section10/Lemma10.3 unused; scoped application/dependency audit only.'
ledger['published_input_pinpoints']['Matomaki2022SparseWarning']='Published JLMS106(2022)1061-1097, DOI10.1112/jlms.12592, Section6 footnote before Lemma6.1 and Remark6.2: first MR II arXiv sparse-sieve calculation used an incorrect version; missing coprimality creates off-diagonal terms. Accepted author v3 pinned and publisher HTML matched online. Selected dense route avoids the affected sparse branch; no full MR II validation claimed.'
ledger.setdefault('internal_proof_extensions',{})['simultaneous_growing_exponents']={
    'canonical_source':PROOF,'claims':ids,
    'provenance':'User-supplied common-offset construction, direct dense parameters and many-incidence argument checked and integrated. Source distinction preserved; no priority claim.',
    'evidence':{name:sha(HERE/name) for name in ('SIMULTANEOUS_EXPONENTS_AND_POWER_SAVING.tex','simultaneous_squarefree_audit.json','simultaneous_squarefree_verification.json')},
    'published_route':'MR2016 gives absolute D(X)<<X/(log X)^(1/50) for one common exceptional set, all m<=n^(1/16), independently of MR II.',
    'additional_preprint_route':'Dense MR II9.2(ii) and ordinary linear sieve, avoiding Section10/Lemma10.3, give D(X)<<X^(1-10^-18) and convergent reciprocal sums. Application/dependency audit, not full-preprint reconstruction.',
    'boundary':'Absolute numerical constants and starting thresholds unextracted. Infinite exceptions and full all-index SFC remain possible/open; no larger exponent range or other Gamma shapes claimed. All assurance flags false.'}
after=(json.dumps(ledger,indent=2)+'\n').encode()
if after!=before:
    assert ledger_path.read_bytes()==before,'Concurrent source edit'
    temp=ledger_path.with_suffix('.simultaneous-tmp');temp.write_bytes(after);os.replace(temp,ledger_path)
print(json.dumps({'registered':ids,'checks':receipt['total'],'simultaneous_range':'m<=n^(1/16)',
    'published_log_fallback':True,'power_input':'dense MR II preprint; qualified',
    'reciprocal_sum_converges_using_power_route':True,'all_index_SFC':'open'}))
