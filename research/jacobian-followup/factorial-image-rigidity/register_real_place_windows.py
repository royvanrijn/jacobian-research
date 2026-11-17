"""Register the elementary real-place refinement without replaying old campaigns."""
from pathlib import Path
import hashlib
import json
import os

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = 'jacobian-followup/factorial-image-rigidity/'
CHECKER = BASE+'verify_mixed_prime_windows.py'
PROOF = BASE+'MIXED_PRIME_RESIDUAL_WINDOWS.tex'
CLAIM = 'FIRREALROOTPARITY20261116'
sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()


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
            assert known[row['id']]==row, ('registered scope changed',row['id'])
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


if any(row['id']=='FIRSLOPECAP20261116' for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']):
    raise SystemExit('Historical real-place registration: use register_slope_coverage.py to preserve the shifted-family scopes and reviews.')

receipt=json.loads((HERE/'mixed_prime_window_verification.json').read_text())
assert receipt['status']=='pass'
for file,expected in receipt['files'].items():
    assert sha(HERE/file)==expected, ('verification receipt stale',file)
registry=json.loads((ROOT/'MATH_STATUS.json').read_text())
known={row['id']:row for row in registry['entries']}
base=known['FIRMULTIPRIMEWINDOW20261116']
new=dict(base,id=CLAIM,
    title='Real-place positivity reduces two-term SFC to squarefreeness and excludes every odd-degree real square divisor',
    scope='For every integer d>=1 and every positive REAL Gamma shape alpha, P_n(t)=E(1+tY^d)^n satisfies P_(n-1)=P_n-(t/n)P_n_prime and P_n(0)=1. Two-term SFC for every nonzero complex A+B*Y^d at every positive starting index is equivalent to squarefreeness of every P_n. Coefficient-zero cases and the already-squarefree linear P1 are treated separately. For even n the polynomial has no real root; for odd n its derivative is strictly positive on the entire real line and it has exactly one negative simple real root. Thus every real H with H^2 dividing P_n has EVEN degree; P1,P2,P3 are squarefree. For rational positive shapes, intersect the p-adic square-divisor degree sets with {0,2,4,...,2*floor(n/4)}. If only0 survives, the actual moment polynomial is squarefree and the adjacent window holds. At d6,n7, prime3 alone leaves local degrees0,1 and the real place excludes1. At n21/23, primes2/3 leave even degree2, so parity does not certify them; the retained prime5 closes those finite inputs. No all-index exponent-six prime-coverage theorem or novelty claim follows.',
    dependencies=['FIRGAMMAADAPTER20261116'], artifact_hash='sha256:'+sha(HERE/'verify_mixed_prime_windows.py'))
new['software_lock']=list(dict.fromkeys(base['software_lock']+[BASE+'register_real_place_windows.py']))
changes={}
for identifier,row in known.items():
    if row['checker']==CHECKER:
        def refresh(row):
            row['artifact_hash']='sha256:'+sha(HERE/'verify_mixed_prime_windows.py')
            if row['id'] in ('FIRMULTIPRIMEWINDOW20261116','FIRMIXEDEXEC20261116'):
                addition=(' FIRREALROOTPARITY20261116 now adds the proved REAL-place restriction for actual positive-shape Gamma moment polynomials: only EVEN positive square-divisor degrees can survive. The exporter records the raw local intersection separately from the real-place-filtered set. At n7,d6, prime3 alone now certifies the adjacent pair despite repeated residual factors. Terminal21/23 retain degree2 with primes2/3, and no all-index prime coverage is claimed. The generic local profile for arbitrary rational polynomials is unchanged.')
                if addition not in row['scope']: row['scope']+=addition
                row['dependencies']=list(dict.fromkeys(row['dependencies']+[CLAIM]))
            return row
        changes[identifier]=refresh

def open_scope(row):
    extra=(' FIRREALROOTPARITY20261116 records the exact all-index squarefreeness reformulation for A+B*Y^6 and all positive Gamma shapes, and proves no repeated real roots. All repeated rational irreducible factors therefore have even degree. Combined Newton certificates exclude only the even positive degree intersection; n7 is certified at prime3 alone, but n21/23 still permit degree2 at primes2/3. The missing all-index exponent-six theorem is prime coverage of every possible repeated even-degree factor, or a stronger simple-root argument. The inactive-variable FC extension remains CLOSED and is not a missing analytic lemma.')
    if extra not in row['scope']: row['scope']+=extra
    row['dependencies']=list(dict.fromkeys(row['dependencies']+[CLAIM]))
    return row
changes['OP-FACTORIAL-GENERAL20261116']=open_scope
update_array(ROOT/'MATH_STATUS.json','entries',changes,[new])
known={row['id']:row for row in json.loads((ROOT/'MATH_STATUS.json').read_text())['entries']}
finger=lambda identifier: hashlib.sha256('\0'.join(known[identifier][key] for key in ('state','title','scope')).encode()).hexdigest()[:16]

def refresh_review(row):
    row['claims']=list(dict.fromkeys(row['claims']+[CLAIM]))
    row['reviewed_claims']={identifier:finger(identifier) for identifier in row['claims']}
    row['sources']=list(dict.fromkeys(row['sources']+[PROOF]))
    return row

def lesson(row):
    row=refresh_review(row)
    extra=(' For actual positive-shape Gamma moments P_n(t)=E(1+tY^d)^n, real positivity proves no real roots at even n and exactly one negative simple real root at odd n. Hence every real square divisor has even degree. Apply this real-place restriction AFTER computing generic local p-adic degree sets, and retain both raw and filtered degrees; prime3 alone then certifies n7,d6. Prove the all-index prime-coverage statement rather than inferring it from a finite table or the lower-density irreducibility theorem.')
    if extra not in row['use']: row['use']+=extra
    avoid=(' Do not apply the even-degree filter to arbitrary rational polynomials: they may have repeated linear factors. A failed first-order real/p-adic test is UNKNOWN, and degree2 at n21/23 survives the parity restriction with primes2/3. The inactive-variable endpoint limit is already proved with full tails; do not list it as an outstanding FC gap.')
    if avoid not in row['avoid']: row['avoid']+=avoid
    return row

def work(row):
    row=refresh_review(row)
    old=('For full two-term SFC, the d6 simple-root/higher-Newton gap remains; the earlier D+1 ODE bound is not a solution.')
    new=('For A+B*Y^6, the exact remaining focused gate is all-index prime coverage excluding every possible repeated EVEN-degree factor, or a stronger simple-root theorem. Combine the real-place restriction with Newton residual multiplicity budgets; degree2 survives primes2/3 at terminal21/23, already certified separately using prime5. A finite table or the earlier D+1 ODE bound is not an all-index solution.')
    row['next_step']=row['next_step'].replace(old,new)
    return row
update_array(ROOT/'knowledge/lessons.json','lessons',{'METHOD-FACTORIAL-RADIAL-PURITY':lesson})
update_array(ROOT/'knowledge/work_items.json','actions',{'WORK-FACTORIAL-GENERAL':work})
audit_path=HERE/'mixed_prime_window_audit.json'
audit=json.loads(audit_path.read_text())
audit['real_place_refinement']={'claim':CLAIM,'all_positive_shapes_and_integer_exponents':True,
    'all_index_equivalence':'two-term SFC iff every actual moment polynomial is squarefree',
    'even_degree_square_divisors_only':True,'one_prime_example':{'terminal_index':7,'exponent':6,'prime':3},
    'two_prime_even_survivors':[21,23],'all_index_exponent_six_closed':False,
    'verified_controls':receipt['total'],'independent_replay':False,'formal_verification':False,'external_review':False}
audit_path.write_text(json.dumps(audit,indent=2)+'\n')
ledger_path=HERE/'source_ledger.json'
ledger=json.loads(ledger_path.read_text())
url='https://www.impan.pl/shop/publication/transaction/download/product/91656'
ledger['additional_primary_references']=list(dict.fromkeys(ledger['additional_primary_references']+[url]))
ledger.setdefault('internal_proof_extensions',{})['real_place_squarefreeness']={
    'canonical_source':PROOF,'claim':CLAIM,
    'provenance':'Elementary derivative/positive-Gamma proof of the user-supplied reduction and parity restriction; no novelty claim. Filaseta--Rocks final published Theorem1.1 verifies lower density>=log2 irreducibility, not an all-index conclusion.',
    'published_density_source':url,
    'boundary':'Parity strengthens actual moment certificates but does not remove the even degree2 survivors at n21/23 for primes2/3. The generic local degree theorem is unchanged.'}
ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
print(json.dumps({'registered':CLAIM,'verification_checks':receipt['total'],'inactive_variable_extension':'already proved','all_index_exponent_six':'open'}))
