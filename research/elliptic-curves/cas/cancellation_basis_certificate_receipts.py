"""Read a bank already certified before use, without inventing a reference run.

Exact minimum exhaustion is performed by the independent, metered construction
checker. This reader binds that record and repeats direct witness arithmetic;
it is not a standalone replacement for the minimum-exhaustion check.
"""
from fractions import Fraction as F
from pathlib import Path

from cancellation_basis_epoch import read, sha, need
from finite_cancellation_corpus import canonical, digest, ROOT
from verify_parity_minimum import ParityMinimumVerifier, integer

CAS = Path(__file__).resolve().parent
CHECKERS = ['verify_parity_minimum.py','verify_cancellation_basis_bank.py','verify_cancellation_basis_bank_v2.py']


def check_bank(packet, case, plan, folder):
    from sage.all import EllipticCurve, QQ, ZZ, matrix
    for name in CHECKERS:
        key = 'elliptic-curves/cas/'+name
        need(plan['bank_certificate_sources'][key] == sha(ROOT/key),'construction checker source differs')
    receipt = read(folder/'verification.json')
    need(receipt['status'] == 'PASS_INDEPENDENT_BANK_CERTIFICATE','unverified certificate bank')
    for path,field in [('bank.json','bank_sha256'),('seed.json','seed_sha256'),
                       ('producer/selection.json','producer_sha256'),('selection-proof.json','selection_proof_sha256')]:
        need(sha(folder/path) == receipt[field],'bank certificate payload changed')
    bank = read(folder/'bank.json'); selection = read(folder/'producer/selection.json'); proof = read(folder/'selection-proof.json')
    n = len(packet['points']); parent_hash = digest(canonical(case['parent_bank']))
    need(read(folder/'seed.json') == bank['seed'] == packet and bank['rank'] == receipt['rank'] == selection['rank'] == n,
         'bank does not use complete current subgroup')
    need(selection['basis'] == packet['points'] and bank['centres'] == selection['centres'][:plan['maximum_centres']],
         'exported bank selection differs')
    need(proof['status'] == 'PASS_INDEPENDENT_BANK_SELECTION_CERTIFICATE' and proof['rank'] == n and
         proof['seed_sha256'] == digest(canonical(packet)) and proof['selection_sha256'] == receipt['producer_sha256'] and
         proof['parent_bank_sha256'] == receipt['parent_bank_sha256'] == parent_hash,
         'selection certificate input binding differs')
    need(proof['complete_centres'] == len(selection['centres']) and receipt['exported_centres'] == len(bank['centres']),
         'bank centre count differs')
    generic = case['parent_bank']['dimension']; extension_count = 1 << (n-generic)
    need(selection['generic_rank'] == generic and 0 <= n-generic <= plan['maximum_extension_dimension'] and
         selection['extensions_per_anchor'] == extension_count and
         selection['full_cosets_scored'] == extension_count*len(selection['anchors']),'extension census differs')
    gram = selection['rounded_gram']; norm = ParityMinimumVerifier(gram).norm
    u = matrix(ZZ,selection['LLL']); need(u.nrows() == u.ncols() == n and abs(u.det()) == 1,'integral LLL transport differs')
    curve = tuple(map(F,packet['curve'])); need(curve[:3] == (0,0,0),'short model required')
    E = EllipticCurve(QQ,list(curve)); basis = [E(list(map(QQ,p))) for p in packet['points']]; cache = {}

    def native(word):
        word = tuple(map(integer,word)); need(len(word) == n,'anchor word dimension differs')
        if word not in cache:
            P = sum((c*p for c,p in zip(word,basis)),E(0)); need(not P.is_zero(),'zero anchor')
            x,y = P.xy(); sign = -1 if y < 0 else 1
            cache[word] = ((F(str(x)),F(str(abs(y)))),[sign*c for c in word])
        return cache[word]

    certificates = {(c['anchor'],c['extension']):c for c in proof['cvp_certificates']}
    need(len(certificates) == len(proof['cvp_certificates']),'duplicate minimum certificate')
    seen = set(); vectors = 0; limit = plan['bank_policy']['exact_cvp_node_limit']
    for ai,anchor in enumerate(selection['anchors']):
        need(anchor['extension_count'] == extension_count,'anchor extension census differs')
        for suffix,field in [('full.npz','full_scores_sha256'),('maps.json','maps_sha256')]:
            need(sha(folder/'producer'/f'anchor-{ai:02d}-{suffix}') == anchor[field],'selection auxiliary bytes differ')
        for row in anchor['refined']:
            key = (ai,row['extension']); need(key in certificates and key not in seen,'minimum certificate census differs'); seen.add(key)
            claimed = certificates[key]; original = row['cvp']; radius = integer(claimed['norm'])
            minima = [tuple(map(integer,v)) for v in claimed['minima']]
            need(claimed['status'] == 'PASS_EXACT_PARITY_MINIMUM' and radius == original['norm'] == row['metric_norm'] and
                 sorted(minima) == sorted(tuple(v) for v in original['minima']) and minima and len(set(minima)) == len(minima),
                 'minimum witness packet differs')
            need(0 < integer(claimed['verification_nodes']) <= integer(original['nodes']) <= limit,'minimum traversal bound differs')
            parity = [c % 2 for c in anchor['anchor']['representative'][:generic]]
            parity += [(row['extension'] >> i) & 1 for i in range(n-generic)]
            choices = []
            for vector in minima:
                need(len(vector) == n,'minimum dimension differs')
                word = [sum(vector[i]*int(u[i,j]) for i in range(n)) for j in range(n)]
                need([w % 2 for w in word] == parity and norm(word) == radius,'direct minimum norm or parity differs')
                choices.append(native(word)); vectors += 1
            point,word = min(choices)
            need(row['point'] == list(map(str,point)) and row['representative'] == word and
                 row['multiplicity'] == len(minima)//2,'native minimum representative differs')
    need(seen == set(certificates),'extra minimum certificate')
    for row in selection['centres']:
        point,word = native(row['representative'])
        need(row['point'] == list(map(str,point)) and row['representative'] == word and
             row['metric_norm'] == norm(word),'native centre identity differs')
    return bank,{'status':'PASS_BOUND_BANK_WITNESSES','rank':n,'minimum_packets':len(certificates),
        'minimum_vectors':vectors,'native_words':len(cache),'verification_sha256':sha(folder/'verification.json'),
        'boundary':'Construction already performs and meters independent exact minimum exhaustion and full selection before use. Final replay binds that certificate, checks all direct integer norm/parity and native representative witnesses, and preserves the saved bank. It does not rerun the closed ellipsoids or claim a new minimum proof from hashes.'}
