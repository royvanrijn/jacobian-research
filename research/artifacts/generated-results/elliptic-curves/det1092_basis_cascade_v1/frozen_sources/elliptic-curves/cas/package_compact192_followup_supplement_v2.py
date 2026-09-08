#!/usr/bin/env python3
"""Standalone fixed five26 follow-up and universal small-prime proof supplement."""
from hashlib import sha256
from pathlib import Path
import zipfile
import certify_compact_r17_candidates as cert
import compact192_specialized_followup as follow
import report_compact192_specialized_followup as report
import audit_r17_scaling_prime_support as support
import classify_r17_other_small_prime_scalings as local
from package_recorded_mod2_audit import dependencies
ROOT=follow.ROOT;ART=follow.ART;CAS=follow.CAS
OUT=ART/'compact192_followup_supplement_evidence_v2.json'

def main():
    archive=OUT.with_suffix('.zip')
    if OUT.exists() or archive.exists():raise FileExistsError('preserve standalone follow-up supplement')
    p=follow.protocol();summary=cert.read(report.OUT)
    if summary['status']!='PASS' or len(p['rows'])!=5:raise ArithmeticError('all five terminal follow-ups required')
    paths={Path(__file__).resolve(),CAS/'verify_compact192_followup_supplement_v2.py',CAS/'replay_retention24_geometry.py',
           CAS/'report_r17_small_prime_minimality.py',CAS/'verify_r17_other_small_prime_scalings.py',
           CAS/'verify_r17_scaling_prime_support.sage',report.OUT,support.OUT,local.OUT,
           ART/'r17_small_prime_minimality_v1.json'}
    for directory in [follow.BATCH,report.CONTROL,support.D,local.D]:
        paths.update(q for q in directory.rglob('*') if q.is_file() and '__pycache__' not in q.parts)
    bindings=[p,cert.read(report.CONTROL/'protocol.json'),cert.read(support.D/'protocol.json'),
              cert.read(local.D/'protocol.json'),summary,cert.read(ART/'r17_small_prime_minimality_v1.json')]
    for data in bindings:
        for key in ('sources','inputs','seed_hashes'):
            for name,h in data.get(key,{}).items():
                q=ROOT/name
                if cert.hashed(q)!=h:raise ArithmeticError('frozen supplement member changed: '+name)
                paths.add(q)
    for row in cert.read(follow.BATCH/'verification-ledger.json')['rows']:
        for key in ('mod2_certificate','modl_certificate'):
            proofpath=ROOT/row[key];paths.add(proofpath)
            for n,h in cert.read(proofpath)['sources'].items():
                q=ROOT/n
                if cert.hashed(q)!=h:raise ArithmeticError('cloud proof source changed')
                paths.add(q)
    for name in ['COMPACT192_SPECIALIZED_FOLLOWUP_2026-09-06.md','R17_SCALING_PRIME_SUPPORT_2026-09-06.md','R17_SMALL_PRIME_MINIMALITY_2026-09-06.md']:
        paths.add(ROOT/'elliptic-curves/notes'/name)
    paths.update([CAS/'audit_recorded_point_mod2_rank_v3.py',CAS/'audit_retained_cloud_modl.py',CAS/'package_compact192_followup_supplement.py',CAS/'verify_compact192_followup_supplement.py',ART/'compact192_followup_supplement_evidence_v1.json',ROOT/'artifacts/local/elliptic-curves/compact192-followup-supplement-portable-v1/controller.supervisor.json'])
    paths=dependencies(paths);members=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED) as z:
        for q in sorted(paths):
            raw=q.read_bytes();name=str(q.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,6,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
            z.writestr(info,raw);members.append({'path':name,'sha256':sha256(raw).hexdigest(),'bytes':len(raw)})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None or any(sha256(z.read(r['path'])).hexdigest()!=r['sha256'] for r in members):raise ArithmeticError('supplement integrity differs')
    cert.write(OUT,{'schema':'elliptic-curves.compact192-followup-supplement-evidence.v2','archive':str(archive.relative_to(ROOT)),
        'archive_sha256':cert.hashed(archive),'archive_bytes':archive.stat().st_size,'files':members,'required_base_archives':[],
        'supersedes_failed_packaging':'compact192_followup_supplement_evidence_v1 omitted a checker invoked by command string; its failed isolated replay is preserved. This version adds the missing cloud checker and every explicit cloud-certificate source. No search input or point proof changes.',
        'claim_boundary':'Standalone replay of all five own26 follow-up point histories, full clouds modulo2,3,5 and exact rational geometry, plus the completed245-box result. Independently replays the six exact resultants, all26 non13 scaling-prime residue exclusions and universal990-pair small-prime minimality theorem. Prior inventory selection and known-control discovery are retained as bound gate context; their full historical campaigns are not rerun. No new point search, full global prime support, rank upper bound, saturation, point absence or universal novelty.'})
    print('STANDALONE FOLLOWUP SUPPLEMENT',len(members),'files',archive.stat().st_size,'bytes',flush=True)
if __name__=='__main__':main()
