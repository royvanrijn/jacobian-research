#!/usr/bin/env python3
"""Package checked MW16 transports and retain portable bounded-run evidence."""
import argparse
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import certify_compact_r17_candidates as cert

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/'artifacts/local/elliptic-curves'
OUTPUT=ROOT/'artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json'


def build(output):
    if output.exists(): raise FileExistsError('preserve immutable atlas')
    families=[]
    for ident in ('01','03','05','06','07'):
        path=LOCAL/'compact-mw16-sections-v1'/('a1-presentation-'+ident+'.json')
        row=cert.read(path)
        if row['status']!='PASS_EXACT_EQUATIONS_AND_16_SECTION_TRANSPORTS':
            raise ArithmeticError('section reconstruction incomplete')
        if len(row['sections'])!=16: raise ArithmeticError('section count changed')
        for name,h in row['checker_sources'].items():
            if cert.hashed(ROOT/name)!=h: raise ArithmeticError('source changed before packaging')
        families.append(row)
    if len({r['fibration_id'] for r in families})!=5: raise ArithmeticError('duplicate family')
    cert.write(output,{'schema':'elliptic-curves.compact-five-mw16-atlas.v1',
        'status':'PASS_EXACT_EQUATIONS_AND_80_SECTION_TRANSPORTS','families':families,
        'parameter_convention':'t is the new compact coordinate, lambda_old=(a*t+b)/(c*t+d); x_old=u^2*x_new/(c*t+d)^4; y_old=u^3*y_new/(c*t+d)^6.',
        'claim_boundary':'Five existing anonymous A1/MW16 fibrations in compact coordinates. Generic independence is inherited from the retained marked bases, reconstructed exactly over Q(t). No new fibration, increased generic rank, specialized rank, global minimality or new curve is claimed.'})
    evidence=output.with_name('compact_five_mw16_evidence_v1.zip')
    if evidence.exists(): raise FileExistsError('preserve evidence archive')
    files=[]
    for name in ('compact-mw16-base-v1','compact-mw16-base-v2','compact-mw16-sections-v1'):
        files.extend(p for p in (LOCAL/name).rglob('*') if p.is_file())
    files.extend(ROOT/'elliptic-curves/cas'/name for name in ('compact_mw16_base.sage','compact_mw16_base_v2.sage','export_compact_mw16_atlas.sage','package_compact_mw16_atlas.py'))
    with ZipFile(evidence,'x',compression=ZIP_DEFLATED) as archive:
        for path in sorted(set(files)): archive.write(path,str(path.relative_to(ROOT)))
    cert.write(output.with_name('compact_five_mw16_evidence_v1.json'),{
        'archive':str(evidence.relative_to(ROOT)),'archive_sha256':cert.hashed(evidence),
        'files':{str(p.relative_to(ROOT)):cert.hashed(p) for p in sorted(set(files))},
        'scope':'Protocols, failed first factor gate, successful bounded continuation, full section generation and retained logs. Generation sources are retained by content; archive presence does not replace equation/marking replay.'})
    print('PACKAGED FIVE MW16 FAMILIES, 80 SECTIONS',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT)
    build(p.parse_args().output)
