#!/usr/bin/env python3
"""Source-guarded replay of canonical Airy tails, cuts and connector support.

The three written proofs supply the analytic arguments. This wrapper pins
every proof/source/receipt and replays all literal algebra/rational controls.
It proves no original-period remainder or Strong Factorial threshold.
"""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Unoptimized execution is required')
from pathlib import Path
from hashlib import sha256
import argparse,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent

LOCAL={
 'NEGATIVE_AIRY_CURVE_OUTER_TAIL.md',
 'verify_negative_airy_curve_outer_tail.py',
 'negative_airy_curve_outer_tail_verification.json',
 'POSITIVE_AIRY_THIN_OUTER_TAIL_AND_CUT_GEOMETRY.md',
 'verify_positive_airy_thin_outer_tail.py',
 'positive_airy_thin_outer_tail_verification.json',
 'COMPACT_AIRY_PHYSICAL_CUT_ADDENDUM.md',
 'verify_compact_airy_physical_cut.py',
 'compact_airy_physical_cut_verification.json',
 'verify_airy_outer_tail_and_cut_package.py',
 'airy_outer_tail_and_cut_inputs.json',
}
REPO={
 'research/factorial/image-rigidity/EXPLICIT_AIRY_CONTOUR_ENVELOPE.md',
 'research/factorial/image-rigidity/verify_explicit_airy_contour_envelope.py',
 'research/factorial/image-rigidity/explicit_airy_contour_envelope_source_ledger.json',
 'research/factorial/image-rigidity/explicit_airy_contour_envelope_verification.json',
}

def sha(p):return sha256(p.read_bytes()).hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=HERE/'airy_outer_tail_and_cut_verification.json')
    args=ap.parse_args();repo=args.repo_root.resolve()
    lp=HERE/'airy_outer_tail_and_cut_source_ledger.json'
    ledger=json.loads(lp.read_text());assert ledger.get('schema')==1
    required={('local',p) for p in LOCAL}|{('repo',p) for p in REPO}
    assert len(ledger['files'])==len(required)
    assert {(r['root'],r['path']) for r in ledger['files']}==required
    for row in ledger['files']:
        assert row['root'] in {'local','repo'}
        p=(HERE if row['root']=='local' else repo)/row['path']
        assert p.is_file() and sha(p)==row['sha256'],row['path']
    inp=json.loads((HERE/'airy_outer_tail_and_cut_inputs.json').read_text())
    assert inp=={
     'schema':1,
     'common_domain':{'m_min':'1','zeta_modulus_max':'1/8','outer_radius':'3/8','connector_radius':'41/100'},
     'positive_domain':{'arg_abs_max':'5*pi/6','imaginary_modulus_max':'1/50000'},
     'negative_domain':{'arg_minus_zeta_abs_max':'pi/6'},
     'compact_domain':{'m_min':'10000000000000000','scaled_argument_modulus_max':'4'},
     'scope':'canonical weighted Airy tails, certified cuts, conditional circle connectors and physical-disk Cauchy replacement only',
     'coverage':['negative canonical tails with weights0,1,2','thin positive canonical tails with weights0,1','compact canonical tails with weights0,1','right and left radius-0.41 cut and connector inequalities','physical-disk Cauchy replacement under certified endpoint matching'],
     'excludes':['original-period nonlocal estimate','actual complete-period remainder','numerical Strong Factorial cutoff'],
    }
    children=[
     ('verify_negative_airy_curve_outer_tail.py','negative_airy_curve_outer_tail_verification.json'),
     ('verify_positive_airy_thin_outer_tail.py','positive_airy_thin_outer_tail_verification.json'),
     ('verify_compact_airy_physical_cut.py','compact_airy_physical_cut_verification.json'),
    ]
    results=[]
    with tempfile.TemporaryDirectory(prefix='airy_outer_tail_replay_') as temp:
        for script,receipt in children:
            output=Path(temp)/receipt
            run=subprocess.run([sys.executable,str(HERE/script),'--output',str(output)],capture_output=True,text=True)
            assert run.returncode==0,(script,run.stdout,run.stderr)
            assert output.read_bytes()==(HERE/receipt).read_bytes(),script
            result=json.loads(output.read_text())
            assert all(result[k] is False for k in ['formal_verification','external_review','independent_replay'])
            results.append({'script':script,'receipt':receipt,'sha256':sha(output),'data':result})
        # Replay the imported Airy package as well; its own mandatory ledger
        # validates the actual cached DLMF primary statements and endpoint proof.
        imported_checker=repo/'research/factorial/image-rigidity/verify_explicit_airy_contour_envelope.py'
        imported_receipt=repo/'research/factorial/image-rigidity/explicit_airy_contour_envelope_verification.json'
        imported_output=Path(temp)/'imported_airy.json'
        run=subprocess.run([sys.executable,str(imported_checker),'--repo-root',str(repo),'--output',str(imported_output)],capture_output=True,text=True)
        assert run.returncode==0,(run.stdout,run.stderr)
        assert imported_output.read_bytes()==imported_receipt.read_bytes()
        imported_replay={'receipt_sha256':sha(imported_output),'data':json.loads(imported_output.read_text())}
    neg,pos,compact=[r['data'] for r in results]
    assert neg['certificates']['phase_loss_identity'] is True
    assert len(neg['rational_guards'])==13 and neg['negative_controls_rejected']==3
    assert pos['universal_identities']==2 and len(pos['rational_guards'])==38 and pos['negative_controls']==3
    assert compact['universal_identities']==1 and len(compact['rational_guards'])==11 and compact['negative_controls']==2
    result={'schema':1,'status':'PASS_CANONICAL_AIRY_TAILS_AND_CUTS',
     'scope':inp['scope'],'coverage':inp['coverage'],'excluded_conclusions':inp['excludes'],
     'source_ledger_sha256':sha(lp),'pins':ledger['files'],'replays':results,'imported_airy_replay':imported_replay,
     'formal_verification':False,'external_review':False,'independent_replay':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','replayed_sources':3,'rational_guards':62,'negative_controls':8,'mandatory_pins':len(required)}))
if __name__=='__main__':main()
