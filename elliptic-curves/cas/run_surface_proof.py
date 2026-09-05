#!/usr/bin/env python3
"""Import or replay retained finite-field proofs, then apply the regulator gate.

Requests name an exact surface and retained Frobenius proofs. Discovery of a
new Frobenius polynomial is a separate explicitly budgeted command. Legacy
ToricControlledReduction certificates can be imported without rerunning CAS.
"""
import argparse
from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
import shlex
import sys

from research_runtime.regulator import Surface, pre_search_gate
from research_runtime.store import checkpoint, default_store, digest
from research_runtime.supervisor import Limits, capture, preserve_previous

ROOT = Path(__file__).resolve().parents[2]


def import_certificate(path):
    from research_runtime.section_gate import surface_from_export
    from research_runtime.toric_surface import PROOF_KIND
    certificate = json.loads(path.read_text())
    if certificate.get('status') != 'PASS_COMPLETE_FROBENIUS_PICARD_BOUND' or certificate.get('good_reduction',{}).get('status') != 'PASS':
        raise ArithmeticError('legacy Frobenius certificate is incomplete')
    for name, expected in certificate['inputs'].items():
        source = (ROOT/name).resolve()
        if not source.is_relative_to(ROOT) or sha256(source.read_bytes()).hexdigest() != expected:
            raise ArithmeticError('legacy Frobenius input changed')
    key = certificate['pair_key']
    surface = surface_from_export({'inputs':certificate['inputs'],
        'candidate':{'kind':'product' if ':' in key else 'singleton','key':key}}, ROOT)
    if surface is None:
        raise ValueError('legacy certificate has no exact surface transport; supply a generic request')
    software = certificate['software']
    raw_path = Path(shlex.split(software['frobenius_invocation'])[-1])
    raw = raw_path.read_text()
    if sha256(raw.encode()).hexdigest() != software['ToricControlledReduction_output_sha256']:
        raise ArithmeticError('retained toric output differs from its certificate')
    proof = {'proof_kind':PROOF_KIND, 'surface_key':surface.key, 'prime':certificate['prime'],
        'coefficients':certificate['elliptic_L']['frobenius_characteristic_coefficients_low_to_high'],
        'moments':certificate['elliptic_L']['power_sums_n1_n2'], 'raw_output':raw,
        'provenance':{'source_commit':software['ToricControlledReduction_commit'],
            'executable_sha256':software['ToricControlledReduction_executable_sha256'],
            'raw_output_sha256':software['ToricControlledReduction_output_sha256'],
            'nondegenerate_driver_completed':certificate['toric_nondegeneracy_jacobian_reduction']=='PASS',
            'legacy_certificate':str(path.resolve()),'legacy_certificate_sha256':sha256(path.read_bytes()).hexdigest()}}
    return surface, proof


def worker(args):
    from research_runtime.surface_repository import SurfaceProofRepository
    from research_runtime.toric_surface import replay_toric
    repository = SurfaceProofRepository()
    if args.import_certificate:
        packets = []
        for path in args.import_certificate:
            surface, proof = import_certificate(path)
            packets.append({'surface':asdict(surface),'proof':proof})
        request = {'packets':packets}
    else:
        request = json.loads((args.replay_request or args.request).read_text())
        if 'packets' in request:
            packets = request['packets']
        else:
            packets = [{'surface':request['surface'],'proof':p} for p in request['reductions']]
    surfaces = {}
    for packet in packets:
        surface = Surface(**packet['surface']); surfaces[surface.key] = surface
        repository.retain(surface,packet['proof'],verify=replay_toric)
    reductions, gates = [], []
    for surface in surfaces.values():
        rows = repository.replay(surface)
        reductions.extend(asdict(r) for r in rows)
        gates.append(pre_search_gate(surface,rows,candidate_rank=request.get('candidate_rank',1),
                                     candidate_regulator=request.get('candidate_regulator')))
    result = {'schema':'elliptic-curves.finite-field-surface-proof.v1','status':'PASS',
              'request_hash':digest(request),'reductions':reductions,'gates':gates,
              'regenerated_frobenius':False,'retained_packets':packets}
    preserve_previous(args.output)
    checkpoint(args.output,result)
    print('SURFACE_PROOF|status=PASS|surfaces='+str(len(surfaces))+'|reductions='+str(len(reductions)),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--request',type=Path)
    group.add_argument('--replay-request',type=Path)
    group.add_argument('--import-certificate',type=Path,nargs='+')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--wall-seconds',type=int,default=45)
    parser.add_argument('--rss-bytes',type=int,default=1_073_741_824)
    parser.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker or args.replay_request:
        worker(args)
    else:
        result = capture(['sage','-python',str(Path(__file__).resolve()),*sys.argv[1:],'--worker'],
            limits=Limits(args.wall_seconds,args.rss_bytes),log_path=args.output.with_suffix('.log'))
        print(result.stdout,end='')


if __name__ == '__main__':
    main()
