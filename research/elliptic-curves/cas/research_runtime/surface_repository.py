"""Portable finite-reduction proof packets, queried before section search."""
from dataclasses import asdict
from pathlib import Path
import json
import os
import shutil

from .regulator import Surface, VerifiedReduction
from .store import FactStore, atomic_write, default_store, digest, encoded


class SurfaceProofRepository:
    def __init__(self, store=None):
        self.store = store or default_store()

    def retain(self, surface, proof, *, verify):
        if proof.get('surface_key') != surface.key or verify(surface, proof) is not True:
            raise ArithmeticError('surface proof failed before publication')
        packet = {'surface': asdict(surface), 'proof': proof}
        key = digest(packet)
        self.store.discover('surface-proof-packet', key, lambda: packet)
        path = self.store.root/'surface-proofs'/f'{surface.key}.json'
        with self.store.lock({'surface-proof-index': surface.key}):
            entries = json.loads(path.read_text()) if path.exists() else []
            if key not in entries:
                atomic_write(path, encoded(sorted([*entries, key])))
        return key

    def packets(self, surface):
        path = self.store.root/'surface-proofs'/f'{surface.key}.json'
        if not path.exists():
            return []
        packets = [self.store.require('surface-proof-packet', key) for key in json.loads(path.read_text())]
        if any(Surface(**p['surface']) != surface for p in packets):
            raise ArithmeticError('surface proof index is misbound')
        return packets

    def replay(self, surface):
        from .sage_surface import SurfaceProofEngine
        from .toric_surface import PROOF_KIND, replay_toric
        engine = SurfaceProofEngine(self.store); reductions = {}
        for packet in self.packets(surface):
            proof = packet['proof']
            if proof.get('proof_kind') != PROOF_KIND:
                raise ValueError('unregistered surface Frobenius proof kind')
            row = engine.reduction(surface, proof, verify_frobenius=replay_toric, discover=True)
            previous = reductions.get(row.prime)
            if previous is not None and (previous.arithmetic_rank_upper, previous.geometric_rank_upper,
                previous.regulator_if_rank_equality) != (row.arithmetic_rank_upper,row.geometric_rank_upper,row.regulator_if_rank_equality):
                raise ArithmeticError('inconsistent proofs at the same good prime')
            reductions[row.prime] = row
        return tuple(reductions[p] for p in sorted(reductions))


_replayed = {}


def available_reductions(surface, *, store=None):
    """Replay retained proofs in one bounded Sage worker; never start a census.

    Empty stores are an explicit missing-data result. Replay failure remains
    UNKNOWN and records its reason; it cannot become a rank exclusion.
    """
    from .supervisor import Limits, capture
    from .store import checkpoint
    repository = SurfaceProofRepository(store)
    packets = repository.packets(surface)
    if not packets:
        return (), {'status':'NO_RETAINED_REDUCTIONS', 'attempted_before_search':True}
    implementation = {}
    from hashlib import sha256
    for name in ('toric_surface.py','sage_surface.py','regulator.py','surface_repository.py'):
        implementation[name] = sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
    key = digest({'packets':packets, 'replayer':implementation})
    if key in _replayed:
        return _replayed[key]
    directory = repository.store.root/'surface-replay'/key
    request_path, output_path = directory/'request.json', directory/'result.json'
    checkpoint(request_path, {'surface':asdict(surface), 'packets':packets})
    command = Path(__file__).resolve().parents[1]/'run_surface_proof.py'
    sage = os.environ.get('SAGE', 'sage')
    if not shutil.which(sage):
        return (), {'status':'REPLAY_UNAVAILABLE', 'reason':'Sage is unavailable', 'attempted_before_search':True}
    try:
        result = capture([sage,'-python',str(command),'--replay-request',str(request_path),'--output',str(output_path)],
            limits=Limits(45,1_073_741_824),log_path=directory/'replay.log')
        data = json.loads(output_path.read_text())
        if data['request_hash'] != digest(json.loads(request_path.read_text())) or data['status'] != 'PASS':
            raise ArithmeticError('misbound finite-reduction replay result')
        rows = tuple(VerifiedReduction(**r) for r in data['reductions'])
        outcome = rows, {'status':'REPLAYED_RETAINED_REDUCTIONS','primes':[r.prime for r in rows],
                         'attempted_before_search':True,'log':str(directory/'replay.log')}
        _replayed[key] = outcome
        return outcome
    except Exception as error:
        return (), {'status':'REPLAY_INCOMPLETE','reason':type(error).__name__+': '+str(error),
                    'attempted_before_search':True,'log':str(directory/'replay.log')}
