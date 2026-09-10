"""Prospective halving admission using an explicitly supplied, sealed footprint.

No prime selection, parameter/orbit enumeration, candidate lookup, or V3 reads.
The existing split_seed_descent engine supplies exact halving/cycle arithmetic.
Only three mathematical statuses leave this API; exhausted bounds stay UNKNOWN.
"""
from split_seed_descent import Classifier, decode_point, point_record, require


def frame_from_packet(packet):
    """Translate an existing finite rank packet, preserving every place and row.

    This is a format adapter, not a new footprint calculation. Classifier
    independently reconstructs and validates the supplied finite groups.
    """
    basis = packet['points']
    proof = packet['proof']
    records, rows = [], []
    for signature in proof['signatures']:
        dim = signature['quotient_dimension']
        local = signature['rows']
        require(len(local) == dim and all(len(row) == len(basis) for row in local),
                'signature shape')
        require(all(type(v) is int and v in (0, 1) for row in local for v in row),
                'binary signature')
        require(signature['doubled_subgroup_order'] * 2**dim == signature['group_order'],
                'quotient order')
        records.append(dict(prime=signature['prime'], order=signature['group_order'],
                            dimension=dim, codes=[sum(local[j][i] << j for j in range(dim))
                                                  for i in range(len(basis))]))
        rows.extend(local)
    require(len({r['prime'] for r in records}) == len(records), 'duplicate footprint place')
    return dict(curve=packet['curve'], basis=basis, inherited_rank=len(basis),
                records=records, finite_rows=rows,
                no_two_torsion_prime=proof['no_rational_2_torsion_prime'],
                candidate_inputs=0)


class ProspectiveAdmission:
    """Reuse a validated M/2M injection; deficient frames are rejected."""
    def __init__(self, frame):
        require(len({r['prime'] for r in frame['records']}) == len(frame['records']),
                'duplicate footprint place')
        self.engine = Classifier(frame)

    def consider(self, point, *, max_steps=8):
        require(type(max_steps) is int and max_steps >= 0, 'nonnegative halving bound')
        result = self.engine.classify(decode_point(self.engine.E, point), max_steps)

        if result['status'] == 'UNKNOWN_STEP_CAP':
            result['status'] = 'UNKNOWN'
            result['reason'] = 'HALVING_STEP_CAP'
        result['max_successful_halves'] = max_steps
        # Even at a zero-step cap retain the uniquely forced subtraction.
        if result['status'] == 'UNKNOWN':
            terminal = decode_point(self.engine.E, point)
            for step in result['history']:
                terminal = decode_point(self.engine.E, step['next'])
            code = self.engine.code(terminal)
            bits = self.engine.M.solve_right(code)
            word = list(map(int, bits))
            result['pending'] = dict(point=point_record(terminal), code=list(map(int, code)),
                                     parity_word=word,
                                     target=point_record(terminal-self.engine.word_point(word)))
        return result
