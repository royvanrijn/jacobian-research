#!/usr/bin/env python3
"""Independent complete branch partitions and Horner polynomial composition."""
import classify_mw16_score_prime_scalings as audit
import certify_compact_r17_candidates as cert


def value(a, x, m):
    return sum(c*pow(x, i, m) for i, c in enumerate(a)) % m


def compose(a, r, step):
    out = [a[-1]]
    for coefficient in a[-2::-1]:
        nxt = [0]*(len(out)+1)
        for i, c in enumerate(out):
            nxt[i] += r*c
            nxt[i+1] += step*c
        nxt[0] += coefficient
        out = nxt
    return out


def main():
    p = audit.protocol()
    d = cert.read(audit.OUT)
    if (d['sources'] != audit.sources() or d['protocol_sha256'] != cert.hashed(audit.D/'protocol.json')
            or len(d['rows']) != 33):
        raise ArithmeticError('all33 bound classification rows required')
    families = {f['fibration_id']: f for f in cert.read(audit.support.INPUT)['families']}
    for pair, row in zip(p['pairs'], d['rows']):
        if any(row[k] != v for k, v in pair.items()) or [c['chart'] for c in row['charts']] != ['affine', 'infinity']:
            raise ArithmeticError('pair or chart ordering differs')
        f = families[row['family']]
        prime = row['prime']
        total = 0
        for chart in row['charts']:
            # A complete claim is accepted only when every branch is discharged.
            if chart['status'] != 'COMPLETE_RESIDUE_CLASSIFICATION':
                raise ArithmeticError('incomplete chart remains UNKNOWN')
            a = list(map(int, f['A_coefficients_low_to_high']))
            b = list(map(int, f['B_coefficients_low_to_high']))
            if chart['chart'] == 'infinity':
                a, b = a[::-1], b[::-1]
            candidates = [0] if chart['chart'] == 'infinity' else list(range(prime))
            balls = []
            if not 1 <= len(chart['levels']) <= 6:
                raise ArithmeticError('fixed depth bound differs')
            for depth, level in enumerate(chart['levels'], 1):
                step = prime**depth
                excluded, admitted, live = [], [], []
                for r in candidates:
                    if value(a, r, prime**min(depth, 4)) or value(b, r, prime**min(depth, 6)):
                        excluded.append(r)
                        continue
                    aa, bb = compose(a, r, step), compose(b, r, step)
                    if all(c % prime**4 == 0 for c in aa) and all(c % prime**6 == 0 for c in bb):
                        admitted.append(r)
                        balls.append({'depth': depth, 'residue': r, 'modulus': step,
                                      'A_divided_coefficients': [str(c//prime**4) for c in aa],
                                      'B_divided_coefficients': [str(c//prime**6) for c in bb]})
                    else:
                        live.append(r)
                if level != {'depth': depth, 'modulus': step, 'excluded_residues': excluded,
                             'admitted_residues': admitted, 'unresolved_residues': live}:
                    raise ArithmeticError('complete independent branch partition differs')
                candidates = sorted(r+step*j for r in live for j in range(prime))
            if candidates or chart['scale_balls'] != balls:
                raise ArithmeticError('unresolved branches or polynomial witnesses differ')
            total += len(balls)
        status = 'CLASSIFIED_SCALE_BALLS' if total else 'NO_REMOVABLE_SCALE'
        if row['scale_balls'] != total or row['status'] != status:
            raise ArithmeticError('pair classification differs')
    print('INDEPENDENT33 COMPLETE MW16 SCALING CLASSIFICATIONS PASS', flush=True)


if __name__ == '__main__':
    main()
