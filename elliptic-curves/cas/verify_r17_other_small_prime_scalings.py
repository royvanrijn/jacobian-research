#!/usr/bin/env python3
"""Independently enumerate complete small residue rings for the26 exclusions."""
import certify_compact_r17_candidates as cert
import classify_r17_other_small_prime_scalings as audit


def value(coefficients,x,modulus):
    return sum(c*pow(x,i,modulus) for i,c in enumerate(coefficients))%modulus


def main():
    p=audit.protocol();d=cert.read(audit.OUT)
    if d['sources']!=audit.sources() or d['protocol_sha256']!=cert.hashed(audit.D/'protocol.json') or len(d['rows'])!=26:raise ArithmeticError('complete26 source binding required')
    families={r['family']:r for r in cert.read(audit.support.INPUT)['families']}
    for selected,row in zip(p['pairs'],d['rows']):
        if any(row[k]!=v for k,v in selected.items()) or row['status']!='NO_REMOVABLE_SCALE' or row['scale_balls']!=0 or [c['chart'] for c in row['charts']]!=['affine','infinity']:raise ArithmeticError('fixed exact26 exclusions required')
        f=families[row['family']];prime=row['prime']
        for chart in row['charts']:
            a=list(map(int,f['A_coefficients_low_to_high']));b=list(map(int,f['B_coefficients_low_to_high']))
            if chart['chart']=='infinity':a,b=a[::-1],b[::-1]
            if chart['status']!='COMPLETE_RESIDUE_CLASSIFICATION' or chart['scale_balls'] or not 1<=len(chart['levels'])<=3:raise ArithmeticError('empty classification by depth3 required')
            previous=None
            for depth,level in enumerate(chart['levels'],1):
                modulus=prime**depth;ma=prime**min(depth,4);mb=prime**min(depth,6)
                if modulus>200000:raise ArithmeticError('fixed exhaustive residue verification bound')
                # Enumerate the whole residue ring, independently of the saved lifting tree.
                domain=range(0,modulus,prime) if chart['chart']=='infinity' else range(modulus)
                actual=[r for r in domain if value(a,r,ma)==0 and value(b,r,mb)==0]
                candidates=([0] if chart['chart']=='infinity' else list(range(prime))) if previous is None else sorted(r+prime**(depth-1)*digit for r in previous for digit in range(prime))
                live=set(actual);excluded=[r for r in candidates if r not in live]
                if level!={'depth':depth,'modulus':modulus,'excluded_residues':excluded,'admitted_residues':[],'unresolved_residues':actual}:raise ArithmeticError('independent complete residue set differs')
                previous=actual
            if previous:raise ArithmeticError('nonempty residue obstruction')
    print('INDEPENDENT COMPLETE RESIDUE RINGS: ALL26 NON13 SCALING PAIRS EXCLUDED',flush=True)


if __name__=='__main__':main()
