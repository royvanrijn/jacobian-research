"""Exact necessary source-lattice data for PF103; no realization classifier."""
import math
import sympy as S


def canonical_packet(gram, indices):
    """gram=-Q in the compact basis dual to ALL omitted primes.

    Return necessary marked homology/canonical data only. Passing does not
    construct a surface, finite map, or contact structure on a proposed end.
    """
    B=S.Matrix(gram)
    if B.rows==0 or B.rows!=B.cols or B!=B.T:
        raise ValueError("nonempty symmetric square source Gram matrix required")
    if any(not x.is_Integer for x in B):
        raise ValueError("integral Gram matrix required")
    r=B.rows
    if len(indices)!=r or any(not S.sympify(e).is_Integer or e<1 for e in indices):
        raise ValueError("one positive integral index per omitted prime required")
    if any(B[:j,:j].det()<=0 for j in range(1,r+1)):
        raise ValueError("negative source intersection form required")
    k=S.Matrix([e-1 for e in indices])
    if any((B[i,i]-k[i])%2 for i in range(r)):
        raise ValueError("canonical characteristic parity fails")
    inv=B.inv();energy=(k.T*inv*k)[0]
    return dict(rank=r,canonical_energy=str(energy),theta=str(r-2-energy),
                H1_order=int(B.det()),
                core_orders=[math.lcm(*(int(inv[j,i].q) for j in range(r)))
                             for i in range(r)],
                canonical_nonzero=any(k),
                necessary_only=True)


def a_type_comparison(packet,p):
    """Compare with ORIENTED negative A_(p-1) lens end, not merely H1."""
    if not isinstance(p,int) or p<2:
        raise ValueError("A-type lens parameter must be an integer >=2")
    return dict(H1_order_matches=packet["H1_order"]==p,
                minimal_filling_rank_bound=packet["rank"]<=p-1,
                contact_theta_matches=S.Rational(packet["theta"])==p-3,
                actual_nonautomorphic_requires_positive_energy=
                    S.Rational(packet["canonical_energy"])>0,
                necessary_only=True)
