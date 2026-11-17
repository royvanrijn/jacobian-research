"""Exact ADE comparison data; oriented end recognition is an external input."""
from fractions import Fraction as F


def ade_data(kind, rank=None):
    if kind=='A':
        if not isinstance(rank,int) or rank<1:raise ValueError('A rank >=1 required')
        return dict(kind=kind,rank=rank,arm_lengths=None,H1_order=rank+1,
                    canonical_theta=rank-2,unique_positive_tight=True,
                    minimal_filling_rank=rank,orientation='negative ADE link')
    if kind=='D':
        if not isinstance(rank,int) or rank<4:raise ValueError('D rank >=4 required')
        arms=[1,1,rank-3]
    elif kind in ['E6','E7','E8']:
        want=int(kind[1:])
        if rank is not None and rank!=want:raise ValueError('exceptional rank mismatch')
        rank=want;arms={'E6':[1,2,2],'E7':[1,2,3],'E8':[1,2,4]}[kind]
    else:raise ValueError('only A, D, E6, E7, E8 are certified')
    schur=2-sum((F(l,l+1) for l in arms),F(0))
    order=schur
    for l in arms:order*=l+1
    assert order.denominator==1 and schur>0
    return dict(kind=kind,rank=rank,arm_lengths=arms,H1_order=int(order),
                Seifert_e0=-2,Seifert_fractions=[str(F(l,l+1)) for l in arms],
                orbifold_Euler=str(schur),central_Schur=str(schur),
                rational_Seifert_Euler=str(-schur),canonical_theta=rank-2,
                unique_positive_tight=True,minimal_filling_rank=rank,
                orientation='negative ADE link')


def cartan(kind, rank=None):
    import sympy as S
    data=ade_data(kind,rank);n=data['rank'];C=2*S.eye(n)
    if kind=='A':edges=[(i,i+1) for i in range(n-1)]
    else:
        edges=[];offset=1
        for length in data['arm_lengths']:
            prev=0
            for i in range(offset,offset+length):edges.append((prev,i));prev=i
            offset+=length
    for i,j in edges:C[i,j]=C[j,i]=-1
    return C


def compare(packet,kind,rank=None,*,oriented_end_identified=False):
    import sympy as S
    data=ade_data(kind,rank)
    return dict(eligible_for_obstruction=oriented_end_identified,
                H1_order_matches=packet['H1_order']==data['H1_order'],
                minimal_rank_matches=packet['rank']==data['minimal_filling_rank'],
                contact_theta_matches=S.Rational(packet['theta'])==data['canonical_theta'],
                positive_source_energy=S.Rational(packet['canonical_energy'])>0,
                necessary_only=True,realization_claimed=False,
                orientation=data['orientation'])
