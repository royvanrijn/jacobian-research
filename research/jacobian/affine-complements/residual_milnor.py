"""Necessary all-branch Milnor/excess filters; passing never realizes a cover."""
from math import gcd


def residual_gate(e,q,ell,n,c,mu,*,axis_source_certified=False,
                  sole_critical_branch_certified=False,
                  chosen_singular_image_branch_certified=False):
    if any(type(x) is not int for x in (e,q,ell,n,c,mu)) or min(e,q,ell,n)<1 or e<2 or c<2 or c%2 or mu<0 or ell<e*q:
        raise ValueError('exact ramified local data and an even positive branch conductor required')
    eta=ell-e*q;predicted=n*(eta-q)*(c-1)
    cyclic=(n==1 or (n<=e-1 and gcd(n,e)==1))
    minimal=(eta!=q or ((e,q,n) in {(2,1,1),(3,2,1),(3,1,2)} and c==2))
    passes=(cyclic and predicted==mu and minimal)
    certified=axis_source_certified and sole_critical_branch_certified and chosen_singular_image_branch_certified
    return dict(excess=eta,predicted_residual_Milnor=predicted,
                necessary_arithmetic_passes=passes,hypotheses_contradicted=certified and not passes,
                source_hypotheses_certified=certified,minimum_excess=q+(e>=4),
                residual_smooth=(mu==0),minimal_excess=(eta==q),
                one_Puiseux_pair_assumed=False,scalar_Morse_assumed=False,
                realization_claimed=False,arbitrary_finite_germ_excluded=False)


def nonprimitive_residue_degree_floor(f):
    if type(f) is not int or f<2:
        raise ValueError('a nontrivial omitted residue degree f>=2 is required')
    return 5*f+(7*f+4)//5+3


def nonimmersion_cost(index,residue,parameter_count,preimage_count=None):
    if any(type(v) is not int for v in (index,residue,parameter_count)) or index<4 or residue<1 or parameter_count<0:
        raise ValueError('index>=4, positive residue and a nonnegative parameter count required')
    count=parameter_count if preimage_count is None else preimage_count
    if type(count) is not int or not parameter_count<=count<=parameter_count*residue:
        raise ValueError('retain between one and f preimages per normalization parameter')
    return parameter_count*residue+count
