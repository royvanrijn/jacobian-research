"""Preload a fixed finite prime bank without changing the rational subgroup."""
from dataclasses import replace
from mod2_reduction_independence import _is_prime


def preload(state,cache,bound=997):
    if type(bound) is not int or not 3<=bound<=1000:raise ValueError('bounded bank required')
    reductions=state.reductions;added=[];excluded=[]
    for prime in range(3,bound+1):
        if not _is_prime(prime) or prime in reductions.primes:continue
        try:reductions=reductions.escalate(prime,cache)
        except ValueError as e:
            excluded.append({'prime':prime,'reason':str(e)});continue
        added.append(prime)
    if reductions.points!=state.reductions.points or not reductions.independent_images:raise ArithmeticError('prime preload changed subgroup')
    result=replace(state,reductions=reductions,parent_state=state.key)
    return result,{'bound':bound,'original_primes':list(state.reductions.primes),'added_primes':added,'excluded_primes':excluded,'final_primes':list(reductions.primes),'basis_unchanged':True}
