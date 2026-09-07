#!/usr/bin/env python3
"""Check hypotheses of the simple-collision 4-division compositum lemma."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc

OUTPUT=r.OUT/'rank_jump_four_division_separation_v1.json'


def build(check=False):
    inp=r.read(lc.INPUT);model=inp['anchor']['short_model_ainvariants'];A,B=map(r.F,model[3:])
    gal=r.galois(model);assert gal['irreducibility_prime'] is not None
    delta=-4*A**3-27*B**2;rows=[]
    for row in inp['rows']:
        u=r.F(row['parameter_u'])
        if not u:continue
        D=1+A*u*u+B*u**3
        def valid(p):
            if p<3 or any(p%d==0 for d in range(2,r.isqrt(p)+1)):return False
            if not u.numerator%p or not u.denominator%p:return False
            if A.denominator%p==0 or B.denominator%p==0 or r.mod(delta,p)==0:return False
            return D.denominator%p!=0 and D.numerator%p==0 and D.numerator%(p*p)!=0
        p=next(p for p in sorted(row['newly_bad_primes_relative_to_anchor']) if valid(p))
        root=pow(r.mod(u,p),-1,p)
        assert (root**3+r.mod(A,p)*root+r.mod(B,p))%p==0
        derivative=(3*root*root+r.mod(A,p))%p;assert derivative
        rows.append({'u':str(u),'simple_good_collision_prime':p,'p_is_prime_by_trial_division':True,
          'D':str(D),'root_mod_p':root,'f_derivative_mod_p':derivative,
          'v_p_D':1,'anchor_good_reduction':True,
          'valuation_witness_matrix_at_three_conjugate_primes':[[1,0,0],[0,1,0],[0,0,1]],
          'matrix_justification':'Good odd-prime reduction makes Q(E0[4])/Q unramified here. Exactly one 1-u*theta_i has valuation 1; transitivity of the roots supplies three conjugate prime ideals isolating the three factors.',
          'relative_degree_over_anchor_four_division_field':8,
          'relative_galois_group':'(Z/2)^3'})
    out={'schema':'rank-jump.four-division-separation.v1','input_sha256':r.digest(lc.INPUT.read_bytes()),
      'script_sha256':r.digest(Path(__file__).read_bytes()),'anchor_galois':gal,'rows':rows,
      'field_formula':'F0*Q(Eu[4]) = F0(sqrt(1-u*theta_1),sqrt(1-u*theta_2),sqrt(1-u*theta_3)), F0=Q(E0[4])',
      'proof':'elliptic-curves/rank-jump/LOCAL_COLLISION_AND_RECIPROCITY.md#the-preserved-cubic-hides-a-full-change-in-four-division-data',
      'source_formula':'https://arxiv.org/pdf/1704.06190, Theorem 1(a)',
      'claim_boundary':'Exact consequences of the written lemma and checked prime hypotheses, not explicit construction of the division fields. The three quadratic field directions are not Mordell-Weil directions; no rank or CT pairing value follows.'}
    if check:
        if r.read(OUTPUT)!=out:raise ValueError('four-division certificate mismatch')
        print('PASS four-division separation replay')
    else:r.write_new(OUTPUT,out)
    for row in rows:print('u',row['u'],'prime',row['simple_good_collision_prime'],'compositum degree',row['relative_degree_over_anchor_four_division_field'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('build','check'));a=p.parse_args();build(a.mode=='check')
