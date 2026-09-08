#!/usr/bin/env sage-python
"""Exact field square-class and valuation certificate proving class-2-rank >=16."""
import argparse
from math import prod
from pathlib import Path
import runpy
from sage.all import QQ, GF, matrix, pari
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint

ROOT,ART,cert = batch.ROOT,batch.ART,batch.cert
SOURCE = Path(__file__).resolve()
OUT = ART/'small_conductor_class_lower16_v1.json'
OLD_SOURCE = ROOT/'elliptic-curves/cas/target_small_conductor_small_base.sage'
PROOF = batch.forms.target.PROOF


def expected():
    old = runpy.run_path(str(OLD_SOURCE))
    data,poly,nf,theta,w = old['setup']()
    proof = cert.read(PROOF)
    # The original pure-rational checker proves the points, discriminant factors,
    # primes and finite-quotient rank; the new character matrix also verifies
    # independence of their field square classes directly.
    curve = runpy.run_path(str(ROOT/'elliptic-curves/cas/certify_small_conductor_curve.py'))
    curve['verify'](proof)
    if int(pari.poldisc(poly))!=256*int(proof['discriminant']):
        raise ArithmeticError('division polynomial discriminant identity failed')
    points = [(QQ(x),QQ(y)) for x,y in proof['integral_points']]
    betas = [pari(4*x)-theta for x,y in points]
    norm_roots = []
    for (x,y),beta in zip(points,betas):
        z = 8*y+4*x
        if pari.nfeltnorm(nf,beta)!=pari(z*z) or z==0:
            raise ArithmeticError('point norm square failed')
        norm_roots.append(str(z))
    valuations = []
    for q,e in proof['discriminant_factorization']:
        q = int(q)
        for P in pari.idealprimedec(nf,q):
            values = [int(pari.idealval(nf,b,P)) for b in betas]
            valuations.append({'p':str(q),'hnf':str(pari.idealhnf(nf,P)),'valuations':values})
    V = matrix(GF(2),[[v%2 for v in r['valuations']] for r in valuations])
    kernel = V.right_kernel().basis_matrix()
    if V.rank()!=4 or kernel.nrows()!=18 or kernel.rank()!=18 or V*kernel.transpose()!=0:
        raise ArithmeticError('valuation kernel differs')
    evaluations = []
    # Fixed finite probe set: first128 rational primes >=1000003. All primality
    # tests are proved, not probable; discard a probe root if any beta has zero
    # or undefined residue. A residue character cannot see a global square.
    q = 1000003
    primes = []
    while len(primes)<128:
        if pari.isprime(q): primes.append(q)
        q += 1
    for q in primes:
        for root in pari.polrootsmod(poly,q):
            root = int(pari.lift(root)); residues = []
            for x,y in points:
                z = 4*x
                if int(z.denominator())%q==0:
                    break
                residues.append((int(z.numerator())*pow(int(z.denominator()),-1,q)-root)%q)
            if len(residues)!=22 or 0 in residues:
                continue
            values = [0 if pow(v,(q-1)//2,q)==1 else 1 for v in residues]
            evaluations.append({'p':q,'theta_root':root,'residues':residues,'characters':values})
    C = matrix(GF(2),[r['characters'] for r in evaluations])
    if C.rank()!=22 or (C*kernel.transpose()).rank()!=18:
        raise ArithmeticError('field square-class independence not certified')
    return {'schema':'elliptic-curves.small-conductor-class-lower16.v1','status':'PASS',
        'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in [SOURCE,OLD_SOURCE,PROOF,batch.forms.OUT]},
        'field_discriminant':str(nf.disc()),'beta_formula':'beta_i=4*x(P_i)-theta',
        'norm_square_roots':norm_roots,'bad_prime_valuations':valuations,
        'valuation_matrix_rank':4,'even_valuation_kernel_basis':[[int(v) for v in row] for row in kernel.rows()],
        'kernel_dimension':18,'probe_primes':primes,'residue_character_evaluations':evaluations,
        'field_square_class_rank':22,'kernel_square_class_rank':18,
        'away_from_bad_primes_argument':'For p not dividing2*Delta(E), F is separable mod p. If x(P) is p-integral, every prime ideal dividing4*x(P)-theta has residue degree1 and is the unique prime corresponding to that rational residue root; its valuation is even because the norm is a rational square. If x(P) has a pole, its valuation is even on an integral Weierstrass model, and theta is integral, so every valuation of4*x(P)-theta above p is even. All primes dividing2*Delta(E) are explicitly checked.',
        'class_rank_argument':'The18 independent kernel products have even valuations at every finite prime. Thus they lie in Sel_2(K). Their norms are positive rational squares, whereas Norm(-1)=-1 in this cubic field, so adjoining -1 gives19 independent field Selmer classes. The unit square-class dimension is r1+r2=3. The exact sequence O_K^*/O_K^{*2} -> Sel_2(K) -> Cl(K)[2] therefore proves dim Cl(K)[2]>=19-3=16.',
        'unconditional_class_two_rank_lower_bound':16,
        'references':['https://arxiv.org/html/1606.07178#S4.SS3'],
        'claim_boundary':'Unconditional lower bound16 for the same cubic field. This is not a class-group upper bound or an exact curve rank. It makes the sufficient upper-bound target16 sharp.'}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');a=p.parse_args()
    result=expected()
    if a.check:
        if cert.read(OUT)!=result:raise ArithmeticError('lower bound certificate differs')
    else:
        if OUT.exists():raise FileExistsError('preserve certificate')
        checkpoint(OUT,result)
    print('CLASS TWO-RANK >=16 UNCONDITIONALLY: PASS',flush=True)
