"""Independent polynomial checks for the q-first residue partitions."""
from fractions import Fraction as F
from math import gcd, isqrt, lcm

from half_lattice_pointed_sieve import binary_transform
from verify_finite_cancellation import constant, valuation


def verify_tree(tree,n,d,q,with_g=False):
    p=tree['prime'];tries={'affine':{},'infinity':{}};mass=F(0)
    for leaf in tree['leaves']:
        k=leaf['depth'];r=leaf['residue'];mass+=F(leaf['mass'])
        assert F(leaf['mass'])==F(1,(p+1)*p**(k-1)) and 0<=r<p**k
        trie=tries[leaf['kind']]
        for j in range(k):
            assert 'terminal' not in trie
            trie=trie.setdefault((r//p**j)%p,{})
        assert not trie;trie['terminal']=True
        matrix=(p**k,r,0,1) if leaf['kind']=='affine' else (0,1,p**k,r)
        if leaf['kind']=='infinity':assert r%p==0
        qq=binary_transform(q,matrix)
        if leaf['status']=='square':
            vq=constant(qq,p,3 if p==2 else 1)
            assert vq is not None and vq%2==0 and vq==leaf['vq']
            unit=qq[0]//p**vq
            assert unit%8==1 if p==2 else pow(unit%p,(p-1)//2,p)==1
        elif leaf['status']=='nonsquare':
            vq=constant(qq,p);assert vq is not None
            if vq%2==0:
                assert constant(qq,p,3 if p==2 else 1)==vq
                unit=qq[0]//p**vq
                assert unit%8!=1 if p==2 else pow(unit%p,(p-1)//2,p)!=1
        else:assert leaf['status']=='unknown'
        if with_g and leaf['vg'] is not None:
            nn=binary_transform(n,matrix);dd=binary_transform(d,matrix)
            v=min(valuation(x,p) for x in nn+dd)
            exact=constant(nn,p)==v or constant(dd,p)==v
            assert leaf['vg']==v and leaf['vg_exact']==exact
            if leaf['status']=='square':assert exact
        elif with_g:assert leaf['status']!='square'
    assert mass==1
    return len(tree['leaves'])


def verify_prepared(prepared):
    A,B=map(F,prepared['curve'][3:]);a,b=map(F,prepared['anchor'])
    assert b*b==a**3+A*a+B
    leaves=0
    for m in prepared['models']:
        n,d,q=([int(v) for v in m[k]] for k in ('N','D','q'))
        matrix=tuple(map(F,m['mapping']['matrix']))
        nn=binary_transform((a**3+4*B,4*a*b,6*a*a+4*A,4*b,a),matrix)
        dd=binary_transform((-3*a*a-4*A,-8*b,-6*a,0,1),matrix)
        den=lcm(*(v.denominator for v in nn+dd));z=[int(v*den) for v in nn+dd];content=gcd(*z)
        assert [v//content for v in z]==n+d and gcd(*(n+d))==1
        ratio=F(m['mapping']['square_ratio']);assert ratio>0 and list(dd)==[ratio*v for v in q]
        assert isqrt(ratio.numerator)**2==ratio.numerator and isqrt(ratio.denominator)**2==ratio.denominator
        for t in m['local']['q']:leaves+=verify_tree(t,n,d,q)
        for t in m['local'].get('g',[]):leaves+=verify_tree(t,n,d,q,True)
    return leaves
