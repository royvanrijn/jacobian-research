#!/usr/bin/env sage-python
"""Portable exact universal section and parent-path coverage check."""
import argparse,json
from pathlib import Path
from sage.all import QQ,PolynomialRing,prod
def verify(path):
    found=json.loads(path.read_text());V=PolynomialRing(QQ,'v');v=V.gen();F=V.fraction_field()
    R=PolynomialRing(F,'T');T=R.gen();S=PolynomialRing(R,'x');x=S.gen()
    centres=[0,(2*v*v+v+2)**2,2*(v+1)**2*(2*v*v+v+1),4*v*v-v+4,
        v*(2*v-1)*(2*v*v+4*v+5),4*v**4+8*v**3+9*v*v-2*v+2]
    product=prod((x-R(c))**2-T**2 for c in centres)
    coeff={6:R(1)}
    for j in range(5,-1,-1):
        known=sum(coeff[i]*coeff[k] for i in coeff for k in coeff if i+k==6+j)
        coeff[j]=(product[6+j]-known)/2
    square=S([coeff[i] for i in range(7)]);residual=square**2-product
    assert residual.degree()==4
    quartic=[]
    for c in residual.list():
        q,rem=c.quo_rem(T*T);assert rem==0;quartic.append(q)
    e,d,c,b,a=quartic
    A=-27*(12*a*e-3*b*d+c*c)
    B=-27*(72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3)
    X=R(list(map(F,found['anti_X'])));Y=R(list(map(F,found['anti_Y_div_sqrt_minus3'])))
    assert -3*Y**2==X**3+A*X+B and Y!=0 and 4*A**3+27*B**2!=0
    assert (X.degree(),Y.degree())==(4,6) and X(-T)==X and Y(-T)==Y
    # A section over sqrt(-3), conjugated to its group inverse.
    # No rational-rank increase is asserted, including after specialization.
    parameters=PolynomialRing(QQ,'z');z=parameters.gen();L=parameters.fraction_field()
    ratio=-z*(8+3*z)/(6*(z+2)*(z+4))
    assert ratio+QQ(1)/2==(5*z+12)/(3*(z+2)*(z+4))
    W=PolynomialRing(F,'w');w=W.gen();inverse=(6*v+3)*w*w+(36*v+8)*w+48*v
    assert inverse.discriminant()==16*(9*v*v+4)
    denominators=prod(c.denominator() for c in list(X)+list(Y))
    return {'schema':'kihara-universal-anti-replay.v1','status':'PASS',
        'parent_ratio':'v=p/q with q=1','section_field':'Q(v)(sqrt(-3))',
        'anti_X_degree':4,'anti_Y_degree':6,'section_denominator_factors':str(denominators.factor()),
        'section_Y_leading_factors':str(Y.leading_coefficient().factor()),
        'rank14_path_ratio':str(ratio),'path_variable':'z=t^2',
        'path_ratio_plus_half':str(ratio+QQ(1)/2),'inverse_discriminant':'16*(9*v^2+4)',
        'path_coverage':'For every nonzero rational path t, z=t^2>0 and -1/2<v<0. Positive parent-ratio coordinates are absent from that path; this does not prove inequivalence under unmarked parent isomorphisms.',
        'scope':'Exact anti-section over the unrestricted rational parent-ratio function field; the parent formula is reconstructed from six centres without the producer model, normalization, elimination or Groebner basis. The section remains nonrational wherever its specialized Y is nonzero. No full generic rank, full constant field, new rational direction, specialization yield or exhaustive parent classification is claimed.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path);a=p.parse_args();r=verify(a.input)
    if a.output:
        with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print('PASS unrestricted Kihara quadratic section and strict rank14-path coordinate coverage',flush=True)
