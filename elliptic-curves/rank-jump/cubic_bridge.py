#!/usr/bin/env python3
"""Explicit norm-conic bridge over the odd-degree cubic descent field."""
import argparse
from itertools import permutations
from pathlib import Path
import retrospective as r
import local_collision as lc
import blocks

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'CUBIC_BRIDGE_PROTOCOL.json'
OUTPUT=r.OUT/'rank_jump_cubic_bridge_v1.json'


class Cubic:
    def __init__(self,A,B):
        self.A,self.B=map(r.F,(A,B))
        self.one=(r.F(1),r.F(0),r.F(0))
        self.theta=(r.F(0),r.F(1),r.F(0))

    def scalar(self,a):return (r.F(a),r.F(0),r.F(0))
    def add(self,a,b):return tuple(x+y for x,y in zip(a,b))
    def neg(self,a):return tuple(-x for x in a)
    def sub(self,a,b):return self.add(a,self.neg(b))
    def scale(self,a,c):return tuple(x*c for x in a)

    def mul(self,a,b):
        out=[r.F(0)]*5
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]+=x*y
        for i in (4,3):
            out[i-2]-=self.A*out[i]
            out[i-3]-=self.B*out[i]
        return tuple(out[:3])

    def square(self,a):return self.mul(a,a)

    def matrix(self,a):
        columns=[self.mul(a,tuple(r.F(i==j) for i in range(3))) for j in range(3)]
        return [[c[i] for c in columns] for i in range(3)]

    def inverse(self,a):
        inv=blocks.exact_inverse(self.matrix(a));out=tuple(row[0] for row in inv)
        assert self.mul(a,out)==self.one
        return out

    def norm(self,a):
        M=self.matrix(a);result=0
        for p in permutations(range(3)):
            sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
            result+=sign*M[0][p[0]]*M[1][p[1]]*M[2][p[2]]
        return result

    def extension_mul(self,z,w,delta):
        s,t=z;v,w=w
        return (self.add(self.mul(s,v),self.mul(delta,self.mul(t,w))),
                self.add(self.mul(s,w),self.mul(t,v)))

    def extension_norm(self,z,delta):
        return self.sub(self.square(z[0]),self.mul(delta,self.square(z[1])))


def encode(x):return [str(c) for c in x]


def build(check=False):
    inp=r.read(lc.INPUT);protocol=r.read(PROTOCOL)
    A,B=map(r.F,inp['anchor']['short_model_ainvariants'][3:])
    K=Cubic(A,B);th=K.theta;th2=K.square(th)
    delta=K.sub(K.scale(th2,-3),K.scalar(4*A))
    derivative=K.add(K.scale(th2,3),K.scalar(A))
    disc=-4*A**3-27*B**2
    assert K.mul(delta,K.square(derivative))==K.scalar(disc)
    points=inp['anchor']['known_points_on_short_model'];witnesses=[];basis=[]
    for i,P in enumerate(points):
        x,y=map(r.F,P);assert y*y==x**3+A*x+B and y
        beta=K.sub(K.scalar(x),th)
        s=K.scale(K.mul(beta,K.add(K.scalar(2*x),th)),1/(2*y))
        t=K.scale(beta,-1/(2*y))
        assert K.norm(beta)==y*y and K.extension_norm((s,t),delta)==beta
        basis.append((beta,(s,t)))
        witnesses.append({'mask':1<<i,'kind':'anchor_basis','beta':encode(beta),'s':encode(s),'t':encode(t)})
    for label,mask in zip(protocol['chain_labels'],protocol['fixed_chain_masks']):
        beta=K.one;z=(K.one,K.scalar(0))
        for i,(b,w) in enumerate(basis):
            if mask>>i&1:beta=K.mul(beta,b);z=K.extension_mul(z,w,delta)
        assert K.extension_norm(z,delta)==beta and z[1]!=K.scalar(0)
        witnesses.append({'mask':mask,'kind':label,'beta':encode(beta),'s':encode(z[0]),'t':encode(z[1])})
    rows=[]
    for u in protocol['parameters']:
        gamma=K.sub(K.one,K.scale(th,u));alpha=K.add(th,K.scale(th2,u))
        D=1+A*u*u+B*u**3
        a=K.add(K.scale(alpha,3),K.scalar(2*A*u))
        factor=K.add(K.add(K.one,K.scale(th,u)),K.scale(K.add(th2,K.scalar(A)),u*u))
        b=K.mul(derivative,factor)
        delta_u=K.sub(K.square(a),K.scale(b,4))
        assert delta_u==K.mul(delta,K.square(gamma))
        assert K.mul(b,gamma)==K.scale(derivative,D)
        assert K.norm(gamma)==D and K.norm(b)==-disc*D*D
        c2=2*A*u;c1=A+3*B*u+A*A*u*u;c0=B+A*B*u*u-B*B*u**3
        assert K.add(K.add(K.mul(K.square(alpha),alpha),K.scale(K.square(alpha),c2)),
                     K.add(K.scale(alpha,c1),K.scalar(c0)))==K.scalar(0)
        assert b==K.add(K.add(K.scale(K.square(alpha),3),K.scale(alpha,2*c2)),K.scalar(c1))
        transported=[]
        for w in witnesses:
            beta,s,t=[tuple(map(r.F,w[k])) for k in ('beta','s','t')]
            t_u=K.mul(t,K.inverse(gamma))
            assert K.extension_norm((s,t_u),delta_u)==beta
            # Explicit K-point on X^2-beta*Y^2=delta_u, no conic solver.
            X=K.mul(K.neg(s),K.inverse(t_u));Y=K.inverse(t_u)
            assert K.sub(K.square(X),K.mul(beta,K.square(Y)))==delta_u
            transported.append(w['mask'])
        coordinate_tests=[]
        for i,P in enumerate(points):
            x=r.F(P[0]);beta,(s,t)=basis[i]
            X=K.mul(gamma,K.add(K.scalar(2*x),th))
            square_gap=K.mul(K.sub(X,a),K.inverse(K.scale(beta,2)))
            candidate_x=K.add(alpha,K.mul(beta,square_gap))
            predicted=K.sub(K.scalar(x),K.scale(K.add(K.add(K.scale(th,x),th2),K.scalar(A)),u))
            assert candidate_x==predicted
            assert (candidate_x[1:]==(0,0))==(u==0)
            if u==0:assert square_gap==K.one
            coordinate_tests.append({'basis_index':i,'candidate_x':encode(candidate_x),
               'rational_x_compatibility':u==0})
        rows.append({'u':u,'alpha':encode(alpha),'gamma':encode(gamma),'a':encode(a),'b':encode(b),
          'isogenous_a':encode(K.scale(a,-2)),'isogenous_b':encode(delta_u),
          'quadratic_discriminant':encode(delta_u),'verified_norm_masks':transported,
          'fixed_norm_point_transport':coordinate_tests})
    ct=[]
    for u in (-1,1):
        base=next(x['W_u_basis'] for x in inp['rows'] if int(x['parameter_u'])==u)
        form=next(x['matrix'] for x in inp['ct'] if x['u']==u)
        e3,f0,e4=[lc.coordinates(m,base) for m in protocol['fixed_chain_masks']]
        ct.append({'u':u,'e3_f0':lc.pairing(e3,f0,form),'e4_f0':lc.pairing(e4,f0,form)})
    assert ct==[{'u':-1,'e3_f0':1,'e4_f0':0},{'u':1,'e3_f0':1,'e4_f0':1}]
    out={'schema':'rank-jump.cubic-bridge.v1','input_sha256':r.digest(lc.INPUT.read_bytes()),
      'protocol_sha256':r.digest(PROTOCOL.read_bytes()),'script_sha256':r.digest(Path(__file__).read_bytes()),
      'constant_quadratic_discriminant':encode(delta),'rational_discriminant':str(disc),
      'norm_witnesses':witnesses,'models':rows,'retained_chain_pairings':ct,
      'summary':{'norm_witness_count':len(witnesses),'models':len(rows),'norm_conic_checks':len(witnesses)*len(rows),
        'rational_x_transport_checks':len(points)*len(rows),'nonzero_parameter_rational_x_transports':0,
        'first_norm_conic_discriminator':'REFUTED_UNIVERSALLY_ON_THE_INHERITED_SPACE'},
      'claim_boundary':'Exact explicit cubic algebra and inherited CT entries. No full CT computation over the cubic field, new rational or cubic point, Selmer complement, or rank claim. The norm conics have points; the isogeny quartics and original rational two-covers are not asserted soluble.'}
    if check:
        if r.read(OUTPUT)!=out:raise ValueError('cubic bridge replay mismatch')
        print('PASS cubic bridge replay')
    else:r.write_new(OUTPUT,out)
    print(out['summary'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('build','check'));a=p.parse_args();build(a.mode=='check')
