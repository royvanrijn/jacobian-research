#!/usr/bin/env python3
"""Independent NumPy finite-field character sums, with no elliptic point counter."""
import argparse
import base64
import struct
from pathlib import Path
import numpy as np
import retrospective as r
import mixed_character as ex


class Field:
    def __init__(self,p,d,modulus):
        self.p,self.d,self.q=p,d,p**d
        self.c0=modulus[0];self.c1=modulus[1] if d==3 else 0
    def digits(self,v):
        return [(v//self.p**i)%self.p for i in range(self.d)]
    def encode(self,v):
        return sum((x%self.p)*self.p**i for i,x in enumerate(v))
    def add(self,v,w):
        return self.encode([a+b for a,b in zip(self.digits(v),self.digits(w))])
    def scale(self,v,n):
        return self.encode([a*n for a in self.digits(v)])
    def mul(self,v,w):
        a=self.digits(v);b=self.digits(w)
        if self.d==1:return a[0]*b[0]%self.p
        if self.d==2:
            return self.encode([a[0]*b[0]-self.c0*a[1]*b[1],a[0]*b[1]+a[1]*b[0]])
        t3=a[1]*b[2]+a[2]*b[1];t4=a[2]*b[2]
        return self.encode([a[0]*b[0]-self.c0*t3,
                           a[0]*b[1]+a[1]*b[0]-self.c1*t3-self.c0*t4,
                           a[0]*b[2]+a[1]*b[1]+a[2]*b[0]-self.c1*t4])
    def power(self,v,n):
        ans=np.ones_like(v)
        while n:
            if n&1:ans=self.mul(ans,v)
            v=self.mul(v,v);n//=2
        return ans


def verify(case,p):
    raw=next(x for x in r.read(ex.COUNTS)["reductions"] if x["case"]==case and x["p"]==p)
    assert raw["status"]=="COUNTED"
    row=r.read(ex.INPUT)["cases"][case]
    A,B=map(r.F,row["model"][3:]);a,b=[r.F(P[0]) for P in row["generic_points"]]
    A,B,a,b=[r.mod(v,p) for v in (A,B,a,b)];records=[]
    for retained in raw["fields"]:
        d=retained["degree"];F=Field(p,d,retained["modulus_ascending"]);q=F.q
        x=np.arange(q,dtype=np.int64);characters=F.power(x,(q-1)//2)
        assert set(map(int,np.unique(characters)))=={0,1,p-1}
        characters=np.where(characters==p-1,-1,characters)
        xx=F.mul(x,x);xxx=F.mul(xx,x)
        y0=F.add(F.add(xxx,F.scale(x,A)),B)
        y1=F.add(F.scale(xx,2*A),F.scale(x,3*B))
        y2=F.add(F.scale(x,A*A),A*B%p)
        expected=struct.unpack("<"+"h"*q,base64.b64decode(retained["trace_values_i16_le_base64"]))
        weights=characters[F.add(F.add(1,F.scale(x,-a-b)),F.scale(xx,a*b))]
        total=0;traces=[]
        for u in range(q):
            uu=F.mul(u,u);uuu=F.mul(uu,u)
            f=F.add(F.add(F.add(y0,F.mul(u,y1)),F.mul(uu,y2)),F.scale(uuu,-B*B))
            value=int(characters[f].sum())
            assert value==expected[u],(case,p,d,u,value,expected[u])
            traces.append(value);total+=int(weights[u])*value
            if u and u%1024==0:print("checkpoint",case,p,d,u,flush=True)
        assert total==retained["residual_H2_trace"]
        records.append({"degree":d,"base_parameters":q,"residual_H2_trace":total,
                        "trace_sha256":r.digest(struct.pack("<"+"h"*q,*traces))})
        print("PASS field",case,p,d,total,flush=True)
    return {"schema":"rank-jump.mixed-character-verification.v1","status":"PASS","case":case,"p":p,
            "counts_sha256":r.digest(ex.COUNTS.read_bytes()),"verifier_sha256":r.digest(Path(__file__).read_bytes()),
            "numpy":np.__version__,"fields":records}


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--case",type=int,required=True)
    p.add_argument("--prime",type=int,required=True);p.add_argument("--destination",type=Path)
    args=p.parse_args();result=verify(args.case,args.prime)
    if args.destination:r.write_new(args.destination,result)
