#!/usr/bin/env python3
"""Verify the two explicit native factorial pair certificates. Read-only."""
import argparse,hashlib,json,math
from pathlib import Path
def need(x,message):
    if not x: raise ValueError(message)
def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
def verify(row):
    q=row["q"];m=q*(4*q+3)**2;D=-3*q*(q+1);p=row["prime"]
    need(q in (1,2),"only the two stated cases")
    need(row["middle_m"]==m and row["D"]==D,"native parameter mismatch")
    need(prime(p) and p>2*m,"prime or full-leading-degree guard")
    need(pow(D%p,(p-1)//2,p)==p-1,"quadratic residue field is not a field")
    zero=(0,0);one=(1,0)
    def add(a,b):return((a[0]+b[0])%p,(a[1]+b[1])%p)
    def mul(a,b):return((a[0]*b[0]+D*a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
    def trim(a):
        while a and a[-1]==zero:a.pop()
        return a
    def pa(a,b):
        out=[zero]*max(len(a),len(b))
        for i,x in enumerate(a):out[i]=add(out[i],x)
        for i,x in enumerate(b):out[i]=add(out[i],x)
        return trim(out)
    def pm(a,b):
        if not a or not b:return []
        out=[zero]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):out[i+j]=add(out[i+j],mul(x,y))
        return trim(out)
    F=[1]
    for j in range(1,2*m+1):F.append(F[-1]*j%p)
    I=[pow(x,-1,p) for x in F]
    z=((-8*q*(q+1)*(2*q+1))%p,1);zp=[one]
    for j in range(m):zp.append(mul(zp[-1],z))
    def H(j):
        # Coefficients of Lambda((r*x^2+x+z)^j), literal multinomial.
        out=[]
        for a in range(j+1):
            v=zero
            for d in range(j-a+1):
                c=F[j]*I[a]*F[j+a-d]*I[j-a-d]*I[d]%p
                v=add(v,(c*zp[d][0]%p,c*zp[d][1]%p))
            out.append(v)
        need(len(out)==j+1 and out[-1]==(F[2*j],0),"native leading coefficient")
        return out
    def T(j):return pa(H(j),[(j*x%p,j*y%p) for x,y in H(j-1)])
    A=T(m-1);B=T(m)
    need(len(A)==m and len(B)==m+1,"native degrees dropped")
    def read_poly(key):
        a=row[key];need(isinstance(a,list),"coefficient vector")
        need(all(isinstance(c,list) and len(c)==2 and all(isinstance(x,int) and 0<=x<p for x in c) for c in a),"field coefficient format")
        return trim([tuple(c) for c in a])
    U=read_poly("Bezout_U_coefficients_ascending")
    V=read_poly("Bezout_V_coefficients_ascending")
    need(pa(pm(U,A),pm(V,B))==[one],"native Bezout product is not one")
    return {"q":q,"middle_m":m,"prime":p,"native_degrees":[len(A)-1,len(B)-1],"native_Bezout_product":"1","pass":True}
def main():
    ap=argparse.ArgumentParser();ap.add_argument("certificate",type=Path)
    a=ap.parse_args();raw=a.certificate.read_bytes();d=json.loads(raw)
    rows=d["cases"];need([r["q"] for r in rows]==[1,2],"exact two-case set")
    result=[verify(r) for r in rows]
    print(json.dumps({"status":"PASS","certificate_sha256":hashlib.sha256(raw).hexdigest(),"cases":result,"scope":"Only q1,q2; no all-q inference."},sort_keys=True))
if __name__=="__main__":main()
