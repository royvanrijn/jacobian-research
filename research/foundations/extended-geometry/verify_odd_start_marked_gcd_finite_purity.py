#!/usr/bin/env python3
"""Standard-library exact sparse replay of marked gcd and prefix controls."""
import argparse,json,sys,hashlib
from fractions import Fraction as Q
from math import comb
from pathlib import Path
if not __debug__:raise SystemExit("Run without -O: exact assertions are required.")

def trim(p):
    p = tuple(Q(x) for x in p)
    while p and p[-1] == 0:
        p = p[:-1]
    return p


Z = ()
ONE = (Q(1),)
XX = (Q(0), Q(1))


def add(p, q):
    return trim(tuple((p[i] if i < len(p) else 0) +
                      (q[i] if i < len(q) else 0)
                      for i in range(max(len(p), len(q)))))


def scale(p, c):
    return trim(tuple(c*x for x in p))


def sub(p, q):
    return add(p, scale(q, -1))


def mul(p, q):
    if not p or not q:
        return Z
    out = [Q(0)]*(len(p)+len(q)-1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            out[i+j] += pi*qj
    return trim(out)


def power(p, n):
    assert n >= 0
    out = ONE
    while n:
        if n & 1:
            out = mul(out, p)
        p = mul(p, p)
        n //= 2
    return out


def deriv(p):
    return trim(tuple(i*p[i] for i in range(1, len(p))))


def divmodp(p, q):
    assert q
    rem = list(p)
    quotient = [Q(0)]*max(0, len(p)-len(q)+1)
    while rem and len(rem) >= len(q):
        shift = len(rem)-len(q)
        c = rem[-1]/q[-1]
        quotient[shift] = c
        for j, qj in enumerate(q):
            rem[shift+j] -= c*qj
        rem = list(trim(rem))
    return trim(quotient), trim(rem)


def exactdiv(p, q):
    quotient, rem = divmodp(p, q)
    assert not rem
    return quotient


def monic(p):
    assert p
    return scale(p, 1/p[-1])


def gcdp(p, q):
    while q:
        p, q = q, divmodp(p, q)[1]
    return monic(p) if p else Z


def degree(p):
    return len(p)-1


def radical(p):
    return monic(exactdiv(p, gcdp(p, deriv(p))))


def max_multiplicity(p):
    count = 0
    while degree(p) > 0:
        count += 1
        p = gcdp(p, deriv(p))
    return count


class RF:
    def __init__(self, num=ONE, den=ONE):
        num, den = trim(num), trim(den)
        assert den
        g = gcdp(num, den)
        num, den = exactdiv(num, g), exactdiv(den, g)
        self.num, self.den = scale(num, 1/den[-1]), scale(den, 1/den[-1])

    def __mul__(self, other):
        return RF(mul(self.num, other.num), mul(self.den, other.den))

    def __truediv__(self, other):
        assert other.num
        return RF(mul(self.num, other.den), mul(self.den, other.num))

    def __add__(self, other):
        return RF(add(mul(self.num, other.den), mul(other.num, self.den)),
                  mul(self.den, other.den))

    def __pow__(self, n):
        if n < 0:
            return RF(power(self.den, -n), power(self.num, -n))
        return RF(power(self.num, n), power(self.den, n))

    def derivative(self):
        return RF(sub(mul(deriv(self.num), self.den),
                      mul(self.num, deriv(self.den))), power(self.den, 2))

    def at_zero(self):
        assert self.den[0] != 0
        return self.num[0]/self.den[0] if self.num else Q(0)

    def polynomial(self):
        return exactdiv(self.num, self.den)

    def __eq__(self, other):
        return self.num == other.num and self.den == other.den



def spclean(p):return {e:c for e,c in p.items()if c}
def spadd(p,q,sign=1):
 out=dict(p)
 for e,c in q.items():out[e]=out.get(e,Q(0))+sign*c
 return spclean(out)
def spscale(p,c):return spclean({e:v*c for e,v in p.items()})
def spmul(p,q):
 out={}
 for e,c in p.items():
  for f,d in q.items():
   k=(e[0]+f[0],e[1]+f[1]);out[k]=out.get(k,Q(0))+c*d
 return spclean(out)
def sppower(p,n):
 out={(0,0):Q(1)}
 while n:
  if n&1:out=spmul(out,p)
  p=spmul(p,p);n//=2
 return out
def xpoly(p):return {(i,0):v for i,v in enumerate(p)if v}
def ypoly(p):return {(0,i):v for i,v in enumerate(p)if v}
def spdeg(p,i):return max((e[i]for e in p),default=-1)
def decode(terms):
 out={}
 for e,c in terms:
  assert len(e)==2 and all(isinstance(v,int)and v>=0 for v in e)
  k=tuple(e);assert k not in out
  out[k]=Q(c)
 return spclean(out)
def spdiag(p):
 out=[Q(0)]*(max((sum(e)for e in p),default=-1)+1)
 for e,c in p.items():out[sum(e)]+=c
 return trim(out)
def divide_diagonal(F):
 # Monic exact division by Y-x, entirely in Q[x,Y].
 rem=dict(F);out={}
 while rem and spdeg(rem,1)>=1:
  m=spdeg(rem,1);lead={(e[0],m-1):c for e,c in rem.items()if e[1]==m}
  out=spadd(out,lead);rem=spadd(rem,spmul(lead,{(0,1):Q(1),(1,0):Q(-1)}),-1)
 assert not rem
 return out
def primitive(p):
 coeffs={}
 for (i,j),v in p.items():
  a=coeffs.setdefault(j,[])
  while len(a)<=i:a.append(Q(0))
  a[i]=v
 content=Z
 for a in coeffs.values():content=gcdp(content,trim(a))
 return degree(content)==0
def sp_order_series(p,j,L):
 out=Z
 for (i,k),v in p.items():
  term=scale(cap_power(j,k,L),v)
  term=(Q(0),)*i+term
  out=add(out,term[:L+1])
 return trim(out[:L+1])
def marked_branch(P,L):
 c=P.get((0,1),Q(0));assert c and P.get((0,0),0)==0
 j=Z
 for k in range(1,L+1):
  j=tuple(j)+(Q(0),)*(k+1-len(j))
  current=sp_order_series(P,j,k)
  v=-(current[k]if len(current)>k else Q(0))/c
  j=list(j);j[k]=v;j=trim(j)
 assert len(j)>1 and j[0]==0 and j[1]==-1
 assert not sp_order_series(P,j,L)
 return j
def order(p):return next((i for i,c in enumerate(p)if c),None)
def cap_mul(p,q,L):return trim(mul(p,q)[:L+1])
def cap_power(p,k,L):
 out=ONE
 while k:
  if k&1:out=cap_mul(out,p,L)
  p=cap_mul(p,p,L);k//=2
 return out
def raw(C,D,a,m):
 if a*m%2:return Q(0)
 out=Q(0)
 for k in range(m+1):
  coeff=cap_mul(cap_power(C,m-k,m),cap_power(D,k,m),m)
  out+=Q(comb(m,k),a*m+2*k+1)*(coeff[m]if len(coeff)>m else Q(0))
 return out
def source_functions(source):
 if source["kind"]=="phase":
  a=source["a"];assert a>=1 and a%2
  C=trim(source["C"]);D=trim(source["D"]);assert C[0]==1 and D[0]==-1
  f=(RF(XX)**2)*(RF(scale(D,-1))**a)/(RF(C)**(a+2))
  W=(RF(scale(D,-1))/RF(C))**((a-1)//2)/RF(C)
  H=RF((Q(2),))*W/f.derivative()
  return f,H,(a,C,D)
 return RF(trim(source["f_num"]),trim(source["f_den"])),RF(trim(source["H_num"]),trim(source["H_den"])),None
def fixture(entry,L=12,maxm=8):
 f,H,phase=source_functions(entry["source"])
 A=trim(entry["f_A"]);B=trim(entry["f_B"]);M=trim(entry["H_M"]);N=trim(entry["H_N"])
 assert RF(A,B)==f and RF(M,N)==H
 assert gcdp(A,B)==ONE and gcdp(M,N)==ONE
 assert B[0]!=0 and A[0]==A[1]==0 and A[2]/B[0]==1
 assert N[0]==0 and N[1]!=0 and M[0]/N[1]==1
 n=max(degree(A),degree(B));h=max(degree(M),degree(N));D0=(n-1)*(n+4*h)
 assert n==entry["n"]and h==entry["h"]and D0==entry["D0"]
 Ax,Bx,Mx,Nx=map(xpoly,[A,B,M,N]);Ay,By,My,Ny=map(ypoly,[A,B,M,N])
 P=divide_diagonal(spadd(spmul(Ay,Bx),spmul(Ax,By),-1))
 cross=spadd(spmul(Mx,Ny),spmul(My,Nx),-1)
 QQ=spadd(spmul(Ax,sppower(cross,2)),spscale(spmul(spmul(Bx,sppower(Nx,2)),sppower(Ny,2)),4),-1)
 assert QQ and spdiag(QQ)==scale(mul(B,power(N,4)),-4)
 assert P.get((0,1))==A[2]*B[0]!=0 and P.get((0,0),0)==0
 assert spdeg(P,0)<=n-1 and spdeg(P,1)<=n-1
 assert spdeg(QQ,0)<=n+2*h and spdeg(QQ,1)<=2*h
 G,P1,Q1,BA,BB,R=[decode(entry[k])for k in ["G","P1","Q1","bezout_A","bezout_B","R"]]
 assert G and primitive(G)
 assert spmul(G,P1)==P and spmul(G,Q1)==QQ
 assert R and all(e[1]==0 for e in R)
 assert spadd(spmul(BA,P1),spmul(BB,Q1))==R
 assert spdeg(R,0)<=D0
 # This polynomial identity makes the residuals coprime over Q(x).
 pure=G.get((0,0),Q(0))==0
 assert pure==entry["expected_pure"]and str(G.get((0,0),Q(0)))==entry["G00"]
 j=marked_branch(P,L)
 branch=sp_order_series(QQ,j,L);qorder=order(branch)
 Gseries=sp_order_series(G,j,L)
 if pure:
  assert not branch and not Gseries
 else:
  assert Gseries and Gseries[0]!=0 and qorder is not None
  assert qorder<=spdeg(R,0)<=D0
 raw_values=[];first=None
 if phase:
  a,C,D=phase
  for m in range(1,maxm+1):
   value=raw(C,D,a,m);raw_values.append(str(value))
   if value and first is None:first=m
  if pure:assert first is None
  else:assert first is not None and first%2==0 and first+4==qorder and first<=D0-4
  if entry["id"]=="nonpure_second_zero_a1":assert first==4
 nonsquare=False
 if "nonsquare_X"in entry["id"]:
  assert phase[1]==ONE and phase[2]==(Q(-1),Q(-1))
  nonsquare=True # -D/C=1+x has a simple geometric zero.
 T1=(f*(H**2)==RF(ONE))
 if entry["id"]in ["pure_trivial_a1","pure_sharp_a1","pure_two_root_a3"]:assert T1 and QQ
 return {"id":entry["id"],"pure":pure,"G_Y_degree":spdeg(G,1),"G00":entry["G00"],
         "n":n,"h":h,"prefix_bound":D0-4,"coarse_prefix":D0,
         "eliminant_degree":spdeg(R,0),"literal_Bezout_terms":len(BA)+len(BB),
         "Q_marked_branch_order":qorder,"first_nonzero_moment":first,
         "raw_moments":raw_values,"nonrational_X_control":nonsquare,
         "T_equals_one_control":T1,"Q_is_nonzero":True}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 parser=argparse.ArgumentParser()
 parser.add_argument("--repo-root",type=Path,default=Path(__file__).resolve().parents[3])
 parser.add_argument("--output",type=Path,default=Path(__file__).with_name("odd_start_marked_gcd_verification.json"))
 args=parser.parse_args();here=Path(__file__).resolve().parent
 cp=here/"odd_start_marked_gcd_certificates.json";lp=here/"odd_start_marked_gcd_source_ledger.json"
 data=json.loads(cp.read_text());ledger=json.loads(lp.read_text())
 pins=[]
 for src in ledger["sources"]:
  p=(args.repo_root if src.get("location","repository")=="repository"else here)/src["path"]
  assert digest(p)==src["sha256"];pins.append({"path":src["path"],"sha256":digest(p),"role":src["role"]})
 results=[fixture(row)for row in data["fixtures"]]
 off=next(r for r in results if r["id"]=="off_base_common_component")
 assert not off["pure"]and off["G_Y_degree"]>0 and Q(off["G00"])!=0
 # Wrong unprimitive multiplication by x would fabricate a base zero.
 G=decode(next(r for r in data["fixtures"]if r["id"]=="off_base_common_component")["G"])
 wrong=spmul({(1,0):Q(1)},G)
 assert wrong.get((0,0),0)==0 and not primitive(wrong)
 # A constant-in-Y second polynomial is a separate resultant algebra control.
 P0={(0,2):Q(1),(0,1):Q(1),(1,0):Q(1)}
 Q0={(3,0):Q(1)};R0=sppower(Q0,2)
 assert R0=={(6,0):Q(1)}
 assert spmul(Q0,Q0)==R0 and spdeg(P0,1)*spdeg(Q0,0)==spdeg(R0,0)
 # Corrupt literal eliminant must fail.
 mutant=dict(data["fixtures"][-1]);mutant["R"]=list(mutant["R"])+[[[0,0],"1"]]
 try:fixture(mutant)
 except AssertionError:pass
 else:raise AssertionError("Corrupt eliminant accepted")
 names=["ODD_START_MARKED_GCD_AND_FINITE_PURITY.md",Path(__file__).name,cp.name,lp.name]
 out={"status":"PASS","scope":"Exact gcd factors and polynomial Bezout eliminants certify finite inputs; the written all-field/all-order proof gives the universal marked criterion and prefix bound.",
      "fixtures":results,"raw_moment_controls":sum(len(r["raw_moments"])for r in results),
      "literal_Bezout_terms":sum(r["literal_Bezout_terms"]for r in results),
      "negative_controls":{"off_base_nonconstant_gcd_rejected":True,"unprimitive_x_multiple_rejected":True,
                           "corrupt_eliminant_rejected":True,"constant_Y_resultant_case":True},
      "discovery_engine":data["discovery_engine"],
      "runtime":{"python":sys.version.split()[0],"standard_library_only":True},
      "source_pins":pins,"packet_sha256":{n:digest(here/n)for n in names},
      "assurance":"No runtime SymPy, numerical root finder, finite-prefix assumption, formal verification or external review."}
 args.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
 print(json.dumps({"status":"PASS","fixtures":len(results),"raw_moments":out["raw_moment_controls"],"literal_Bezout_terms":out["literal_Bezout_terms"],"output":str(args.output)},sort_keys=True))
if __name__=="__main__":main()

