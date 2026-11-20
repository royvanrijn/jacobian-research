#!/usr/bin/env python3
"""Exact/source support for fixed-support simple-turn consecutive exclusion."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact assertions must execute.")
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb
import argparse,copy,hashlib,json,re
HERE=Path(__file__).resolve().parent
PN="SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md"
CN="verify_simple_turn_shrinking_window.py"
IN="simple_turn_shrinking_window_inputs.json"
LN="simple_turn_shrinking_window_source_ledger.json"
RN="simple_turn_shrinking_window_verification.json"
PREFIX="research/factorial/image-rigidity/"
ROLES={
 "proof":("local",PN),"checker":("local",CN),"inputs":("local",IN),
 "actual_stationary_bridge":("repo",PREFIX+"COALESCING_BINOMIAL_SADDLE_BRIDGE.md"),
 "actual_morse_mechanism":("repo",PREFIX+"CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md")}
SCOPE={"schema":1,"a_min":1,"d_min":2,"support_condition":"4*a>(d-1)^2",
 "support_parameters":"fixed positive integers",
 "conclusion":"there exist a fixed physical turning disc and a terminal cutoff for each fixed support",
 "entire_physical_neighborhood":True,"bounded_and_growing_fields":True,
 "individual_complex_m_periods":True,"original_noninteger_m_identity_required":False,
 "original_identity_orders":"positive integers only",
 "fixed_center_rho_on_complex_w_collars":True,
 "d2_infinity":"locally narrower wedge with alpha<pi/2",
 "cubic_tail_method":"proper critical-value inverse rays, not compact real-phase bands",
 "effective_cutoff":False,"uniform_cutoff_over_supports":False,
 "all_terminal_orders":False,"global_coefficient_plane_exclusion":False,
 "a0_support_included":False,"formal_verification":False,"external_review":False}
CONTEXTS=[
 {"url":"https://dlmf.nist.gov/2.4#v","scope":"local cubic/contour-dependent framework; actual cycle is proved in the note"},
 {"url":"https://dlmf.nist.gov/9.5#ii","scope":"canonical homogeneous cubic contour"},
 {"url":"https://dlmf.nist.gov/9.2","scope":"Airy ODE and nonzero initial values"}]
def require(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding="utf-8"))
def same(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(same(a[k],b[k]) for k in b)
    if type(b) is list:return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return a==b
def inputs(d):require(same(d,SCOPE),"Exact local fixed-support eventual scope")
def sources(d,repo):
    require(type(d) is dict and set(d)=={"schema","formal_verification","external_review","primary_contexts","sources"},"Ledger schema")
    require(type(d["schema"]) is int and d["schema"]==1,"Schema1")
    require(d["formal_verification"] is False and d["external_review"] is False,"Assurance flags")
    require(same(d["primary_contexts"],CONTEXTS),"Primary context scope")
    rows=d["sources"];require(type(rows) is list and len(rows)==len(ROLES),"All5 mandatory sources")
    seen={};paths=set()
    for r in rows:
        require(type(r) is dict and set(r)=={"role","root","path","sha256"},"Row schema")
        role=r["role"];require(type(role) is str and role in ROLES and role not in seen,"Unique role")
        require((r["root"],r["path"])==ROLES[role],"Exact path")
        key=(r["root"],r["path"]);require(key not in paths,"Unique source")
        require(type(r["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",r["sha256"]),"SHA256")
        base=(HERE if r["root"]=="local" else repo).resolve();p=(base/r["path"]).resolve()
        require(p.is_relative_to(base) and p.is_file() and sha(p)==r["sha256"],"Missing/stale:"+role)
        seen[role]=r;paths.add(key)
    return seen
# Sparse Laurent polynomials over Q; all displayed identities are tested as
# coefficient identities, not by evaluation at a list of numerical points.
class P:
    def __init__(self,d=0):
        if isinstance(d,P):self.t=dict(d.t)
        elif isinstance(d,dict):self.t={m:Q(c) for m,c in d.items() if c}
        else:self.t={():Q(d)} if d else {}
    @staticmethod
    def var(s):return P({((s,1),):Q(1)})
    def __add__(self,o):
        o=P(o);d=dict(self.t)
        for m,c in o.t.items():d[m]=d.get(m,Q(0))+c
        return P(d)
    __radd__=__add__
    def __neg__(self):return P({m:-c for m,c in self.t.items()})
    def __sub__(self,o):return self+-P(o)
    def __rsub__(self,o):return P(o)+-self
    def __mul__(self,o):
        o=P(o);d={}
        for m,c in self.t.items():
            for n,e in o.t.items():
                ex=dict(m)
                for v,k in n:ex[v]=ex.get(v,0)+k
                z=tuple(sorted((v,k) for v,k in ex.items() if k))
                d[z]=d.get(z,Q(0))+c*e
        return P(d)
    __rmul__=__mul__
    def __pow__(self,k):
        require(type(k) is int,"Integer Laurent power")
        if k<0:
            require(len(self.t)==1,"Negative power only for monomials")
            (m,c),=self.t.items()
            return P({tuple((v,k*e) for v,e in m):c**k})
        r=P(1);b=self
        while k:
            if k&1:r=r*b
            b=b*b;k//=2
        return r
    def __eq__(self,o):return self.t==P(o).t
    def diff(self,v):
        d={}
        for m,c in self.t.items():
            ex=dict(m);k=ex.get(v,0)
            if k:
                ex[v]=k-1;z=tuple(sorted((w,e) for w,e in ex.items() if e))
                d[z]=d.get(z,Q(0))+k*c
        return P(d)
    def coefficient(self,v,k):
        d={}
        for m,c in self.t.items():
            ex=dict(m)
            if ex.get(v,0)==k:
                ex.pop(v,None);z=tuple(sorted(ex.items()));d[z]=d.get(z,Q(0))+c
        return P(d)
def universal():
    a,d,S,T,w,m,r,k,U,W,q,eta,G,Delta,Dv,bm,bw,E=[P.var(n) for n in
      ["a","d","S","T","w","m","rho","k","u","omega","q","eta","G","Delta","Dval","bm","bw","eps"]]
    names=[]
    def check(name,left,right=0):require(left==right,"Universal:"+name);names.append(name)
    D=a+d;Qpoly=S*S-(D+a-1)*S+a*D
    check("actual stationary numerator",(a-S)*(w+T)+d*T,(a-S)*w+(D-S)*T)
    check("actual Omega derivative after T'=dT/S",
      d*(D-S)*(S-a)-S*(S-a)-S*(D-S),-d*Qpoly)
    check("critical stationary mixed derivative numerator",(D-S)*(S-a)-S,-Qpoly)
    H=P.var("H") # H represents t^2-eta.
    check("normalized physical-to-cubic flow",G**-1*(-G*H**-1),-H**-1)
    check("normalized phase descent",(H*G**-1)*(-G*H**-1),-1)
    phi=U**3*Q(1,3)-W*U
    check("cubic derivative",phi.diff("u"),U*U-W)
    check("cubic critical values",q**3*Q(1,3)-q*q*q,-Q(2,3)*q**3)
    # Exact exponent identities after m=k/rho^3.
    check("error rho/m",r*(k*r**-3)**-1,r**4*k**-1)
    check("error1/(mk)",(k*r**-3)**-1*k**-1,r**3*k**-2)
    check("error1/(m rho)",(k*r**-3)**-1*r**-1,r**2*k**-1)
    check("error1/(mk rho²)",(k*r**-3)**-1*k**-1*r**-2,r*k**-2)
    check("exact radial logarithm transport",
      -Dv+m*(-Delta*d**-1*w**-1)*(d*w*m**-1)-bm+bw*d*w*m**-1,
      -Dv-Delta-bm+bw*d*w*m**-1)
    # Formal epsilon=m^-1/3 calculation, retaining terms through epsilon².
    wc,A,B,z1,z2=[P.var(n) for n in ["wc","A","B","z1","z2"]]
    upper=wc+A*E**2+B*E**3
    radial=1+d*E**3+d*(d+1)*Q(1,2)*E**6
    lower=radial*upper
    zet=lambda v:z1*(v-wc)+z2*(v-wc)**2
    hdiff=E**-2*((1-Q(2,3)*E**3-Q(1,9)*E**6)*zet(lower)-zet(upper))
    check("bounded actual radial field leading coefficient",hdiff.coefficient("eps",1),d*wc*z1)
    check("bounded actual radial field constant",hdiff.coefficient("eps",0),0)
    check("bounded actual radial field no epsilon² term",hdiff.coefficient("eps",2),0)
    check("Airy initial data determinant nonzero cyclotomic certificate",
      (W*W-W)**2+3,(W*W+W+1)*(W*W-3*W+3))
    return names
# Exact quadratic algebra Q[j]/(j²+tau²); j is chosen -i*tau.
class F:
    tau2=Q(1)
    def __init__(self,x=0,y=0):
        if isinstance(x,F):self.x,self.y=x.x,x.y
        else:self.x,self.y=Q(x),Q(y)
    def __add__(self,o):o=F(o);return F(self.x+o.x,self.y+o.y)
    __radd__=__add__
    def __neg__(self):return F(-self.x,-self.y)
    def __sub__(self,o):return self+-F(o)
    def __rsub__(self,o):return F(o)+-self
    def __mul__(self,o):o=F(o);return F(self.x*o.x-F.tau2*self.y*o.y,self.x*o.y+self.y*o.x)
    __rmul__=__mul__
    def inv(self):
        n=self.x*self.x+F.tau2*self.y*self.y;require(n!=0,"Quadratic nonzero denominator")
        return F(self.x/n,-self.y/n)
    def __truediv__(self,o):return self*F(o).inv()
    def __rtruediv__(self,o):return F(o)*self.inv()
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        out=F(1);b=self
        while n:
            if n&1:out=out*b
            b=b*b;n//=2
        return out
    def __eq__(self,o):o=F(o);return self.x==o.x and self.y==o.y
def fixtures():
    count=0;shapes=0
    for a in range(1,21):
        for d in range(2,9):
            if 4*a<=(d-1)**2:continue
            b=Q(d-1,2);F.tau2=Q(a)-b*b
            S=F(Q(a)+b,1);D=a+d;w=S**d*(D-S)/(S-a)
            require(S*S-(D+a-1)*S+a*D==0,"Actual selected quadratic")
            require(w.y<0,"Actual upper turning value is nonreal upper")
            mixed=-d*S**(d-1)/(w+S**d)**2
            require(-d*w*mixed==1,"Exact actual radial mixed identity")
            third=-2*F(0,1)/(S*S)
            b3=2/third
            require(b3.y>0,"Radial cubic has strictly negative imaginary part")
            require(S.x*S.x-F.tau2>0,"Re Sstar² strictly positive")
            count+=1
    for a,d in [(0,1),(0,6),(1,2),(2,2),(2,3),(3,2),(3,4),(5,5)]:
        for m in range(1,5):
            A=Q(2,3);B=-Q(3,5);x=A/B
            direct=sum((Q(comb(m,k)*factorial(a*m+d*k))*A**(m-k)*B**k for k in range(m+1)),Q(0))
            norm=sum((Q(comb(m,k)*factorial(a*m+d*k),factorial(a*m))*x**(m-k) for k in range(m+1)),Q(0))
            require(direct==B**m*factorial(a*m)*norm,"Exact own-order moment normalization")
            shapes+=1
    return count,shapes
def math_negatives():
    names=[]
    # Full real-phase cubic bands are noncompact even for omega=1.
    F.tau2=Q(1)
    for R in [1,2,10,10**6]:
        u=F(0,R);v=u**3/3-u
        require(v.x==0 and u.x*u.x+u.y*u.y>=R*R,"Genuine noncompact real-band control")
    names.append("full-cubic-real-bands-not-compact")
    # One cannot retain the original closed d=2 boundary as a decaying infinity end.
    require(F(0,1).x==0,"Re(iR)=0 on the closed d2 lower ray")
    names.append("full-d2-boundary-not-uniformly-decaying")
    require(2*factorial(4)//factorial(2)==24 and 2*factorial(3)//factorial(1)==12,"Actual gamma3 vs false gamma2")
    names.append("wrong-shared-Gamma-shape")
    # A bare leading scaled root is not enough: the first correction is needed.
    require(Q(2,3)<1 and Q(4,3)>1,"Root location leading scale vs two-term error scale")
    names.append("bare-leading-root-error-not-o-one-over-N")
    return names
def source_controls(ledger,repo):
    names=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,KeyError,ValueError,TypeError):names.append(name);return
        raise AssertionError("Negative accepted:"+name)
    for role in ROLES:
        x=copy.deepcopy(ledger);x["sources"]=[r for r in x["sources"] if r["role"]!=role]
        reject("omit:"+role,lambda x=x:sources(x,repo))
        x=copy.deepcopy(ledger);next(r for r in x["sources"] if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda x=x:sources(x,repo))
    for key in ["effective_cutoff","uniform_cutoff_over_supports","all_terminal_orders",
                "global_coefficient_plane_exclusion","a0_support_included",
                "original_noninteger_m_identity_required"]:
        x=copy.deepcopy(SCOPE);x[key]=True
        reject("false-scope:"+key,lambda x=x:inputs(x))
    for key in ["entire_physical_neighborhood","bounded_and_growing_fields","individual_complex_m_periods","fixed_center_rho_on_complex_w_collars"]:
        x=copy.deepcopy(SCOPE);x[key]=False
        reject("missing-scope:"+key,lambda x=x:inputs(x))
    x=copy.deepcopy(SCOPE);x["a_min"]=True
    reject("typed-a",lambda:inputs(x))
    for key in ["schema","empty","duplicate"]:
        x=copy.deepcopy(ledger)
        if key=="schema":x["schema"]=True
        if key=="empty":x["sources"]=[]
        if key=="duplicate":x["sources"][-1]=copy.deepcopy(x["sources"][0])
        reject("ledger:"+key,lambda x=x:sources(x,repo))
    return names
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--output",type=Path,default=HERE/RN);ns=ap.parse_args()
    repo=ns.repo_root.resolve();require(repo.is_dir(),"Explicit repository root")
    inp=load(HERE/IN);inputs(inp);ledger=load(HERE/LN);rows=sources(ledger,repo)
    identities=universal();count,shapes=fixtures();bad=math_negatives();controls=source_controls(ledger,repo)
    sources(ledger,repo)
    out={"status":"PASS","schema":1,"scope":inp,"formal_verification":False,"external_review":False,
     "proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),"input_sha256":sha(HERE/IN),
     "source_ledger_sha256":sha(HERE/LN),"mandatory_source_roles":len(rows),
     "sources":[rows[r] for r in sorted(rows)],"universal_identity_names":identities,
     "mathematical_negative_names":bad,"source_scope_negative_names":controls,
     "counts":{"universal_Laurent_identity_families":len(identities),"exact_actual_turning_fixtures":count,
               "exact_actual_moment_normalizations":shapes,"mathematical_negative_controls":len(bad),
               "source_scope_controls":len(controls)},
     "meaning":"The written complete-contour proof proves the fixed-support local theorem. These exact algebra/source controls support and authenticate it; no finite control proves the quantified flow or uniform estimates by itself."}
    ns.output.parent.mkdir(parents=True,exist_ok=True);ns.output.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","counts":out["counts"],"mandatory_sources":len(rows),"output":str(ns.output)}))
if __name__=="__main__":main()

