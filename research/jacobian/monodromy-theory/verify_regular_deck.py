"""Boundary controls for the universal regular-deck rigidity proof."""
from pathlib import Path
import hashlib,json,platform,time
import sympy as S
ROOT=Path(__file__).resolve().parent
started=time.monotonic();checks=[]
def check(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
x=S.symbols("x")
for p in [2,3,5,7]:
 F=x**p-x
 check("positive-characteristic Keller derivative p="+str(p),S.Poly(S.diff(F,x)+1,x,modulus=p).is_zero)
 for a in range(p):
  check(f"positive-characteristic regular translation p={p},a={a}",S.Poly(F.subs(x,x+a)-F,x,modulus=p).is_zero)
 iterate=x
 for i in range(p):iterate=iterate+1
 check("positive-characteristic finite order p="+str(p),S.Poly(iterate-x,x,modulus=p).is_zero)
for n in [2,3,4,5]:
 # Use a formal root of unity in Q[zeta]/Phi_n.
 zeta=S.symbols("zeta");minimum=S.cyclotomic_poly(n,zeta)
 check("ramified affine regular symmetry n="+str(n),S.rem((zeta*x)**n-x**n,minimum,zeta)==0)
 check("ramified affine fixed origin n="+str(n),(zeta*x).subs(x,0)==0)
 check("affine power not Keller n="+str(n),S.degree(S.diff(x**n,x),x)>0)
 check("punctured source allows etale deck group n="+str(n),S.diff(x**n,x)/x**(n-1)==n)
manifest=json.loads((ROOT/"sources/regular-deck-manifest.json").read_text())
for src in manifest["sources"]:
 check("Euler primary-source pin",hashlib.sha256((ROOT/src["path"]).read_bytes()).hexdigest()==src["sha256"])
files=["GENERAL_DECK_RIGIDITY.md","verify_regular_deck.py","sources/regular-deck-manifest.json"]
rec=dict(result="pass",check_count=len(checks),checks=checks,python=platform.python_version(),sympy=S.__version__,
 runtime_seconds=round(time.monotonic()-started,2),sha256={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in files},
 assurance="Written all-dimensional endomorphism, free-action and Euler proof. Exact hypothesis-boundary controls are not formal verification or independent whole-proof review.",
 monodromy_spectrum_complete=False,goal_complete=False)
(ROOT/"regular_deck_receipt.json").write_text(json.dumps(rec,indent=2)+"\n")
print(json.dumps(dict(result="pass",checks=len(checks),seconds=rec["runtime_seconds"])))
