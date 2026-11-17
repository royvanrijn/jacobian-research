#!/usr/bin/env python3
"""Finite quotient certificates for arbitrary-degree raw-output target corrections."""
from pathlib import Path
from itertools import product
import hashlib,json,time
import sympy as s
from graph36 import solve_or_separate
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def multiply(a,b,order_ideal):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(q+r for q,r in zip(e,f))
            if g in order_ideal:out[g]=out.get(g,s.Integer(0))+c*d
    return {e:c for e,c in out.items() if c}
def algebra_columns(generators,T):
    maxdegree=max(map(sum,T));zero=(0,)*12;columns=[];words=[];seen=set();visited=0
    def visit(p,word,last):
        nonlocal visited
        if word:
            visited+=1;high={e:c for e,c in p.items() if sum(e)>3}
            key=tuple(sorted((e,str(c)) for e,c in high.items()))
            if high and key not in seen:seen.add(key);columns.append(high);words.append(word)
        if len(word)==maxdegree:return
        for j in range(last,len(generators)):
            child=multiply(p,generators[j],T)
            if child:visit(child,word+(j,),j)
    visit({zero:s.Integer(1)},(),0)
    return columns,words,visited
def main():
    source=HERE/"cubic12-integer-witness.json";d=json.loads(source.read_text());X=s.symbols(d["variables"])
    K=[s.expand(sum(s.Rational(t["coefficient"])*s.prod(v**n for v,n in zip(X,t["exponent"])) for t in row)) for row in d["coefficients"]]
    folder=HERE/"all-degree36";folder.mkdir(exist_ok=True);cases=[]
    for j in range(3,12):
        start=time.monotonic()
        raw_path=HERE/"graph36"/("pivot-%02d.json"%(j+1));raw=json.loads(raw_path.read_text())
        assert raw["source_sha256"]==sha(HERE/"graph36.py") and raw["input_sha256"]==sha(source)
        old=next(t for t in raw["tests"] if t["status"]=="separated");i=old["component"]-1
        T=set()
        for row in old["functional"]:T.update(product(*(range(e+1) for e in row["exponent"])))
        graph=s.expand(X[j]-K[j]);kept={a:s.expand(f.subs(X[j],graph)) for a,f in enumerate(K) if a!=j}
        allowed=[a for a in kept if a!=i]
        generators=[]
        for a in allowed:
            polynomial={e:c for e,c in s.Poly(kept[a],*X).terms() if e in T}
            assert not polynomial.get((0,)*12,0);generators.append(polynomial)
        columns,words,visited=algebra_columns(generators,T)
        target={e:c for e,c in s.Poly(kept[i],*X).terms() if e in T and sum(e)>3}
        result,nrows=solve_or_separate(columns,target)
        row=dict(pivot=j+1,component=i+1,allowed_outputs=[a+1 for a in allowed],
            quotient_monomials=[list(e) for e in sorted(T)],maximum_jet=max(map(sum,T)),
            columns=len(columns),nonzero_words_visited=visited,coefficient_rows=nrows,**result)
        if result["status"]=="span-member":
            solution=[]
            correction=s.Integer(0)
            for word,c in zip(words,result["solution"]):
                if not s.Rational(c):continue
                ix=[allowed[q] for q in word]
                correction+=s.Rational(c)*s.prod(kept[q] for q in ix)
                solution.append(dict(indices=[q+1 for q in ix],coefficient=c))
            residual=s.expand(kept[i]-correction)
            row["solution"]=solution;row["residual_degree"]=s.Poly(residual,*X).total_degree()
            row["residual_high_terms"]=[dict(exponent=list(e),coefficient=str(c)) for e,c in s.Poly(residual,*X).terms() if sum(e)>3]
            row["status"]="candidate-cubic-completion" if row["residual_degree"]<=3 else "finite-quotient-inconclusive"
        else:row["status"]="arbitrary-degree-obstructed"
        row.update(source_sha256=sha(Path(__file__)),linear_solver_sha256=sha(HERE/"graph36.py"),input_sha256=sha(source),raw_certificate_sha256=sha(raw_path),seconds=round(time.monotonic()-start,3))
        (folder/("pivot-%02d.json"%(j+1))).write_text(json.dumps(row,indent=2,sort_keys=True)+"\n")
        cases.append(dict(pivot=j+1,status=row["status"],file="pivot-%02d.json"%(j+1)))
        checkpoint=dict(status="running",source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases,
            scope="Functional support is finite. All target monomials through its maximum source degree are tested in a downward-closed monomial quotient; higher-degree target monomials vanish there because generators have no constants.")
        (folder/"checkpoint.json").write_text(json.dumps(checkpoint,indent=2,sort_keys=True)+"\n")
        print(json.dumps({key:row[key] for key in ["pivot","component","status","maximum_jet","columns","nonzero_words_visited","seconds"]}),flush=True)
    checkpoint["status"]="complete";(folder/"checkpoint.json").write_text(json.dumps(checkpoint,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

