#!/usr/bin/env python3
"""Adaptive finite-jet test for arbitrary-degree elementary target completion."""
from pathlib import Path
from itertools import product
import hashlib,json,time
import sympy as s
from graph36 import solve_or_separate
from graph36_all_degree import algebra_columns
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def closure(exponents):
    out=set()
    for e in exponents:out.update(product(*(range(q+1) for q in e)))
    return out
def times(a,b,cap):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(q+r for q,r in zip(e,f))
            if sum(g)<=cap:out[g]=out.get(g,s.Integer(0))+c*d
    return {e:c for e,c in out.items() if c}
def trial_residual(polys,allowed,words,solution,target,cap):
    cache={():{(0,)*12:s.Integer(1)}}
    def word_poly(word):
        if word not in cache:cache[word]=times(word_poly(word[:-1]),polys[allowed[word[-1]]],cap)
        return cache[word]
    out={e:c for e,c in target.items() if 3<sum(e)<=cap}
    terms=[]
    for word,coefficient in zip(words,solution):
        c=s.Rational(coefficient)
        if not c:continue
        terms.append(dict(indices=[allowed[j]+1 for j in word],coefficient=str(c)))
        for e,d in word_poly(word).items():
            if sum(e)>3:out[e]=out.get(e,s.Integer(0))-c*d
    return {e:c for e,c in out.items() if c},terms
def main():
    source=HERE/"cubic12-integer-witness.json";d=json.loads(source.read_text());X=s.symbols(d["variables"])
    K=[s.expand(sum(s.Rational(t["coefficient"])*s.prod(v**n for v,n in zip(X,t["exponent"])) for t in row)) for row in d["coefficients"]]
    folder=HERE/"adaptive36-wide";folder.mkdir(exist_ok=True);cases=[]
    for j in range(3,12):
        start=time.monotonic();raw_path=HERE/"graph36"/("pivot-%02d.json"%(j+1));raw=json.loads(raw_path.read_text())
        assert raw["source_sha256"]==sha(HERE/"graph36.py") and raw["input_sha256"]==sha(source)
        old=next(t for t in raw["tests"] if t["status"]=="separated");i=old["component"]-1
        seeds=[row["exponent"] for row in old["functional"]]
        jet_path=HERE/"jet36"/("pivot-%02d.json"%(j+1));jet=json.loads(jet_path.read_text())
        for test in jet["tests"]:
            if test["component"]==i+1 and test["status"]=="separated":seeds += [row["exponent"] for row in test["functional"]]
        T=closure(seeds)
        graph=s.expand(X[j]-K[j]);kept={a:s.expand(f.subs(X[j],graph)) for a,f in enumerate(K) if a!=j}
        polys={a:dict(s.Poly(f,*X).terms()) for a,f in kept.items()};allowed=[a for a in kept if a!=i]
        history=[];row=None
        for iteration in range(100):
            D=max(map(sum,T))
            if len(T)>5000 or D>7:
                row=dict(status="support-cap-inconclusive");break
            generators=[{e:c for e,c in polys[a].items() if e in T} for a in allowed]
            columns,words,visited=algebra_columns(generators,T)
            target={e:c for e,c in polys[i].items() if e in T and sum(e)>3}
            result,nrows=solve_or_separate(columns,target)
            stage=dict(iteration=iteration,quotient_dimension=len(T),maximum_jet=D,
                columns=len(columns),visited_words=visited,coefficient_rows=nrows,status=result["status"],
                rank=result["rank"],augmented_rank=result["augmented_rank"])
            history.append(stage)
            print(json.dumps(dict(pivot=j+1,component=i+1,**stage)),flush=True)
            if result["status"]=="separated":
                row=dict(status="arbitrary-degree-obstructed",functional=result["functional"],
                    certificate_maximum_jet=max(sum(t["exponent"]) for t in result["functional"]),
                    certificate_quotient_dimension=len(closure([t["exponent"] for t in result["functional"]])))
                break
            residual,terms=trial_residual(polys,allowed,words,result["solution"],polys[i],D+1)
            assert not any(e in T for e in residual)
            if not residual:
                # A positive jet check is never called a full solution.
                row=dict(status="higher-jet-needed",trial_correction=terms,checked_through_degree=D+1);break
            minimum=min(map(sum,residual))
            choices=[e for e in residual if sum(e)==minimum]
            choices.sort(key=lambda e:(not bool(e[i]),len(closure([e])-T),e))
            new=choices[:32];T.update(closure(new));stage["added_monomials"]=[list(e) for e in new]
            partial=dict(status="running",pivot=j+1,component=i+1,history=history,source_sha256=sha(Path(__file__)),input_sha256=sha(source))
            (folder/("pivot-%02d.json"%(j+1))).write_text(json.dumps(partial,indent=2,sort_keys=True)+"\n")
        if row is None:row=dict(status="iteration-cap-inconclusive")
        row.update(pivot=j+1,component=i+1,allowed_outputs=[a+1 for a in allowed],history=history,
           source_sha256=sha(Path(__file__)),algebra_sha256=sha(HERE/"graph36_all_degree.py"),
           linear_solver_sha256=sha(HERE/"graph36.py"),input_sha256=sha(source),seconds=round(time.monotonic()-start,3))
        (folder/("pivot-%02d.json"%(j+1))).write_text(json.dumps(row,indent=2,sort_keys=True)+"\n")
        cases.append(dict(pivot=j+1,status=row["status"],file="pivot-%02d.json"%(j+1)))
        checkpoint=dict(status="running",source_sha256=sha(Path(__file__)),input_sha256=sha(source),cases=cases,
            scope="Exact finite-jet dual certificates exclude one elementary raw-output target correction of any polynomial degree when status is obstructed. Other statuses are inconclusive; no general dimension theorem.")
        (folder/"checkpoint.json").write_text(json.dumps(checkpoint,indent=2,sort_keys=True)+"\n")
    checkpoint["status"]="complete";(folder/"checkpoint.json").write_text(json.dumps(checkpoint,indent=2,sort_keys=True)+"\n")
if __name__=="__main__":main()

