#!/usr/bin/env python3
"""Replay exact cubic cube-zero shuffle certificates; standard library only."""
if not __debug__:
    raise SystemExit("Refusing optimized Python: exact certificate checks are required.")
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction
from functools import cache
from hashlib import sha256
from heapq import heapify, heappop, heappush
from importlib.util import module_from_spec, spec_from_file_location
from math import comb, prod, isqrt
from pathlib import Path
import gzip,json,sys,time

P=1000000009
HERE=Path(__file__).resolve().parent
def require(test,msg):
    if not test:raise ValueError(msg)
def fp(c):
    c=Fraction(c)
    require(c.denominator%P!=0,"coefficient denominator vanishes modulo proof prime")
    return c.numerator%P*pow(c.denominator%P,-1,P)%P
def nf_mod(row,pivots):
    pending={i:int(c)%P for i,c in row.items() if int(c)%P}
    heap=[-i for i in pending];heapify(heap);out={}
    while heap:
        i=-heappop(heap)
        if i not in pending:continue
        a=pending.pop(i)
        if i not in pivots:out[i]=a;continue
        for j,b in pivots[i].items():
            if j==i:continue
            prior=pending.get(j,0);value=(prior-a*b)%P
            if value:
                if not prior:heappush(heap,-j)
                pending[j]=value
            else:pending.pop(j,None)
    return out
def insert_mod(row,pivots):
    reduced=nf_mod(row,pivots)
    if not reduced:return False
    i=max(reduced);inv=pow(reduced[i],-1,P)
    pivots[i]={j:c*inv%P for j,c in reduced.items()}
    return True
def reduce_rref(row,pivots):
    out=dict(row)
    for i in tuple(row):
        if i not in pivots:continue
        a=out.pop(i)
        for j,b in pivots[i].items():
            if j==i:continue
            value=out.get(j,0)-a*b
            if value:out[j]=value
            else:out.pop(j,None)
    return out
def validate_rref(pp,ids):
    universe=set(ids);piv=set(pp);free=universe-piv
    require(piv<=universe,"pivot outside actual prequotient")
    for i,row in pp.items():
        require(row.get(i)==1,"pivot coefficient is not one")
        require(all(c and (j==i or (j in free and j<i)) for j,c in row.items()),"nonreduced or zero pivot term")
    return free
def make_children(q,k):
    children={}
    for cc in q.forests(k-1,3):
        i=len(q.grade);q.grade.append(k);q.tree_ids[cc]=i;children[i]=cc
    return children
def update_count(q,k,mod):
    count=sum(prod(comb(q.unordered_counts[size]+amount-1,amount)
                   for size,amount in Counter(ss).items())
              for ss in mod.partitions(k-1,3))
    q.unordered_counts.append(count)
def verify_exact(q,cert,mod):
    k=cert["grade"];require(k==len(q.basis),"exact grades not consecutive")
    children=make_children(q,k);ids=list(children)
    require(ids==cert["ids"],"prequotient reconstruction mismatch")
    pp={int(i):{int(j):Fraction(c) for j,c in row.items()}
        for i,row in cert["pivot_rows"].items()}
    free=validate_rref(pp,ids)
    # Every rational coefficient used by the field prefix must specialize.
    pp_mod={i:{j:fp(c) for j,c in row.items()} for i,row in pp.items()}
    selections=cert["selected_rows"]
    require(selections==sorted(set(selections)),"repeated or unsorted selected row index")
    selected=set(selections);bp={};ordinal=0
    for bg in range(k-2):
        for side in q.forests(k-3-bg,6):
            for bottom in q.basis[bg]:
                row=q.root_shuffles(side,bottom)
                require(not reduce_rref(row,pp),"actual root relation not in certified rational span")
                if ordinal in selected:
                    block={i:fp(c) for i,c in row.items() if i in pp and fp(c)}
                    require(insert_mod(block,bp),"selected pivot block is singular modulo proof prime")
                ordinal+=1
    require(all(i<ordinal for i in selected),"selected row does not exist")
    require(len(bp)==len(pp)==len(selections),"selected pivot block has wrong rank")
    q.pivots.append(pp);q.basis.append(tuple(i for i in ids if i in free))
    coefficient={}
    for tree,cc in children.items():
        a=mod.multinomial(cc)*prod(q.inverse[q.grade[c]].get(c,0) for c in cc)
        if a:coefficient[tree]=a
    remainder=reduce_rref(coefficient,pp);q.inverse.append(remainder);update_count(q,k,mod)
    require(len(remainder)==cert["expected_inverse_support"],"inverse support mismatch")
    return {"grade":k,"columns":len(ids),"root_relations":ordinal,
            "rank":len(pp),"dimension":len(free),"inverse_support":len(remainder)},pp_mod

def validate_substitution(q,k,side,bottom):
    require(len(side)==6 and tuple(side)==tuple(sorted(side)),"malformed off-path multiset")
    all_ids=(*side,bottom)
    require(all(isinstance(i,int) and 0<=i<len(q.grade) for i in all_ids),"bad substitution identifier")
    require(all(i in q.basis[q.grade[i]] for i in all_ids),"substitution is not in a previous quotient basis")
    require(q.grade[bottom]+sum(q.grade[i] for i in side)+3==k,"root relation has wrong total grade")

def controls():
    require(P>720 and all(P%d for d in range(2,isqrt(P)+1)),"proof field modulus is not prime")
    pp={};require(insert_mod({9:1,4:3},pp),"control insertion failed")
    require(not insert_mod({9:1,4:3},pp),"duplicate relation was accepted as independent")
    rejected=0
    for bad in [{7:{7:2}},{7:{7:1,4:1},4:{4:1}},{8:{8:1,99:1}},{7:{7:1,8:1}}]:
        try:validate_rref(bad,[4,7,8])
        except ValueError:rejected+=1
        else:raise ValueError("malformed rational pivot certificate accepted")
    # The propagation interval is quadratic, not the older arity-dependent one.
    require(all(i>=15 or 29-i>=15 for i in range(1,29)),"gap endpoint arithmetic")
    return {"duplicate_relation_rejected":True,"malformed_pivot_controls":rejected,
            "quadratic_gap_first_inductive_grade":29}
def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path)
    ap.add_argument("--ledger",type=Path,default=HERE/"cubic_cube_zero_shuffle_source_ledger.json")
    ap.add_argument("--input",type=Path,default=HERE/"cubic_cube_zero_shuffle_certificates.json.gz")
    ap.add_argument("--output",type=Path)
    args=ap.parse_args();start=time.time()
    if args.repo_root is None:
        roots=[Path.cwd(),*Path.cwd().parents,HERE,*HERE.parents]
        args.repo_root=next((p for p in roots if (p/"research/jacobian/subtree-shuffling/verify_subtree_shuffle.py").is_file()),None)
        require(args.repo_root is not None,"cannot locate repository; use --repo-root")
    ledger=json.loads(args.ledger.read_text())
    for entry in ledger["inputs"]:
        base=args.repo_root if entry["scope"]=="repository" else HERE
        path=base/entry["path"]
        require(path.is_file(),"missing pinned source: "+entry["path"])
        require(sha256(path.read_bytes()).hexdigest()==entry["sha256"],"changed pinned source: "+entry["path"])
    input_entry=next(e for e in ledger["inputs"] if e["role"]=="certificate")
    require(sha256(args.input.read_bytes()).hexdigest()==input_entry["sha256"],"selected certificate hash mismatch")
    source=args.repo_root/"research/jacobian/subtree-shuffling/verify_subtree_shuffle.py"
    spec=spec_from_file_location("retained_shuffle_certificate_kernel",source)
    mod=module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    packet=json.loads(gzip.decompress(args.input.read_bytes()))
    require(packet["schema_version"]==1 and packet["arity"]==3 and packet["path_length"]==3,"wrong certificate problem")
    require([c["grade"] for c in packet["exact_grades"]]==list(range(1,19)),"incomplete exact initial interval")
    require([c["grade"] for c in packet["modular_zero_grades"]]==list(range(19,29)),"incomplete modular zero interval")
    q=mod.ShuffleQuotient(3,3);exact=[];mod_pivots=[{}]
    for cert in packet["exact_grades"]:
        row,pp=verify_exact(q,cert,mod);exact.append(row);mod_pivots.append(pp)
        print("EXACT",json.dumps(row),flush=True)
    require([r["dimension"] for r in exact[13:]]==[35,16,13,9,2],"wrong final exact dimensions")
    require(all(q.inverse[k] for k in range(1,15)),"missing lower nonzero inverse coefficient")
    require(all(not q.inverse[k] for k in range(15,19)),"missing exact vanishing coefficient")
    require(sum(len(b) for b in q.basis[:15])==195,"wrong sharp lower-example dimension")
    dual_index=min(q.inverse[14]);dual_value=q.inverse[14][dual_index]
    require(dual_value!=0,"grade14 lower witness vanished")
    class ModularQuotient(mod.ShuffleQuotient):
        @cache
        def lower_shuffle(self,side,bottom):
            if not side:return {bottom:1}
            grade=self.grade[bottom]+sum(self.grade[i] for i in side)+len(side)//2
            require(grade<len(self.basis),"using an unproved lower grade")
            if grade>=19:
                require(not self.basis[grade],"zero shortcut used for a nonzero lower grade")
                return {}
            return nf_mod(self.root_shuffles(side,bottom),self.pivots[grade])
    mq=ModularQuotient(3,3);mq.grade=list(q.grade);mq.tree_ids=dict(q.tree_ids)
    mq.basis=list(q.basis);mq.pivots=mod_pivots
    mq.inverse=[{i:fp(c) for i,c in rr.items()} for rr in q.inverse]
    mq.unordered_counts=list(q.unordered_counts)
    zeros=[]
    for cert in packet["modular_zero_grades"]:
        k=cert["grade"];require(cert["prime"]==P,"wrong zero-piece field")
        require(k==len(mq.basis),"zero-piece grades not consecutive")
        cc=make_children(mq,k);pp={}
        require(len(cert["relations"])==len(cc),"wrong number of selected zero-piece relations")
        for side,bottom in cert["relations"]:
            side=tuple(side);validate_substitution(mq,k,side,bottom)
            require(insert_mod(mq.root_shuffles(side,bottom),pp),"dependent selected zero-piece relation")
        require(len(pp)==len(cc),"zero-piece rank is not full")
        mq.pivots.append({i:{i:1} for i in pp});mq.basis.append(());mq.inverse.append({})
        update_count(mq,k,mod)
        row={"grade":k,"columns":len(cc),"independent_valid_relations":len(pp),"dimension":0}
        zeros.append(row);print("ZERO",json.dumps(row),flush=True)
    claim=packet["conclusion"]
    require(claim=={"first_zero_inverse_grade":15,"quadratic_gap_last_grade":28,
                   "sharp_inverse_degree":29,"lower_witness_truncation_dimension":195},"conclusion metadata mismatch")
    result={"status":"PASS","arithmetic":"QQ and F_1000000009; standard library only",
            "exact_grades":exact,"modular_zero_grades":zeros,
            "grade14_dual_coordinate":dual_index,"grade14_nonzero_dual_value":str(dual_value),
            "sharp_inverse_degree":29,"sharp_example_dimension":195,"controls":controls(),
            "input_sha256":sha256(args.input.read_bytes()).hexdigest(),
            "kernel_sha256":sha256(source.read_bytes()).hexdigest(),
            "checker_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
            "source_ledger_sha256":sha256(args.ledger.read_bytes()).hexdigest(),
            "source_pins_checked":len(ledger["inputs"]),
            "seconds":time.time()-start}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":"PASS","sharp_inverse_degree":29,"seconds":result["seconds"]}))
if __name__=="__main__":main()

