#!/usr/bin/env python3
"""Coordinate-gated, checkpointed sparse Keller collision cover through nine terms.

This enumerates only supports passing necessary coordinate gates. It does not
count all raw determinant-balanced supports. Branching at a failed gate adds
every possible repair atom, so no hypothetical Keller collision is lost.
"""
from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
SEVEN=HERE.parent/"sparse-seven"
sys.path.insert(0,str(SEVEN))
import search as old
SCHEMA="coordinate-gated-keller-shard.v1"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def indices(mask):
    while mask:
        bit=mask&-mask
        yield bit.bit_length()-1
        mask^=bit

class Shard(old.Shard):
    def __init__(self,degree,bound,p,q,r,compiler=None):
        if degree>10:
            raise ValueError("The elementary-shear gate imports Moh only through degree ten")
        super().__init__(degree,bound,p,q,r,compiler)
        self.row_masks=tuple(sum(bit for bit,(row,e) in zip(self.bits,self.atoms) if row==i)
                             for i in range(3))
        self.own_masks=tuple(sum(bit for bit,(row,e) in zip(self.bits,self.atoms)
                                if row==i and e[i]>0) for i in range(3))
        self.cross_masks=tuple(sum(bit for bit,(row,e) in zip(self.bits,self.atoms)
                                  if row==i and any(e[j]>0 for j in range(3) if j!=i))
                               for i in range(3))
        self.first_own_masks=tuple(sum(bit for bit,(row,e) in zip(self.bits,self.atoms)
                                      if row==i and e[i]==1) for i in range(3))
        self.cut_masks=tuple(
            (cut,sum(bit for bit,(row,e) in zip(self.bits,self.atoms)
                     if cut&(1<<row) and any(e[j]>0 and not cut&(1<<j) for j in range(3))))
            for cut in range(1,7)
        )

    def requirements(self,mask):
        """Return missing necessary existential support conditions.

        Each entry is (name, eligible allowed repair atoms). All gates are
        necessities only for a completed exact-support Keller collision.
        """
        result=[]
        for i in range(3):
            row=mask&self.row_masks[i]
            if not mask&self.own_masks[i]:
                result.append((f"own-{i}",self.own_masks[i]&self.allowed&~mask))
            if not mask&self.cross_masks[i]:
                result.append((f"cross-{i}",self.cross_masks[i]&self.allowed&~mask))
            if row and not row&~self.own_masks[i] and row&self.first_own_masks[i]:
                result.append((f"normal-{i}",self.row_masks[i]&~self.own_masks[i]&self.allowed&~mask))
        for cut,atoms in self.cut_masks:
            if not mask&atoms:
                result.append((f"cut-{cut}",atoms&self.allowed&~mask))
        return result

    def repair_lower_bound(self,mask):
        # Requirements confined to different output rows need distinct atoms.
        # In a transverse row an axis-count-one repair is disjoint from its
        # own-variable repair, so those two requirements need two atoms.
        needed=0
        for i in range(3):
            row=mask&self.row_masks[i]
            structural=(not mask&self.own_masks[i] or not mask&self.cross_masks[i] or
                        (row and not row&~self.own_masks[i] and row&self.first_own_masks[i]))
            collision=i in (1,2) and (mask&self.axis_masks[i]).bit_count()==1
            if collision:
                needed+=1
                # The axis atom does not contain this transverse own variable.
                needed+=int(not bool(mask&self.own_masks[i]))
            else:
                needed+=int(bool(structural))
        return needed

    def next_nodes(self,mask):
        size=mask.bit_count()
        if size+self.repair_lower_bound(mask)>self.bound:
            return False,(),"repair-lower-bound"
        balanced,children=super().children(mask)
        if balanced:
            requirements=self.requirements(mask)
            if requirements:
                name,atoms=min(requirements,key=lambda item:(item[1].bit_count(),item[0]))
                if size==self.bound:
                    return False,(),name
                children=tuple(mask|self.bits[i] for i in indices(atoms))
                return False,children,name
            return True,children,"retained"
        return False,children,"determinant-or-collision"

    def run_gated(self,path,seconds=120,node_limit=500000,progress_every=30):
        parameters=dict(degree=self.degree,bound=self.bound,seed=list(self.seed))
        hashes=dict(source_sha256=sha(Path(__file__)),predecessor_sha256=sha(SEVEN/"search.py"),
                    legacy_sha256=sha(old.LEGACY/"jcsearch/low_degree_census.py"))
        if path.exists():
            state=json.loads(path.read_text())
            assert state["schema"]==SCHEMA and state["parameters"]==parameters
            assert all(state[k]==v for k,v in hashes.items()),"Source changed: use a new checkpoint directory"
            if state["status"]=="complete":
                return state
            stack=[int(x,16) for x in state["stack"]]
            seen={int(x,16) for x in state["seen"]}
            retained={int(x,16) for x in state["retained_masks"]}
            counts=Counter(state["gate_visits"])
        else:
            state=dict(schema=SCHEMA,parameters=parameters,elapsed_seconds=0.,**hashes)
            stack=[self.seed_mask];seen=set();retained=set();counts=Counter()
        started=last=time.monotonic()
        before=state["elapsed_seconds"]
        initial=len(seen)
        def save():
            state.update(
                status="complete" if not stack else "incomplete",visited=len(seen),pending=len(stack),
                retained_masks=[hex(x) for x in sorted(retained)],
                supports=[self.support(x).to_json() for x in sorted(retained)],
                stack=[hex(x) for x in stack],seen=[hex(x) for x in sorted(seen)] if stack else [],
                gate_visits=dict(counts),elapsed_seconds=before+time.monotonic()-started)
            old.write_json(path,state)
            print(json.dumps(dict(parameters=parameters,status=state["status"],visited=len(seen),
                                  pending=len(stack),retained=len(retained),
                                  elapsed_seconds=state["elapsed_seconds"])),flush=True)
        while stack and time.monotonic()-started<seconds and len(seen)-initial<node_limit:
            mask=stack.pop()
            if mask in seen:
                continue
            seen.add(mask)
            keep,children,gate=self.next_nodes(mask)
            counts[gate]+=1
            if keep:
                retained.add(mask)
            stack.extend(c for c in reversed(children) if c not in seen)
            if time.monotonic()-last>progress_every:
                save();last=time.monotonic()
        save()
        return state

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--degree",type=int,default=7)
    p.add_argument("--bound",type=int,default=9)
    p.add_argument("--seed",default="2,2,3")
    p.add_argument("--seconds",type=int,default=120)
    p.add_argument("--nodes",type=int,default=500000)
    p.add_argument("--output",type=Path,default=HERE/"results/p2-q2-r3.json")
    a=p.parse_args()
    seed=tuple(map(int,a.seed.split(",")))
    shard=Shard(a.degree,a.bound,*seed)
    shard.run_gated(a.output,a.seconds,a.nodes)
if __name__=="__main__":
    main()

