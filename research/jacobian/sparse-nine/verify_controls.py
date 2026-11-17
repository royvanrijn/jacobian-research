#!/usr/bin/env python3
"""Exact positive witness, one-parameter frame minimum, and search regressions."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import tempfile
import sys
import sympy as sp
import gated_search as gated

HERE=Path(__file__).resolve().parent
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    x,y,z,w=sp.symbols("x y z w")
    vars=(x,y,z)
    u=1+x*y
    F=sp.Matrix([
        u**3*z+y**2*u*(4+3*x*y),
        y+3*x*u**2*z+3*x*y**2*(4+3*x*y),
        2*x-3*x*x*y-x**3*z])
    assert sp.expand(F.jacobian(vars).det())==-2
    lam=sp.symbols("lambda",nonzero=True)
    torus=dict(zip(vars,(lam*x,y/lam,z/lam**2)))
    weighted=sp.diag(lam**-2,lam**-1,lam)*F
    assert (F.subs(torus,simultaneous=True)-weighted).applyfunc(sp.expand)==sp.zeros(3,1)
    a=sp.Matrix([0,4*w+6,-63*w*w-189*w-142])
    b=sp.Matrix([1,w,2-3*w])
    suba=dict(zip(vars,a));subb=dict(zip(vars,b))
    assert (F.subs(suba)-F.subs(subb)).applyfunc(sp.expand)==sp.zeros(3,1)
    L=sp.Matrix.hstack(b-a,sp.Matrix([0,1,0]),sp.Matrix([0,0,1]))
    A=F.jacobian(vars).subs(suba)
    assert sp.expand(L.det())==1 and sp.expand(A.det())==-2
    sub=dict(zip(vars,a+L*sp.Matrix(vars)))
    H=((A*L).inv()*(F.subs(sub,simultaneous=True)-F.subs(suba))).applyfunc(sp.expand)
    assert H.subs({x:0,y:0,z:0})==sp.zeros(3,1)
    assert H.subs({x:1,y:0,z:0}).applyfunc(sp.expand)==sp.zeros(3,1)
    assert H.jacobian(vars).subs({x:0,y:0,z:0}).applyfunc(sp.expand)==sp.eye(3)
    coefficients=[]
    coefficient_rows=[]
    for i,entry in enumerate(H):
        for exponent,c in sp.Poly(entry,*vars).terms():
            if sum(exponent)>=2:
                p=sp.Poly(c,w,domain=sp.QQ)
                coefficients.append(p)
                coefficient_rows.append(dict(row=i,exponent=list(exponent),coefficient=str(p.as_expr())))
    assert len(coefficients)==67
    factors={}
    for p in coefficients:
        prefactor,parts=sp.factor_list(p)
        reconstructed=sp.Poly(prefactor,w,domain=sp.QQ)
        for f,multiplicity in parts:
            reconstructed*=f**multiplicity
            f=f.monic()
            factors[str(f.as_expr())]=f
        assert reconstructed==p
    strata=[]
    for f in factors.values():
        assert f.is_irreducible
        zeros=[i for i,p in enumerate(coefficients) if p.rem(f).is_zero]
        strata.append(dict(factor=str(f.as_expr()),degree=f.degree(),
                           support=67-len(zeros),zero_coefficient_indices=zeros))
    strata.sort(key=lambda row:(row["support"],row["factor"]))
    assert len(strata)==30 and strata[0]["support"]==49 and strata[0]["factor"]=="w + 2"
    assert strata[1]["support"]==54
    K=H.subs(w,-2).applyfunc(sp.expand)
    assert sp.expand(K.jacobian(vars).det())==1
    assert K.subs({x:0,y:0,z:0})==K.subs({x:1,y:0,z:0})==sp.zeros(3,1)
    assert K.jacobian(vars).subs({x:0,y:0,z:0})==sp.eye(3)
    rows=tuple(tuple(sorted((e for e,c in sp.Poly(h,*vars).terms() if sum(e)>=2),
                           key=lambda e:(sum(e),e))) for h in K)
    support=gated.old.Support(rows)
    assert [len(row) for row in rows]==[5,18,26] and support.nonlinear_size==49
    assert [sp.Poly(h,*vars).total_degree() for h in K]==[4,6,7]
    witness=dict(
        schema="normalized-keller-49.v1",coordinate_degrees=[4,6,7],
        nonlinear_counts=[5,18,26],nonlinear_total=49,
        source_points=[[0,0,0],[1,0,0]],common_target=[0,0,0],
        parameter=-2,old_source_points=[[0,-2,-16],[1,-2,8]],
        normalization_matrix_L=[[1,0,0],[0,1,0],[24,0,1]],
        coefficients=[[
            dict(exponent=list(e),coefficient=str(c)) for e,c in sp.Poly(h,*vars).terms()
        ] for h in K],
        formula=[str(h) for h in K],
        support=support.to_json())
    gated.old.write_json(HERE/"normalized-49-witness.json",witness)
    gated.old.write_json(HERE/"frame-family-certificate.json",dict(
        status="pass",source_sha256=sha(Path(__file__)),
        scope="All foundational collision frames based on x=0 with transverse columns e_y,e_z, up to the verified diagonal torus equivalence. This is not all affine frames.",
        generic_support=67,minimum=49,unique_minimizing_parameter=-2,
        coefficient_rows=coefficient_rows,strata=strata))

    # Positive control: the actual noninjective degree-seven map must survive
    # every new gate, and a branch chain inside its support must reach it.
    axis=support.collision_axis_counts()
    assert axis[0]>=1 and axis[1]>=2
    p=min(e[0] for e in rows[0] if e[1:]==(0,0))
    q,r=sorted(e[0] for e in rows[1] if e[1:]==(0,0))[:2]
    compiler=gated.old.SupportEnumerator(7,49)
    shard=gated.Shard(7,49,p,q,r,compiler)
    target=sum(shard.bits[shard.compiler.atom_index[atom]] for atom in support.atoms)
    assert shard.requirements(target)==[] and shard.repair_lower_bound(target)==0
    keep,children,gate=shard.next_nodes(target)
    assert keep and not children
    mask=shard.seed_mask
    trace=[]
    while mask!=target:
        keep,children,gate=shard.next_nodes(mask)
        choices=[child for child in children if child&target==child]
        assert choices,(mask.bit_count(),gate)
        child=max(choices,key=lambda x:(x.bit_count(),-x))
        assert child.bit_count()>mask.bit_count()
        trace.append(dict(size=mask.bit_count(),gate=gate,next_size=child.bit_count()))
        mask=child
    assert len(trace)<=46

    # Whole small-degree legacy comparison, using the old unseeded closure.
    legacy=gated.old.SupportEnumerator(3,8)
    full,statistics=legacy.enumerate()
    expected=set()
    for size,items in full.items():
        for s in items:
            if not s.collision_axis_unimodularity_possible():continue
            counts=s.collision_axis_counts()
            if counts[1]<2:s=s.swapped_23()
            ps=min(e[0] for e in s.rows[0] if e[1:]==(0,0))
            qs,rs=sorted(e[0] for e in s.rows[1] if e[1:]==(0,0))[:2]
            trial=gated.Shard(3,8,ps,qs,rs,legacy)
            m=sum(trial.bits[trial.compiler.atom_index[atom]] for atom in s.atoms)
            if not trial.requirements(m):expected.add(s.canonical_under_collision_stabilizer())
    actual=set()
    with tempfile.TemporaryDirectory() as folder:
        for seed in ((2,2,3),(3,2,3)):
            trial=gated.Shard(3,8,*seed,legacy)
            path=Path(folder)/("%s.json"%(seed[0]))
            partial=trial.run_gated(path,60,1)
            assert partial["status"]=="incomplete" and partial["pending"]>0
            state=trial.run_gated(path,60,500000)
            assert state["status"]=="complete"
            actual.update(gated.old.read_support(s).canonical_under_collision_stabilizer()
                          for s in state["supports"])
    assert actual==expected
    receipt=dict(
        status="pass",source_sha256=sha(Path(__file__)),
        witness_sha256=sha(HERE/"normalized-49-witness.json"),
        frame_family_sha256=sha(HERE/"frame-family-certificate.json"),
        witness_nonlinear_support=49,witness_degrees=[4,6,7],
        witness_all_gates_pass=True,positive_branch_trace=trace,torus_equivariance=True,
        family_minimum=49,generic_family_support=67,exceptional_factors=30,
        small_degree=3,small_bound=8,legacy_balanced=sum(map(len,full.values())),
        gated_candidates=len(actual),old_new_set_equality=True,checkpoint_resume=True,
        boundary="The sharp minimum 49 covers natural boundary-based frames (one source point on x=0, transverse columns e_y,e_z), using diagonal equivariance. It is not a global support minimum.")
    gated.old.write_json(HERE/"controls-validation.json",receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k!="positive_branch_trace"}),flush=True)
if __name__=="__main__":main()

