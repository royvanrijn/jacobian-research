#!/usr/bin/env sage -python
"""Export physical effective-zero markings for the corrected A3+2A2 -> 5A1 q8 edge."""

import hashlib, json
from pathlib import Path
from sage.all import ZZ, block_diagonal_matrix, matrix, pari, vector

ROOT = Path(__file__).resolve().parents[2]
GENDIR = ROOT / "artifacts/generated-results"
SOURCE = GENDIR / "elkies-k3-h3-q4o208-corrected-a3-2a2-old_a11_component_10-marking.json"
CERT = GENDIR / "elkies-k3-h3-q4o208-corrected-a3-2a2-to-5a1-physical-q8-c10-certificate.json"
OUTPUT = GENDIR / "elkies-k3-h3-corrected-a3-2a2-q8-5a1-effective-zero-markings.json"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

def load(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.startswith("#")])
def ent(v): return [int(x) for x in vector(ZZ, v)]
def rows(m): return [[int(x) for x in r] for r in m.rows()]
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

source=json.loads(SOURCE.read_text()); cert=json.loads(CERT.read_text())
assert source["status"] == "PASS_EXACT_CORRECTED_A3_2A2_PHYSICAL_EFFECTIVE_ZERO_MARKING"
assert cert["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
frame=load(ROOT/source["frame_output"]); gram=block_diagonal_matrix(U2,-frame)
fibre=vector(ZZ,cert["source_to_child_basis"][0]); explicit={n:vector(ZZ,v) for n,v in source["equation_explicit_curves_in_child"].items()}
zeros=sorted(n for n,c in explicit.items() if c*gram*fibre==1); assert len(zeros)==4
out={}
for zero_name in zeros:
    zero=explicit[zero_name]; u=matrix(ZZ,[list(fibre),list(zero+fibre)])
    ker=matrix(ZZ,[list(fibre*gram),list((zero+fibre)*gram)]).right_kernel_matrix(); pre=u.stack(ker)
    assert abs(pre.det())==1; invpre=pre.inverse().change_ring(ZZ); comp=-(pre*gram*pre.transpose())[2:,2:]
    half=[vector(ZZ,c) for c in matrix(ZZ,pari(comp).qfminim(2)[2]).columns()]; roots=half+[-r for r in half]
    assert len(roots)==10
    known=[]
    for curve in {tuple(c):c for c in explicit.values()}.values():
        if curve*gram*fibre: continue
        incidence=ZZ(zero*gram*curve); assert incidence in (0,1)
        pos=curve if incidence==0 else fibre-curve
        tail=vector(ZZ,(pos*invpre)[2:])
        if tuple(tail) not in {tuple(x) for x in known}: known.append(tail)
    chosen=[]; seen=set(); unique_explicit={tuple(c):c for c in explicit.values()}.values()
    for root in roots:
        if tuple(root) in seen: continue
        seen.add(tuple(root)); seen.add(tuple(-root))
        pair=(root,-root); pick=[r for r in pair if tuple(r) in {tuple(k) for k in known}]
        if not pick:
            old_fibre=vector(ZZ,[1,0]+[0]*17)
            by_old_degree=[r for r in pair if old_fibre*gram*(vector(ZZ,[0,0]+list(r))*pre)>0]
            if len(by_old_degree)==1:
                pick=by_old_degree
                chosen.append(pick[0])
                continue
            pinned=vector(ZZ,source["target_fibres_in_root_adapted_hub"]["pinned_R17"])
            by_pinned_degree=[r for r in pair if pinned*gram*(vector(ZZ,[0,0]+list(r))*pre)>0]
            if len(by_pinned_degree)==1:
                pick=by_pinned_degree
                chosen.append(pick[0])
                continue
            viable=[]
            for r in pair:
                curve=vector(ZZ,[0,0]+list(r))*pre
                if all(curve*gram*c>=0 for c in unique_explicit): viable.append(r)
            if len(viable) != 1:
                print("AMBIGUOUS_ROOT", zero_name, len(viable), ent(root), [
                    [(name, int((vector(ZZ,[0,0]+list(r))*pre)*gram*c))
                     for name,c in explicit.items() if (vector(ZZ,[0,0]+list(r))*pre)*gram*c < 0]
                    for r in pair
                ], flush=True)
            assert len(viable)==1; pick=viable
        assert len(pick)==1; chosen.append(pick[0])
    assert len(chosen)==5
    simple_source=[vector(ZZ,[0,0]+list(r))*pre for r in chosen]
    partial=u.stack(matrix(ZZ,[list(-c) for c in simple_source])); smith,left,right=partial.smith_form()
    assert smith[:,:7]==matrix.identity(ZZ,7) and not any(smith[:,7:].list())
    basis=block_diagonal_matrix(left.inverse().change_ring(ZZ),matrix.identity(ZZ,12))*right.inverse().change_ring(ZZ)
    for i in range(7,19):
        row=basis.row(i); row-=ZZ(row*gram*(zero+fibre))*fibre; row-=ZZ(row*gram*fibre)*(zero+fibre); basis[i]=row
    assert abs(basis.det())==1; inverse=basis.inverse().change_ring(ZZ); childgram=basis*gram*basis.transpose(); child=-childgram[2:,2:]
    assert child[:5,:5]==2*matrix.identity(ZZ,5)
    slug=zero_name.lower(); fp=GENDIR/f"elkies-k3-h3-physical-q8-5a1-{slug}-frame.txt"; mp=GENDIR/f"elkies-k3-h3-physical-q8-5a1-{slug}-marking.json"
    fp.write_text(f"# physical q8 5A1; zero={zero_name}\n"+"\n".join(" ".join(map(str,r)) for r in child.rows())+"\n")
    eq=matrix(ZZ,source["equation_A11_to_root_adapted_hub_basis"])
    payload={"schema":"elkies-k3.h3-physical-q8-5a1-effective-zero-marking.v1","status":"PASS_EXACT_PHYSICAL_Q8_5A1_EFFECTIVE_ZERO_MARKING","hub":f"physical_q8_5A1_{zero_name}_zero","zero":zero_name,"root_data":[5,10,32],"ade":"5A1","frame_output":str(fp.relative_to(ROOT)),"frame_sha256":sha(fp),"basis_in_source":rows(basis),"source_in_basis":rows(inverse),"equation_A11_to_root_adapted_hub_basis":rows(basis*eq),"target_fibres_in_root_adapted_hub":{n:ent(vector(ZZ,v)*inverse) for n,v in source["target_fibres_in_root_adapted_hub"].items()},"equation_explicit_curves_in_child":{n:ent(c*inverse) for n,c in explicit.items()},"prefix_operational_score":None,"proof_boundary":"Exact physical 5A1 chamber, effective zero, inherited curves, targets, and unimodular transports.","inputs":{"paths":[str(SOURCE.relative_to(ROOT)),str(CERT.relative_to(ROOT))],"sha256":{str(p.relative_to(ROOT)):sha(p) for p in (SOURCE,CERT)}}}
    mp.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n"); out[zero_name]={"frame":str(fp.relative_to(ROOT)),"marking":str(mp.relative_to(ROOT)),"frame_max_abs":max(abs(int(x)) for x in basis.list())}
summary={"schema":"elkies-k3.h3-physical-q8-5a1-effective-zero-markings.v1","status":"PASS_EXACT_PHYSICAL_Q8_5A1_ALL_EFFECTIVE_ZERO_MARKINGS","outputs":out,"proof_boundary":"All outputs use actual degree-one inherited curves as zero."}
OUTPUT.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n"); print(f"PHYSQ85A1ZEROS|count={len(out)}|status={summary['status']}|output={OUTPUT}")
