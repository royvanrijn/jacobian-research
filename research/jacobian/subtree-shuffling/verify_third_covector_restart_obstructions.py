#!/usr/bin/env python3
"""Complete finite cubic Gamma image and two all-helper polynomial obstructions."""
if not __debug__:raise SystemExit("Refusing optimized Python: exact assertions required")
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
import argparse,json,hashlib,importlib.util,copy
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def kernel(columns):
    """Exact column elimination; returned combinations form the FULL kernel."""
    pivots={};null=[]
    for j,col in enumerate(columns):
        v={k:F(c) for k,c in col.items() if c};combo={j:F(1)}
        while v:
            key=max(v);a=v[key]
            if key not in pivots:
                v={k:c/a for k,c in v.items()};combo={k:c/a for k,c in combo.items()}
                pivots[key]=(v,combo);break
            p,q=pivots[key]
            for k,c in p.items():
                v[k]=v.get(k,F(0))-a*c
                if not v[k]:del v[k]
            for k,c in q.items():
                combo[k]=combo.get(k,F(0))-a*c
                if not combo[k]:del combo[k]
        else:
            assert combo.get(j)==1 and all(k<=j for k in combo)
            null.append(combo)
    assert len(pivots)+len(null)==len(columns)
    return len(pivots),null
def readpoly(rows):return {tuple(e):F(c) for e,c in rows}
def pack(p):return [[list(e),str(c)] for e,c in sorted(p.items())]
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path,required=True)
    ap.add_argument("--input",type=Path,default=HERE/"third_covector_restart_inputs.json")
    ap.add_argument("--ledger",type=Path,default=HERE/"third_covector_restart_source_ledger.json")
    ap.add_argument("--output",type=Path,default=HERE/"third_covector_restart_verification.json")
    a=ap.parse_args();repo=a.repo_root.resolve()
    ledger=json.loads(a.ledger.read_text())
    assert ledger["schema_version"]==1
    assert {x["role"] for x in ledger["sources"]}=={"proof","checker","input","joint_proof","joint_checker","literal_proof","literal_data"}
    for item in ledger["sources"]:
        p=((repo/"research") if item["scope"]=="repo" else a.ledger.parent)/item["path"]
        assert p.is_file() and sha(p)==item["sha256"],"Source pin mismatch: "+str(p)
        if item["role"]=="checker":assert p.resolve()==Path(__file__).resolve()
        if item["role"]=="input":assert p.resolve()==a.input.resolve()
    source=repo/"research/jacobian/subtree-shuffling/verify_commuting_syzygy_slices.py"
    spec=importlib.util.spec_from_file_location("pinned_slices",source)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    inp=json.loads(a.input.read_text());assert inp["schema_version"]==1
    dr,du=m.rows()
    ga={8:m.mon({0:1,1:1},9),9:m.mon({1:2},-36),
        10:m.plus(m.mon({1:1,3:1},-54),m.mon({0:1,4:1},36))}
    for i in range(11):
        z=m.mon({i:1})
        assert m.apply(m.apply(z,dr),ga)==m.apply(m.apply(z,ga),dr)
        assert m.apply(m.apply(z,du),ga)==m.apply(m.apply(z,ga),du)
    t=m.mon({0:-2,2:1})
    norm={i:m.projection(m.mon({i:1}),t,dr) for i in range(3,11)}
    kc=m.plus(norm[9],m.times(m.mon({0:-1,1:1},F(8,5)),norm[8]))
    kr=m.plus(norm[10],m.scale(m.times(m.mon({0:-2,1:-4},F(1,15)),m.times(du[10],norm[8])),-1))
    assert m.apply(kc,ga)==m.mon({1:2},F(-108,5))
    assert m.apply(kr,ga)==m.scale(m.times(m.mon({1:1}),norm[3]),F(108,5))
    triple_six={i:norm[i] for i in (3,4,5,6,7)}
    triple_six[10]=m.plus(kr,m.times(m.mon({1:-1}),m.times(norm[3],kc)))
    for q in triple_six.values():
        assert not m.apply(q,dr) and not m.apply(q,du) and not m.apply(q,ga)
    mons=[]
    for ids in combinations_with_replacement(range(11),3):
        e=[0]*11
        for i in ids:e[i]+=1
        mons.append(tuple(e))
    assert len(mons)==286 and len(set(mons))==286
    def columns(operators):
        return [{(tag,e):c for tag,op in operators for e,c in m.apply({mon:F(1)},op).items()} for mon in mons]
    rank,jcomb=kernel(columns([("d",dr),("e",du)]))
    trank,tcomb=kernel(columns([("d",dr),("e",du),("g",ga)]))
    joint=[{mons[i]:c for i,c in q.items()} for q in jcomb]
    triple=[{mons[i]:c for i,c in q.items()} for q in tcomb]
    assert rank==265 and len(joint)==inp["expected_joint_dimension"]==21
    assert trank==268 and len(triple)==inp["expected_triple_dimension"]==18
    for q in joint:assert not m.apply(q,dr) and not m.apply(q,du)
    for q in triple:
        assert not m.apply(q,dr) and not m.apply(q,du) and not m.apply(q,ga)
        assert all(not any(e[i] for i in (8,9,10)) and sum(e[i] for i in (5,6,7))<=1 for e in q)
    T=m.plus(m.mon({0:1,3:1}),m.mon({1:1,2:1},-2))
    targets=[m.mon({0:2,1:2}),m.mon({0:1,1:3}),m.times(m.mon({1:2}),T)]
    images=[m.apply(q,ga) for q in joint]
    for p in images:
        coeff=[p.get(next(iter(targets[0])),F(0)),p.get(next(iter(targets[1])),F(0)),
               p.get((1,2,0,1,0,0,0,0,0,0,0),F(0))]
        assert p==m.plus(*(m.scale(q,c) for q,c in zip(targets,coeff)))
    irank,_=kernel(images);assert irank==inp["expected_gamma_rank"]==3
    Y=[readpoly(inp[k]) for k in ("Y0","Y1")]
    assert m.apply(Y[0],ga)==m.mon({0:2,1:2},-108)
    assert m.apply(Y[1],ga)==m.mon({0:1,1:3},-108)
    for q in Y:assert not m.apply(q,dr) and not m.apply(q,du)
    literal=json.loads((repo/"research/jacobian/subtree-shuffling/cubic_index_four_syzygy58.json").read_text())
    mp={0:0,47:1,**{i:i-46 for i in range(48,57)}}
    H={}
    for tr in literal["H_terms"]:
        if 48<=tr["output"]<=57:
            e={}
            for i in tr["inputs"]:e[mp[i]]=e.get(mp[i],0)+1
            j=tr["output"]-46;H[j]=m.plus(H.get(j,{}),m.mon(e,F(tr["coefficient"])))
    for i in range(2,12):
        assert not m.apply(H[i],dr)
        assert not m.apply(H[i],du)
        if i<11:assert dr[i]==m.diff(H[i],1)
    eta={i:m.plus(*(m.times(m.diff(H[i],j),m.diff(H[j],0)) for j in range(2,i))) for i in range(2,11)}
    for i in range(11):
        z=m.mon({i:1})
        theta=m.plus(m.apply(m.apply(z,eta),dr),m.scale(m.apply(m.apply(z,dr),eta),-1),m.diff(du.get(i,{}),1))
        assert not theta
    assert m.apply(H[11],ga)==m.scale(m.times(m.mon({1:2}),T),-108)
    aa=[m.mon({0:2}),m.mon({0:1,1:1})]
    alpha=[m.plus(m.times(c,m.diff(H[11],1)),m.scale(m.times(T,m.diff(y,1)),-1)) for c,y in zip(aa,Y)]
    beta=[m.plus(m.times(c,m.apply(H[11],eta)),m.scale(m.times(T,m.apply(y,eta)),-1)) for c,y in zip(aa,Y)]
    for x,y in zip(alpha,beta):assert m.apply(y,dr)==m.apply(x,du)
    target=tuple(inp["first_persistent_monomial"])
    assert beta[0][target]==F(inp["first_value"])==63
    assert all(all(e[1]>=4 or e[4] for e in q) for q in du.values())
    modI=lambda p:{e:c for e,c in p.items() if e[1]<4 and not e[4]}
    assert modI(beta[0])==m.plus(m.mon({0:4,1:3,6:1},48),m.mon({0:3,1:3,3:2},63))
    variations=[m.times(T,m.apply(q,eta)) for q in triple]
    assert all(not q.get(target) for q in variations)
    projected=[modI(q) for q in variations]
    prank,_=kernel(projected);arank,_=kernel(projected+[modI(beta[0])])
    assert prank==1 and arank==2
    functional=inp["functional"];assert len(functional)==6 and all(tag=="d" for tag,e,c in functional)
    def lam(q):return sum((F(c)*q.get(tuple(e),F(0)) for tag,e,c in functional),F(0))
    assert all(not lam(m.apply({e:F(1)},dr)) for e in mons)
    assert all(not lam(m.times(T,m.diff(q,1))) for q in triple)
    assert lam(alpha[1])==F(inp["second_value"])==1312
    # Compare the complete actual polynomial rows for both chosen helpers.
    gm={0:0,1:47,**{i:i+46 for i in range(2,11)}}
    def ordinary(q):
        out={}
        for e,c in q.items():
            assert all(k>=0 for k in e)
            mon=tuple(sorted(gm[i] for i,k in enumerate(e) for _ in range(k)))
            out[mon]=out.get(mon,F(0))+c
        return {e:c for e,c in out.items() if c}
    def product(p,q):return m.multiply({(0,0):p},{(0,0):q}).get((0,0),{})
    cb=m.jacobian([tr for tr in literal["H_terms"] if tr["output"]<48]);cb2=m.multiply(cb,cb)
    old_j=m.jacobian(literal["H_terms"]);old_j2=m.multiply(old_j,old_j)
    omega_fibre={48:m.mon({1:2}),49:m.mon({0:1,1:1},-2),50:m.mon({0:2})}
    fibre_checks=0
    full_row_checks=0
    second_row_checks=0
    a0=m.plus(m.mon({1:1,3:1},6),m.mon({0:1,4:1},-4))
    b0=m.plus(m.mon({0:1,3:1},-12),m.mon({1:1,2:1},24))
    fibres={
      8:{48:m.mon({0:1,1:3},9),49:m.mon({0:2,1:2},-3),50:m.mon({0:3,1:1},-1)},
      9:{48:m.mon({1:4},-36),49:m.mon({0:1,1:3},48),50:m.mon({0:2,1:2},-20)},
      10:{48:m.scale(m.times(m.mon({1:2}),a0),-9),
          49:m.plus(m.scale(m.times(m.mon({0:1,1:1}),a0),3),m.scale(m.times(m.mon({1:2}),b0),3)),
          50:m.plus(m.times(m.mon({0:2}),a0),m.scale(m.times(m.mon({0:1,1:1}),b0),-2))}}
    for i in range(2,11):
        for col in range(58):
            want=product(ordinary(dr[i]),cb.get((47,col),{}))
            if col==0:
                for e,c in ordinary(eta[i]).items():want[e]=want.get(e,F(0))+c
            for e,c in ordinary(fibres.get(i,{}).get(col,{})).items():want[e]=want.get(e,F(0))+c
            assert old_j2.get((gm[i],col),{})=={e:c for e,c in want.items() if c}
            second_row_checks+=1
    assert second_row_checks==522
    for q in joint:
        grad={(0,gm[i]):ordinary(m.diff(q,i)) for i in range(11) if m.diff(q,i)}
        actual_row=m.multiply(grad,old_j2)
        gq=m.apply(q,ga)
        for col in range(58):
            want=product(ordinary(m.diff(q,1)),cb2.get((47,col),{}))
            if col==0:
                for e,c in ordinary(m.apply(q,eta)).items():want[e]=want.get(e,F(0))+c
            for e,c in ordinary(m.times(gq,omega_fibre.get(col,{}))).items():want[e]=want.get(e,F(0))+c
            assert actual_row.get((0,col),{})=={e:c for e,c in want.items() if c}
            full_row_checks+=1
        for col,factor in omega_fibre.items():
            assert actual_row.get((0,col),{})==ordinary(m.times(gq,factor))
            fibre_checks+=1
    assert fibre_checks==63 and full_row_checks==1218
    actual=[]
    for case,y in enumerate(Y):
        ts=list(literal["H_terms"])+[{"output":58,"inputs":list(e),"coefficient":str(c)} for e,c in ordinary(y).items()]
        j=m.jacobian(ts);j2=m.multiply(j,j);j3=m.multiply(j2,j);assert not m.multiply(j3,j)
        ts+=[{"output":59,"inputs":[0,0 if case==0 else 47,57],"coefficient":"1"},
             {"output":59,"inputs":[0,49,58],"coefficient":"-1"},{"output":59,"inputs":[47,48,58],"coefficient":"2"}]
        j=m.jacobian(ts);j2=m.multiply(j,j);j3=m.multiply(j2,j);j4=m.multiply(j3,j)
        assert j4 and all(i==59 for i,col in j4)
        ap,bp=ordinary(alpha[case]),ordinary(beta[case])
        for col in range(60):
            want=product(ap,cb2.get((47,col),{}))
            if col==0:
                for e,c in bp.items():want[e]=want.get(e,F(0))+c
            assert j4.get((59,col),{})=={e:c for e,c in want.items() if c}
        actual.append({"helper_case":case,"valid59block_fourth_entries":0,"candidate_fourth_nonzero_entries":len(j4),"row_columns_compared":60})
    # A malformed third derivation must fail its explicit helper-image identity.
    broken=copy.deepcopy(ga);broken[8]=m.scale(broken[8],2)
    assert m.apply(Y[0],broken)!=m.mon({0:2,1:2},-108)
    out={"status":"PASS","source_ledger_sha256":sha(a.ledger),"ordinary_cubic_domain":286,
      "joint_matrix_rank":rank,"joint_kernel_dimension":len(joint),"Gamma_image_rank":irank,
      "triple_matrix_rank":trank,"triple_kernel_dimension":len(triple),
      "first_projected_variation_rank":prank,"first_augmented_rank":arank,
      "first_persistent_value":"63","second_functional_value":"1312",
      "all_joint_basis_fibre_entries_checked":fibre_checks,"all_joint_basis_full_entries_checked":full_row_checks,"all_preS_second_row_entries_checked":second_row_checks,"all_components_joint_invariance_checked":10,"actual_block_checks":actual,"joint_kernel_basis":[pack(q) for q in joint],"triple_kernel_basis":[pack(q) for q in triple],
      "scope":"Exact finite cubic image and both ALL-helper ordinary-polynomial correction obstructions; no universal parallel-block impossibility or infinite-family claim."}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+"\n")
    print("PASS complete21/18cubicbases,Gammaimage3,ALLhelperobstructions63/1312,1218fulljointrows,522preSsecondrows,120actualresidualcolumns")
    print(a.output)
if __name__=="__main__":main()
