#!/usr/bin/env python3
"""Audit the completed gated cover and independently verify the 49-term certificate.

The default performs exact standard-library rational polynomial checks.
--replay-controls also reruns the symbolic construction and regression suite.
It does not relabel a cached search comparison as a fresh enumeration.
"""
from fractions import Fraction as Q
from itertools import permutations
from pathlib import Path
import argparse,ast,hashlib,json,subprocess,sys

HERE=Path(__file__).resolve().parent
SEVEN=HERE.parent/"sparse-seven"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clean(p):return {e:c for e,c in p.items() if c}
def add(a,b,scale=Q(1)):
    r=dict(a)
    for e,c in b.items():r[e]=r.get(e,Q(0))+scale*c
    return clean(r)
def mul(a,b):
    r={}
    for e,c in a.items():
        for f,d in b.items():
            g=tuple(x+y for x,y in zip(e,f))
            r[g]=r.get(g,Q(0))+c*d
    return clean(r)
def derivative(p,j):
    r={}
    for e,c in p.items():
        if e[j]:
            f=list(e);f[j]-=1
            r[tuple(f)]=c*e[j]
    return clean(r)

# Separate one-variable rational polynomial arithmetic for the frame strata.
def uadd(a,b,scale=Q(1)):
    r=dict(a)
    for e,c in b.items():r[e]=r.get(e,Q(0))+scale*c
    return clean(r)
def umul(a,b):
    r={}
    for e,c in a.items():
        for f,d in b.items():r[e+f]=r.get(e+f,Q(0))+c*d
    return clean(r)
def udiv(a,b):
    assert b
    r=dict(a);q={}
    n=max(b)
    while r and max(r)>=n:
        exponent=max(r)-n;coefficient=r[max(r)]/b[n]
        q[exponent]=q.get(exponent,Q(0))+coefficient
        r=uadd(r,{e+exponent:c*coefficient for e,c in b.items()},Q(-1))
    return clean(q),r
def ugcd(a,b):
    while b:
        _,r=udiv(a,b);a,b=b,r
    return a
def parse_polynomial(text):
    def visit(node):
        if isinstance(node,ast.Constant) and type(node.value) is int:
            return clean({0:Q(node.value)})
        if isinstance(node,ast.Name) and node.id=="w":return {1:Q(1)}
        if isinstance(node,ast.UnaryOp):
            value=visit(node.operand)
            if isinstance(node.op,ast.USub):return {e:-c for e,c in value.items()}
            if isinstance(node.op,ast.UAdd):return value
        if isinstance(node,ast.BinOp):
            left=visit(node.left)
            if isinstance(node.op,ast.Pow):
                assert isinstance(node.right,ast.Constant) and type(node.right.value) is int
                exponent=node.right.value
                assert 0<=exponent<=20
                result={0:Q(1)}
                for _ in range(exponent):result=umul(result,left)
                return result
            right=visit(node.right)
            if isinstance(node.op,ast.Add):return uadd(left,right)
            if isinstance(node.op,ast.Sub):return uadd(left,right,Q(-1))
            if isinstance(node.op,ast.Mult):return umul(left,right)
            if isinstance(node.op,ast.Div):
                assert set(right)=={0} and right[0]
                return {e:c/right[0] for e,c in left.items()}
        raise ValueError("Non-polynomial certificate expression")
    return visit(ast.parse(text,mode="eval").body)

def independent_witness():
    path=HERE/"normalized-49-witness.json"
    witness=json.loads(path.read_text())
    outputs=[]
    counts=[]
    degrees=[]
    for row in witness["coefficients"]:
        p={tuple(term["exponent"]):Q(term["coefficient"]) for term in row}
        assert len(p)==len(row) and all(p.values())
        assert all(len(e)==3 and min(e)>=0 for e in p)
        outputs.append(p)
        counts.append(sum(sum(e)>=2 for e in p))
        degrees.append(max(map(sum,p)))
    assert counts==[5,18,26] and degrees==[4,6,7]
    assert all(p.get((0,0,0),Q(0))==0 for p in outputs)
    assert all(sum(c for e,c in p.items() if e[1:]==(0,0))==0 for p in outputs)
    jac=[[derivative(p,j) for j in range(3)] for p in outputs]
    assert [[entry.get((0,0,0),Q(0)) for entry in row] for row in jac]==[
        [1,0,0],[0,1,0],[0,0,1]]
    determinant={}
    for perm in permutations(range(3)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        contribution={(0,0,0):Q(sign)}
        for i in range(3):contribution=mul(contribution,jac[i][perm[i]])
        determinant=add(determinant,contribution)
    assert determinant=={(0,0,0):Q(1)}
    return dict(determinant_one=True,collision_0_e1=True,tangent_identity=True,
                nonlinear_counts=counts,coordinate_degrees=degrees)

def independent_family_inputs(coefficient_rows):
    # Independent construction in Q[x,y,z,w], using an explicit inverse
    # normalization rather than a symbolic matrix inverse or substitution.
    zero=(0,0,0,0)
    def constant(c):return clean({zero:Q(c)})
    def variable(i):
        e=[0]*4;e[i]=1
        return {tuple(e):Q(1)}
    def scale(p,c):return clean({e:v*Q(c) for e,v in p.items()})
    def power(p,n):
        r=constant(1)
        for _ in range(n):r=mul(r,p)
        return r
    X,y,z,w=[variable(i) for i in range(4)]
    one=constant(1)
    m=add(scale(w,4),constant(6))
    n=add(add(scale(power(w,2),-63),scale(w,-189)),constant(-142))
    t=add(constant(2),scale(w,-3))
    dy=add(w,m,Q(-1));dz=add(t,n,Q(-1))
    Y=add(add(m,mul(dy,X)),y)
    Z=add(add(n,mul(dz,X)),z)
    U=add(one,mul(X,Y))
    bracket=add(constant(4),scale(mul(X,Y),3))
    P=add(mul(power(U,3),Z),mul(mul(power(Y,2),U),bracket))
    B=add(add(Y,scale(mul(mul(X,power(U,2)),Z),3)),
          scale(mul(mul(X,power(Y,2)),bracket),3))
    C=add(add(scale(X,2),scale(mul(power(X,2),Y),-3)),
          scale(mul(power(X,3),Z),-1))
    k=add(add(scale(n,3),scale(power(m,2),12)),dy)
    h=add(add(add(scale(mul(m,n),3),scale(power(m,3),7)),
              scale(mul(m,dy),8)),dz)
    p0=mul(add(one,w),add(constant(2),w))
    assert add(add(n,scale(power(m,2),4)),p0,Q(-1))=={}
    H1=scale(C,Q(1,2))
    H2=add(add(B,m,Q(-1)),mul(k,H1),Q(-1))
    H3=add(add(add(P,p0,Q(-1)),mul(h,H1),Q(-1)),scale(mul(m,H2),8),Q(-1))
    expected={}
    for i,entry in enumerate((H1,H2,H3)):
        linear={}
        for exponent,c in entry.items():
            xyz=exponent[:3]
            if sum(xyz)>=2:
                key=(i,xyz)
                expected.setdefault(key,{})[exponent[3]]=c
            else:
                assert sum(xyz)==1 and exponent[3]==0
                linear[xyz]=c
        unit=tuple(1 if j==i else 0 for j in range(3))
        assert linear=={unit:Q(1)}
    actual={(row["row"],tuple(row["exponent"])):parse_polynomial(row["coefficient"])
            for row in coefficient_rows}
    assert len(actual)==len(coefficient_rows)==67
    assert actual==expected,"Independent frame coefficient compilation differs"

def independent_family():
    certificate=json.loads((HERE/"frame-family-certificate.json").read_text())
    assert certificate["generic_support"]==67 and len(certificate["coefficient_rows"])==67
    independent_family_inputs(certificate["coefficient_rows"])
    polynomials=[parse_polynomial(row["coefficient"]) for row in certificate["coefficient_rows"]]
    factors=[parse_polynomial(row["factor"]) for row in certificate["strata"]]
    assert len(factors)==30 and all(f and max(f)>=1 and f[max(f)]==1 for f in factors)
    for i,f in enumerate(factors):
        for g in factors[:i]:
            gcd=ugcd(f,g)
            assert gcd and max(gcd)==0,"Factors must have disjoint roots"
    for polynomial in polynomials:
        assert polynomial
        remainder=polynomial
        for factor in factors:
            while remainder and max(remainder)>=max(factor):
                quotient,r=udiv(remainder,factor)
                if r:break
                remainder=quotient
        assert remainder and max(remainder)==0,"Incomplete factor basis"
    sizes=[]
    for factor,row in zip(factors,certificate["strata"]):
        zeros=[i for i,p in enumerate(polynomials) if not udiv(p,factor)[1]]
        assert zeros==row["zero_coefficient_indices"]
        size=67-len(zeros)
        assert row["support"]==size
        sizes.append(size)
    assert min(sizes)==49 and sizes.count(49)==1
    best=certificate["strata"][sizes.index(49)]
    assert best["factor"]=="w + 2" and certificate["unique_minimizing_parameter"]==-2
    return dict(pairwise_coprime_factors=30,complete_coefficient_factorization=True,
                generic_support=67,minimum=49,unique_minimizing_parameter=-2,
                independent_coefficient_inputs=True,proof_uses_irreducibility_labels=False)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-controls",action="store_true")
    a=parser.parse_args()
    manifest_path=HERE/"results/manifest.json"
    manifest=json.loads(manifest_path.read_text())
    assert manifest["status"]==manifest["enumeration_status"]=="complete"
    assert (manifest["degree"],manifest["bound"],manifest["shard_count"])==(7,9,90)
    assert manifest["source_sha256"]==sha(HERE/"gated_search.py")
    assert manifest["runner_sha256"]==sha(HERE/"run_all.py")
    assert manifest["predecessor_sha256"]==sha(SEVEN/"search.py")
    assert manifest["legacy_sha256"]==sha(HERE.parents[1]/"foundations/jcsearch/low_degree_census.py")
    seeds={(p,q,r) for p in range(2,8) for q in range(2,8) for r in range(q+1,8)}
    assert len(manifest["shards"])==90 and {tuple(r["seed"]) for r in manifest["shards"]}==seeds
    for row in manifest["shards"]:
        path=HERE/"results"/row["path"]
        assert sha(path)==row["sha256"]
        shard=json.loads(path.read_text())
        assert shard["source_sha256"]==manifest["source_sha256"]
        assert shard["parameters"]==dict(degree=7,bound=9,seed=row["seed"])
        assert shard["status"]=="complete" and shard["pending"]==0 and not shard["stack"]
        assert shard["supports"]==shard["retained_masks"]==[]
    assert sum(r["visited"] for r in manifest["shards"])==manifest["visited"]==6040767
    assert manifest["oriented_retained"]==manifest["orbit_count"]==0
    assert manifest["algebra_batches"]==[] and manifest["algebra_backend_count"]==0
    assert sha(HERE/"results/supports.json")==manifest["supports_sha256"]
    coverage=json.loads((HERE/"coverage-validation.json").read_text())
    assert coverage["status"]=="pass" and coverage["retained"]==0 and coverage["supports"]==[]
    assert (coverage["degree"],coverage["bound"],coverage["visited"])==(7,9,7228226)
    assert coverage["manifest_sha256"]==sha(manifest_path)
    assert coverage["source_sha256"]==sha(HERE/"verify_coverage.py")
    if a.replay_controls:
        subprocess.run([sys.executable,str(HERE/"verify_controls.py")],check=True)
    controls=json.loads((HERE/"controls-validation.json").read_text())
    assert controls["status"]=="pass" and controls["witness_all_gates_pass"]
    assert controls["source_sha256"]==sha(HERE/"verify_controls.py")
    assert controls["witness_sha256"]==sha(HERE/"normalized-49-witness.json")
    assert controls["frame_family_sha256"]==sha(HERE/"frame-family-certificate.json")
    assert controls["old_new_set_equality"] and controls["checkpoint_resume"] and controls["torus_equivariance"]
    witness=independent_witness()
    family=independent_family()
    receipt=dict(
        status="pass",claims=["GLDC4","KFRAME1"],source_sha256=sha(Path(__file__)),
        manifest_sha256=sha(manifest_path),coverage_sha256=sha(HERE/"coverage-validation.json"),
        controls_sha256=sha(HERE/"controls-validation.json"),
        production_visits=6040767,independent_visits=7228226,retained_candidates=0,
        lower_bound=10,upper_bound=49,degree_bound=7,arity=3,
        normalized_witness_independent=witness,frame_certificate_independent=family,
        controls_rerun=a.replay_controls,
        excluded_nine_term_ideals=0,
        proof_dependency="KPAIR1 and Moh degree-100 plane theorem",
        boundary="The interval [10,49] concerns the normalized raw-degree-at-most-seven frame; 49 is sharp only among the specified boundary frames.")
    (HERE/"validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps(receipt),flush=True)
if __name__=="__main__":main()

