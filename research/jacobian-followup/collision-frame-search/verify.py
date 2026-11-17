#!/usr/bin/env python3
"""Independent rational construction and certificate audit for the 48-term map."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
from collections import Counter
import ast,gzip,hashlib,importlib.util,json,re

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/"sparse-nine"
spec=importlib.util.spec_from_file_location("retained_rational_checker",OLD/"verify.py")
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
add,mul,clean=base.add,base.mul,base.clean
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n")
def variable(n,i):
    e=[0]*n;e[i]=1
    return {tuple(e):Q(1)}
def constant(n,c):return clean({(0,)*n:Q(c)})
def power(p,n):
    result=constant(len(next(iter(p))),1) if p else {}
    if n==0:return constant(len(next(iter(p))),1)
    for _ in range(n):result=mul(result,p)
    return result
def parse2(text):
    def visit(node):
        if isinstance(node,ast.Constant) and type(node.value)is int:return constant(2,node.value)
        if isinstance(node,ast.Name) and node.id in ("w","t"):return variable(2,("w","t").index(node.id))
        if isinstance(node,ast.UnaryOp):
            p=visit(node.operand)
            if isinstance(node.op,ast.USub):return {e:-c for e,c in p.items()}
            if isinstance(node.op,ast.UAdd):return p
        if isinstance(node,ast.BinOp):
            left=visit(node.left)
            if isinstance(node.op,ast.Pow):
                assert isinstance(node.right,ast.Constant) and type(node.right.value)is int
                assert 0<=node.right.value<=30
                result=constant(2,1)
                for _ in range(node.right.value):result=mul(result,left)
                return result
            right=visit(node.right)
            if isinstance(node.op,ast.Add):return add(left,right)
            if isinstance(node.op,ast.Sub):return add(left,right,Q(-1))
            if isinstance(node.op,ast.Mult):return mul(left,right)
            if isinstance(node.op,ast.Div):
                assert set(right)=={(0,0)}
                return {e:c/right[0,0] for e,c in left.items()}
        raise ValueError("Invalid polynomial expression")
    return visit(ast.parse(text,mode="eval").body)
def monic_key(p):
    lead=p[max(p)]
    return tuple(sorted((e,c/lead) for e,c in p.items()))
def divide2(p,f):
    p=dict(p);q={};r={};lead=max(f);lc=f[lead]
    while p:
        e=max(p);c=p[e]
        difference=tuple(a-b for a,b in zip(e,lead))
        if min(difference)>=0:
            q[difference]=q.get(difference,Q(0))+c/lc
            p=add(p,{tuple(a+b for a,b in zip(j,difference)):v*c/lc for j,v in f.items()},Q(-1))
        else:r[e]=c;del p[e]
    return clean(q),r
def compose(outputs,images):
    n=len(next(iter(images[0])))
    cache={(i,0):constant(n,1) for i in range(len(images))}
    def image_power(i,e):
        if (i,e) not in cache:cache[i,e]=mul(image_power(i,e-1),images[i])
        return cache[i,e]
    result=[]
    for row in outputs:
        out={}
        for exponent,c in row.items():
            term=constant(n,c)
            for i,e in enumerate(exponent):
                if e:term=mul(term,image_power(i,e))
            out=add(out,term)
        result.append(out)
    return result
def shear5(outputs,i,j,k):
    images=[variable(5,h) for h in range(5)]
    images[i]=add(images[i],mul(images[k],images[j]))
    transformed=compose(outputs,images)
    transformed[i]=add(transformed[i],mul(images[k],transformed[j]),Q(-1))
    return transformed
def embed49():
    data=json.loads((OLD/"normalized-49-witness.json").read_text())
    return [{tuple(v["exponent"])+(0,0):Q(v["coefficient"]) for v in row} for row in data["coefficients"]]
def boundary_family():
    data=json.loads((OLD/"frame-family-certificate.json").read_text())
    base.independent_family_inputs(data["coefficient_rows"])
    outputs=[variable(5,i) for i in range(3)]
    for row in data["coefficient_rows"]:
        p=base.parse_polynomial(row["coefficient"])
        for degree,c in p.items():
            e=tuple(row["exponent"])+(degree,0)
            outputs[row["row"]][e]=c
    return outputs
def coefficient_inputs(outputs):
    result={}
    for i,row in enumerate(outputs):
        linear={}
        for e,c in row.items():
            if sum(e[:3])<2:
                linear[e]=c
            else:result.setdefault((i,e[:3]),{})[e[3:]]=c
        assert linear=={tuple(int(j==i) for j in range(5)):Q(1)}
    return result
def check_map(outputs):
    assert [sum(sum(e)>=2 for e in p) for p in outputs]==[6,19,23]
    assert [max(map(sum,p)) for p in outputs]==[4,6,7]
    assert all(p.get((0,0,0),0)==0 for p in outputs)
    assert all(sum(c for e,c in p.items() if e[1:]==(0,0))==0 for p in outputs)
    J=[[base.derivative(p,j) for j in range(3)] for p in outputs]
    assert [[p.get((0,0,0),0) for p in row] for row in J]==[[1,0,0],[0,1,0],[0,0,1]]
    determinant={}
    for perm in permutations(range(3)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
        term=constant(3,sign)
        for i in range(3):term=mul(term,J[i][perm[i]])
        determinant=add(determinant,term)
    assert determinant=={(0,0,0):Q(1)}
def extract_witness():
    maps=shear5(embed49(),2,1,4)
    outputs=[]
    for p in maps:
        out={}
        for e,c in p.items():
            assert e[3]==0
            out[e[:3]]=out.get(e[:3],Q(0))+c*Q(14)**e[4]
        outputs.append(clean(out))
    check_map(outputs)
    witness=dict(status="pass",nonlinear_total=48,nonlinear_counts=[6,19,23],
                 coordinate_degrees=[4,6,7],source_points=[[0,0,0],[1,0,0]],
                 derivative_at_origin="identity",jacobian_determinant=1,
                 source_frame=[["1","0","0"],["0","1","0"],["24","14","1"]],
                 construction="S^(-1) H49(S X), S=I+14 E_32",
                 coefficients=[[dict(exponent=list(e),coefficient=str(c))
                                for e,c in sorted(row.items())] for row in outputs])
    save(HERE/"normalized-48-witness.json",witness)
    return outputs

def verify_family(folder,outputs,expected_count,expected_cases,expected_generic):
    family=json.loads((folder/"two-parameter-family.json").read_text())
    expected=coefficient_inputs(outputs)
    actual={(r["row"],tuple(r["exponent"])):parse2(r["coefficient"]) for r in family["coefficients"]}
    assert actual==expected and len(actual)==expected_count
    factors=[parse2(row["polynomial"]) for row in family["factors"]]
    coefficients=list(actual.values())
    for p in coefficients:
        residual=p
        for f in factors:
            while True:
                q,r=divide2(residual,f)
                if r:break
                residual=q
        assert residual and set(residual)=={(0,0)}
    for f,row in zip(factors,family["factors"]):
        zeros=[i for i,c in enumerate(coefficients) if not divide2(c,f)[1]]
        assert zeros==row["vanishing_coefficients"]
    moving=[i for i,f in enumerate(factors) if max(e[1] for e in f)>0]
    vertical=[i for i in range(len(factors)) if i not in moving]
    generic=expected_count-max(len(family["factors"][i]["vanishing_coefficients"]) for i in moving)
    assert generic==expected_generic
    projection=json.loads((folder/"projection-factors.json").read_text())
    records=json.loads((folder/"projection-checkpoint.json").read_text())["resultants"]
    products=json.loads((folder/"projection-factor-products.json").read_text())["products"]
    assert len(records)==len(products)==len(factors)+len(moving)*(len(moving)-1)//2
    by_key={monic_key(f):i for i,f in enumerate(factors)}
    pairs=[];vert=[];leading=[]
    for record,product_record in zip(records,products):
        result=base.parse_polynomial(record["polynomial"])
        reconstructed={0:Q(product_record["constant"])}
        for index,exponent in product_record["factors"]:
            f=base.parse_polynomial(projection["factors"][index]["polynomial"])
            for _ in range(exponent):reconstructed=base.umul(reconstructed,f)
        assert reconstructed==result
        source=record["source"]
        if source[0]=="pair":
            pairs.append(tuple(sorted(by_key[monic_key(parse2(p))] for p in source[1:])))
        elif source[0]=="vertical":vert.append(monic_key(parse2(source[1])))
        elif source[0]=="leading":leading.append(tuple(sorted(base.parse_polynomial(source[1]).items())))
        else:raise AssertionError
    assert len(pairs)==len(set(pairs))
    assert set(pairs)=={tuple(sorted((i,j))) for k,i in enumerate(moving) for j in moving[:k]}
    assert Counter(vert)==Counter(monic_key(factors[i]) for i in vertical)
    expected_leading=[]
    for i in moving:
        degree=max(e[1] for e in factors[i])
        expected_leading.append(tuple(sorted((e[0],c) for e,c in factors[i].items() if e[1]==degree)))
    assert Counter(leading)==Counter(expected_leading)
    projection_log=(folder/"singular-replay/projection.log").read_text()
    assert "PROJECTION|PASS" in projection_log and "FAIL|" not in projection_log
    projection_script=(folder/"singular-replay/projection.sing").read_text()
    assert projection_script.count("independent projection ")==len(records)
    assert projection_script.count("projection irreducibility ")==expected_cases
    manifest=json.loads((folder/"specialization-manifest-v2.json").read_text())
    assert manifest["status"]=="complete" and manifest["case_count"]==expected_cases and manifest["minimum"]==48
    assert manifest["source_sha256"]==sha(HERE/"specialize.py")
    all_cases=[]
    for path in sorted((folder/"singular-replay").glob("cases-*.json")):
        receipt=json.loads(path.read_text())
        assert receipt["status"]=="pass" and receipt["script_sha256"]==sha(path.with_suffix(".sing"))
        actual_cases=[list(map(int,line.split("|")[1:])) for line in path.with_suffix(".log").read_text().splitlines()
                      if line.startswith("CASE|")]
        assert actual_cases==receipt["cases"]
        all_cases.extend(actual_cases)
    assert [v[0] for v in all_cases]==list(range(expected_cases))
    assert [v[1] for v in all_cases]==[c["minimum"] for c in manifest["cases"]]
    for case in manifest["cases"]:
        assert sha(folder/case["path"])==case["sha256"]
        with gzip.open(folder/case["path"],"rt") as f:data=json.load(f)
        assert data["status"]=="complete" and data["minimum"]==case["minimum"]
    invalid=json.loads((folder/"specialization-manifest.json").read_text())
    assert invalid["status"]=="invalid"
    return dict(coefficients=expected_count,curves=len(factors),projection_cases=expected_cases,
                generic_minimum=expected_generic,minimum=48,
                independent_inputs=True,independent_projection_products=True,
                singular_projection=True,singular_all_specializations=True,
                manifest_sha256=sha(folder/"specialization-manifest-v2.json"))

def main():
    assert base.independent_witness()["nonlinear_counts"]==[5,18,26]
    extract_witness()
    lower=shear5(boundary_family(),2,1,4)
    full=shear5(shear5(embed49(),2,1,3),1,2,4)
    lower_result=verify_family(HERE,lower,70,469,65)
    full_result=verify_family(HERE/"full-transverse",full,111,824,52)
    singular=json.loads((HERE/"singular-validation.json").read_text())
    assert singular["status"]=="pass" and singular["source_sha256"]==sha(HERE/"verify_singular.py")
    receipt=dict(status="pass",claims=["KFRAME2"],source_sha256=sha(Path(__file__)),
                 previous_checker_sha256=sha(OLD/"verify.py"),
                 witness_sha256=sha(HERE/"normalized-48-witness.json"),
                 witness_independent_determinant=True,witness_independent_collision=True,
                 normalized_support_interval=[10,48],
                 boundary_shear_family=lower_result,full_transverse_at_fixed_pair=full_result,
                 singular_sha256=sha(HERE/"singular-validation.json"),
                 invalid_v1_retained=True,
                 boundary="48 is sharp only in the two specified affine frame families. Global optimality remains open.")
    save(HERE/"validation.json",receipt)
    print(json.dumps(receipt),flush=True)
if __name__=="__main__":main()

