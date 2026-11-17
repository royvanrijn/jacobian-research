#!/usr/bin/env python3
"""Independent Singular replay of factor, projection and finite-fibre certificates."""
from pathlib import Path
import gzip,hashlib,json,re,subprocess,time
import sympy as sp
HERE=Path(__file__).resolve().parent
w,t=sp.symbols("w t")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def text(expr):
    result=re.sub(r"([A-Za-z_]\w*)\*\*(\d+)",r"(\1^\2)",str(expr))
    assert "**" not in result
    return result
def save(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+".tmp")
    tmp.write_text(json.dumps(v,indent=2,sort_keys=True)+"\n");tmp.replace(p)
def load_case(folder,row):
    path=folder/row["path"];assert sha(path)==row["sha256"]
    with gzip.open(path,"rt") as f:return json.load(f)
def field_poly(record):
    parts=[]
    for term in record:
        coefficient="+".join("("+c+")*(a^"+str(i)+")" for i,c in enumerate(term["w_coefficients"])) or "0"
        parts.append("("+coefficient+")*(t^"+str(term["t_degree"])+")")
    return "+".join(parts) or "0"

def checked_run(script,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(script)
    result=subprocess.run(["Singular","-q",str(path)],capture_output=True,text=True,timeout=300)
    (path.with_suffix(".log")).write_text(result.stdout+result.stderr)
    assert result.returncode==0 and "FAIL|" not in result.stdout and "error" not in result.stdout.lower(),result.stdout[-3000:]
    return result.stdout

def projection_check(folder):
    family=json.loads((folder/"two-parameter-family.json").read_text())
    projection=json.loads((folder/"projection-factors.json").read_text())
    checkpoint=json.loads((folder/"projection-checkpoint.json").read_text())
    factors=[row["polynomial"] for row in projection["factors"]]
    index={str(sp.Poly(f,w,domain=sp.QQ).monic().as_expr()):i for i,f in enumerate(factors)}
    lines=["ring R=0,(w,t),dp;","option(redSB);",
           'proc require(int ok,string label){if(ok==0){print("FAIL|"+label);exit;}}',
           "ideal Projection="+",".join(text(f) for f in factors)+";"]
    for i in range(len(factors)):
        lines.extend([f"list fact{i}=factorize(Projection[{i+1}]);",
                      f"ideal ip{i}=fact{i}[1];",
                      "int count=0;int k;",
                      f'for(k=1;k<=size(ip{i});k++){{if(deg(ip{i}[k])>0){{count++;require(deg(ip{i}[k])==deg(Projection[{i+1}]),"projection squarefree irreducibility {i}");}}}}',
                      f'require(count==1,"projection irreducibility {i}");',
                      "kill count,k;"])
    curves=[row["polynomial"] for row in family["factors"]]
    coefficients=[row["coefficient"] for row in family["coefficients"]]
    lines.extend(["ideal Curves="+",".join(text(f) for f in curves)+";",
                  "ideal Coefficients="+",".join(text(f) for f in coefficients)+";"])
    for i in range(len(curves)):
        for j in range(i):
            lines.append(f'require(deg(gcd(Curves[{i+1}],Curves[{j+1}]))==0,"curve coprimality {i},{j}");')
    lines.extend(["int i;int j;poly p;",
                  "for(i=1;i<=size(Coefficients);i++){",
                  "p=Coefficients[i];",
                  "for(j=1;j<=size(Curves);j++){while(deg(p)>=deg(Curves[j]) && reduce(p,ideal(Curves[j]))==0){p=p/Curves[j];}}",
                  'require(p!=0 && deg(p)==0,"complete coefficient factorization");',
                  "}","kill i,j,p;"])
    products=[]
    for number,row in enumerate(checkpoint["resultants"]):
        expression=sp.Poly(row["polynomial"],w,domain=sp.QQ)
        constant,parts=sp.factor_list(expression)
        selected=[]
        for f,m in parts:
            constant*=f.LC()**m
            selected.append((index[str(f.monic().as_expr())],m))
        product_expr=text(constant)+"".join(f"*(Projection[{i+1}]^{m})" for i,m in selected)
        lines.append(f"poly result{number}={text(expression.as_expr())};")
        lines.append(f'require(result{number}==({product_expr}),"projection product {number}");')
        source=row["source"]
        if source[0]=="pair":
            independent=f"resultant({text(source[1])},{text(source[2])},t)"
        else:independent=text(source[1])
        lines.append(f'require(result{number}==({independent}),"independent projection {number}");')
        products.append(dict(constant=str(constant),factors=selected))
    lines.extend(['print("PROJECTION|PASS");',"exit;"])
    output=checked_run("\n".join(lines)+"\n",folder/"singular-replay/projection.sing")
    assert "PROJECTION|PASS" in output
    save(folder/"projection-factor-products.json",dict(status="pass",products=products,
        source_sha256=sha(Path(__file__)),projection_sha256=sha(folder/"projection-factors.json")))
    return family,projection

def cases_check(folder,family,projection):
    manifest_path=folder/"specialization-manifest-v2.json"
    manifest=json.loads(manifest_path.read_text())
    assert manifest["status"]=="complete" and manifest["minimum"]==48
    records=manifest["cases"]
    results=[]
    for offset in range(0,len(records),25):
        batch=records[offset:offset+25]
        script_lines=['proc require(int ok,string label){if(ok==0){print("FAIL|"+label);exit;}}']
        expected=[]
        for local,row in enumerate(batch):
            case=load_case(folder,row)
            r=sp.Poly(case["projection_polynomial"],w,domain=sp.QQ)
            name=f"K{local}"
            if r.degree()==1:
                value=-r.nth(0)/r.nth(1)
                script_lines.append(f"ring {name}=0,(t),dp;")
                coeff=[text(sp.expand(sp.sympify(c["coefficient"]).subs(w,value))) for c in family["coefficients"]]
                def basis_expr(record):
                    parts=[]
                    for term in record:
                        c=sum(sp.Rational(v)*value**i for i,v in enumerate(term["w_coefficients"]))
                        parts.append(f"({text(c)})*(t^{term['t_degree']})")
                    return "+".join(parts) or "0"
            else:
                script_lines.extend([f"ring {name}=(0,a),(t),dp;",
                                     "minpoly="+text(r.as_expr()).replace("w","a")+";"])
                coeff=[text(c["coefficient"]).replace("w","a") for c in family["coefficients"]]
                basis_expr=field_poly
            basis=[basis_expr(row["polynomial"]) for row in case["bases"]]
            assert basis
            script_lines.extend(["list P="+",".join("("+c+")" for c in coeff)+";",
                                 "list B="+",".join("("+c+")" for c in basis)+";",
                                 "int i;int j;int k;int count;int current;int best=0;poly q;poly b;"])
            script_lines.append("for(i=1;i<=size(P);i++){if(P[i]!=0){best++;}}")
            script_lines.append(f'require(best=={case["generic_specialized_support"]},"generic specialization");')
            script_lines.extend(["for(i=1;i<=size(B);i++){",
                                 'require(deg(B[i])>0 && leadcoef(B[i])==1,"monic positive root polynomial");',
                                 "for(j=1;j<i;j++){",
                                 'require(deg(gcd(B[i],B[j]))==0,"coprime root basis");',"}","}"])
            script_lines.extend(["for(i=1;i<=size(P);i++){q=P[i];if(q!=0){",
                                 "for(j=1;j<=size(B);j++){",
                                 "while(deg(q)>=deg(B[j]) && reduce(q,ideal(B[j]))==0){q=q/B[j];}",
                                 "}",
                                 'require(q!=0 && deg(q)==0,"complete specialized factors");',
                                 "}}"])
            script_lines.extend(["for(i=1;i<=size(B);i++){count=0;",
                                 "for(j=1;j<=size(P);j++){if(reduce(P[j],ideal(B[i]))==0){count++;}}",
                                 f"current={family['generic_support']}-count;",
                                 "if(current<best){best=current;}",
                                 "}",
                                 f'require(best=={case["minimum"]},"minimum {offset+local}");',
                                 f'print("CASE|{offset+local}|"+string(best));',
                                 f"kill {name};"])
            expected.append((offset+local,case["minimum"]))
        script_lines.append("exit;")
        script="\n".join(script_lines)+"\n"
        path=folder/f"singular-replay/cases-{offset:04d}.sing"
        receipt_path=path.with_suffix(".json")
        digest=hashlib.sha256(script.encode()).hexdigest()
        if receipt_path.exists():
            receipt=json.loads(receipt_path.read_text())
            assert receipt["script_sha256"]==digest
        else:
            output=checked_run(script,path)
            actual=[tuple(map(int,line.split("|")[1:])) for line in output.splitlines() if line.startswith("CASE|")]
            assert actual==expected,(actual,expected)
            receipt=dict(status="pass",script_sha256=digest,cases=[list(v) for v in actual])
            save(receipt_path,receipt)
        results.extend(receipt["cases"])
        print(json.dumps(dict(family=folder.name,verified=len(results),total=len(records))),flush=True)
    assert len(results)==len(records) and min(v[1] for v in results)==48
    return dict(case_count=len(results),minimum=48,manifest_sha256=sha(manifest_path))

def main():
    summaries=[]
    for folder in (HERE,HERE/"full-transverse"):
        family,projection=projection_check(folder)
        summary=cases_check(folder,family,projection)
        summaries.append(dict(directory=str(folder.relative_to(HERE)),**summary))
    save(HERE/"singular-validation.json",dict(status="pass",source_sha256=sha(Path(__file__)),
                                             families=summaries))
if __name__=="__main__":main()

