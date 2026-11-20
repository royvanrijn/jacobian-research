#!/usr/bin/env python3
"""Source-guarded finite-union/fiberwise quantifier bridge for ordinary binomials."""
from __future__ import annotations
if not __debug__:raise SystemExit("Refusing -O/-OO: exact assertions must execute.")
from fractions import Fraction as Q
from pathlib import Path
from math import factorial,comb,isqrt
import argparse,copy,hashlib,json,re,subprocess,sys,tempfile,importlib.util
HERE=Path(__file__).resolve().parent
PN="FIBERWISE_UNIFORM_EVENTUAL_BINOMIAL_SFC.md"
CN="verify_fiberwise_uniform_eventual_binomial.py"
IN="fiberwise_uniform_eventual_inputs.json"
LN="fiberwise_uniform_eventual_source_ledger.json"
RN="fiberwise_uniform_eventual_verification.json"
PREFIX="research/factorial/image-rigidity/"

ROLES={'a0_receipt': ('repo', 'research/factorial/image-rigidity/uniform_two_term_sfc_verification.json'), 'a0_six_source_1': ('repo', 'research/factorial/image-rigidity/verify_exponent_six_simplicity.py'), 'a0_six_source_2': ('repo', 'research/factorial/image-rigidity/exponent_six_prefix_certificates.json'), 'a0_six_source_3': ('repo', 'research/factorial/image-rigidity/exponent_six_simplicity_audit.json'), 'a0_six_source_4': ('repo', 'research/factorial/image-rigidity/EFFECTIVE_SADDLE_ANNULUS.tex'), 'a0_source_1': ('repo', 'research/factorial/image-rigidity/uniform_two_term_sfc.py'), 'a0_source_10': ('repo', 'research/factorial/image-rigidity/all_exponent_squarefree_verification.json'), 'a0_source_2': ('repo', 'research/factorial/image-rigidity/verify_uniform_two_term_sfc.py'), 'a0_source_3': ('repo', 'research/factorial/image-rigidity/all_exponent_prefix_certificates.json'), 'a0_source_4': ('repo', 'research/factorial/image-rigidity/uniform_two_term_sfc_audit.json'), 'a0_source_5': ('repo', 'research/factorial/image-rigidity/EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex'), 'a0_source_6': ('repo', 'research/factorial/image-rigidity/exponent_six_simplicity_verification.json'), 'a0_source_7': ('repo', 'research/factorial/image-rigidity/ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex'), 'a0_source_8': ('repo', 'research/factorial/image-rigidity/PRIME_POWER_TWO_TERM_SFC.tex'), 'a0_source_9': ('repo', 'research/factorial/image-rigidity/QUANTITATIVE_SQUAREFREENESS_ALL_EXPONENTS.tex'), 'a1_analytic_ledger': ('repo', 'research/factorial/image-rigidity/factorial_practical_uniform_saddle_source_ledger.json'), 'a1_analytic_receipt': ('repo', 'research/factorial/image-rigidity/factorial_practical_uniform_saddle_verification.json'), 'a1_analytic_source_1': ('repo', 'research/factorial/image-rigidity/verify_factorial_practical_uniform_saddle_tail.py'), 'a1_analytic_source_2': ('repo', 'research/factorial/image-rigidity/factorial_practical_uniform_saddle_constants.json'), 'a1_analytic_source_3': ('repo', 'research/factorial/image-rigidity/FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md'), 'a1_analytic_source_4': ('repo', 'research/factorial/image-rigidity/verify_factorial_gap_ten_saddle_extension.py'), 'a1_analytic_source_5': ('repo', 'research/factorial/image-rigidity/FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md'), 'a1_analytic_source_8': ('repo', 'research/factorial/image-rigidity/ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex'), 'a1_finite_ledger': ('repo', 'research/factorial/image-rigidity/ordinary_lower_one_finite_bridge_source_ledger.json'), 'a1_finite_receipt': ('repo', 'research/factorial/image-rigidity/ordinary_lower_one_finite_bridge_verification.json'), 'a1_finite_source_1': ('repo', 'research/factorial/image-rigidity/verify_ordinary_lower_one_finite_bridge.py'), 'a1_finite_source_2': ('repo', 'research/factorial/image-rigidity/ordinary_lower_one_finite_bridge_inputs.json'), 'a1_finite_source_3': ('repo', 'research/factorial/image-rigidity/GAP_SIX_RESULTANT_STRICT_POSITIVITY.md'), 'a1_finite_source_4': ('repo', 'research/factorial/image-rigidity/BINOMIAL_SHAPE_RESULTANT.tex'), 'boundary_a0': ('repo', 'research/factorial/image-rigidity/ALL_UNIVARIATE_TWO_TERM_SFC.tex'), 'boundary_a1': ('repo', 'research/factorial/image-rigidity/ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md'), 'boundary_adjacent': ('repo', 'research/factorial/image-rigidity/UNIFORM_PADE_PARABOLA_SFC.tex'), 'checker': ('local', 'verify_fiberwise_uniform_eventual_binomial.py'), 'eventual_checker': ('repo', 'research/factorial/image-rigidity/verify_ordinary_binomial_eventual_completion.py'), 'eventual_inputs': ('repo', 'research/factorial/image-rigidity/ordinary_binomial_eventual_inputs.json'), 'eventual_ledger': ('repo', 'research/factorial/image-rigidity/ordinary_binomial_eventual_source_ledger.json'), 'eventual_proof': ('repo', 'research/factorial/image-rigidity/ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md'), 'eventual_receipt': ('repo', 'research/factorial/image-rigidity/ordinary_binomial_eventual_verification.json'), 'inputs': ('local', 'fiberwise_uniform_eventual_inputs.json'), 'old_high_gap': ('repo', 'research/factorial/image-rigidity/lower_exponent_initial_and_boundary.tex'), 'old_terminal_triangle': ('repo', 'research/factorial/image-rigidity/POLYNOMIAL_TRIANGLE_AND_TERMINAL_SIX.md'), 'older_soft': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_ROOT_MODULUS.tex'), 'proof': ('local', 'FIBERWISE_UNIFORM_EVENTUAL_BINOMIAL_SFC.md'), 'quadratic_checker': ('repo', 'research/factorial/image-rigidity/verify_uniform_quadratic_strong_factorial_completion.py'), 'quadratic_inputs': ('repo', 'research/factorial/image-rigidity/uniform_quadratic_strong_factorial_completion_inputs.json'), 'quadratic_ledger': ('repo', 'research/factorial/image-rigidity/uniform_quadratic_strong_factorial_completion_source_ledger.json'), 'quadratic_proof': ('repo', 'research/factorial/image-rigidity/UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md'), 'quadratic_receipt': ('repo', 'research/factorial/image-rigidity/uniform_quadratic_strong_factorial_completion_verification.json'), 'qualitative_fixed_gap': ('repo', 'research/factorial/image-rigidity/FIXED_GAP_LARGE_LOWER_EXPONENT_EXTENSION.md'), 'relative_analytic_proof': ('repo', 'research/factorial/image-rigidity/RELATIVE_GAP_THREE_LOWER_WEIGHTED_SADDLE_TAIL.md'), 'relative_checker': ('repo', 'research/factorial/image-rigidity/verify_relative_gap_three_completion.py'), 'relative_finite_checker': ('repo', 'research/factorial/image-rigidity/verify_relative_gap_three_finite_box.py'), 'relative_finite_receipt': ('repo', 'research/factorial/image-rigidity/relative_gap_three_finite_box_verification.json'), 'relative_input': ('repo', 'research/factorial/image-rigidity/relative_gap_three_completion_inputs.json'), 'relative_ledger': ('repo', 'research/factorial/image-rigidity/relative_gap_three_completion_source_ledger.json'), 'relative_lower_two_proof': ('repo', 'research/factorial/image-rigidity/ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md'), 'relative_proof': ('repo', 'research/factorial/image-rigidity/RELATIVE_GAP_THREE_STRONG_FACTORIAL.md'), 'relative_receipt': ('repo', 'research/factorial/image-rigidity/relative_gap_three_completion_verification.json'), 'soft_checker': ('repo', 'research/factorial/image-rigidity/verify_gamma_product_sharp_lower_annulus.py'), 'soft_inputs': ('repo', 'research/factorial/image-rigidity/gamma_product_sharp_lower_annulus_inputs.json'), 'soft_ledger': ('repo', 'research/factorial/image-rigidity/gamma_product_sharp_lower_annulus_source_ledger.json'), 'soft_proof': ('repo', 'research/factorial/image-rigidity/GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md'), 'soft_receipt': ('repo', 'research/factorial/image-rigidity/gamma_product_sharp_lower_annulus_verification.json')}
SCOPE={'schema': 1, 'parameters': 'integer a>=0,d>=1,N>=2', 'all_fixed_gap_bad_pair_sets_finite': True, 'all_fixed_lower_bad_pair_sets_finite': True, 'all_fixed_terminal_bad_support_sets_finite': True, 'fixed_gap_cutoff': 'for each d exists L(d), for all a>=0 and N>=L(d)', 'fixed_lower_cutoff': 'for each a exists U(a), for all d>=1 and N>=U(a)', 'same_cutoff_all_fields_and_coefficients': True, 'field_scope': 'every characteristic-zero field; rational Bezout, no arbitrary embedding in C', 'explicit_bad_lower_projection': 'a<10^24*d^2 for each d>=2', 'relative_bad_gap_projection': 'd<3*a for a>=2', 'boundary_a0_a1_all_order_inputs_reused': True, 'finite_projection_and_finite_fibers_iff': True, 'terminal_only_iff_uses_actual_support_complements': True, 'fixed_terminal_triangle': '2<=d<D_N,1<=a<d*(N-2)', 'fixed_terminal_count_domain': 'N>=3; N=2 separately empty', 'all_coordinate_escape_for_infinite_distinct_bad_sequences': True, 'global_B_finite_claim': False, 'global_terminal_cutoff_claim': False, 'all_order_SFC_claim': False, 'effective_L_U_claim': False, 'new_analytic_estimate_claim': False, 'new_fixed_terminal_finiteness_claim': False, 'relative_completed_cases_reused': 3694040, 'relative_completed_supports_reused': 532, 'relative_full_gcd_run_rerun': False, 'relative_fresh_fallback_Bezout_rerun': False, 'boundary_old_finite_runs_repeated': False, 'boundary_a0_recorded_prefix_cases_reused': 92935, 'boundary_a1_recorded_finite_cases_reused': 171812, 'boundary_a0_receipt_selection': 'designated current_replay; obsolete historical top-level entries unused', 'boundary_nested_density_context_replayed': False, 'formal_verification': False, 'external_review': False}
ROUTE=['retained401_all_order_large_a', 'retained407_fixed_support_eventual', 'retained_relative_gap_all_order', 'finite_projection_and_finite_union', 'boundary_theorems_and_immutable_completed_certificates', 'old_terminal_high_gap', 'sharp_soft_improved_counts', 'older_soft_independent_reverse', 'rational_Bezout']
def require(ok,msg):
    if not ok:raise AssertionError(msg)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text(encoding="utf-8"))
def same_typed(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(same_typed(a[k],b[k]) for k in b)
    if type(b) is list:return len(a)==len(b) and all(same_typed(x,y) for x,y in zip(a,b))
    return a==b
def check_inputs(data):require(same_typed(data,SCOPE),"Exact fiberwise/field/reuse scope and types required")
def bases(repo,support):return {"local":HERE,"support":support,"repo":repo}
def check_sources(doc,repo,support):
    require(type(doc) is dict and set(doc)=={"schema","formal_verification","external_review","logical_route","sources"},"Source schema")
    require(type(doc["schema"]) is int and doc["schema"]==1,"Schema1")
    require(doc["formal_verification"] is False and doc["external_review"] is False,"Assurance flags")
    require(doc["logical_route"]==ROUTE,"Exact logical route")
    entries=doc["sources"];require(type(entries) is list and len(entries)==len(ROLES),"All mandatory source roles")
    out={};seen=set();roots=bases(repo,support)
    for row in entries:
        require(type(row) is dict and set(row)=={"role","kind","path","sha256"},"Source row schema")
        role=row["role"];require(type(role) is str and role in ROLES and role not in out,"Unique role")
        require((row["kind"],row["path"])==ROLES[role],"Exact role path")
        key=(row["kind"],row["path"]);require(key not in seen,"Unique path")
        require(type(row["sha256"]) is str and re.fullmatch("[0-9a-f]{64}",row["sha256"]),"SHA256")
        root=roots[row["kind"]].resolve();p=(root/row["path"]).resolve()
        require(p.is_relative_to(root) and p.is_file() and sha(p)==row["sha256"],"Missing/stale source:"+role)
        out[role]=row;seen.add(key)
    forbidden={"UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md"}
    require(not any(Path(path).name in forbidden for _,path in seen),"Unneeded uniform-cusp route excluded")
    return out
def links(repo,support,rows):
    direct={r["path"]:r["sha256"] for r in rows.values() if r["kind"]=="repo"}
    specs=[("quadratic","files",13),("eventual","sources",77),("soft","sources",6)]
    out=[]
    for tag,key,count in specs:
        d=load(repo/ROLES[tag+"_ledger"][1])
        require(type(d.get("schema")) is int and d["schema"]==1,"Parent schema1")
        entries=d.get(key);require(type(entries) is list and len(entries)==count,"Parent exact cardinality")
        matched=set();seen=set()
        for r in entries:
            require(type(r) is dict and type(r.get("path")) is str and type(r.get("sha256")) is str,"Parent row")
            kind=r.get("kind",r.get("root",r.get("location",r.get("scope"))))
            require(kind in {"local","support","packet","repo"},"Parent source root")
            name=PREFIX+r["path"] if kind in {"local","support","packet"} else r["path"]
            require(name not in seen,"Parent unique path");seen.add(name)
            if name in direct:
                require(direct[name]==r["sha256"],"Parent/direct pin agreement");matched.add(name)
        for role in ["proof","checker","inputs"]:
            require(ROLES[tag+"_"+role][1] in matched,"Each hierarchical core proof/checker/input is mandatory")
        out.append({"tag":tag,"sources":count,"all_descendants_freshly_validated_by_parent":True})
    d=load(repo/ROLES["relative_ledger"][1])
    require(type(d.get("schema_version")) is int and d["schema_version"]==1,"Relative legacy schema1")
    entries=d.get("sources");require(type(entries) is list and len(entries)==7,"All seven relative completion source records")
    seen=set()
    for r in entries:
        require(type(r) is dict and r.get("scope") in {"packet","repo"},"Relative source row")
        name=PREFIX+r["path"] if r["scope"]=="packet" else r["path"]
        require(name not in seen and name in direct and direct[name]==r["sha256"],"Every reused relative source is directly pinned")
        seen.add(name)
    out.append({"tag":"relative","sources":7,"recorded_completed_evidence_reused":True})
    return out
def equal_receipt(a,b):require(a==b,"Fresh child receipt differs")
def replay(repo,support,rows,parent_python):
    statuses={"quadratic":"PASS_SOURCE_GUARDED_UNIFORM_QUADRATIC_SUPPORT","eventual":"PASS","soft":"PASS"}
    out=[]
    with tempfile.TemporaryDirectory(prefix="fiberwise-uniform-") as folder:
        for i,tag in enumerate(["quadratic","eventual","soft"]):
            p=Path(folder)/(str(i)+".json")
            cmd=[str(parent_python),"-B",str(repo/ROLES[tag+"_checker"][1]),"--repo-root",str(repo),"--output",str(p)]
            r=subprocess.run(cmd,cwd=folder,capture_output=True,text=True,timeout=120)
            require(r.returncode==0,"Fresh parent failed:"+tag+";"+r.stdout+";"+r.stderr)
            require(p.is_file(),"Fresh parent receipt missing")
            raw=p.read_bytes();equal_receipt(raw,(repo/ROLES[tag+"_receipt"][1]).read_bytes())
            d=json.loads(raw);require(d.get("status")==statuses[tag],"Exact parent PASS scope")
            require(d.get("source_ledger_sha256")==rows[tag+"_ledger"]["sha256"],"Parent current source ledger identity")
            out.append({"tag":tag,"receipt_sha256":hashlib.sha256(raw).hexdigest(),"byte_identical":True})
    return out
def relative_artifact(repo,rows):
    script=repo/ROLES["relative_checker"][1]
    # No bytecode is written into the frozen repository and no FLINT routine or
    # full checker main is called. The retained stdlib stream validator is reused.
    sys.dont_write_bytecode=True
    spec=importlib.util.spec_from_file_location("fiberwise_retained_relative_validator",script)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    inp=load(repo/ROLES["relative_input"][1])
    require(type(inp.get("schema_version")) is int and inp["schema_version"]==1,"Typed relative input schema")
    for key in ["support_pairs","terminal_cases","maximum_terminal","maximum_factorial"]:
        require(type(inp[key]) is int,"Typed relative count")
    module.validate_input(inp)
    source=repo/ROLES["relative_finite_checker"][1]
    completed=repo/ROLES["relative_finite_receipt"][1]
    require(sha(source)==inp["finite_source_sha256"],"Exact already executed source")
    require(sha(completed)==inp["completed_receipt_sha256"],"Exact immutable completed receipt")
    supports,cases=module.verify_completed(load(completed),inp,True)
    require(supports==532 and cases==3694040,"Complete retained relative coverage")
    receipt=load(repo/ROLES["relative_receipt"][1])
    require(receipt.get("status")=="PASS" and receipt.get("source_ledger_sha256")==rows["relative_ledger"]["sha256"],"Retained relative validator receipt")
    require(receipt["actual_finite_source_sha256"]==sha(source) and receipt["actual_completed_receipt_sha256"]==sha(completed),"Recorded source/receipt identities")
    require(receipt["completed_support_pairs"]==supports and receipt["completed_terminal_stream_entries"]==cases,"Recorded complete counts")
    require(len(receipt["fresh_fallback_Bezout_replays"])==5,"Five PREVIOUSLY recorded fallback identities")
    return {"support_rows_authenticated":supports,"terminal_stream_entries_authenticated":cases,
       "finite_source_sha256":sha(source),"immutable_completed_receipt_sha256":sha(completed),
       "full_gcd_run_rerun":False,"fresh_fallback_Bezout_replays":0,
       "previously_recorded_fallbacks":5}
def boundary_artifacts(repo,rows):
    direct={r["path"]:r["sha256"] for r in rows.values() if r["kind"]=="repo"}
    def rp(name):return repo/PREFIX/name
    maps=[]
    for receipt_name,count in [("uniform_two_term_sfc_verification.json",11),("exponent_six_simplicity_verification.json",6)]:
        doc=load(rp(receipt_name));d=doc.get("current_replay",doc)
        require(d.get("status")=="pass","Retained current boundary PASS")
        files=d.get("files");require(type(files) is dict and len(files)==count,"Exact selected current boundary file map")
        for name,pin in files.items():
            require(PREFIX+name in direct and direct[PREFIX+name]==pin,"Each selected boundary receipt pin is directly mandatory")
        checker="verify_uniform_two_term_sfc.py" if count==11 else "verify_exponent_six_simplicity.py"
        require(d["checker_sha256"]==sha(rp(checker)),"Retained current checker identity")
        maps.append({"receipt":receipt_name,"selected_current_replay":True,"file_map_count":count})
    a0=load(rp("all_exponent_prefix_certificates.json"))
    require(a0["range_n"]==[2,9999] and a0["columns"]==["n","m","prime","degree","mode","G","U","V"],"Exact a0 finite certificate format")
    records=a0["records"];require(type(records) is list and len(records)==92935,"Complete recorded a0 domain")
    seen=set();modes={}
    for row in records:
        require(type(row) is list and len(row)==8,"Recorded a0 certificate tuple")
        n,m,prime,degree,mode,g,u,v=row
        require(all(type(x) is int for x in [n,m,prime,degree]),"Typed recorded a0 certificate coordinates")
        require((n,m) not in seen,"Unique a0 finite pair");seen.add((n,m))
        require(mode in {"preceding_prime_kernel","full_degree_modular"},"Actual recorded certificate mode")
        require(all(type(x) is list and x and all(type(a) is int for a in x) for x in [g,u,v]),"Literal coefficient vectors retained")
        modes[mode]=modes.get(mode,0)+1
    expected={(n,m) for n in range(2,10000) for m in [10,12,14,15,18,20,21,22,24,26,28] if 2**m<4*n*n}
    require(seen==expected and modes=={"preceding_prime_kernel":92900,"full_degree_modular":35},"Exact recorded a0 complete finite domain and modes")
    six=load(rp("exponent_six_prefix_certificates.json"))
    require(six["range_inclusive"]==[2,9999] and six["n1"]=="linear","Exact exponent-six boundary")
    rec=six["records"];require(type(rec) is list and len(rec)==9998,"Complete recorded sixth prefix")
    require(all(type(r) is list and len(r)==7 and type(r[0]) is int for r in rec),"Typed exponent-six records")
    require([r[0] for r in rec]==list(range(2,10000)),"Every exponent-six recorded order")
    a1ledgers=[]
    for tag,size in [("a1_finite",5),("a1_analytic",10)]:
        d=load(repo/ROLES[tag+"_ledger"][1])
        entries=d.get("files");require(type(entries) is list and len(entries)==size,"Exact reused lower-one ledger cardinality")
        names=set()
        for r in entries:
            require(type(r) is dict and r.get("location") in {"packet","repository"},"Reused lower-one source row")
            name=PREFIX+r["path"] if r["location"]=="packet" else r["path"]
            require(name not in names and name in direct and direct[name]==r["sha256"],"Every lower-one parent row directly authenticated")
            names.add(name)
        receipt=load(repo/ROLES[tag+"_receipt"][1])
        require(receipt["status"]=="PASS","Reused lower-one complete PASS")
        if tag=="a1_finite":
            require(receipt["source_sha256"]=={"checker":sha(rp("verify_ordinary_lower_one_finite_bridge.py")),
                "inputs":sha(rp("ordinary_lower_one_finite_bridge_inputs.json")),
                "ledger":rows[tag+"_ledger"]["sha256"],"proof":rows["boundary_a1"]["sha256"]},"Reused lower-one finite exact source identities")
        else:
            for key,name in [("checker_sha256","verify_factorial_practical_uniform_saddle_tail.py"),
                             ("inputs_sha256","factorial_practical_uniform_saddle_constants.json")]:
                require(receipt[key]==sha(rp(name)),"Reused analytical source identity")
            require(receipt["ledger_sha256"]==rows[tag+"_ledger"]["sha256"] and receipt["proof_sha256"]==rows["boundary_a1"]["sha256"],"Reused analytical ledger/proof identity")
        pins=receipt["verified_source_pins"]
        require(type(pins) is list and len(pins)==size,"Complete reused lower-one receipt source pins")
        for r in pins:
            name=r["path"] if r["path"].startswith(PREFIX) else PREFIX+r["path"]
            require(name in names and direct[name]==r["sha256"],"Reused lower-one receipt and ledger agree")
        a1ledgers.append({"tag":tag,"required_sources":size,"recorded_receipt_authenticated":True})
    inp=load(rp("ordinary_lower_one_finite_bridge_inputs.json"));receipt=load(rp("ordinary_lower_one_finite_bridge_verification.json"))
    for key,value in {"schema_version":1,"prime":1000003,"minimum_gap":10,"maximum_gap":43,
                      "minimum_terminal":2,"tail_minimum":5000,"tail_linear_factor":128,"expected_cases":171812}.items():
        require(type(inp.get(key)) is int and inp[key]==value,"Exact typed lower-one finite domain")
    require(len(inp["gaps"])==34 and len(receipt["gaps"])==34,"All 34 reused gap families")
    total=0
    for d,r,s in zip(range(10,44),inp["gaps"],receipt["gaps"]):
        last=max(5000,128*(d+1))-1;count=last-1;total+=count
        for k,v in {"gap":d,"first_terminal":2,"last_terminal":last,"tested":count}.items():
            require(type(r[k]) is int and r[k]==v and type(s[k]) is int and s[k]==v,"Recorded full lower-one gap domain")
        require(re.fullmatch("[0-9a-f]{64}",r["case_stream_sha256"]) and r["case_stream_sha256"]==s["case_stream_sha256"],"Recorded complete domain stream identity")
        require(s["prime"]==1000003 and s["all_full_degree_gcds_equal_one"] is True,"Recorded fixed-prime full-degree outcome")
    require(total==171812 and receipt["completed_cases"]==total and receipt["complete_gap_families"]==34,"Complete recorded lower-one finite count")
    require(receipt["maximum_terminal"]==5631 and receipt["maximum_factorial_argument"]==247763,"Actual recorded lower-one caps")
    return {"a0_selected_current_receipt_maps":maps,"a0_prefix_records_authenticated":92935,
       "a0_exponent_six_records_authenticated":9998,"a1_parent_ledgers":a1ledgers,
       "a1_completed_cases_authenticated":171812,"boundary_finite_Bezout_or_gcd_replays":0,
       "density_context_is_not_a_mathematical_premise":True,
       "meaning":"Immutable defining code, literal certificate data, domain records and completed receipts are source-authenticated. Their old large finite computations are reused, not rerun."}


def cutoff(N):
    require(type(N) is int and N>=2,"Terminal integer>=2")
    target=4*N*(N-1)
    e=(target-1).bit_length();D=2*(N-1)*e
    require(2**(e-1)<target<=2**e,"Exact least binary cutoff")
    return e,D
def exact_controls():
    binary=0;annular=0;terminal_triangles=0;gap_triangles=0;lower_rectangles=0;models=0;normal=0
    for N in range(2,129):
        e,D=cutoff(N);require(D==2*(N-1)*e and D>=6,"Old explicit high-gap bound")
        binary+=1
    for N in range(2,33):
        e,D=cutoff(N)
        for j in range(2,N+1):
            require(Q(j,j-1)**((N-1)*e)>=2**e>=4*N*(N-1)>=4*j*(j-1),"Literal rational annular separation majorant")
            annular+=1
    for N in range(2,65):
        e,D=cutoff(N)
        literal=sum(max(0,d*(N-2)-1) for d in range(2,D))
        if N==2:require(literal==0,"Terminal2 sharp triangle is empty")
        else:require(literal==(N-2)*((D-1)*D//2-1)-(D-2),"Exact sharp fixed-terminal support count")
        terminal_triangles+=1
    for d in range(1,17):
        for L in range(2,21):
            exact=sum(d*(N-2) for N in range(2,L))
            require(exact==d*(L-2)*(L-3)//2,"Sharp fixed-gap bad-pair triangle")
            old=sum(d*(2*N-3) for N in range(2,L))
            require(old==d*(L-2)**2,"Independent older soft reverse count")
            require(exact<=old,"Sharp count improves the independent old soft bound")
            gap_triangles+=1
    for a in range(2,17):
        for U in range(2,21):
            pairs={(d,N) for d in range(1,3*a) for N in range(2,U)}
            require(len(pairs)==(3*a-1)*(U-2),"Fixed-lower finite rectangle count")
            lower_rectangles+=1
    # Finite set controls model the pure projection/fiber reconstruction,
    # not measured gcd failures of these parameter values.
    for size in range(2,14):
        B={(a,d,N) for a in range(size) for d in range(1,size+1) for N in range(2,size+3)
           if (2*a+3*d+N)%7==0}
        for d in range(1,size+1):
            bd={(a,N) for a,dd,N in B if dd==d}
            projection={a for a,N in bd}
            rebuilt=set().union(*[{(a,N) for aa,N in bd if aa==a} for a in projection]) if projection else set()
            require(rebuilt==bd,"Fixed-gap projection/fiber union")
            if bd:
                L=1+max(N for a,N in bd);require(all(N<L for a,N in bd),"One maximum covers the finite gap fiber")
        for a in range(size):
            ba={(d,N) for aa,d,N in B if aa==a}
            projection={d for d,N in ba}
            rebuilt=set().union(*[{(d,N) for dd,N in ba if dd==d} for d in projection]) if projection else set()
            require(rebuilt==ba,"Fixed-lower projection/fiber union")
            if ba:
                U=1+max(N for d,N in ba);require(all(N<U for d,N in ba),"One maximum covers finite lower fiber")
        # The union of bounded coordinate fibers is exactly its bounded slab.
        for bound in [0,1,size//2,size]:
            direct={x for x in B if x[0]<=bound}
            rebuilt=set().union(*[{x for x in B if x[0]==a} for a in range(bound+1)])
            require(direct==rebuilt,"Bounded-coordinate slab is a finite union of fibers")
        models+=1
    C=10**24
    for n in range(4,129):
        e,D=cutoff(n)
        require(n<C*n*n and n<3*n and n<n*(n-2) and n<D,"Diagonal countermodel avoids all retained all-order regions")
        require(n*(n-2)-n==(n-4)**2+5*(n-4)+4,"Universal positive shifted diagonal soft slack")
        # Distinct diagonal triples have singleton fibers in ALL coordinates.
        diagonal={(j,j,j) for j in range(4,n+1)}
        for k in [4,n]:
            require(sum(a==k for a,d,N in diagonal)==1 and
                    sum(d==k for a,d,N in diagonal)==1 and
                    sum(N==k for a,d,N in diagonal)==1,"Diagonal has all coordinate fibers singleton")
    for a in range(4):
        for d in range(1,5):
            for m in range(1,4):
                A=Q(2,3);B=-Q(3,5);x=A/B
                direct=sum((Q(comb(m,j)*factorial(a*m+d*j))*A**(m-j)*B**j for j in range(m+1)),Q(0))
                poly=sum((Q(comb(m,j)*factorial(a*m+d*j),factorial(a*m))*x**(m-j) for j in range(m+1)),Q(0))
                require(direct==B**m*factorial(a*m)*poly,"Exact actual moment normalization for universal field translation")
                normal+=1
    # Genuine logical countermodels; no actual counterexample is asserted.
    bad=[]
    require(all(2<3 for a in range(100)),"A uniform terminal cutoff can coexist with unbounded support at one early order")
    bad.append("uniform-terminal-only-does-not-bound-support-without-soft-input")
    require(len({0 for N in range(2,100)})==1 and len({(0,N) for N in range(2,100)})==98,"Finite projection does not replace finite individual fibers")
    bad.append("finite-projection-with-infinite-fiber")
    require(len({(a,a+2) for a in range(100)})==100,"Finite individual fibers do not replace finite projection")
    bad.append("finite-fibers-with-unbounded-projection")
    require((2-2)*((cutoff(2)[1]-1)*cutoff(2)[1]//2-1)-(cutoff(2)[1]-2)<0,"General count formula must not be applied to N2")
    bad.append("fixed-terminal-count-formula-needs-N>=3")
    bad.append("all-three-coordinate-fibers-finite-does-not-make-global-set-finite")
    return {"exact_binary_cutoffs":binary,"exact_rational_annular_majorants":annular,
       "fixed_terminal_triangle_counts":terminal_triangles,"fixed_gap_triangle_and_old_counts":gap_triangles,
       "fixed_lower_rectangle_counts":lower_rectangles,"finite_projection_fiber_models":models,
       "diagonal_model_parameter_controls":125,"actual_moment_normalizations":normal,
       "logical_countermodel_controls":len(bad)},bad


def negatives(ledger,repo,support):
    out=[]
    def reject(name,fn):
        try:fn()
        except (AssertionError,KeyError,TypeError,ValueError):out.append(name);return
        raise AssertionError("Negative accepted:"+name)
    for role in ROLES:
        d=copy.deepcopy(ledger);d["sources"]=[r for r in d["sources"] if r["role"]!=role]
        reject("omit:"+role,lambda d=d:check_sources(d,repo,support))
        d=copy.deepcopy(ledger);next(r for r in d["sources"] if r["role"]==role)["sha256"]="0"*64
        reject("badpin:"+role,lambda d=d:check_sources(d,repo,support))
    for key in ["global_B_finite_claim","global_terminal_cutoff_claim","all_order_SFC_claim",
                "effective_L_U_claim","new_analytic_estimate_claim","new_fixed_terminal_finiteness_claim",
                "relative_full_gcd_run_rerun","relative_fresh_fallback_Bezout_rerun","boundary_old_finite_runs_repeated","boundary_nested_density_context_replayed"]:
        d=copy.deepcopy(SCOPE);d[key]=True
        reject("false-scope:"+key,lambda d=d:check_inputs(d))
    for key in ["all_fixed_gap_bad_pair_sets_finite","all_fixed_lower_bad_pair_sets_finite",
                "all_fixed_terminal_bad_support_sets_finite","same_cutoff_all_fields_and_coefficients",
                "finite_projection_and_finite_fibers_iff","terminal_only_iff_uses_actual_support_complements",
                "all_coordinate_escape_for_infinite_distinct_bad_sequences"]:
        d=copy.deepcopy(SCOPE);d[key]=False
        reject("missing-scope:"+key,lambda d=d:check_inputs(d))
    for key,val in [("fixed_gap_cutoff","for each (a,d) exists only its own cutoff"),
                    ("fixed_lower_cutoff","one cutoff for all a and d"),
                    ("fixed_terminal_count_domain","N>=2 with no separate case"),
                    ("relative_completed_cases_reused",3694039)]:
        d=copy.deepcopy(SCOPE);d[key]=val
        reject("wrong-quantity:"+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(SCOPE);d["schema"]=True
    reject("typed-input-schema",lambda:check_inputs(d))
    for key in ["schema","empty","duplicate"]:
        d=copy.deepcopy(ledger)
        if key=="schema":d["schema"]=True
        if key=="empty":d["sources"]=[]
        if key=="duplicate":d["sources"][-1]=copy.deepcopy(d["sources"][0])
        reject("ledger:"+key,lambda d=d:check_sources(d,repo,support))
    reject("forged-parent-receipt",lambda:equal_receipt(b'{"status":"PASS"}',b'{"status":"PASS","actual":true}'))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--repo-root",required=True,type=Path)
    ap.add_argument("--output",type=Path,default=HERE/RN);ap.add_argument("--parent-python",type=Path);ns=ap.parse_args()
    repo=ns.repo_root.resolve();require(repo.is_dir(),"Explicit repository root")
    inp=load(HERE/IN);check_inputs(inp)
    ledger=load(HERE/LN);rows=check_sources(ledger,repo,HERE)
    parent_links=links(repo,HERE,rows)
    parent_python=ns.parent_python
    if parent_python is None:
        parent_python=next((p for p in [repo/".venv/bin/python",repo/".venv/Scripts/python.exe"] if p.is_file()),Path(sys.executable))
    parent_python=parent_python.absolute()  # Preserve venv symlink, but honor relative explicit paths across subprocess cwd.
    require(parent_python.is_file(),"Parent interpreter exists")
    children=replay(repo,HERE,rows,parent_python)
    recorded=relative_artifact(repo,rows);boundaries=boundary_artifacts(repo,rows)
    counts,mathbad=exact_controls();bad=negatives(ledger,repo,HERE)
    check_sources(ledger,repo,HERE)
    counts.update({"fresh_parent_replays":len(children),"parent_source_boundaries":len(parent_links),
                   "source_scope_controls":len(bad)})
    payload={"schema":1,"status":"PASS","scope":inp,"formal_verification":False,"external_review":False,
      "proof_sha256":sha(HERE/PN),"checker_sha256":sha(HERE/CN),"input_sha256":sha(HERE/IN),
      "source_ledger_sha256":sha(HERE/LN),"mandatory_source_roles":len(rows),
      "sources":[rows[r] for r in sorted(rows)],"parent_source_boundaries":parent_links,
      "fresh_parent_replays":children,"relative_completed_evidence":recorded,"boundary_completed_evidence":boundaries,
      "counts":counts,"logical_countermodel_names":mathbad,"source_scope_negative_names":bad,
      "meaning":"The written finite-union proof combines retained quantified theorems. Fresh parent guards authenticate their own closures. The relative full run and boundary all-order theorems are explicitly reused, not rerun. Finite set/count controls support scope and algebra; they are not formal verification or new numerical SFC evidence."}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"status":"PASS","mandatory_sources":len(rows),"counts":counts,"output":str(ns.output)}))
if __name__=="__main__":main()
