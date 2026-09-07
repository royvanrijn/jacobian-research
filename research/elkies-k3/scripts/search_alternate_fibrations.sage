from sage.all import *
from pathlib import Path
import argparse, math

ap=argparse.ArgumentParser(description="Search alternate elliptic fibrations (U embeddings) in NS=U+(-M).")
ap.add_argument("--frame",type=Path,
                default=Path("elkies-k3/data/lattice/rank17_gram.txt"),
                help="positive frame Gram to search from")
ap.add_argument("--min-qnorm",type=int,default=2)
ap.add_argument("--max-qnorm",type=int,default=20,help="enumerate M quadratic norms q(v)=M(v,v)/2 up to this")
ap.add_argument("--max-candidates",type=int,default=200)
ap.add_argument("--report",type=int,default=20)
ap.add_argument("--enum-seed",type=int,default=0)
ap.add_argument("--enum-restarts",type=int,default=1)
ap.add_argument("--enum-cap",type=int,default=750)
ap.add_argument("--enum-baseline-cap",type=int,default=4000)
ap.add_argument("--proper-factors-only",action="store_true",
                help="skip the a=1 or b=1 presentations, whose frame is the original rootless lattice")
ap.add_argument("--per-root-data-cap",type=int,default=0,
                help="retain at most this many frames for each (root rank, count, determinant); zero keeps all")
ap.add_argument("--quiet-candidates",action="store_true",
                help="suppress the per-candidate progress lines")
ap.add_argument("--root-method",choices=("pari","sage"),default="pari",
                help="exact norm-2 enumeration backend")
ap.add_argument("--rank-order", choices=("low", "high"), default="low",
                help="report and retain low-MW or high-MW frames first")
ap.add_argument("--one-factor-order",action="store_true",
                help="search only (a,b) with a<=b; swapping the U coordinates gives an isometric frame")
ap.add_argument("--fixed-coordinate",action="append",default=[],metavar="INDEX:VALUE",
                help="restrict original frame coordinates during norm enumeration; repeatable")
ap.add_argument("--require-outside-root-span",action="store_true",
                help="count sampled norm-shell vectors only when they are outside the source root lattice")
ap.add_argument("--out",type=Path,
                default=Path("artifacts/local/elkies-k3/alternate-fibration-search.txt"))
ap.add_argument("--frames-dir",type=Path,
                help="optionally write the reported frame Grams and a hits.tsv manifest")
a=ap.parse_args()

fixed_coordinates={}
for item in a.fixed_coordinate:
    index_text,value_text=item.split(":",1)
    fixed_coordinates[ZZ(index_text)]=ZZ(value_text)

Mpath=a.frame
M=matrix(ZZ,[[ZZ(x) for x in ln.split()] for ln in Mpath.read_text().splitlines()
             if ln.strip() and not ln.startswith("#")])
assert M.is_square() and M.is_positive_definite()
frame_rank=M.nrows()
ns_rank=frame_rank+2
assert all(0<=index<M.nrows() for index in fixed_coordinates)
U=matrix(ZZ,[[0,1],[1,0]])
NS=block_diagonal_matrix(U,-M)
# Since M is exactly positive definite and U has signature (1,1), the block
# form has signature (1,rank(M)+1).  Avoid a redundant floating eigenvalue
# conversion here: for ill-conditioned but exact Gram presentations Sage can
# return tiny spurious imaginary parts and reject the conversion to RR.
print(f"FIBSEARCH|stage=start|rankM={frame_rank}|detM={M.det()}|max_qnorm={a.max_qnorm}",flush=True)

source_root_module=None
if a.require_outside_root_span:
    source_minimum=pari(M).qfminim(2)
    source_roots=matrix(ZZ,source_minimum[2]).transpose()
    source_root_module=source_roots.row_module()
    print(f"FIBSEARCH|stage=source_roots|rank={source_root_module.rank()}|roots={source_minimum[0]}|filter=outside",flush=True)

# Sage QuadraticForm uses q(v)=1/2 v^T M v for even Hessian M.
co=[]
for i in range(frame_rank):
    for j in range(i,frame_rank):
        co.append(M[i,i]//2 if i==j else M[i,j])
QM=QuadraticForm(ZZ,frame_rank,co)

# Representatives up to sign, grouped by q(v).
def vectors_of_qnorm(n):
    # Sample an exact norm shell from several unimodularly equivalent bases.
    # This remains a DISCOVERY sampler unless the traversal completes.
    target = ZZ(2*n)
    print(f"FIBSEARCH|stage=enumerate|q={n}|target_norm={target}|method=multistart_fincke_pohst|restarts={a.enum_restarts}", flush=True)

    all_found = {}
    exhaustive = False
    for restart in range(a.enum_restarts):
        # Restart 0 reproduces the original deterministic scan so multistart
        # can never be worse than the old search. Later restarts are cheap
        # alternative traversals.
        cap_per_restart = max(1, a.enum_baseline_cap if restart == 0 else a.enum_cap)
        set_random_seed(a.enum_seed + 1000003*restart + 97*n)

        # Random unimodular change of basis. y in transformed coordinates maps
        # to the original row vector v = y*T.
        T = identity_matrix(ZZ, M.nrows())
        perm = list(range(M.nrows()))
        if restart > 0:
            # Permuting coordinates changes the Cholesky/enumeration order
            # without destroying conditioning. Random shears were much slower.
            shuffle(perm)
            T = matrix(ZZ, M.nrows(), M.nrows(),
                       lambda i,j: ZZ(1) if j==perm[i] else ZZ(0))

        # Every transformation used here is a coordinate permutation.  Since
        # v=y*T, an original constraint v[j]=c becomes y[perm^-1(j)]=c.
        fixed_transformed={perm.index(index):value
                           for index,value in fixed_coordinates.items()}

        Ms = T*M*T.transpose()
        MR = Ms.change_ring(RR)
        Rchol = MR.cholesky().transpose()
        dim = Ms.nrows()
        x = [ZZ(0)]*dim
        found_this = 0
        nodes = 0
        eps = RR("1e-10")

        def rec(k, partial):
            nonlocal found_this, nodes
            nodes += 1
            if found_this >= cap_per_restart:
                return
            if k < 0:
                y = vector(ZZ,x)
                if y == 0 or ZZ(y*Ms*y) != target:
                    return
                v = y*T
                if ZZ(v*M*v) != target:
                    raise RuntimeError("basis-map norm mismatch")
                if source_root_module is not None and v in source_root_module:
                    return
                tup=tuple(v)
                neg=tuple(-v)
                canon=min(tup,neg)
                if canon not in all_found:
                    all_found[canon]=vector(ZZ,canon)
                    found_this += 1
                return

            off=sum(Rchol[k,j]*RR(x[j]) for j in range(k+1,dim))
            rem=RR(target)-partial
            if rem < -eps:
                return
            rad=sqrt(max(RR(0),rem))
            diag=Rchol[k,k]
            lo=ceil((-rad-off)/diag-eps)
            hi=floor((rad-off)/diag+eps)

            center=ZZ(round(-off/diag))
            vals=[]
            if k in fixed_transformed:
                fixed_value=fixed_transformed[k]
                if lo<=fixed_value<=hi:
                    vals.append(fixed_value)
            else:
                if lo<=center<=hi:
                    vals.append(center)
                d=1
                while center-d>=lo or center+d<=hi:
                    if center-d>=lo:
                        vals.append(center-d)
                    if center+d<=hi:
                        vals.append(center+d)
                    d+=1

            for z in vals:
                x[k]=ZZ(z)
                yy=diag*RR(z)+off
                rec(k-1,partial+yy*yy)
                if found_this>=cap_per_restart:
                    break
            x[k]=ZZ(0)

        rec(dim-1,RR(0))
        # Reaching the cap may truncate the traversal.  Conversely, returning
        # with fewer hits than the cap certifies that this coordinate-ordered
        # Fincke--Pohst traversal visited the complete constrained shell.
        restart_exhaustive = found_this < cap_per_restart
        exhaustive = exhaustive or restart_exhaustive
        print(f"FIBSEARCH|stage=restart_done|q={n}|restart={restart}|new_pairs={found_this}|union_pairs={len(all_found)}|nodes={nodes}|exhaustive={str(restart_exhaustive).lower()}",flush=True)

    vv=list(all_found.values())
    print(f"FIBSEARCH|stage=enumerate_done|q={n}|pairs={len(vv)}|method=multistart_fincke_pohst|exhaustive={str(exhaustive).lower()}", flush=True)
    return vv

def divisors_nontrivial(n):
    # isotropic f=(aa,bb,v): aa*bb=q(v).
    # IMPORTANT: aa=1 or bb=1 are valid primitive isotropic classes and
    # need not be equivalent to the standard fiber in this one-U lattice.
    out=[]
    for aa in divisors(n):
        bb=n//aa
        if aa>=1 and bb>=1 and (not a.proper_factors_only or (aa>1 and bb>1)):
            out.append((ZZ(aa),ZZ(bb)))
    return out

def primitive_div(v):
    pair=NS*v
    return gcd([abs(ZZ(x)) for x in pair])

def bezout_vector_for_pairing(f):
    p=list(NS*f)
    g=0; coeff=[]
    # iterative xgcd for coefficients c with sum c_i p_i = gcd
    cs=[]
    cur=ZZ(0)
    vec=[ZZ(0)]*ns_rank
    for i,pi in enumerate(p):
        if pi==0: continue
        gg,s,t=xgcd(cur,ZZ(pi))
        vec=[s*x for x in vec]
        vec[i]+=t
        cur=gg
    if abs(cur)!=1:
        return None
    if cur==-1: vec=[-x for x in vec]
    return vector(ZZ,vec)

def build_frame(f):
    g=bezout_vector_for_pairing(f)
    if g is None: return None
    assert f*NS*g==1
    gsq=ZZ(g*NS*g)
    assert gsq%2==0
    # Replace g by isotropic mate.
    g0=g-(gsq//2)*f
    assert g0*NS*g0==0 and f*NS*g0==1
    B=matrix(ZZ,[list(f*NS),list(g0*NS)])
    K=B.right_kernel_matrix()
    assert K.rank()==frame_rank
    G=K*NS*K.transpose()
    # should be negative definite
    P=-G
    if not P.is_positive_definite():
        return None
    return g0,K,P

def root_data(P):
    if a.root_method=="pari":
        # qfminim returns one vector from each +/- pair as columns, while its
        # first entry counts both signs.  This is exact and substantially
        # faster than rebuilding a Sage QuadraticForm for every candidate.
        z=pari(P).qfminim(2)
        count=ZZ(z[0])
        if count==0:
            return 0,0,1
        R=matrix(ZZ,z[2]).transpose()
        rank=R.rank()
        B=R.row_module().basis_matrix()
        Gram=B*P*B.transpose()
        return rank,count,abs(Gram.det())
    co=[]
    for i in range(frame_rank):
        for j in range(i,frame_rank):
            co.append(P[i,i]//2 if i==j else P[i,j])
    Q=QuadraticForm(ZZ,frame_rank,co)
    arr=Q.short_vector_list_up_to_length(2,True)
    roots_half=arr[1] if len(arr)>1 else []
    if not roots_half:
        return 0,0,1
    R=matrix(ZZ,[list(r) for r in roots_half])
    rank=R.rank()
    L=R.row_module()
    B=L.basis_matrix()
    Gram=B*P*B.transpose()
    return rank,2*len(roots_half),abs(Gram.det())

seen=set()
results=[]
root_data_counts={}
tested=0

for qnorm in range(a.min_qnorm,a.max_qnorm+1):
    pairs=divisors_nontrivial(qnorm)
    if not pairs:
        continue
    try:
        qvectors=vectors_of_qnorm(qnorm)
    except Exception as e:
        print(f"FIBSEARCH|stage=enumerate_failed|q={qnorm}|error={e}", flush=True)
        continue
    for vv in qvectors:
        v=vector(ZZ,vv)
        for aa,bb in pairs:
            if a.proper_factors_only and aa>bb:
                continue
            swaps=[False] if a.one_factor_order else ([False,True] if aa!=bb else [False])
            for swap in swaps:
                A,B=(bb,aa) if swap else (aa,bb)
                f=vector(ZZ,[A,B]+list(v))
                if f*NS*f!=0 or primitive_div(f)!=1: continue
                built=build_frame(f)
                if built is None: continue
                g0,K,P=built
                # cheap invariant hash using determinant + diagonal sorted + minimum if available.
                key=(tuple(sorted(P.diagonal())),)
                # Do not rely on this for proof; just avoid obvious repeats.
                tested+=1
                rr,rc,rd=root_data(P)
                mw=frame_rank-rr
                root_key=(rr,rc,rd)
                if a.per_root_data_cap and root_data_counts.get(root_key,0)>=a.per_root_data_cap:
                    continue
                # Deduplicate aggressively. For discovery we care about inequivalent-looking
                # frame/root data, not repeated isotropic representatives.
                try:
                    Pred = P.LLL_gram()
                except Exception:
                    Pred = P
                frame_key = tuple(Pred.list())
                key = (rr,rc,rd,frame_key)
                if key in seen:
                    continue
                seen.add(key)
                root_data_counts[root_key]=root_data_counts.get(root_key,0)+1

                row=(mw,rr,rc,rd,qnorm,A,B,tuple(v),P,g0)
                results.append(row)
                if not a.quiet_candidates:
                    print(f"FIBSEARCH|cand={len(results)}|q={qnorm}|ab={A},{B}|root_rank={rr}|roots={rc}|rootdet={rd}|MW={mw}",flush=True)

if a.rank_order == "low":
    results.sort(key=lambda x:(x[0],-x[2],x[4]))
else:
    results.sort(key=lambda x:(-x[0],x[2],x[4]))
print(f"FIBSEARCH|stage=summary|tested={tested}|unique_frames={len(results)}|sampled_vectors={len(qvectors) if 'qvectors' in locals() else -1}",flush=True)
for i,row in enumerate(results[:a.report],1):
    mw,rr,rc,rd,qnorm,A,B,v,P,g0=row
    print(f"FIBBEST|rank={i}|MW={mw}|root_rank={rr}|roots={rc}|rootdet={rd}|q={qnorm}|ab={A},{B}|v={v}",flush=True)
    print(f"FIBBEST|frame_gram={P}",flush=True)

out=a.out
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    for i,row in enumerate(results[:a.report],1):
        mw,rr,rc,rd,qnorm,A,B,v,P,g0=row
        h.write(f"rank={i} MW={mw} root_rank={rr} roots={rc} rootdet={rd} q={qnorm} ab={A},{B} v={v}\n")
        h.write("frame=\n"+str(P)+"\n\n")
if a.frames_dir is not None:
    a.frames_dir.mkdir(parents=True,exist_ok=True)
    with (a.frames_dir / "hits.tsv").open("w") as h:
        h.write("rank\tMW\troot_rank\troots\trootdet\tq\ta\tb\tv\tframe\tparent_frame\n")
        for i,row in enumerate(results[:a.report],1):
            mw,rr,rc,rd,qnorm,A,B,v,P,g0=row
            frame_name=f"frame-{i:03d}.txt"
            (a.frames_dir / frame_name).write_text(
                "\n".join(" ".join(str(P[r,c]) for c in range(P.ncols()))
                          for r in range(P.nrows()))+"\n"
            )
            h.write(f"{i}\t{mw}\t{rr}\t{rc}\t{rd}\t{qnorm}\t{A}\t{B}\t{v}\t{frame_name}\t{Mpath}\n")
print(f"FIBSEARCH|stage=done|out={out}",flush=True)
