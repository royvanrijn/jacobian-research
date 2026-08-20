from sage.all import *
from pathlib import Path
import argparse, math

ap=argparse.ArgumentParser(description="Search alternate elliptic fibrations (U embeddings) in NS=U+(-M17).")
ap.add_argument("--min-qnorm",type=int,default=2)
ap.add_argument("--max-qnorm",type=int,default=20,help="enumerate M quadratic norms q(v)=M(v,v)/2 up to this")
ap.add_argument("--max-candidates",type=int,default=200)
ap.add_argument("--report",type=int,default=20)
ap.add_argument("--enum-seed",type=int,default=0)
ap.add_argument("--enum-restarts",type=int,default=1)
ap.add_argument("--enum-cap",type=int,default=750)
ap.add_argument("--enum-baseline-cap",type=int,default=4000)
a=ap.parse_args()

Mpath=Path("elkies-k3/data/lattice/rank17_gram.txt")
M=matrix(ZZ,[[ZZ(x) for x in ln.split()] for ln in Mpath.read_text().splitlines()
             if ln.strip() and not ln.startswith("#")])
assert M.nrows()==17 and M.det()==948
U=matrix(ZZ,[[0,1],[1,0]])
NS=block_diagonal_matrix(U,-M)
ev=[RR(x) for x in NS.change_ring(RR).eigenvalues()]; assert sum(x>0 for x in ev)==1 and sum(x<0 for x in ev)==18
print(f"FIBSEARCH|stage=start|detM={M.det()}|max_qnorm={a.max_qnorm}",flush=True)

# Sage QuadraticForm uses q(v)=1/2 v^T M v for even Hessian M.
co=[]
for i in range(17):
    for j in range(i,17):
        co.append(M[i,i]//2 if i==j else M[i,j])
QM=QuadraticForm(ZZ,17,co)

# Representatives up to sign, grouped by q(v).
def vectors_of_qnorm(n):
    # Sample an exact norm shell from several unimodularly equivalent bases.
    # This remains a DISCOVERY sampler unless the traversal completes.
    target = ZZ(2*n)
    print(f"FIBSEARCH|stage=enumerate|q={n}|target_norm={target}|method=multistart_fincke_pohst|restarts={a.enum_restarts}", flush=True)

    all_found = {}
    for restart in range(a.enum_restarts):
        # Restart 0 reproduces the original deterministic scan so multistart
        # can never be worse than the old search. Later restarts are cheap
        # alternative traversals.
        cap_per_restart = max(1, a.enum_baseline_cap if restart == 0 else a.enum_cap)
        set_random_seed(a.enum_seed + 1000003*restart + 97*n)

        # Random unimodular change of basis. y in transformed coordinates maps
        # to the original row vector v = y*T.
        T = identity_matrix(ZZ, M.nrows())
        if restart > 0:
            # Permuting coordinates changes the Cholesky/enumeration order
            # without destroying conditioning. Random shears were much slower.
            perm = list(range(M.nrows()))
            shuffle(perm)
            T = matrix(ZZ, M.nrows(), M.nrows(),
                       lambda i,j: ZZ(1) if j==perm[i] else ZZ(0))

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
        print(f"FIBSEARCH|stage=restart_done|q={n}|restart={restart}|new_pairs={found_this}|union_pairs={len(all_found)}|nodes={nodes}",flush=True)

    vv=list(all_found.values())
    print(f"FIBSEARCH|stage=enumerate_done|q={n}|pairs={len(vv)}|method=multistart_fincke_pohst|exhaustive=false", flush=True)
    return vv

def divisors_nontrivial(n):
    # isotropic f=(aa,bb,v): aa*bb=q(v).
    # IMPORTANT: aa=1 or bb=1 are valid primitive isotropic classes and
    # need not be equivalent to the standard fiber in this one-U lattice.
    out=[]
    for aa in divisors(n):
        bb=n//aa
        if aa>=1 and bb>=1:
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
    vec=[ZZ(0)]*19
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
    assert K.rank()==17
    G=K*NS*K.transpose()
    # should be negative definite
    P=-G
    if not P.is_positive_definite():
        return None
    return g0,K,P

def root_data(P):
    co=[]
    for i in range(17):
        for j in range(i,17):
            co.append(P[i,i]//2 if i==j else P[i,j])
    Q=QuadraticForm(ZZ,17,co)
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
            for swap in [False,True] if aa!=bb else [False]:
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
                mw=17-rr
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

                row=(mw,rr,rc,rd,qnorm,A,B,tuple(v),P,g0)
                results.append(row)
                print(f"FIBSEARCH|cand={len(results)}|q={qnorm}|ab={A},{B}|root_rank={rr}|roots={rc}|rootdet={rd}|MW={mw}",flush=True)

results.sort(key=lambda x:(x[0],-x[2],x[4]))
print(f"FIBSEARCH|stage=summary|tested={tested}|unique_frames={len(results)}|sampled_vectors={len(qvectors) if 'qvectors' in locals() else -1}",flush=True)
for i,row in enumerate(results[:a.report],1):
    mw,rr,rc,rd,qnorm,A,B,v,P,g0=row
    print(f"FIBBEST|rank={i}|MW={mw}|root_rank={rr}|roots={rc}|rootdet={rd}|q={qnorm}|ab={A},{B}|v={v}",flush=True)
    print(f"FIBBEST|frame_gram={P}",flush=True)

out=Path("artifacts/local/elkies-k3/alternate-fibration-search.txt")
out.parent.mkdir(parents=True,exist_ok=True)
with out.open("w") as h:
    for i,row in enumerate(results[:a.report],1):
        mw,rr,rc,rd,qnorm,A,B,v,P,g0=row
        h.write(f"rank={i} MW={mw} root_rank={rr} roots={rc} rootdet={rd} q={qnorm} ab={A},{B} v={v}\n")
        h.write("frame=\n"+str(P)+"\n\n")
print(f"FIBSEARCH|stage=done|out={out}",flush=True)
