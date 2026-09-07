from sage.all import *
from pathlib import Path

FRAME=Path("elkies-k3/data/fibrations/q15_mw3_frame.txt")
F=matrix(ZZ,[[ZZ(x) for x in ln.split()] for ln in FRAME.read_text().splitlines() if ln.strip()])
assert F.nrows()==17 and F.is_symmetric() and F.det()==948
print(f"MW3NODE|stage=input|rank={F.rank()}|det={F.det()}|positive={F.is_positive_definite()}",flush=True)

def qform(H):
    co=[]
    for i in range(H.nrows()):
        for j in range(i,H.ncols()):
            co.append(H[i,i]//2 if i==j else H[i,j])
    return QuadraticForm(ZZ,H.nrows(),co)

Q=qform(F)
arr=Q.short_vector_list_up_to_length(2,True)
half=[vector(ZZ,r) for r in arr[1]]
roots=half+[-r for r in half]
R=matrix(ZZ,[list(r) for r in half])
L=matrix(ZZ,[list(r) for r in roots]).row_module()
B=L.basis_matrix()
G=B*F*B.transpose()
print(f"MW3NODE|stage=roots|half={len(half)}|signed={len(roots)}|rank={R.rank()}|det={abs(G.det())}",flush=True)

def ade_solutions(rank,nroots,det):
    atoms=[]
    for n in range(1,18): atoms.append((f"A{n}",n,n*(n+1),n+1))
    for n in range(4,18): atoms.append((f"D{n}",n,2*n*(n-1),4))
    for n,nr,d in [(6,72,3),(7,126,2),(8,240,1)]:
        atoms.append((f"E{n}",n,nr,d))
    sols=[]
    def rec(start,rleft,nleft,dleft,cur):
        if rleft==0:
            if nleft==0 and dleft==1: sols.append(tuple(cur))
            return
        for i in range(start,len(atoms)):
            nm,r,nr,d=atoms[i]
            if r<=rleft and nr<=nleft and dleft%d==0:
                rec(i,rleft-r,nleft-nr,dleft//d,cur+[nm])
    rec(0,rank,nroots,det,[])
    return sols

sols=ade_solutions(R.rank(),len(roots),abs(G.det()))
print(f"MW3NODE|stage=ADE|solutions={sols}",flush=True)

U=matrix(ZZ,[[0,1],[1,0]])
NS=block_diagonal_matrix(U,-F)
rootdet=abs(G.det())
mw_rank=17-R.rank()
reg_t1=QQ(abs(NS.det()))/QQ(rootdet)
print(f"MW3NODE|stage=shioda|MW_rank={mw_rank}|rootdet={rootdet}|regulator_if_torsion1={reg_t1}",flush=True)

out=Path("artifacts/local/elkies-k3/mw3-node.txt")
out.parent.mkdir(parents=True,exist_ok=True)
out.write_text(
    f"ADE={sols}\nMW_rank={mw_rank}\nrootdet={rootdet}\nregulator_t1={reg_t1}\n"
    f"root_basis=\n{B}\nroot_gram=\n{G}\n"
)
print(f"MW3NODE|stage=done|out={out}",flush=True)
