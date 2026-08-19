from pathlib import Path
import numpy as np

BASE=Path(__file__).resolve().parents[1]
J=BASE/"checkpoints/24A1-JACKPOT-948-r0-n1311"

V=np.loadtxt(J/"all_1311_in_short_basis.txt",dtype=np.int64)
H=np.loadtxt(J/"short_vector_basis_gram.txt",dtype=np.int64)

P=(V@H)@V.T
n=len(V)

for rel in [-2,-1,0,1,2]:
    A=(P==rel).astype(float)
    np.fill_diagonal(A,0)

    degrees=A.sum(axis=1)
    eig=np.linalg.eigvalsh(A)

    # cluster numerical eigenvalues
    vals=[]
    for x in eig:
        if not vals or abs(x-vals[-1][0])>1e-6:
            vals.append([x,1])
        else:
            vals[-1][1]+=1

    print("\nRELATION",rel)
    print("degree values:",sorted(set(map(int,degrees))))
    print("number eigenvalues:",len(vals))
    print("spectrum:")
    for x,m in vals:
        print(f"  {x:.10f} x {m}")
