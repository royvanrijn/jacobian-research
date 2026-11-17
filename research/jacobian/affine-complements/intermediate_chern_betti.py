"""Exact small controls for PF124; geometry and all-degree bounds are written."""


def affine_betti(b,cusps):
    assert b>=4 and b&(b-1)==0 and cusps in (0,1)
    r=b//4
    return dict(degree=b,critical_curves=r,cusps=cusps,chi=r*(3-cusps),b2=r*(3-cusps)-1)


def pullback_matrix(packets):
    rows=[]
    for j,packet in enumerate(packets):
        for k in packet:
            assert k>=1
            row=[0]*len(packets);row[j]=k;rows.append(row)
    return rows
