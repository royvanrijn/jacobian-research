"""Weighted all-scalar obstructions for primitive quadratic wave pairs."""
import sympy as S
def vertical_sector(weight):
    # Weight(w)=2, weight(z)=3.
    basis=[(i,j) for j in range(weight//3+1) for i in range(weight//2+1) if 2*i+3*j==weight]
    basis.sort(key=lambda t:t[1])
    matrix=S.zeros(len(basis))
    for row,(i,j) in enumerate(basis):
        matrix[row,row]=2*(i+2)*(j+1)
        if row:
            matrix[row,row-1]=(i+3)*(i+2)
        if row+1<len(basis):
            matrix[row,row+1]=(j+1)*(j+2)
    return basis,matrix
def curved_sector(weight,k):
    # Weight(w)=1, weight(z)=2, v=(w^2+z,k*w*z).
    last=weight//2
    matrix=S.zeros(last+2,last+1)
    for j in range(last+1):
        i=weight-2*j
        matrix[j,j]=(i+4)*(i+3)+2*k*(i+3)*(j+1)+k*k*(j+1)*(j+2)
        matrix[j+1,j]=2*(i+1)*(i+2+k*(j+2))
        if j<last:matrix[j+2,j]=i*(i-1)
    return matrix
def bezout_certificate(weight,k):
    matrix=curved_sector(weight,k)
    top=S.Poly(matrix[:-1,:].det(method="domain-ge"),k,domain=S.QQ).primitive()[1]
    bottom=S.Poly(matrix[1:,:].det(method="domain-ge"),k,domain=S.QQ).primitive()[1]
    U,V,G=S.gcdex(top,bottom)
    return dict(weight=weight,rows=matrix.rows,columns=matrix.cols,
      top=top,bottom=bottom,U=U,V=V,G=G)
