"""Mixed polynomial logarithmic boundary matrices and two-block gates."""
import sympy as S
from boundary_multiplier import weighted_matrix
def mixed_matrix(Q,ell,xs,e,J):
    row=S.Matrix([[0,*[S.diff(J,v)/J for v in xs]]])
    v=S.Matrix([(S.Matrix([list(ell)])*S.Matrix(xs))[0],*list(Q)])
    return weighted_matrix(Q,ell,xs,e)+v*row
def binary_matrix(R,x,y,e,j):
    R=S.Matrix(R)
    C=S.Matrix([[(e+1)*x,1,0],[e*R[0],S.diff(R[0],x),S.diff(R[0],y)],
                [e*R[1],S.diff(R[1],x),S.diff(R[1],y)]])
    return C+S.Matrix([x,*list(R)])*S.Matrix([[0,S.diff(j,x)/j,S.diff(j,y)/j]])
def two_block(R,x,y,e,t,eta):
    return S.Matrix([[(e+1)*t,t*eta],[e*R,S.diff(R,y)+R*eta]])
