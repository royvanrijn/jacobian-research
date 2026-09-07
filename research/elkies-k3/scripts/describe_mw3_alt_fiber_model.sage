from sage.all import *
import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--branch",choices=["e6","a6"],required=True)
args=ap.parse_args()

if args.branch=="e6":
    print("BRANCHMODEL|name=e6|ADE=E6+A3+A3+A1+A1")
    print("BRANCHMODEL|preferred=IV*+I4+I4+I2+I2+4I1")
    print("BRANCHMODEL|euler=8+4+4+2+2+4=24")
    print("BRANCHMODEL|normalization=IV*@infinity,I4@0,I4@1,I2@lambda,I2@mu")
    print("BRANCHMODEL|benefit=IV* at infinity forces deg(A)<=5,deg(B)<=8 directly; no I11 cancellation staircase")
else:
    print("BRANCHMODEL|name=a6|ADE=A6+A4+A1^4")
    print("BRANCHMODEL|preferred=I7+I5+I2+I2+I2+I2+4I1")
    print("BRANCHMODEL|euler=7+5+2+2+2+2+4=24")
    print("BRANCHMODEL|normalization=I7@infinity,I5@0,I2@1,I2@lambda,I2@mu,I2@nu")
    print("BRANCHMODEL|benefit=all multiplicative; component-label machinery from A10 branch transfers directly")
