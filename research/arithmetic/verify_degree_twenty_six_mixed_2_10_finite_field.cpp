// Independent complete four-coordinate census for the mixed (2,10) system.
// Reads prime, then four integer-polynomial coefficient/exponent tables.
// Mode 0 uses all terms; mode 1 selects weighted-leading terms independently.
// Every c is tested explicitly; this does not use the Python linear solver.
#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
struct Term { int64_t coefficient; array<int,4> exponent; };
int reduce(int64_t x,int p) { int y=int(x%p); return y<0?y+p:y; }
array<int,16> powers(int x,int p) {
    array<int,16> result; result[0]=1;
    for(int j=1;j<16;j++) result[j]=result[j-1]*x%p;
    return result;
}
int weight(const Term& term) {
    return term.exponent[0]+2*term.exponent[1]+term.exponent[2]+term.exponent[3];
}
int main() {
    int p; if(!(cin>>p)||p!=61) return 2;
    vector<vector<Term>> input(4);
    for(auto& polynomial:input) {
        int size; cin>>size;
        for(int j=0;j<size;j++) {
            Term term; cin>>term.coefficient;
            for(int& exponent:term.exponent) {
                cin>>exponent; if(exponent<0||exponent>=16) return 3;
            }
            polynomial.push_back(term);
        }
    }
    const array<int,4> degrees={7,6,11,12};
    for(int mode=0;mode<2;mode++) {
        vector<vector<Term>> polys(4);
        for(int j=0;j<4;j++) {
            int maximum=0;
            for(auto term:input[j]) maximum=max(maximum,weight(term));
            if(maximum!=degrees[j]) return 4;
            for(auto term:input[j])
                if(mode==0||weight(term)==maximum) polys[j].push_back(term);
        }
        int64_t candidates=0;
        vector<array<int,4>> points;
        for(int av=0;av<p;av++) for(int bv=0;bv<p;bv++) for(int rv=0;rv<p;rv++) {
            auto A=powers(av,p),B=powers(bv,p),R=powers(rv,p);
            int coefficient[2][2]={{0,0},{0,0}};
            for(int j=0;j<2;j++) for(auto term:polys[j]) {
                if(term.exponent[2]>1||term.exponent[3]>1) return 5;
                int value=reduce(term.coefficient,p);
                value=value*A[term.exponent[0]]%p;
                value=value*B[term.exponent[1]]%p;
                value=value*R[term.exponent[3]]%p;
                int index=term.exponent[2];
                coefficient[j][index]=(coefficient[j][index]+value)%p;
            }
            // Literal full c-loop supplies an independent completeness control.
            for(int cv=0;cv<p;cv++) {
                if((coefficient[0][1]*cv+coefficient[0][0])%p||
                   (coefficient[1][1]*cv+coefficient[1][0])%p) continue;
                candidates++;
                auto C=powers(cv,p); bool zero=true;
                for(int j=2;j<4;j++) {
                    int sum=0;
                    for(auto term:polys[j]) {
                        int value=reduce(term.coefficient,p);
                        value=value*A[term.exponent[0]]%p;
                        value=value*B[term.exponent[1]]%p;
                        value=value*C[term.exponent[2]]%p;
                        value=value*R[term.exponent[3]]%p;
                        sum=(sum+value)%p;
                    }
                    if(sum) { zero=false; break; }
                }
                if(zero) points.push_back({av,bv,cv,rv});
            }
        }
        cout<<"{\"mode\":"<<mode<<",\"first_two_candidates\":"<<candidates
            <<",\"points\":[";
        for(size_t i=0;i<points.size();i++) {
            if(i) cout<<",";
            cout<<"["<<points[i][0]<<","<<points[i][1]<<","
                <<points[i][2]<<","<<points[i][3]<<"]";
        }
        cout<<"]}"<<endl;
    }
}
