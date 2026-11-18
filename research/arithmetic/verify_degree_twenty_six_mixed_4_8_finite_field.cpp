// Complete F29 census for the normalized mixed (4,8) parameter equations.
// The first four equations are affine linear in c,d,r,cr,dr at fixed V.
// Every affine solution of those linear equations is enumerated by RREF;
// the two product identities recover exactly the original triples c,d,r.
#include <array>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <vector>
#include <utility>
using namespace std;
constexpr int prime=29;
struct Term { int64_t coefficient; array<unsigned,7> exponent; };
using Polynomial=vector<Term>;
int mod(int64_t x) { int a=int(x%prime); return a<0?a+prime:a; }
void points_json(vector<array<int,7>> const& points) {
    cout<<"[";
    for(size_t i=0;i<points.size();i++) {
        if(i) cout<<",";
        cout<<"[";
        for(unsigned j=0;j<7;j++) { if(j) cout<<","; cout<<points[i][j]; }
        cout<<"]";
    }
    cout<<"]";
}
int main() {
    int supplied_prime=0; cin>>supplied_prime;
    if(supplied_prime!=prime) return 2;
    vector<Polynomial> raw(7),strong(4);
    for(auto collection:{&raw,&strong}) for(auto& polynomial:*collection) {
        unsigned size; if(!(cin>>size)) return 3;
        for(unsigned i=0;i<size;i++) {
            Term t; if(!(cin>>t.coefficient)) return 4;
            for(auto& e:t.exponent) if(!(cin>>e)||e>=32) return 5;
            t.coefficient=mod(t.coefficient); polynomial.push_back(t);
        }
    }
    for(unsigned i=0;i<4;i++) for(auto const& t:raw[i])
        if(t.exponent[4]+t.exponent[5]>1||t.exponent[6]>1) return 6;
    array<int,prime> inverse{};
    for(int i=1;i<prime;i++) for(int j=1;j<prime;j++) if(i*j%prime==1) inverse[i]=j;
    array<array<int,32>,prime> power{};
    for(int x=0;x<prime;x++) {
        power[x][0]=1;
        for(unsigned j=1;j<32;j++) power[x][j]=mod(power[x][j-1]*x);
    }
    int64_t outer_count=0,linear_candidates=0,raw_count=0,strong_count=0;
    array<int64_t,5> consistent_rank_count{};
    vector<array<int,7>> raw_points,strong_points;
    array<int,7> x{};
    auto evaluate=[&](Polynomial const& f) {
        int sum=0;
        for(auto const& t:f) {
            int64_t value=t.coefficient;
            for(unsigned j=0;j<7;j++) if(t.exponent[j])
                value=mod(value*power[x[j]][t.exponent[j]]);
            sum=mod(sum+value);
        }
        return sum;
    };
    for(x[0]=0;x[0]<prime;x[0]++) for(x[1]=0;x[1]<prime;x[1]++)
    for(x[2]=0;x[2]<prime;x[2]++) for(x[3]=0;x[3]<prime;x[3]++) {
        ++outer_count;
        array<array<int,6>,4> matrix{};
        for(unsigned row=0;row<4;row++) for(auto const& t:raw[row]) {
            unsigned column;
            if(t.exponent[4]) column=t.exponent[6]?3:0;
            else if(t.exponent[5]) column=t.exponent[6]?4:1;
            else column=t.exponent[6]?2:5;
            int64_t value=t.coefficient;
            for(unsigned j=0;j<4;j++) if(t.exponent[j])
                value=mod(value*power[x[j]][t.exponent[j]]);
            matrix[row][column]=mod(matrix[row][column]+value);
        }
        // Unknowns are (c,d,r,cr,dr); the last column is the right side.
        for(auto& row:matrix) row[5]=mod(-row[5]);
        array<int,4> pivot{}; unsigned rank=0;
        for(unsigned column=0;column<5&&rank<4;column++) {
            unsigned found=rank;
            while(found<4&&!matrix[found][column]) ++found;
            if(found==4) continue;
            swap(matrix[rank],matrix[found]);
            int scale=inverse[matrix[rank][column]];
            for(unsigned j=column;j<6;j++) matrix[rank][j]=mod(matrix[rank][j]*scale);
            for(unsigned i=0;i<4;i++) if(i!=rank&&matrix[i][column]) {
                int multiple=matrix[i][column];
                for(unsigned j=column;j<6;j++) matrix[i][j]=mod(matrix[i][j]-multiple*matrix[rank][j]);
            }
            pivot[rank++]=column;
        }
        bool consistent=true;
        for(unsigned i=rank;i<4;i++) if(matrix[i][5]) consistent=false;
        if(!consistent) continue;
        ++consistent_rank_count[rank];
        array<bool,5> is_pivot{};
        for(unsigned i=0;i<rank;i++) is_pivot[pivot[i]]=true;
        vector<unsigned> free_columns;
        for(unsigned j=0;j<5;j++) if(!is_pivot[j]) free_columns.push_back(j);
        array<int,5> y{};
        function<void(unsigned)> enumerate=[&](unsigned depth) {
            if(depth<free_columns.size()) {
                for(int value=0;value<prime;value++) {
                    y[free_columns[depth]]=value; enumerate(depth+1);
                }
                return;
            }
            for(unsigned i=0;i<rank;i++) {
                int value=matrix[i][5];
                for(unsigned j:free_columns) value=mod(value-matrix[i][j]*y[j]);
                y[pivot[i]]=value;
            }
            if(y[3]!=mod(y[0]*y[2])||y[4]!=mod(y[1]*y[2])) return;
            x[4]=y[0]; x[5]=y[1]; x[6]=y[2];
            ++linear_candidates;
            // Redundant direct evaluation guards the affine-linear encoding.
            for(unsigned j=0;j<4;j++) if(evaluate(raw[j])) { cerr<<"linear encoding mismatch\n"; exit(7); }
            for(unsigned j=4;j<7;j++) if(evaluate(raw[j])) return;
            ++raw_count; raw_points.push_back(x);
            for(auto const& f:strong) if(evaluate(f)) return;
            ++strong_count; strong_points.push_back(x);
        };
        enumerate(0);
    }
    cout<<"{\"prime\":"<<prime<<",\"outer_count\":"<<outer_count
        <<",\"first_four_candidates\":"<<linear_candidates
        <<",\"raw_count\":"<<raw_count<<",\"strong_count\":"<<strong_count
        <<",\"consistent_rank_counts\":[";
    for(unsigned i=0;i<5;i++) { if(i) cout<<","; cout<<consistent_rank_count[i]; }
    cout<<"],\"raw_points\":"; points_json(raw_points);
    cout<<",\"strong_points\":"; points_json(strong_points); cout<<"}\n";
}
