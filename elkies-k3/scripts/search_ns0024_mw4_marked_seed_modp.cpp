// Deterministic finite-field seed search for the marked NS0024 MW4 chart.
// Compile with: c++ -O3 -std=c++17 -o /tmp/ns0024-mw4-seed this_file.cpp

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using Poly = std::vector<int>;

static int p;
static bool partial_mode=false, partial_all_mode=false;
static long long fibre_tests = 0, fibre_hits = 0;
static long long p1_tests = 0, p1_hits = 0, p2_tests = 0, p2_hits = 0;
static long long p3_tests = 0, p3_hits = 0, p4_tests = 0;

static inline int mod(long long value) {
    value %= p;
    return value < 0 ? int(value + p) : int(value);
}

static int power(int base, int exponent) {
    int answer = 1;
    while (exponent) {
        if (exponent & 1) answer = mod(1LL * answer * base);
        base = mod(1LL * base * base);
        exponent >>= 1;
    }
    return answer;
}

static int integer_power(int base, int exponent) {
    int answer = 1;
    while (exponent--) answer *= base;
    return answer;
}

static inline int inverse(int value) { return power(value, p - 2); }

static int square_root(int value) {
    for (int candidate = 0; candidate < p; ++candidate)
        if (mod(1LL * candidate * candidate) == value) return candidate;
    return -1;
}

static Poly multiply(const Poly& left, const Poly& right, int maximum_degree) {
    Poly answer(maximum_degree + 1, 0);
    for (int i = 0; i < int(left.size()); ++i)
        for (int j = 0; j < int(right.size()) && i + j <= maximum_degree; ++j)
            answer[i + j] = mod(answer[i + j] + 1LL * left[i] * right[j]);
    return answer;
}

static Poly branch_jet(const Poly& a, int root, int precision) {
    const int inv3 = inverse(3), inv2root = inverse(mod(2 * root));
    Poly h(precision, 0);
    h[0] = root;
    for (int degree = 1; degree < precision; ++degree) {
        int known = 0;
        for (int left = 1; left < degree; ++left)
            known = mod(known + 1LL * h[left] * h[degree - left]);
        h[degree] = mod((-1LL * a[degree] * inv3 - known) * inv2root);
    }
    Poly h2 = multiply(h, h, precision - 1);
    Poly h3 = multiply(h2, h, precision - 1);
    for (int& value : h3) value = mod(2 * value);
    return h3;
}

static Poly curve_rhs(const Poly& x, const Poly& a, const Poly& b, int maximum_degree) {
    Poly x2 = multiply(x, x, maximum_degree);
    Poly x3 = multiply(x2, x, maximum_degree);
    Poly ax = multiply(a, x, maximum_degree);
    Poly answer(maximum_degree + 1, 0);
    for (int index = 0; index <= maximum_degree; ++index) {
        int value = x3[index] + ax[index];
        if (index < int(b.size())) value += b[index];
        answer[index] = mod(value);
    }
    return answer;
}

static bool polynomial_square_root(const Poly& value, int maximum_y_degree, Poly& root) {
    int valuation = -1;
    for (int index = 0; index < int(value.size()); ++index)
        if (value[index]) { valuation = index; break; }
    if (valuation < 0 || (valuation & 1)) return false;
    const int initial = valuation / 2;
    if (initial > maximum_y_degree) return false;
    const int first = square_root(value[valuation]);
    if (first < 0 || first == 0) return false;
    root.assign(maximum_y_degree + 1, 0);
    root[initial] = first;
    const int denominator_inverse = inverse(mod(2 * first));
    for (int offset = 1; initial + offset <= maximum_y_degree; ++offset) {
        const int degree = 2 * initial + offset;
        int known = 0;
        for (int left = initial + 1; left < initial + offset; ++left)
            known = mod(known + 1LL * root[left] * root[degree - left]);
        const int target = degree < int(value.size()) ? value[degree] : 0;
        root[initial + offset] = mod(1LL * (target - known) * denominator_inverse);
    }
    Poly check = multiply(root, root, 2 * maximum_y_degree);
    for (int index = 0; index <= 2 * maximum_y_degree; ++index) {
        const int target = index < int(value.size()) ? value[index] : 0;
        if (check[index] != target) return false;
    }
    return true;
}

// Allocation-free hot path for Y^2 = X^3 + A*X*afactor + B*bfactor.
// The factors are 1 for polynomial sections and h^4,h^6 for a simple pole.
static bool section_square_root(const Poly& x, const Poly& a, const Poly& b,
                                int maximum_y_degree, Poly& root,
                                const Poly& afactor = Poly{1},
                                const Poly& bfactor = Poly{1}) {
    const int maximum_degree = 2 * maximum_y_degree;
    std::array<int,19> value{};
    std::array<int,19> first_product{};
    for (int i=0;i<int(x.size());++i)
      for (int j=0;j<int(x.size()) && i+j<=maximum_degree;++j)
        first_product[i+j]=mod(first_product[i+j]+1LL*x[i]*x[j]);
    for (int i=0;i<=maximum_degree;++i) if (first_product[i])
      for (int j=0;j<int(x.size()) && i+j<=maximum_degree;++j)
        value[i+j]=mod(value[i+j]+1LL*first_product[i]*x[j]);
    first_product.fill(0);
    for (int i=0;i<int(a.size());++i)
      for (int j=0;j<int(x.size()) && i+j<=maximum_degree;++j)
        first_product[i+j]=mod(first_product[i+j]+1LL*a[i]*x[j]);
    for (int i=0;i<=maximum_degree;++i) if (first_product[i])
      for (int j=0;j<int(afactor.size()) && i+j<=maximum_degree;++j)
        value[i+j]=mod(value[i+j]+1LL*first_product[i]*afactor[j]);
    for (int i=0;i<int(b.size());++i)
      for (int j=0;j<int(bfactor.size()) && i+j<=maximum_degree;++j)
        value[i+j]=mod(value[i+j]+1LL*b[i]*bfactor[j]);
    int valuation=-1;
    for (int index=0;index<=maximum_degree;++index)
        if (value[index]) { valuation=index; break; }
    if (valuation<0 || (valuation&1)) return false;
    const int initial=valuation/2;
    if (initial>maximum_y_degree) return false;
    const int first=square_root(value[valuation]);
    if (first<=0) return false;
    std::array<int,10> answer{};
    answer[initial]=first;
    const int denominator_inverse=inverse(mod(2*first));
    for (int offset=1;initial+offset<=maximum_y_degree;++offset) {
        const int degree=2*initial+offset;
        int known=0;
        for (int left=initial+1;left<initial+offset;++left)
            known=mod(known+1LL*answer[left]*answer[degree-left]);
        answer[initial+offset]=mod(1LL*(value[degree]-known)*denominator_inverse);
    }
    for (int degree=0;degree<=maximum_degree;++degree) {
        int check=0;
        for (int left=0;left<=maximum_y_degree;++left) {
            int right=degree-left;
            if (right>=0 && right<=maximum_y_degree)
                check=mod(check+1LL*answer[left]*answer[right]);
        }
        if (check!=value[degree]) return false;
    }
    root.assign(answer.begin(),answer.begin()+maximum_y_degree+1);
    return true;
}

static int evaluate(const Poly& polynomial, int point) {
    int answer = 0;
    for (int index = int(polynomial.size()) - 1; index >= 0; --index)
        answer = mod(1LL * answer * point + polynomial[index]);
    return answer;
}

static int binomial_mod(int n, int k) {
    if (k < 0 || k > n) return 0;
    long long answer = 1;
    for (int index = 1; index <= k; ++index)
        answer = mod(answer * (n - k + index) * inverse(index));
    return int(answer);
}

static bool exact_fibre_orders(const Poly& a, const Poly& b) {
    Poly a2 = multiply(a, a, 24), a3 = multiply(a2, a, 24);
    Poly b2 = multiply(b, b, 24), delta(25, 0);
    for (int index = 0; index <= 24; ++index)
        delta[index] = mod(4LL * a3[index] + 27LL * b2[index]);
    for (int index = 0; index < 7; ++index) if (delta[index]) return false;
    if (!delta[7]) return false;
    for (int jet = 0; jet < 5; ++jet) {
        int value = 0;
        for (int index = jet; index <= 24; ++index)
            value = mod(value + 1LL * delta[index] * binomial_mod(index, jet));
        if (value) return false;
    }
    int fifth = 0;
    for (int index = 5; index <= 24; ++index)
        fifth = mod(fifth + 1LL * delta[index] * binomial_mod(index, 5));
    if (!fifth) return false;
    for (int index = 21; index <= 24; ++index) if (delta[index]) return false;
    return delta[20] != 0;
}

static bool infinity_node(const Poly& x, const Poly& y, int ri, int x_weight, int y_weight) {
    const int xv = int(x.size()) - 1 == x_weight ? x.back() : 0;
    const int yv = int(y.size()) - 1 == y_weight ? y.back() : 0;
    return xv == ri && yv == 0;
}

static void print_poly(const Poly& value) {
    for (int index = 0; index < int(value.size()); ++index) {
        if (index) std::cout << ',';
        std::cout << value[index];
    }
}

static void trim(Poly& value) {
    while (value.size() > 1 && value.back() == 0) value.pop_back();
}

static bool poly_zero(const Poly& value) {
    for (int entry : value) if (entry) return false;
    return true;
}

static Poly poly_add(const Poly& left, const Poly& right, int sign = 1) {
    Poly answer(std::max(left.size(), right.size()), 0);
    for (int index = 0; index < int(answer.size()); ++index)
        answer[index] = mod((index < int(left.size()) ? left[index] : 0)
                            + sign * (index < int(right.size()) ? right[index] : 0));
    trim(answer); return answer;
}

static std::pair<Poly,Poly> poly_divmod(Poly dividend, Poly divisor) {
    trim(dividend); trim(divisor);
    if (poly_zero(divisor)) throw std::runtime_error("polynomial division by zero");
    Poly quotient(std::max(1, int(dividend.size()) - int(divisor.size()) + 1), 0);
    const int leading_inverse = inverse(divisor.back());
    while (!poly_zero(dividend) && dividend.size() >= divisor.size()) {
        int shift = int(dividend.size() - divisor.size());
        int coefficient = mod(1LL * dividend.back() * leading_inverse);
        quotient[shift] = coefficient;
        for (int index = 0; index < int(divisor.size()); ++index)
            dividend[index + shift] = mod(dividend[index + shift] - 1LL * coefficient * divisor[index]);
        trim(dividend);
    }
    trim(quotient); return {quotient, dividend};
}

static Poly poly_gcd(Poly left, Poly right) {
    while (!poly_zero(right)) {
        Poly remainder = poly_divmod(left, right).second;
        left = right; right = remainder;
    }
    if (poly_zero(left)) return Poly{1};
    int scale = inverse(left.back());
    for (int& entry : left) entry = mod(1LL * entry * scale);
    return left;
}

struct Rat {
    Poly numerator{0}, denominator{1};
    Rat() = default;
    Rat(Poly value): numerator(std::move(value)), denominator{1} { normalize(); }
    Rat(Poly n, Poly d): numerator(std::move(n)), denominator(std::move(d)) { normalize(); }
    void normalize() {
        trim(numerator); trim(denominator);
        if (poly_zero(numerator)) { numerator = {0}; denominator = {1}; return; }
        Poly common = poly_gcd(numerator, denominator);
        numerator = poly_divmod(numerator, common).first;
        denominator = poly_divmod(denominator, common).first;
        int scale = inverse(denominator.back());
        for (int& entry : numerator) entry = mod(1LL * entry * scale);
        for (int& entry : denominator) entry = mod(1LL * entry * scale);
    }
};

static bool operator==(const Rat& left, const Rat& right) {
    return poly_add(multiply(left.numerator, right.denominator,
                             int(left.numerator.size()+right.denominator.size()-2)),
                    multiply(right.numerator, left.denominator,
                             int(right.numerator.size()+left.denominator.size()-2)), -1)
           == Poly{0};
}

static Rat rat_add(const Rat& left, const Rat& right, int sign = 1) {
    int degree = int(std::max(left.numerator.size()+right.denominator.size(),
                              right.numerator.size()+left.denominator.size())) - 2;
    Poly numerator = poly_add(multiply(left.numerator,right.denominator,degree),
                              multiply(right.numerator,left.denominator,degree),sign);
    Poly denominator = multiply(left.denominator,right.denominator,
                                int(left.denominator.size()+right.denominator.size()-2));
    return Rat(numerator, denominator);
}

static Rat rat_mul(const Rat& left, const Rat& right) {
    return Rat(multiply(left.numerator,right.numerator,
                        int(left.numerator.size()+right.numerator.size()-2)),
               multiply(left.denominator,right.denominator,
                        int(left.denominator.size()+right.denominator.size()-2)));
}

static Rat rat_div(const Rat& left, const Rat& right) {
    return Rat(multiply(left.numerator,right.denominator,
                        int(left.numerator.size()+right.denominator.size()-2)),
               multiply(left.denominator,right.numerator,
                        int(left.denominator.size()+right.numerator.size()-2)));
}

struct ECPoint {
    bool zero = true;
    Rat x, y;
};

static ECPoint point_neg(ECPoint point) {
    if (!point.zero) for (int& entry : point.y.numerator) entry = mod(-entry);
    return point;
}

static ECPoint point_add(const ECPoint& left, const ECPoint& right, const Poly& a) {
    if (left.zero) return right;
    if (right.zero) return left;
    if (left.x == right.x && rat_add(left.y,right.y).numerator == Poly{0}) return ECPoint{};
    Rat slope;
    if (left.x == right.x) {
        Rat three(Poly{3}), two(Poly{2}), ar(a);
        slope = rat_div(rat_add(rat_mul(three,rat_mul(left.x,left.x)),ar),rat_mul(two,left.y));
    } else {
        slope = rat_div(rat_add(right.y,left.y,-1),rat_add(right.x,left.x,-1));
    }
    Rat x3 = rat_add(rat_add(rat_mul(slope,slope),left.x,-1),right.x,-1);
    Rat y3 = rat_add(rat_mul(slope,rat_add(left.x,x3,-1)),left.y,-1);
    return ECPoint{false,x3,y3};
}

static ECPoint point_multiple(int multiplier, ECPoint point, const Poly& a) {
    if (multiplier < 0) return point_multiple(-multiplier,point_neg(point),a);
    ECPoint answer;
    while (multiplier) {
        if (multiplier & 1) answer = point_add(answer,point,a);
        point = point_add(point,point,a); multiplier >>= 1;
    }
    return answer;
}

static int rat_finite_value(const Rat& value, int support, bool& finite) {
    int denominator = evaluate(value.denominator,support);
    finite = denominator != 0;
    return finite ? mod(1LL * evaluate(value.numerator,support) * inverse(denominator)) : 0;
}

static bool point_hits_node(const ECPoint& point, int fibre, int r1, int ri) {
    if (point.zero) return false;
    if (fibre < 2) {
        bool xf, yf; int support = fibre;
        int x = rat_finite_value(point.x,support,xf), y = rat_finite_value(point.y,support,yf);
        return xf && yf && x == (fibre ? r1 : 1) && y == 0;
    }
    int x_excess = int(point.x.numerator.size()-point.x.denominator.size());
    int y_excess = int(point.y.numerator.size()-point.y.denominator.size());
    int x = x_excess < 4 ? 0 : x_excess == 4
        ? mod(1LL*point.x.numerator.back()*inverse(point.x.denominator.back())) : -1;
    int y = y_excess < 6 ? 0 : y_excess == 6
        ? mod(1LL*point.y.numerator.back()*inverse(point.y.denominator.back())) : -1;
    return x == ri && y == 0;
}

static int component_label(const ECPoint& point, const ECPoint& reference,
                           int order, int fibre, int r1, int ri, const Poly& a) {
    int answer = -1;
    for (int multiplier = 0; multiplier < order; ++multiplier) {
        ECPoint difference = point_add(point,point_neg(point_multiple(multiplier,reference,a)),a);
        if (!point_hits_node(difference,fibre,r1,ri)) {
            if (answer >= 0) return -1;
            answer = multiplier;
        }
    }
    return answer;
}

static int section_intersection(const ECPoint& left, const ECPoint& right, const Poly& a) {
    ECPoint difference = point_add(left,point_neg(right),a);
    if (difference.zero) return -2;
    int denominator_degree = int(difference.x.denominator.size())-1;
    int infinity_excess = std::max(0,int(difference.x.numerator.size()-difference.x.denominator.size())-4);
    if ((denominator_degree+infinity_excess)&1) return -99;
    return (denominator_degree+infinity_excess)/2;
}

struct Section {
    Poly x, y, h;
};

static std::vector<std::array<Section,3>> all_partial_triples;

static ECPoint ec_point(const Section& section) {
    Poly h2 = multiply(section.h,section.h,2*int(section.h.size()-1));
    Poly h3 = multiply(h2,section.h,3*int(section.h.size()-1));
    return ECPoint{false,Rat(section.x,h2),Rat(section.y,h3)};
}

static bool orient_first_three(std::array<Section,3>& sections,
                               const Poly& a, int r1, int ri) {
    const int expected_profiles[3][3]={{4,0,0},{1,1,3},{1,1,1}};
    const int expected_intersections[3][3]={{-2,1,2},{1,-2,0},{2,0,-2}};
    const int orders[3]={7,5,4};
    std::array<ECPoint,3> raw;
    for (int index=0;index<3;++index) raw[index]=ec_point(sections[index]);
    for (int mask=0;mask<8;++mask) {
        std::array<ECPoint,3> point=raw;
        for (int index=0;index<3;++index)
            if ((mask>>index)&1) point[index]=point_neg(point[index]);
        bool pass=true;
        for (int index=0;index<3 && pass;++index)
          for (int fibre=0;fibre<3;++fibre)
            if (component_label(point[index],point[2],orders[fibre],fibre,r1,ri,a)
                != expected_profiles[index][fibre]) pass=false;
        for (int i=0;i<3 && pass;++i) for (int j=0;j<i;++j)
            if (section_intersection(point[i],point[j],a)
                != expected_intersections[i][j]) pass=false;
        if (!pass) continue;
        for (int index=0;index<3;++index) if ((mask>>index)&1)
            for (int& entry:sections[index].y) entry=mod(-entry);
        return true;
    }
    return false;
}

static bool orient_exact_marked_basis(std::array<Section,4>& sections,
                                      const Poly& a, int r1, int ri) {
    const int expected_profiles[4][3] = {{1,0,0},{2,1,3},{2,1,1},{1,1,1}};
    const int expected_intersections[4][4] = {
        {-2,1,2,1},{1,-2,0,1},{2,0,-2,1},{1,1,1,-2}
    };
    std::array<ECPoint,4> raw;
    for (int index=0;index<4;++index) raw[index]=ec_point(sections[index]);
    for (int mask=0;mask<16;++mask) {
        std::array<ECPoint,4> point=raw;
        for (int index=0;index<4;++index) if ((mask>>index)&1) point[index]=point_neg(point[index]);
        const ECPoint references[3]={point[3],point[3],point[3]};
        bool pass=true;
        const int orders[3]={7,5,4};
        for (int index=0;index<4 && pass;++index)
            for (int fibre=0;fibre<3;++fibre)
                if (component_label(point[index],references[fibre],orders[fibre],fibre,r1,ri,a)
                    != expected_profiles[index][fibre]) pass=false;
        for (int i=0;i<4 && pass;++i) for (int j=0;j<i;++j)
            if (section_intersection(point[i],point[j],a)!=expected_intersections[i][j]) pass=false;
        if (!pass) continue;
        for (int index=0;index<4;++index) if ((mask>>index)&1)
            for (int& entry:sections[index].y) entry=mod(-entry);
        return true;
    }
    return false;
}

static std::vector<Section> collect_p4(
    const Poly& a, const Poly& b, int r1, int ri
) {
    std::vector<Section> result;
    const int maximum = integer_power(p, 4);
    for (int c = 2; c < p; ++c) {
        Poly h = {mod(-c), 1};
        Poly h2 = multiply(h, h, 18), h4 = multiply(h2, h2, 18);
        Poly h6 = multiply(h4, h2, 18);
        for (int code = 0; code < maximum; ++code) {
            ++p4_tests;
            int cursor = code;
            Poly x(7, 0);
            x[0] = mod(1LL * c * c);
            x[6] = ri;
            for (int index = 1; index <= 4; ++index) { x[index] = cursor % p; cursor /= p; }
            int x_at_one = mod(1LL * r1 * mod(1LL * (1-c) * (1-c)));
            int partial_sum = 0;
            for (int index = 0; index <= 4; ++index) partial_sum = mod(partial_sum + x[index]);
            x[5] = mod(x_at_one - partial_sum - x[6]);
            if (!evaluate(x, c)) continue;
            Poly y;
            if (!section_square_root(x,a,b,9,y,h4,h6) || y[0] != 0
                || evaluate(y,1) != 0 || y[9] != 0 || !evaluate(y, c)) continue;
            result.push_back(Section{x,y,h});
        }
    }
    return result;
}

static bool scan_sections(const Poly& a, const Poly& b, int r1, int ri, std::array<Section,4>& answer) {
    std::vector<Section> first_sections, all_node_sections;
    const int p4code = integer_power(p, 4);
    for (int code1 = 0; code1 < p4code; ++code1) {
        ++p1_tests;
        int cursor1 = code1;
        Poly x1(5,0); x1[0]=1;
        for (int index=1;index<=4;++index) { x1[index]=cursor1%p; cursor1/=p; }
        Poly y1;
        if (!section_square_root(x1,a,b,6,y1)) continue;
        if (y1[0]) continue;
        if (evaluate(x1,1)==r1 && evaluate(y1,1)==0) continue;
        if (x1[4]==ri && y1[6]==0) continue;
        ++p1_hits;
        first_sections.push_back(Section{x1,y1,{1}});
    }
    if (first_sections.empty()) return false;

    // Components +/-2 at the I7 fibre are the strict-transform chart
    // X=center mod t^2.  Since center^2=-A/3 and center(0)=1, its linear
    // coefficient is -a1/6.  Substituting it here removes the singular node
    // tangent and one full finite-field loop from both all-node sections.
    const int p2code = p;
    for (int code = 0; code < p2code; ++code) {
        ++p2_tests;
        int cursor=code;
        Poly x(5,0); x[0]=1; x[4]=ri;
        x[1]=mod(-1LL*a[1]*inverse(6)); x[2]=cursor%p;
        x[3]=mod(r1-1-ri-x[1]-x[2]);
        Poly y;
        if (!section_square_root(x,a,b,6,y) || y[0] || evaluate(y,1) || y[6]) continue;
        ++p2_hits;
        all_node_sections.push_back(Section{x,y,{1}});
    }
    if (all_node_sections.size()<2) return false;

    std::vector<std::array<Section,3>> triples;
    for (const Section& P1:first_sections)
      for (int second=0;second<int(all_node_sections.size());++second)
       for (int third=0;third<int(all_node_sections.size());++third) {
        if (second==third) continue;
        ++p3_tests;
        const Section& P2=all_node_sections[second];
        const Section& P3=all_node_sections[third];
        std::array<Section,3> oriented={P1,P2,P3};
        if (!orient_first_three(oriented,a,r1,ri)) continue;
        ++p3_hits;
        triples.push_back(oriented);
    }
    if (triples.empty()) return false;
    if (partial_mode) {
        if (partial_all_mode) all_partial_triples=triples;
        answer={triples[0][0],triples[0][1],triples[0][2],Section{}};
        return true;
    }
    std::vector<Section> fourth_sections=collect_p4(a,b,r1,ri);
    for (const auto& triple:triples) for (const Section& P4:fourth_sections) {
        std::array<Section,4> candidate={triple[0],triple[1],triple[2],P4};
        if (orient_exact_marked_basis(candidate,a,r1,ri)) {
            answer=candidate;
            return true;
        }
    }
    return false;
}

int main(int argc, char** argv) {
    if (argc != 2 && argc != 4 && argc != 5) {
        std::cerr << "usage: search_ns0024_mw4_marked_seed_modp PRIME [R1 RI [partial|partial-all]]\n";
        return 2;
    }
    p = std::stoi(argv[1]);
    if (p != 11 && p != 13 && p != 17) { std::cerr << "supported primes: 11,13,17\n"; return 2; }
    const int total = integer_power(p, 6);
    const int choose[13][6] = {
        {1,0,0,0,0,0},{1,1,0,0,0,0},{1,2,1,0,0,0},{1,3,3,1,0,0},
        {1,4,6,4,1,0},{1,5,10,10,5,1},{1,6,15,20,15,6},
        {1,7,21,35,35,21},{1,8,28,56,70,56},{1,9,36,84,126,126},
        {1,10,45,120,210,252},{1,11,55,165,330,462},{1,12,66,220,495,792}
    };
    partial_mode = argc == 5 && (
        std::string(argv[4]) == "partial" || std::string(argv[4]) == "partial-all"
    );
    partial_all_mode = argc == 5 && std::string(argv[4]) == "partial-all";
    if (argc == 5 && !partial_mode) {
        std::cerr << "the fourth argument must be 'partial' or 'partial-all'\n";
        return 2;
    }
    const int r1_begin = argc >= 4 ? std::stoi(argv[2]) : 1;
    const int r1_end = argc >= 4 ? r1_begin + 1 : p;
    const int ri_begin = argc >= 4 ? std::stoi(argv[3]) : 1;
    const int ri_end = argc >= 4 ? ri_begin + 1 : p;
    if (r1_begin < 1 || r1_end > p || ri_begin < 1 || ri_end > p) {
        std::cerr << "R1 and RI must be nonzero residues modulo PRIME\n";
        return 2;
    }
    bool found_any = false;
    for (int r1 = r1_begin; r1 < r1_end; ++r1)
      for (int ri = ri_begin; ri < ri_end; ++ri) {
        if (square_root(mod(3 * r1)) < 0 || square_root(mod(3 * ri)) < 0) continue;
        for (int code = 0; code < total; ++code) {
            ++fibre_tests;
            int cursor = code;
            Poly a(9, 0); a[0] = mod(-3); a[8] = mod(-3LL * ri * ri);
            for (int index = 1; index <= 6; ++index) { a[index] = cursor % p; cursor /= p; }
            int sum = 0; for (int index = 0; index <= 6; ++index) sum = mod(sum + a[index]);
            a[7] = mod(-3LL * r1 * r1 - sum - a[8]);
            Poly zero(a.begin(), a.begin() + 7);
            Poly ar = {a[8],a[7],a[6],a[5]};
            Poly bz = branch_jet(zero, 1, 7), bi = branch_jet(ar, ri, 4);
            Poly b(13, 0);
            for (int i = 0; i < 7; ++i) b[i] = bz[i];
            for (int i = 0; i < 4; ++i) b[12-i] = bi[i];
            Poly at1(5, 0);
            for (int jet = 0; jet < 5; ++jet)
                for (int i = jet; i <= 8; ++i) at1[jet] = mod(at1[jet] + 1LL * a[i] * choose[i][jet]);
            Poly target = branch_jet(at1, r1, 5);
            int known0 = 0, known1 = 0;
            for (int i = 0; i <= 12; ++i) if (i != 7 && i != 8) {
                known0 = mod(known0 + b[i]);
                known1 = mod(known1 + 1LL * i * b[i]);
            }
            const int rhs0 = mod(target[0] - known0), rhs1 = mod(target[1] - known1);
            b[8] = mod(rhs1 - 7 * rhs0); b[7] = mod(rhs0 - b[8]);
            bool compatible = true;
            for (int jet = 0; jet < 5; ++jet) {
                int actual = 0;
                for (int i = jet; i <= 12; ++i) actual = mod(actual + 1LL * b[i] * choose[i][jet]);
                if (actual != target[jet]) { compatible = false; break; }
            }
            if (!compatible) continue;
            if (!exact_fibre_orders(a, b)) continue;
            ++fibre_hits;
            std::array<Section,4> sections;
            if (!scan_sections(a, b, r1, ri, sections)) continue;
            found_any = true;
            const int output_count = partial_all_mode ? int(all_partial_triples.size()) : 1;
            for (int output_index = 0; output_index < output_count; ++output_index) {
            if (partial_all_mode) {
                for (int index=0;index<3;++index) sections[index]=all_partial_triples[output_index][index];
            }
            std::cout << (partial_mode ? "NS0024MW3SEED|p=" : "NS0024MW4SEED|p=")
                      << p << "|r1=" << r1 << "|ri=" << ri << "|A=";
            print_poly(a); std::cout << "|B="; print_poly(b);
            for (int index = 0; index < (partial_mode ? 3 : 4); ++index) {
                std::cout << "|P" << index+1 << "X="; print_poly(sections[index].x);
                std::cout << "|P" << index+1 << "Y="; print_poly(sections[index].y);
                if (index == 3) { std::cout << "|P4H="; print_poly(sections[index].h); }
            }
            std::cout << "|fibre_tests=" << fibre_tests << "|fibre_hits=" << fibre_hits
                      << "|p1_tests=" << p1_tests << "|p1_hits=" << p1_hits
                      << "|p2_tests=" << p2_tests << "|p2_hits=" << p2_hits
                      << "|p3_tests=" << p3_tests << "|p3_hits=" << p3_hits
                      << "|p4_tests=" << p4_tests << "\n";
            }
            if (!partial_all_mode) return 0;
        }
    }
    if (found_any) return 0;
    std::cout << (partial_mode ? "NS0024MW3SEED|p=" : "NS0024MW4SEED|p=") << p
              << "|status=NO_SEED|fibre_tests=" << fibre_tests
              << "|fibre_hits=" << fibre_hits << "|p4_tests=" << p4_tests << "\n";
    return 1;
}
