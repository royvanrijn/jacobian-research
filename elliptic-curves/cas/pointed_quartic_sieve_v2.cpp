// Complete primitive rational-slope boxes for integral binary quartics.
// Build: g++ -O3 -std=c++17 pointed_quartic_sieve.cpp -lgmpxx -lgmp -o sieve
// stdin: H denominator_start denominator_end seconds; then 5 ascending coefficients.
// Modular square and primitive-pair filters; empty denominator rows skip whole scans.
// Equivalent primitive hits, separate version: original worker stays frozen.
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <vector>

struct Point {
    mpq_class x = 0, y = 0;
    bool infinity = true;
};

Point add(const Point &p, const Point &q, const mpq_class &a) {
    if (p.infinity) return q;
    if (q.infinity) return p;
    mpq_class slope;
    if (p.x == q.x) {
        if (p.y == -q.y) return {};
        slope = (3*p.x*p.x+a)/(2*p.y);
    } else slope = (q.y-p.y)/(q.x-p.x);
    Point result;
    result.infinity = false;
    result.x = slope*slope-p.x-q.x;
    result.y = slope*(p.x-result.x)-p.y;
    return result;
}

int base_point() {
    int count;
    mpq_class a,b;
    if (!(std::cin >> count >> a >> b) || count < 1 || count > 1024) return 2;
    a.canonicalize(); b.canonicalize();
    if (4*a*a*a+27*b*b == 0) return 2;
    std::vector<Point> points(count);
    std::vector<mpz_class> coefficients(count);
    size_t bits = 0;
    for (int i = 0; i < count; ++i) {
        auto &p = points[i];
        if (!(std::cin >> p.x >> p.y >> coefficients[i])) return 2;
        p.x.canonicalize(); p.y.canonicalize(); p.infinity = false;
        if (p.y*p.y != p.x*p.x*p.x+a*p.x+b) return 2;
        if (coefficients[i] < 0) { coefficients[i] = -coefficients[i]; p.y = -p.y; }
        bits = std::max(bits, mpz_sizeinbase(coefficients[i].get_mpz_t(),2));
    }
    Point result;
    // Simultaneous signed binary multiplication avoids forming each large
    // multiple separately. All arithmetic is exact GMP rational arithmetic.
    for (size_t bit = bits; bit-- > 0;) {
        result = add(result,result,a);
        for (int i = 0; i < count; ++i)
            if (mpz_tstbit(coefficients[i].get_mpz_t(),bit)) result = add(result,points[i],a);
    }
    if (result.infinity) std::cout << "INFINITY\n";
    else {
        if (result.y*result.y != result.x*result.x*result.x+a*result.x+b) return 3;
        std::cout << "BASE " << result.x << ' ' << result.y << '\n';
    }
    return 0;
}

struct Filter {
    int p, count = 0;
    std::vector<unsigned char> allowed;
    std::vector<uint64_t> words;
    int word_step = 0;
    std::vector<unsigned char> empty_rows;
};

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--base-point") return base_point();
    if (argc == 2 && std::string(argv[1]) == "--version") {
        std::cout << "GMP " << gmp_version << "; compiler " << __VERSION__ << '\n';
        return 0;
    }
    int height, first, last;
    double seconds;
    std::array<mpz_class, 5> f;
    if (!(std::cin >> height >> first >> last >> seconds) || height < 1 ||
        height > 1000000 || first < 1 || last < first || last > height || seconds <= 0)
        return 2;
    for (auto &c : f) if (!(std::cin >> c)) return 2;
    const auto start = std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count(); };
    std::vector<Filter> filters;
    for (int p = 3; p <= 151; p += 2) {
        bool prime = true;
        for (int q = 3; q*q <= p; q += 2) if (p % q == 0) prime = false;
        if (!prime) continue;
        Filter filter{p, 0, std::vector<unsigned char>(p*p), {}};
        std::vector<bool> square(p, false);
        for (int i = 0; i < p; ++i) square[i*i % p] = true;
        std::array<int,5> c;
        for (int i = 0; i < 5; ++i) c[i] = mpz_fdiv_ui(f[i].get_mpz_t(), p);
        for (int d = 0; d < p; ++d) for (int n = 0; n < p; ++n) {
            int value = c[4], dpow = d;
            for (int i = 3; i >= 0; --i) {
                value = (value*n + c[i]*dpow) % p;
                dpow = dpow*d % p;
            }
            filter.allowed[d*p+n] = square[value] && (d != 0 || n != 0);
            filter.count += filter.allowed[d*p+n];
        }
        filter.empty_rows.resize(p,1);
        for (int d=0; d<p; ++d) for (int n=0; n<p; ++n)
            if (filter.allowed[d*p+n]) filter.empty_rows[d]=0;
        if (filter.count < p*p) filters.push_back(std::move(filter));
    }
    std::sort(filters.begin(), filters.end(), [](const Filter &a, const Filter &b) {
        const int lhs = a.count*b.p*b.p, rhs = b.count*a.p*a.p;
        return lhs == rhs ? a.p < b.p : lhs < rhs;
    });
    const int word_filters = std::min(8, int(filters.size()));
    for (int i = 0; i < word_filters; ++i) {
        auto &s = filters[i];
        s.words.resize(s.p*s.p);
        s.word_step = 64%s.p;
        for (int d = 0; d < s.p; ++d) for (int n = 0; n < s.p; ++n) {
            uint64_t word = 0;
            for (int j = 0; j < 64; ++j)
                word |= uint64_t(s.allowed[d*s.p+(n+j)%s.p]) << j;
            s.words[d*s.p+n] = word;
        }
    }
    std::cout << "PRIMES";
    for (const auto &s : filters) std::cout << ' ' << s.p;
    std::cout << '\n';
    uint64_t word_survivors = 0, modular_survivors = 0, exact_tests = 0, hits = 0;
    int completed = first-1;
    for (int d = first; d <= last; ++d) {
        if (elapsed() >= seconds) break; // Finish each denominator atomically.
        std::vector<int> residues(word_filters), offsets(filters.size());
        bool possible = true;
        for (size_t i = 0; i < filters.size(); ++i) {
            int residue = d%filters[i].p;
            offsets[i] = residue*filters[i].p;
            if (filters[i].empty_rows[residue]) possible=false;
        }
        if (!possible) { completed=d; continue; }
        for (int i = 0; i < word_filters; ++i)
            residues[i] = ((-height)%filters[i].p+filters[i].p)%filters[i].p;
        for (int begin = -height; begin <= height; begin += 64) {
            uint64_t word = ~uint64_t(0);
            int remaining = height-begin+1;
            if (remaining < 64) word = (uint64_t(1) << remaining)-1;
            for (int i = 0; i < word_filters; ++i) {
                const auto &s = filters[i];
                word &= s.words[offsets[i]+residues[i]];
                residues[i] += s.word_step;
                if (residues[i] >= s.p) residues[i] -= s.p;
            }
            word_survivors += __builtin_popcountll(word);
            while (word) {
                const int n = begin+__builtin_ctzll(word);
                word &= word-1;
                bool pass = true;
                for (size_t i = word_filters; i < filters.size(); ++i) {
                    const auto &s = filters[i];
                    if (!s.allowed[offsets[i]+(n%s.p+s.p)%s.p]) { pass = false; break; }
                }
                if (!pass) continue;
                ++modular_survivors;
                if (std::gcd(std::abs(n),d) != 1) continue;
                ++exact_tests;
                mpz_class value = f[4], dpow = d;
                for (int i = 3; i >= 0; --i) { value = value*n+f[i]*dpow; dpow *= d; }
                if (value < 0 || !mpz_perfect_square_p(value.get_mpz_t())) continue;
                mpz_class root;
                mpz_sqrt(root.get_mpz_t(), value.get_mpz_t());
                std::cout << "POINT " << n << ' ' << d << ' ' << root << '\n';
                ++hits;
            }
        }
        completed = d;
    }
    std::cout << "DONE " << completed << ' ' << word_survivors << ' '
              << modular_survivors << ' ' << exact_tests << ' ' << hits << ' '
              << elapsed() << '\n';
    return 0;
}
