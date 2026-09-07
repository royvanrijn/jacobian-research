#include "antipode/base/exact.h"
#include "antipode/cli/logging.h"
#include "antipode/linalg/linalg.h"

#include "absl/log/check.h"
#include "absl/log/log.h"
#include "gflags/gflags.h"

#include <array>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

DEFINE_string(ambient_file, "", "Exported Niemeier ambient.txt.");
DEFINE_string(repair_file, "", "Saved search state containing basis_coords.");
DEFINE_string(output_dir, "/tmp/elkies-niemeier-repair", "Output directory.");

DEFINE_int64(target_det, 948, "Target determinant.");
DEFINE_int64(target_norm4, 1311, "Target orthogonal norm-4 representatives.");

DEFINE_bool(
    scan_replacement, false,
    "Scan row-replacement neighborhood.");

DEFINE_bool(
    scan_additive, true,
    "Scan k_i -> k_i +/- v neighborhood.");

DEFINE_bool(
    scan_two_replacement, false,
    "Scan simultaneous replacement of two rows by norm-4 ambient vectors.");

DEFINE_bool(
    scan_two_additive, false,
    "Scan simultaneous k_i += s*u, k_j += t*v using norm-4 ambient vectors.");


DEFINE_bool(
    two_add_rootwalk, false,
    "Root-walk mode: save only a small directed frontier near norm4=1311.");

DEFINE_int32(
    two_add_rootwalk_per_mask, 4,
    "In root-walk mode, maximum saved (roots=1,norm4=1311) "
    "states per resulting root mask.");

DEFINE_int32(
    two_add_rootwalk_side_per_mask, 1,
    "In root-walk mode, maximum saved roots=1 norm4=1310/1312 "
    "states per resulting root mask.");

DEFINE_int32(
    two_add_rootwalk_two_per_mask, 1,
    "In root-walk mode, maximum saved (roots=2,norm4=1311) "
    "states per resulting root mask.");

DEFINE_int64(
    two_add_expand_r1, -1,
    "First changed row for concrete two-row additive expansion.");

DEFINE_int64(
    two_add_expand_r2, -1,
    "Second changed row for concrete two-row additive expansion.");

DEFINE_int64(
    two_add_expand_h, 99,
    "Required <u,v> for concrete two-row additive expansion.");

DEFINE_int64(
    two_add_expand_s1, 0,
    "Required sign (+/-1) for first additive row during concrete expansion.");

DEFINE_int64(
    two_add_expand_s2, 0,
    "Required sign (+/-1) for second additive row during concrete expansion.");

DEFINE_bool(
    two_add_signature_only, false,
    "Only write simultaneous two-row additive signature survivors.");


DEFINE_bool(
    two_add_require_rootless, true,
    "Reject simultaneous two-row additive candidates having roots before "
    "saturation. Set false when harvesting low-root intermediate states.");

DEFINE_bool(
    two_add_index1_only, false,
    "During concrete two-row additive expansion, only consider raw_det "
    "equal to target_det (predicted saturation index 1).");

DEFINE_int64(
    two_add_exact_index, 0,
    "If positive, during concrete two-row additive expansion only consider "
    "this predicted saturation index.");



DEFINE_int32(
    two_add_only_r1, -1,
    "If nonnegative, restrict simultaneous two-row additive scanning "
    "to row pairs whose first row equals this value.");

DEFINE_int32(
    two_add_only_r2, -1,
    "If nonnegative, restrict simultaneous two-row additive scanning "
    "to row pairs whose second row equals this value.");

DEFINE_bool(
    two_add_expand_all_index1, false,
    "Expand every realizable simultaneous two-row additive signature cell "
    "with raw_det == target_det in one census pass.");

DEFINE_int64(
    two_expand_r1, -1,
    "First row for concrete simultaneous two-row replacement expansion.");

DEFINE_int64(
    two_expand_r2, -1,
    "Second row for concrete simultaneous two-row replacement expansion.");

DEFINE_int64(
    two_expand_h, 99,
    "Required inner product <u,v> for concrete two-row expansion; rooted "
    "Niemeier norm-4 vectors require checking h=-3..3.");

DEFINE_bool(
    two_signature_only, false,
    "Only write simultaneous two-row signature survivors; do not expand "
    "concrete vector pairs.");

DEFINE_int64(
    norm4_max_roots, 3,
    "Compute exact norm-4 count for det-948 candidates with at most this many roots.");

DEFINE_int64(
    save_max_roots, 3,
    "Save candidate embeddings with at most this many roots.");

DEFINE_int64(
    save_n4_band, 15,
    "Save candidate embeddings with |norm4-1311| within this band.");



namespace antipode {
namespace {

constexpr int kRank = 7;
constexpr int kDim = 24;

struct AmbientData {
  IntMat gram;
  IntMat roots;
  IntMat norm4;
  IntMat roots_dual;
  IntMat norm4_dual;
};

struct State {
  IntMat basis;
  IntMat gram;
  Int det;
};


IntMat ReadMatrix(
    std::ifstream& in,
    long rows,
    long cols) {

  IntMat M(rows, cols);

  for (long i = 0; i < rows; ++i) {
    for (long j = 0; j < cols; ++j) {

      long x;
      CHECK(in >> x);

      M.Set(i, j, x);
    }
  }

  return M;
}


IntMat MetricDual(
    const IntMat& V,
    const IntMat& G) {

  IntMat VG(
      V.Rows(),
      G.Cols());

  fmpz_mat_mul(
      VG.raw(),
      V.raw(),
      G.raw());

  return VG;
}


AmbientData LoadAmbient(
    const std::string& path) {

  std::ifstream in(path);
  CHECK(in.good());

  std::string tag;
  long n;

  CHECK(in >> tag >> n);
  CHECK_EQ(tag, "GRAM");
  CHECK_EQ(n, kDim);

  AmbientData a;

  a.gram =
      ReadMatrix(
          in,
          kDim,
          kDim);

  long count;

  CHECK(in >> tag >> count);
  CHECK_EQ(tag, "ROOT_REPS");

  a.roots =
      ReadMatrix(
          in,
          count,
          kDim);

  CHECK(in >> tag >> count);
  CHECK_EQ(tag, "NORM4_REPS");

  a.norm4 =
      ReadMatrix(
          in,
          count,
          kDim);

  CHECK_EQ(
      Det(a.gram),
      Int(1));

  a.roots_dual =
      MetricDual(
          a.roots,
          a.gram);

  a.norm4_dual =
      MetricDual(
          a.norm4,
          a.gram);

  LOG(INFO)
      << "ambient roots="
      << a.roots.Rows()
      << " norm4="
      << a.norm4.Rows();

  return a;
}


IntMat ReadSavedBasis(
    const std::string& path) {

  std::ifstream in(path);
  CHECK(in.good());

  const std::string text(
      (std::istreambuf_iterator<char>(in)),
      std::istreambuf_iterator<char>());

  const std::string marker =
      "basis_coords =";

  size_t p =
      text.find(marker);

  CHECK_NE(
      p,
      std::string::npos);

  p += marker.size();

  std::vector<long> xs;
  xs.reserve(kRank * kDim);

  while (p < text.size() &&
         xs.size() <
             static_cast<size_t>(
                 kRank * kDim)) {

    while (p < text.size() &&
           text[p] != '-' &&
           !std::isdigit(
               static_cast<unsigned char>(
                   text[p]))) {

      ++p;
    }

    CHECK_LT(p, text.size());

    bool negative = false;

    if (text[p] == '-') {
      negative = true;
      ++p;
    }

    CHECK(
        std::isdigit(
            static_cast<unsigned char>(
                text[p])));

    long x = 0;

    while (p < text.size() &&
           std::isdigit(
               static_cast<unsigned char>(
                   text[p]))) {

      x =
          10 * x +
          (text[p] - '0');

      ++p;
    }

    xs.push_back(
        negative ? -x : x);
  }

  CHECK_EQ(
      xs.size(),
      static_cast<size_t>(
          kRank * kDim));

  IntMat B(
      kRank,
      kDim);

  size_t q = 0;

  for (int i = 0; i < kRank; ++i)
    for (int j = 0; j < kDim; ++j)
      B.Set(
          i,
          j,
          xs[q++]);

  return B;
}


IntMat MetricGram(
    const IntMat& B,
    const IntMat& G) {

  IntMat BG(
      B.Rows(),
      kDim);

  fmpz_mat_mul(
      BG.raw(),
      B.raw(),
      G.raw());

  IntMat BT(
      kDim,
      B.Rows());

  fmpz_mat_transpose(
      BT.raw(),
      B.raw());

  IntMat out(
      B.Rows(),
      B.Rows());

  fmpz_mat_mul(
      out.raw(),
      BG.raw(),
      BT.raw());

  return out;
}


State SaturatedState(
    const AmbientData& a,
    const IntMat& tuple) {

  State s;

  s.basis =
      SaturateRowSpan(
          tuple);

  CHECK_EQ(
      s.basis.Rows(),
      kRank);

  s.gram =
      MetricGram(
          s.basis,
          a.gram);

  s.det =
      Det(s.gram);

  return s;
}


long OrthogonalCount(
    const IntMat& basis,
    const IntMat& shell_dual) {

  long count = 0;

  for (long i = 0;
       i < shell_dual.Rows();
       ++i) {

    bool orth = true;

    for (int r = 0;
         r < kRank;
         ++r) {

      if (Dot(
              basis.RowPtr(r),
              shell_dual.RowPtr(i),
              kDim) !=
          Int(0)) {

        orth = false;
        break;
      }
    }

    if (orth)
      ++count;
  }

  return count;
}


std::string GramKey(
    const IntMat& G) {

  std::string key;

  for (long i = 0;
       i < G.Rows();
       ++i) {

    for (long j = 0;
         j <= i;
         ++j) {

      key +=
          G(i, j).get_str();

      key.push_back(',');
    }
  }

  return key;
}



std::string EmbeddedKey(
    const IntMat& basis) {

  const IntMat h =
      HermiteBasis(basis);

  std::string key;

  for (long i = 0; i < h.Rows(); ++i) {
    for (long j = 0; j < h.Cols(); ++j) {
      key += h(i,j).get_str();
      key.push_back(',');
    }
  }

  return key;
}


void SaveHit(
    const State& s,
    long roots,
    long norm4,
    int replace,
    const char* shell_name,
    long shell_index,
    const Int& raw_det,
    const Int& raw_index,
    uint64_t hit) {

  std::filesystem::create_directories(
      FLAGS_output_dir);

  const std::string path =
      FLAGS_output_dir +
      "/hit-" +
      std::to_string(hit) +
      "-row" +
      std::to_string(replace) +
      "-" +
      shell_name +
      std::to_string(shell_index) +
      "-r" +
      std::to_string(roots) +
      "-n" +
      std::to_string(norm4) +
      ".txt";

  std::ofstream out(path);

  out << "replaced_row = "
      << replace << "\n";

  out << "shell = "
      << shell_name << "\n";

  out << "shell_index = "
      << shell_index << "\n";

  out << "raw_det = "
      << raw_det << "\n";

  out << "raw_index = "
      << raw_index << "\n";

  out << "determinant = "
      << s.det << "\n";

  out << "root_reps_orthogonal = "
      << roots << "\n";

  out << "norm4_reps_orthogonal = "
      << norm4 << "\n";

  out << "\ngram =\n"
      << s.gram << "\n";

  out << "\nbasis_coords =\n"
      << s.basis << "\n";

  out << "\nbasis_hnf =\n"
      << HermiteBasis(s.basis)
      << "\n";

  LOG(INFO)
      << "saved "
      << path;
}




using EligibleByPair =
    std::array<std::array<std::vector<long>, kRank>, kRank>;

EligibleByPair BuildEligibleByChangedPair(
    const IntMat& basis,
    const IntMat& shell_dual) {

  EligibleByPair out;

  for (long i = 0; i < shell_dual.Rows(); ++i) {

    std::array<int, kRank> bad{};
    int nbad = 0;

    for (int r = 0; r < kRank; ++r) {
      if (Dot(
              basis.RowPtr(r),
              shell_dual.RowPtr(i),
              kDim) != Int(0)) {
        bad[nbad++] = r;
        if (nbad > 2)
          break;
      }
    }

    if (nbad > 2)
      continue;

    for (int r1 = 0; r1 < kRank; ++r1) {
      for (int r2 = r1 + 1; r2 < kRank; ++r2) {

        bool eligible = true;

        for (int q = 0; q < nbad; ++q) {
          if (bad[q] != r1 && bad[q] != r2) {
            eligible = false;
            break;
          }
        }

        if (eligible)
          out[r1][r2].push_back(i);
      }
    }
  }

  return out;
}


long OrthogonalCountTwoChangedRows(
    const IntVec& row1,
    const IntVec& row2,
    const IntMat& shell_dual,
    const std::vector<long>& eligible) {

  long count = 0;

  for (long i : eligible) {
    if (Dot(
            row1.data(),
            shell_dual.RowPtr(i),
            kDim) != Int(0)) {
      continue;
    }

    if (Dot(
            row2.data(),
            shell_dual.RowPtr(i),
            kDim) != Int(0)) {
      continue;
    }

    ++count;
  }

  return count;
}


uint64_t RootMaskTwoChangedRows(
    const IntVec& changed1,
    const IntVec& changed2,
    const IntMat& roots_dual,
    const std::vector<long>& eligible) {
  uint64_t mask = 0;

  for (long ri : eligible) {
    if (Dot(
            changed1.data(),
            roots_dual.RowPtr(ri),
            kDim) != Int(0)) {
      continue;
    }

    if (Dot(
            changed2.data(),
            roots_dual.RowPtr(ri),
            kDim) != Int(0)) {
      continue;
    }

    mask |= (uint64_t{1} << ri);
  }

  return mask;
}

uint64_t CurrentRootMask(
    const IntMat& basis,
    const IntMat& roots_dual) {
  uint64_t mask = 0;

  CHECK_LE(roots_dual.Rows(), 64);

  for (long ri = 0; ri < roots_dual.Rows(); ++ri) {
    bool orthogonal = true;

    for (int r = 0; r < kRank; ++r) {
      if (Dot(
              basis.RowPtr(r),
              roots_dual.RowPtr(ri),
              kDim) != Int(0)) {
        orthogonal = false;
        break;
      }
    }

    if (orthogonal) {
      mask |= (uint64_t{1} << ri);
    }
  }

  return mask;
}




void SaveTwoRowHit(
    const State& state,
    long roots,
    long norm4,
    int r1,
    int r2,
    long ui,
    long vi,
    long h,
    const Int& raw_det,
    const Int& raw_index,
    uint64_t hit,
    uint64_t rootmask = 0) {

  std::filesystem::create_directories(
      FLAGS_output_dir);

  const std::string path =
      FLAGS_output_dir +
      "/two-hit-" +
      std::to_string(hit) +
      "-rows" +
      std::to_string(r1) +
      "-" +
      std::to_string(r2) +
      "-h" +
      std::to_string(h) +
      "-u" +
      std::to_string(ui) +
      "-v" +
      std::to_string(vi) +
      "-r" +
      std::to_string(roots) +
      "-n" +
      std::to_string(norm4) +
      ".txt";

  std::ofstream out(path);

  out << "rootmask = " << rootmask << "\n";

  out << "rows = "
      << r1 << "," << r2 << "\n";
  out << "h = " << h << "\n";
  out << "u = " << ui << "\n";
  out << "v = " << vi << "\n";
  out << "raw_det = " << raw_det << "\n";
  out << "raw_index = " << raw_index << "\n";
  out << "determinant = " << state.det << "\n";
  out << "root_reps_orthogonal = " << roots << "\n";
  out << "norm4_reps_orthogonal = " << norm4 << "\n";
  out << "\ngram =\n" << state.gram << "\n";
  out << "\nbasis_coords =\n" << state.basis << "\n";
  out << "\nbasis_hnf =\n" << HermiteBasis(state.basis) << "\n";

  LOG(INFO) << "saved " << path;
}


using EligibleByRow =
    std::array<std::vector<long>, kRank>;

/*
 * For an additive move changing only row r, a shell vector can become
 * orthogonal to the new section only if it is already orthogonal to
 * all six unchanged rows.
 *
 * Build those candidate lists once.
 *
 * Efficiently:
 *
 *   bad == 0:
 *       shell vector is currently orthogonal to all seven rows;
 *       it is eligible for every changed row.
 *
 *   bad == 1:
 *       it is eligible only when we change that one offending row.
 *
 *   bad >= 2:
 *       changing one row can never make it orthogonal.
 */
EligibleByRow BuildEligibleByChangedRow(
    const IntMat& basis,
    const IntMat& shell_dual) {

  EligibleByRow out;

  for (long i = 0;
       i < shell_dual.Rows();
       ++i) {

    int bad = 0;
    int bad_row = -1;

    for (int r = 0;
         r < kRank;
         ++r) {

      if (Dot(
              basis.RowPtr(r),
              shell_dual.RowPtr(i),
              kDim) != Int(0)) {

        ++bad;
        bad_row = r;

        if (bad >= 2)
          break;
      }
    }

    if (bad == 0) {

      for (int r = 0;
           r < kRank;
           ++r) {
        out[r].push_back(i);
      }

    } else if (bad == 1) {

      out[bad_row].push_back(i);
    }
  }

  return out;
}


long OrthogonalCountOneChangedRow(
    const IntVec& changed_row,
    const IntMat& shell_dual,
    const std::vector<long>& eligible) {

  long count = 0;

  for (const long i : eligible) {

    if (Dot(
            changed_row.data(),
            shell_dual.RowPtr(i),
            kDim) == Int(0)) {

      ++count;
    }
  }

  return count;
}

int Run() {

  CHECK(!FLAGS_ambient_file.empty());
  CHECK(!FLAGS_repair_file.empty());

  const AmbientData a =
      LoadAmbient(
          FLAGS_ambient_file);

  const IntMat saved =
      ReadSavedBasis(
          FLAGS_repair_file);

  const State start =
      SaturatedState(
          a,
          saved);

  const long start_roots =
      OrthogonalCount(
          start.basis,
          a.roots_dual);

  const long start_norm4 =
      OrthogonalCount(
          start.basis,
          a.norm4_dual);

  LOG(INFO)
      << "START det="
      << start.det
      << " roots="
      << start_roots
      << " norm4="
      << start_norm4;



  /*
   * This catches accidentally pointing the repair at an older best.
   */
  // Rootful exact-det states are deliberately allowed here.
  // In particular we want to repair (948, roots=1, norm4=1312).

  const Int target(
      FLAGS_target_det);

  uint64_t tested = 0;
  uint64_t square_survivors = 0;
  uint64_t det948_hits = 0;

  std::set<std::string>
      distinct_grams;

  std::set<std::pair<long,long>>
      fingerprints;

  std::set<std::string>
      emitted_embeddings;

  /*
   * +/- representatives are sufficient:
   *
   * replacing a row by v or -v gives exactly the same row lattice
   * together with the six unchanged rows.
   */
  auto scan_shell =
      [&](const char* shell_name,
          const IntMat& shell,
          const IntMat& shell_dual,
          long shell_norm) {

        for (int replace = 0;
             replace < kRank;
             ++replace) {

          LOG(INFO)
              << "===== "
              << shell_name
              << " replacing row "
              << replace
              << " =====";

          for (long vi = 0;
               vi < shell.Rows();
               ++vi) {

            ++tested;

            /*
             * Construct the raw 7x7 Gram cheaply from the existing
             * det-946 Gram.
             */
            IntMat raw =
                start.gram;

            raw.Set(
                replace,
                replace,
                shell_norm);

            for (int j = 0;
                 j < kRank;
                 ++j) {

              if (j == replace)
                continue;

              const Int ip =
                  Dot(
                      start.basis.RowPtr(j),
                      shell_dual.RowPtr(vi),
                      kDim);

              raw.Set(
                  replace,
                  j,
                  ip);

              raw.Set(
                  j,
                  replace,
                  ip);
            }

            const Int raw_det =
                Det(raw);

            if (raw_det <= Int(0))
              continue;

            /*
             * Necessary condition for saturation to determinant 948:
             *
             *   raw_det = 948 * index^2.
             */
            if (raw_det % target !=
                Int(0)) {
              continue;
            }

            const Int q =
                raw_det /
                target;

            if (q <= Int(0))
              continue;

            if (!fmpz_is_square(
                    q.raw())) {
              continue;
            }

            ++square_survivors;

            Int predicted_index;

            fmpz_sqrt(
                predicted_index.raw(),
                q.raw());

            IntMat tuple =
                start.basis;

            tuple.SetRow(
                replace,
                shell.RowVec(vi));

            const State candidate =
                SaturatedState(
                    a,
                    tuple);

            if (candidate.det !=
                target) {
              continue;
            }

            ++det948_hits;

            const long roots =
                OrthogonalCount(
                    candidate.basis,
                    a.roots_dual);

            /*
             * The Elkies complement must be rootless.
             *
             * Root shell is tiny; norm-4 shell is huge.  Do not scan
             * ~100k norm-4 representatives for a candidate that is
             * already disqualified by an orthogonal root.
             */
            long norm4 = -1;

            if (roots == 0) {
              norm4 =
                  OrthogonalCount(
                      candidate.basis,
                      a.norm4_dual);
            }

            const bool new_gram =
                distinct_grams
                    .insert(
                        GramKey(
                            candidate.gram))
                    .second;

            const bool new_fp =
                fingerprints
                    .insert(
                        {roots, norm4})
                    .second;

            LOG(INFO)
                << "TARGET DET HIT"
                << " row="
                << replace
                << " shell="
                << shell_name
                << " index="
                << vi
                << " raw_det="
                << raw_det
                << " raw_index="
                << predicted_index
                << " roots="
                << roots
                << " norm4="
                << norm4
                << " new_fp="
                << new_fp
                << " new_gram="
                << new_gram;

            if (new_fp ||
                (roots == 0 &&
                 std::abs(
                     norm4 -
                     FLAGS_target_norm4) <= 10)) {

              SaveHit(
                  candidate,
                  roots,
                  norm4,
                  replace,
                  shell_name,
                  vi,
                  raw_det,
                  predicted_index,
                  det948_hits);
            }

            if (roots == 0 &&
                norm4 ==
                    FLAGS_target_norm4) {

              LOG(INFO)
                  << "************************************************";

              LOG(INFO)
                  << "*** JACKPOT det=948 roots=0 norm4=1311 ***";

              LOG(INFO)
                  << "************************************************";

              SaveHit(
                  candidate,
                  roots,
                  norm4,
                  replace,
                  shell_name,
                  vi,
                  raw_det,
                  predicted_index,
                  det948_hits);

              return true;
            }
          }
        }

        return false;
      };



  /*
   * Precompute the only shell vectors that can possibly be orthogonal
   * after changing each individual row.
   */
  const EligibleByRow root_eligible =
      BuildEligibleByChangedRow(
          start.basis,
          a.roots_dual);

  const EligibleByRow norm4_eligible =
      BuildEligibleByChangedRow(
          start.basis,
          a.norm4_dual);

  for (int r = 0; r < kRank; ++r) {

    LOG(INFO)
        << "ELIGIBLE row="
        << r
        << " roots="
        << root_eligible[r].size()
        << "/"
        << a.roots.Rows()
        << " norm4="
        << norm4_eligible[r].size()
        << "/"
        << a.norm4.Rows();
  }

  auto scan_additive_shell =
      [&](const char* shell_name,
          const IntMat& shell,
          const IntMat& shell_dual,
          long shell_norm) {

        for (int replace = 0;
             replace < kRank;
             ++replace) {

          LOG(INFO)
              << "===== additive "
              << shell_name
              << " row "
              << replace
              << " =====";

          for (long vi = 0;
               vi < shell.Rows();
               ++vi) {

            for (int sign : {-1, +1}) {

              ++tested;

              /*
               * k_r' = k_r + sign*v.
               *
               * Cross terms:
               *
               *   <k_r',k_j>
               *     = <k_r,k_j> + sign*<v,k_j>
               *
               * Diagonal:
               *
               *   |k_r'|^2
               *     = |k_r|^2
               *       + 2 sign <k_r,v>
               *       + |v|^2.
               */
              IntMat raw =
                  start.gram;

              const Int self_ip =
                  Dot(
                      start.basis.RowPtr(replace),
                      shell_dual.RowPtr(vi),
                      kDim);

              raw.Set(
                  replace,
                  replace,
                  start.gram(replace, replace)
                      + Int(2 * sign) * self_ip
                      + Int(shell_norm));

              for (int j = 0;
                   j < kRank;
                   ++j) {

                if (j == replace)
                  continue;

                const Int ip =
                    Dot(
                        start.basis.RowPtr(j),
                        shell_dual.RowPtr(vi),
                        kDim);

                const Int x =
                    start.gram(replace, j)
                    + Int(sign) * ip;

                raw.Set(
                    replace,
                    j,
                    x);

                raw.Set(
                    j,
                    replace,
                    x);
              }

              const Int raw_det =
                  Det(raw);

              if (raw_det <= Int(0))
                continue;

              if (raw_det %
                      Int(FLAGS_target_det) !=
                  Int(0)) {
                continue;
              }

              const Int q =
                  raw_det /
                  Int(FLAGS_target_det);

              if (q <= Int(0))
                continue;

              if (!fmpz_is_square(
                      q.raw())) {
                continue;
              }

              ++square_survivors;

              Int predicted_index;

              fmpz_sqrt(
                  predicted_index.raw(),
                  q.raw());

              IntMat tuple =
                  start.basis;

              IntVec changed =
                  start.basis.RowVec(
                      replace);

              changed.AddScaled(
                  shell.RowVec(vi),
                  Int(sign));

              tuple.SetRow(
                  replace,
                  changed);

              const State candidate =
                  SaturatedState(
                      a,
                      tuple);

              if (candidate.det !=
                  Int(FLAGS_target_det)) {
                continue;
              }

              ++det948_hits;

              /*
               * Cheap disqualifier first.
               */
              /*
               * Exact incremental orthogonality count.
               *
               * Saturation changes the Z-basis but not span_Q(tuple), so
               * orthogonality can be tested against the raw tuple:
               * six unchanged rows + this one changed row.
               */
              const long roots =
                  OrthogonalCountOneChangedRow(
                      changed,
                      a.roots_dual,
                      root_eligible[replace]);

              long norm4 = -1;

              if (roots <= FLAGS_norm4_max_roots) {

                norm4 =
                    OrthogonalCountOneChangedRow(
                        changed,
                        a.norm4_dual,
                        norm4_eligible[replace]);
              }

              const bool new_gram =
                  distinct_grams
                      .insert(
                          GramKey(
                              candidate.gram))
                      .second;

              const bool new_fp =
                  fingerprints
                      .insert(
                          {roots, norm4})
                      .second;

              if (new_fp ||
                  roots == 0) {

                LOG(INFO)
                    << "ADDITIVE TARGET DET HIT"
                    << " row="
                    << replace
                    << " shell="
                    << shell_name
                    << " index="
                    << vi
                    << " sign="
                    << sign
                    << " raw_det="
                    << raw_det
                    << " raw_index="
                    << predicted_index
                    << " roots="
                    << roots
                    << " norm4="
                    << norm4
                    << " new_fp="
                    << new_fp
                    << " new_gram="
                    << new_gram;
              }


              /*
               * Preserve one representative of every newly discovered
               * rootless exact-det fingerprint.  These become nodes for
               * the next local-search generation.
               */
              const bool promising =
                  roots <= FLAGS_save_max_roots &&
                  norm4 >= 0 &&
                  std::abs(norm4 - FLAGS_target_norm4)
                      <= FLAGS_save_n4_band;

              const bool new_embedding =
                  promising &&
                  emitted_embeddings
                      .insert(
                          EmbeddedKey(candidate.basis))
                      .second;

              if (new_embedding) {

                const std::string label =
                    std::string(shell_name) +
                    (sign > 0 ? "_add" : "_sub");

                SaveHit(
                    candidate,
                    roots,
                    norm4,
                    replace,
                    label.c_str(),
                    vi,
                    raw_det,
                    predicted_index,
                    det948_hits);
              }

              if (roots == 0 &&
                  norm4 ==
                      FLAGS_target_norm4) {

                LOG(INFO)
                    << "************************************************";

                LOG(INFO)
                    << "*** JACKPOT det=948 roots=0 norm4=1311 ***";

                LOG(INFO)
                    << "************************************************";

                const std::string label =
                    std::string(shell_name) +
                    (sign > 0 ? "_add" : "_sub");

                SaveHit(
                    candidate,
                    roots,
                    norm4,
                    replace,
                    label.c_str(),
                    vi,
                    raw_det,
                    predicted_index,
                    det948_hits);

                return true;
              }
            }
          }
        }

        return false;
      };

  if (FLAGS_scan_replacement) {

    if (scan_shell(
            "root",
            a.roots,
            a.roots_dual,
            2)) {
      return 0;
    }

    if (scan_shell(
            "norm4",
            a.norm4,
            a.norm4_dual,
            4)) {
      return 0;
    }
  }

  if (FLAGS_scan_additive) {

    if (scan_additive_shell(
            "root",
            a.roots,
            a.roots_dual,
            2)) {
      return 0;
    }

    if (scan_additive_shell(
            "norm4",
            a.norm4,
            a.norm4_dual,
            4)) {
      return 0;
    }
  }

  
  if (FLAGS_scan_two_replacement) {

    using TwoSig5 = std::array<long, 5>;

    const EligibleByPair root_pair_eligible =
        BuildEligibleByChangedPair(
            start.basis,
            a.roots_dual);

    const EligibleByPair norm4_pair_eligible =
        BuildEligibleByChangedPair(
            start.basis,
            a.norm4_dual);

    std::filesystem::create_directories(
        FLAGS_output_dir);

    const std::string survivor_path =
        FLAGS_output_dir +
        "/two-row-signature-survivors.tsv";

    std::ofstream two_out(survivor_path);
    CHECK(two_out.good());

    two_out
        << "r1\tr2\th\tbi\tbj\traw_det\tindex"
        << "\tbucket_i\tbucket_j\tpair_upper_bound";

    for (int q = 0; q < 5; ++q)
      two_out << "\tb" << q;
    for (int q = 0; q < 5; ++q)
      two_out << "\tc" << q;
    two_out << "\n";

    uint64_t two_sig_pairs = 0;
    uint64_t two_tests = 0;
    uint64_t two_survivors = 0;
    uint64_t two_pair_upper = 0;
    uint64_t two_concrete_pairs = 0;
    uint64_t two_det_hits = 0;

    for (int r1 = 0; r1 < kRank; ++r1) {
      for (int r2 = r1 + 1; r2 < kRank; ++r2) {

        /*
         * Concrete-expansion mode is intentionally one selected cell.
         * Signature-only mode still scans the entire census.
         */
        if (!FLAGS_two_add_signature_only &&
            FLAGS_two_add_expand_r1 >= 0 &&
            FLAGS_two_add_expand_r2 >= 0 &&
            (r1 != FLAGS_two_add_expand_r1 ||
             r2 != FLAGS_two_add_expand_r2)) {
          continue;
        }

        std::array<int, 5> fixed_rows{};
        {
          int q = 0;
          for (int r = 0; r < kRank; ++r) {
            if (r != r1 && r != r2)
              fixed_rows[q++] = r;
          }
          CHECK_EQ(q, 5);
        }

        IntMat A(5, 5);
        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j < 5; ++j) {
            A.Set(
                i, j,
                start.gram(
                    fixed_rows[i],
                    fixed_rows[j]));
          }
        }

        const Int d_exact = Det(A);
        CHECK(d_exact.fits_slong_p());
        const long d = d_exact.get_si();
        CHECK_GT(d, 0);

        const IntMat adj_exact = Adjugate(A);

        std::array<std::array<long, 5>, 5> adj{};
        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j < 5; ++j) {
            CHECK(adj_exact(i,j).fits_slong_p());
            adj[i][j] = adj_exact(i,j).get_si();
          }
        }

        struct Bucket {
          TwoSig5 sig{};
          uint64_t count = 0;
          __int128 quad = 0;
          std::vector<long> members;
        };

        std::map<TwoSig5, std::vector<long>> member_map;

        for (long vi = 0; vi < a.norm4.Rows(); ++vi) {
          TwoSig5 sig{};

          for (int q = 0; q < 5; ++q) {
            const Int ip =
                Dot(
                    start.basis.RowPtr(fixed_rows[q]),
                    a.norm4_dual.RowPtr(vi),
                    kDim);

            CHECK(ip.fits_slong_p());
            sig[q] = ip.get_si();
          }

          member_map[sig].push_back(vi);
        }

        std::vector<Bucket> buckets;
        buckets.reserve(member_map.size());

        for (const auto& [sig, members] : member_map) {
          Bucket b;
          b.sig = sig;
          b.count = members.size();
          b.members = members;

          __int128 x = 0;
          for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
              x +=
                  (__int128)sig[i] *
                  adj[i][j] *
                  sig[j];
            }
          }

          b.quad = x;
          buckets.push_back(std::move(b));
        }

        const bool expand_pair =
            !FLAGS_two_signature_only &&
            r1 == FLAGS_two_expand_r1 &&
            r2 == FLAGS_two_expand_r2;

        const uint64_t n = buckets.size();
        two_sig_pairs += n * (n + 1) / 2;

        uint64_t row_survivors = 0;
        uint64_t row_pair_upper = 0;

        std::array<uint64_t, 7> h_survivors{};
        std::array<uint64_t, 7> h_pair_upper{};

        for (size_t bi = 0; bi < buckets.size(); ++bi) {

          const Bucket& B = buckets[bi];

          const __int128 left =
              (__int128)4 * d -
              B.quad;

          for (size_t bj = bi; bj < buckets.size(); ++bj) {

            const Bucket& C = buckets[bj];

            const __int128 right =
                (__int128)4 * d -
                C.quad;

            __int128 z = 0;
            for (int i = 0; i < 5; ++i) {
              for (int j = 0; j < 5; ++j) {
                z +=
                    (__int128)B.sig[i] *
                    adj[i][j] *
                    C.sig[j];
              }
            }

            uint64_t concrete_upper = 0;
            if (bi == bj) {
              concrete_upper =
                  B.count * (B.count - 1) / 2;
            } else {
              concrete_upper =
                  B.count * C.count;
            }

            for (int h = -3; h <= 3; ++h) {

              ++two_tests;

              const __int128 cross =
                  (__int128)h * d -
                  z;

              const __int128 numerator =
                  left * right -
                  cross * cross;

              if (numerator <= 0)
                continue;

              if (numerator % d != 0)
                continue;

              const __int128 raw_det_128 =
                  numerator / d;

              if (raw_det_128 >
                  static_cast<__int128>(
                      std::numeric_limits<long>::max())) {
                continue;
              }

              const long raw_det_long =
                  static_cast<long>(raw_det_128);

              if (raw_det_long <= 0)
                continue;

              const Int raw_det(raw_det_long);

              if (raw_det %
                      Int(FLAGS_target_det) !=
                  Int(0)) {
                continue;
              }

              const Int q =
                  raw_det /
                  Int(FLAGS_target_det);

              if (q <= Int(0) ||
                  !fmpz_is_square(q.raw())) {
                continue;
              }

              Int predicted_index;
              fmpz_sqrt(
                  predicted_index.raw(),
                  q.raw());

              ++row_survivors;
              ++two_survivors;
              row_pair_upper += concrete_upper;
              two_pair_upper += concrete_upper;

              const int hi = h + 3;
              ++h_survivors[hi];
              h_pair_upper[hi] += concrete_upper;

              two_out
                  << r1 << "\t"
                  << r2 << "\t"
                  << h << "\t"
                  << bi << "\t"
                  << bj << "\t"
                  << raw_det_long << "\t"
                  << predicted_index << "\t"
                  << B.count << "\t"
                  << C.count << "\t"
                  << concrete_upper;

              for (int zq = 0; zq < 5; ++zq)
                two_out << "\t" << B.sig[zq];
              for (int zq = 0; zq < 5; ++zq)
                two_out << "\t" << C.sig[zq];
              two_out << "\n";

              if (!expand_pair ||
                  h != FLAGS_two_expand_h) {
                continue;
              }

              const Int wanted_h(h);

              for (size_t ai = 0;
                   ai < B.members.size();
                   ++ai) {

                const size_t aj_start =
                    bi == bj ? ai + 1 : 0;

                for (size_t aj = aj_start;
                     aj < C.members.size();
                     ++aj) {

                  const long ui = B.members[ai];
                  const long vi = C.members[aj];

                  if (Dot(
                          a.norm4.RowPtr(ui),
                          a.norm4_dual.RowPtr(vi),
                          kDim) != wanted_h) {
                    continue;
                  }

                  ++two_concrete_pairs;

                  IntMat tuple =
                      start.basis;

                  tuple.SetRow(
                      r1,
                      a.norm4.RowVec(ui));

                  tuple.SetRow(
                      r2,
                      a.norm4.RowVec(vi));

                  const IntMat raw_gram =
                      MetricGram(
                          tuple,
                          a.gram);

                  const Int concrete_raw_det =
                      Det(raw_gram);

                  CHECK_EQ(
                      concrete_raw_det,
                      raw_det);

                  const State candidate =
                      SaturatedState(
                          a,
                          tuple);

                  if (candidate.det !=
                      Int(FLAGS_target_det)) {
                    continue;
                  }

                  ++two_det_hits;

                  const IntVec changed1 =
                      a.norm4.RowVec(ui);
                  const IntVec changed2 =
                      a.norm4.RowVec(vi);

                  const long roots =
                      OrthogonalCountTwoChangedRows(
                          changed1,
                          changed2,
                          a.roots_dual,
                          root_pair_eligible[r1][r2]);

                  long norm4 = -1;

                  if (roots <=
                      FLAGS_norm4_max_roots) {
                    norm4 =
                        OrthogonalCountTwoChangedRows(
                            changed1,
                            changed2,
                            a.norm4_dual,
                            norm4_pair_eligible[r1][r2]);
                  }

                  LOG(INFO)
                      << "TWO_ROW TARGET DET HIT"
                      << " rows=" << r1 << "," << r2
                      << " h=" << h
                      << " u=" << ui
                      << " v=" << vi
                      << " raw_det=" << raw_det
                      << " index=" << predicted_index
                      << " roots=" << roots
                      << " norm4=" << norm4;

                  const bool promising =
                      roots <= FLAGS_save_max_roots &&
                      norm4 >= 0 &&
                      std::abs(
                          norm4 -
                          FLAGS_target_norm4) <=
                          FLAGS_save_n4_band;

                  if (promising) {
                    SaveTwoRowHit(
                        candidate,
                        roots,
                        norm4,
                        r1,
                        r2,
                        ui,
                        vi,
                        h,
                        raw_det,
                        predicted_index,
                        two_det_hits);
                  }

                  if (roots == 0 &&
                      norm4 ==
                          FLAGS_target_norm4) {

                    LOG(INFO)
                        << "**********************************************";
                    LOG(INFO)
                        << "*** JACKPOT det=948 roots=0 norm4=1311 ***";
                    LOG(INFO)
                        << "**********************************************";

                    SaveTwoRowHit(
                        candidate,
                        roots,
                        norm4,
                        r1,
                        r2,
                        ui,
                        vi,
                        h,
                        raw_det,
                        predicted_index,
                        two_det_hits);

                    return 0;
                  }
                }
              }
            }
          }
        }

        LOG(INFO)
            << "TWO_ROW_FILTER rows="
            << r1 << "," << r2
            << " signatures="
            << buckets.size()
            << " survivors="
            << row_survivors
            << " pair_upper_bound="
            << row_pair_upper;

        for (int h = -3; h <= 3; ++h) {
          const int hi = h + 3;

          if (h_survivors[hi] == 0)
            continue;

          LOG(INFO)
              << "TWO_ROW_H rows="
              << r1 << "," << r2
              << " h=" << h
              << " survivors="
              << h_survivors[hi]
              << " pair_upper_bound="
              << h_pair_upper[hi];
        }
      }
    }

    two_out.close();

    LOG(INFO)
        << "TWO_ROW_FILTER SUMMARY"
        << " signature_pairs="
        << two_sig_pairs
        << " tests="
        << two_tests
        << " survivors="
        << two_survivors
        << " pair_upper_bound="
        << two_pair_upper
        << " concrete_pairs="
        << two_concrete_pairs
        << " det948_hits="
        << two_det_hits
        << " file="
        << survivor_path;
  }


  /*
   * TRUE simultaneous two-row additive repair:
   *
   *   x = k_r1 + s1*u
   *   y = k_r2 + s2*v
   *
   * with u,v ambient norm-4 vectors.
   *
   * Keep the other five rows fixed.  A candidate vector is bucketed by:
   *   - its five inner products with the fixed rows;
   *   - its inner product with k_r1;
   *   - its inner product with k_r2.
   *
   * For a bucket pair, signs s1,s2 and h=<u,v>, the complete raw 7x7
   * determinant is determined by the Schur complement, so concrete u,v
   * pairs are expanded only for determinant-square survivors.
   */
  if (FLAGS_scan_two_additive) {

    // Root-walk frontier caps for this parent state.
    std::unordered_map<uint64_t, int> rootwalk_main_saved;
    std::unordered_map<uint64_t, int> rootwalk_side_saved;
    std::unordered_map<uint64_t, int> rootwalk_two_saved;

    uint64_t rootwalk_considered = 0;
    uint64_t rootwalk_saved = 0;


    using AddSig8 = std::array<long, 8>;

    /*
     * Track one old root as an optional signature heuristic.
     *
     * Rootless parents are allowed (unique_old_root == -1).
     * Correctness for rooted parents is enforced later using the complete
     * parent_root_mask, not this single-root heuristic.
     */
    long unique_old_root = -1;
    uint64_t parent_root_mask = 0;

    for (long ri = 0; ri < a.roots.Rows(); ++ri) {

      bool orth = true;

      for (int r = 0; r < kRank; ++r) {
        if (Dot(
                start.basis.RowPtr(r),
                a.roots_dual.RowPtr(ri),
                kDim) != Int(0)) {
          orth = false;
          break;
        }
      }

      if (orth) {
        
        if (unique_old_root < 0) {
          unique_old_root = ri;
        }
        parent_root_mask |= (uint64_t{1} << ri);
      }
    }

    parent_root_mask =
        CurrentRootMask(start.basis, a.roots_dual);

    CHECK_EQ(
        __builtin_popcountll(parent_root_mask),
        start_roots)
        << "parent root mask/count mismatch";

    LOG(INFO)
        << "TWO_ADD unique old root index="
        << unique_old_root
        << " parent_root_mask=0x"
        << std::hex << parent_root_mask << std::dec;

    const EligibleByPair root_pair_eligible =
        BuildEligibleByChangedPair(
            start.basis,
            a.roots_dual);

    const EligibleByPair norm4_pair_eligible =
        BuildEligibleByChangedPair(
            start.basis,
            a.norm4_dual);

    std::filesystem::create_directories(
        FLAGS_output_dir);

    const std::string survivor_path =
        FLAGS_output_dir +
        "/two-add-signature-survivors.tsv";

    std::ofstream out(survivor_path);
    CHECK(out.good());

    out << "r1\tr2\ts1\ts2\th\tbi\tbj\traw_det\tindex"
        << "\tbucket_i\tbucket_j\tpair_upper_bound";
    for (int q = 0; q < 8; ++q)
      out << "\tb" << q;
    for (int q = 0; q < 8; ++q)
      out << "\tc" << q;
    out << "\n";

    uint64_t total_sig_tests = 0;
    uint64_t total_survivors = 0;
    uint64_t total_pair_upper = 0;
    uint64_t actual_h_pairs = 0;
    uint64_t old_root_killing_pairs = 0;
    uint64_t pre_rootless_pairs = 0;
    uint64_t concrete_pairs = 0;
    uint64_t det_hits = 0;

    for (int r1 = 0; r1 < kRank; ++r1) {
      for (int r2 = r1 + 1; r2 < kRank; ++r2) {

        if (!FLAGS_two_add_signature_only &&
            !FLAGS_two_add_expand_all_index1 &&
            FLAGS_two_add_expand_r1 >= 0 &&
            FLAGS_two_add_expand_r2 >= 0 &&
            (r1 != FLAGS_two_add_expand_r1 ||
             r2 != FLAGS_two_add_expand_r2)) {
          continue;
        }

        if (FLAGS_two_add_only_r1 >= 0 &&
            r1 != FLAGS_two_add_only_r1) {
          continue;
        }

        if (FLAGS_two_add_only_r2 >= 0 &&
            r2 != FLAGS_two_add_only_r2) {
          continue;
        }

        std::array<int, 5> fixed_rows{};
        {
          int q = 0;
          for (int r = 0; r < kRank; ++r) {
            if (r != r1 && r != r2)
              fixed_rows[q++] = r;
          }
          CHECK_EQ(q, 5);
        }

        IntMat A(5, 5);
        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j < 5; ++j) {
            A.Set(
                i, j,
                start.gram(
                    fixed_rows[i],
                    fixed_rows[j]));
          }
        }

        const Int d_exact = Det(A);
        CHECK(d_exact.fits_slong_p());
        const long d = d_exact.get_si();
        CHECK_GT(d, 0);

        const IntMat adj_exact = Adjugate(A);
        std::array<std::array<long, 5>, 5> adj{};

        for (int i = 0; i < 5; ++i) {
          for (int j = 0; j < 5; ++j) {
            CHECK(adj_exact(i,j).fits_slong_p());
            adj[i][j] = adj_exact(i,j).get_si();
          }
        }

        /*
         * Role-specific buckets.  L represents u added to row r1:
         *   sig[5] = <k_r1,u>, sig[6] = <k_r2,u>.
         *
         * R represents v added to row r2:
         *   sig[5] = <k_r2,v>, sig[6] = <k_r1,v>.
         */
        std::map<AddSig8, std::vector<long>> left_map;
        std::map<AddSig8, std::vector<long>> right_map;

        for (long vi = 0; vi < a.norm4.Rows(); ++vi) {

          AddSig8 L{};
          AddSig8 R{};

          for (int q = 0; q < 5; ++q) {
            const Int ip =
                Dot(
                    start.basis.RowPtr(fixed_rows[q]),
                    a.norm4_dual.RowPtr(vi),
                    kDim);
            CHECK(ip.fits_slong_p());
            L[q] = ip.get_si();
            R[q] = ip.get_si();
          }

          {
            const Int self1 =
                Dot(
                    start.basis.RowPtr(r1),
                    a.norm4_dual.RowPtr(vi),
                    kDim);
            const Int other2 =
                Dot(
                    start.basis.RowPtr(r2),
                    a.norm4_dual.RowPtr(vi),
                    kDim);
            CHECK(self1.fits_slong_p());
            CHECK(other2.fits_slong_p());
            L[5] = self1.get_si();
            L[6] = other2.get_si();
          }

          {
            const Int self2 =
                Dot(
                    start.basis.RowPtr(r2),
                    a.norm4_dual.RowPtr(vi),
                    kDim);
            const Int other1 =
                Dot(
                    start.basis.RowPtr(r1),
                    a.norm4_dual.RowPtr(vi),
                    kDim);
            CHECK(self2.fits_slong_p());
            CHECK(other1.fits_slong_p());
            R[5] = self2.get_si();
            R[6] = other1.get_si();
          }


          const Int root_ip =
              unique_old_root >= 0
                  ? Dot(
                        a.norm4.RowPtr(vi),
                        a.roots_dual.RowPtr(unique_old_root),
                        kDim)
                  : Int(0);

          CHECK(root_ip.fits_slong_p());



          L[7] = root_ip.get_si();


          R[7] = root_ip.get_si();



          left_map[L].push_back(vi);


          right_map[R].push_back(vi);
        }

        struct AddBucket {
          AddSig8 sig{};
          uint64_t count = 0;
          std::vector<long> members;
        };

        std::vector<AddBucket> left_buckets;
        std::vector<AddBucket> right_buckets;

        for (const auto& [sig, members] : left_map) {
          AddBucket b;
          b.sig = sig;
          b.count = members.size();
          b.members = members;
          left_buckets.push_back(std::move(b));
        }

        for (const auto& [sig, members] : right_map) {
          AddBucket b;
          b.sig = sig;
          b.count = members.size();
          b.members = members;
          right_buckets.push_back(std::move(b));
        }

        const bool expand_pair =
            !FLAGS_two_add_signature_only &&
            (FLAGS_two_add_expand_all_index1 ||
             (r1 == FLAGS_two_add_expand_r1 &&
              r2 == FLAGS_two_add_expand_r2));

        uint64_t row_survivors = 0;
        uint64_t row_pair_upper = 0;

        std::array<uint64_t, 7> h_survivors{};
        std::array<uint64_t, 7> h_pair_upper{};

        for (size_t bi = 0; bi < left_buckets.size(); ++bi) {

          const AddBucket& B = left_buckets[bi];

          for (size_t bj = 0; bj < right_buckets.size(); ++bj) {

            const AddBucket& C = right_buckets[bj];

            /*
             * Necessary condition for destroying the unique old root.
             * Initially every section row is orthogonal to that root.
             * If both additive vectors are also orthogonal to it, the
             * root certainly survives.
             */
            if (unique_old_root >= 0 &&
                          B.sig[7] == 0 &&
                          C.sig[7] == 0) {
                        continue;
                      }


            /*
             * Concrete upper bound.  Unlike replacement, u and v are
             * role-specific and ordered; only exclude the exact same vector
             * when it would make the final rows dependent later (exact rank
             * and determinant checks still certify every concrete hit).
             */
            const uint64_t pair_upper =
                B.count * C.count;

            for (int s1 : {-1, +1}) {
              for (int s2 : {-1, +1}) {

                if (!FLAGS_two_add_signature_only &&
                    !FLAGS_two_add_expand_all_index1 &&
                    (s1 != FLAGS_two_add_expand_s1 ||
                     s2 != FLAGS_two_add_expand_s2)) {
                  continue;
                }

                std::array<long,5> b{};
                std::array<long,5> c{};

                for (int q = 0; q < 5; ++q) {
                  b[q] =
                      start.gram(r1, fixed_rows[q]).get_si()
                      + s1 * B.sig[q];

                  c[q] =
                      start.gram(r2, fixed_rows[q]).get_si()
                      + s2 * C.sig[q];
                }

                __int128 xq = 0;
                __int128 yq = 0;
                __int128 zq = 0;

                for (int i = 0; i < 5; ++i) {
                  for (int j = 0; j < 5; ++j) {
                    xq +=
                        (__int128)b[i] *
                        adj[i][j] *
                        b[j];
                    yq +=
                        (__int128)c[i] *
                        adj[i][j] *
                        c[j];
                    zq +=
                        (__int128)b[i] *
                        adj[i][j] *
                        c[j];
                  }
                }

                const long alpha =
                    start.gram(r1,r1).get_si()
                    + 2 * s1 * B.sig[5]
                    + 4;

                const long beta =
                    start.gram(r2,r2).get_si()
                    + 2 * s2 * C.sig[5]
                    + 4;

                const __int128 left =
                    (__int128)alpha * d -
                    xq;

                const __int128 right =
                    (__int128)beta * d -
                    yq;

                for (int h = -3; h <= 3; ++h) {

                  if (!FLAGS_two_add_signature_only &&
                      !FLAGS_two_add_expand_all_index1 &&
                      h != FLAGS_two_add_expand_h) {
                    continue;
                  }

                  ++total_sig_tests;

                  const long gamma =
                      start.gram(r1,r2).get_si()
                      + s1 * B.sig[6]
                      + s2 * C.sig[6]
                      + s1 * s2 * h;

                  const __int128 cross =
                      (__int128)gamma * d -
                      zq;

                  const __int128 numerator =
                      left * right -
                      cross * cross;

                  if (numerator <= 0)
                    continue;

                  if (numerator % d != 0)
                    continue;

                  const __int128 raw_det_128 =
                      numerator / d;

                  if (raw_det_128 >
                      static_cast<__int128>(
                          std::numeric_limits<long>::max())) {
                    continue;
                  }

                  const long raw_det_long =
                      static_cast<long>(raw_det_128);

                  if (raw_det_long <= 0)
                    continue;

                  const Int raw_det(raw_det_long);

                  if (raw_det %
                          Int(FLAGS_target_det) !=
                      Int(0)) {
                    continue;
                  }

                  const Int quotient =
                      raw_det /
                      Int(FLAGS_target_det);

                  if (quotient <= Int(0) ||
                      !fmpz_is_square(quotient.raw())) {
                    continue;
                  }

                  Int predicted_index;
                  fmpz_sqrt(
                      predicted_index.raw(),
                      quotient.raw());

                  if (!FLAGS_two_add_signature_only &&
                      FLAGS_two_add_index1_only &&
                      predicted_index != Int(1)) {
                    continue;
                  }

                  if (!FLAGS_two_add_signature_only &&
                      FLAGS_two_add_exact_index > 0 &&
                      predicted_index !=
                          Int(FLAGS_two_add_exact_index)) {
                    continue;
                  }

                  ++total_survivors;
                  ++row_survivors;
                  total_pair_upper += pair_upper;
                  row_pair_upper += pair_upper;

                  const int hi = h + 3;
                  ++h_survivors[hi];
                  h_pair_upper[hi] += pair_upper;

                  out << r1 << "\t"
                      << r2 << "\t"
                      << s1 << "\t"
                      << s2 << "\t"
                      << h << "\t"
                      << bi << "\t"
                      << bj << "\t"
                      << raw_det_long << "\t"
                      << predicted_index << "\t"
                      << B.count << "\t"
                      << C.count << "\t"
                      << pair_upper;

                  for (int q = 0; q < 8; ++q)
                    out << "\t" << B.sig[q];
                  for (int q = 0; q < 8; ++q)
                    out << "\t" << C.sig[q];
                  out << "\n";

                  if (!expand_pair) {
                    continue;
                  }

                  if (!FLAGS_two_add_expand_all_index1 &&
                      (h != FLAGS_two_add_expand_h ||
                       s1 != FLAGS_two_add_expand_s1 ||
                       s2 != FLAGS_two_add_expand_s2)) {
                    continue;
                  }

                  const Int wanted_h(h);

                  for (long ui : B.members) {
                    for (long vi : C.members) {

                      if (Dot(
                              a.norm4.RowPtr(ui),
                              a.norm4_dual.RowPtr(vi),
                              kDim) != wanted_h) {
                        continue;
                      }

                      ++actual_h_pairs;

                      /*
                       * Necessary condition for becoming rootless:
                       * at least one changed row must cease to be
                       * orthogonal to the unique starting root.
                       *
                       * Since <k_r,root>=0 initially, this is exactly
                       * <u,root> != 0 or <v,root> != 0.
                       */
                      if (unique_old_root >= 0) {
                        const Int u_old_root =
                            Dot(
                                a.norm4.RowPtr(ui),
                                a.roots_dual.RowPtr(unique_old_root),
                                kDim);

                        const Int v_old_root =
                            Dot(
                                a.norm4.RowPtr(vi),
                                a.roots_dual.RowPtr(unique_old_root),
                                kDim);

                        if (u_old_root == Int(0) &&
                            v_old_root == Int(0)) {
                          continue;
                        }

                        ++old_root_killing_pairs;
                      }
                      ++concrete_pairs;

                      IntMat tuple =
                          start.basis;

                      IntVec changed1 =
                          start.basis.RowVec(r1);

                      changed1.AddScaled(
                          a.norm4.RowVec(ui),
                          Int(s1));

                      IntVec changed2 =
                          start.basis.RowVec(r2);

                      changed2.AddScaled(
                          a.norm4.RowVec(vi),
                          Int(s2));

                      tuple.SetRow(r1, changed1);
                      tuple.SetRow(r2, changed2);

                      /*
                       * Exact rootlessness depends only on the rational
                       * row span.  Test it before Gram construction,
                       * saturation and determinant certification.
                       *
                       * root_pair_eligible already contains every ambient
                       * root that could possibly become orthogonal after
                       * changing these two rows.
                       */
                      /*
                       * Directed requirement: every root already present
                       * in the parent must be destroyed by this move.
                       *
                       * For the current one-root search this means the
                       * old root may not remain orthogonal to both changed
                       * rows.
                       */
                      bool parent_roots_destroyed = true;

                      if (parent_root_mask != 0) {
                        for (long pri = 0; pri < 60; ++pri) {
                          if ((parent_root_mask &
                               (uint64_t{1} << pri)) == 0) {
                            continue;
                          }

                          const bool survives =
                              Dot(
                                  changed1.data(),
                                  a.roots_dual.RowPtr(pri),
                                  kDim) == Int(0) &&
                              Dot(
                                  changed2.data(),
                                  a.roots_dual.RowPtr(pri),
                                  kDim) == Int(0);

                          if (survives) {
                            parent_roots_destroyed = false;
                            break;
                          }
                        }
                      }

                      if (!parent_roots_destroyed) {
                        continue;
                      }

                      const uint64_t new_root_mask =
                          RootMaskTwoChangedRows(
                              changed1,
                              changed2,
                              a.roots_dual,
                              root_pair_eligible[r1][r2]);

                      const long pre_roots =
                          __builtin_popcountll(new_root_mask);

                      if (FLAGS_two_add_require_rootless &&
                          pre_roots != 0) {
                        continue;
                      }

                      if (pre_roots == 0) {
                        ++pre_rootless_pairs;
                      }

                      const IntMat raw_gram =
                          MetricGram(
                              tuple,
                              a.gram);

                      const Int concrete_raw_det =
                          Det(raw_gram);

                      CHECK_EQ(
                          concrete_raw_det,
                          raw_det);

                      const State candidate =
                          SaturatedState(
                              a,
                              tuple);

                      if (candidate.det !=
                          Int(FLAGS_target_det)) {
                        continue;
                      }

                      ++det_hits;

                      const long roots = pre_roots;

                      long norm4 = -1;

                      if (roots <=
                          FLAGS_norm4_max_roots) {
                        norm4 =
                            OrthogonalCountTwoChangedRows(
                                changed1,
                                changed2,
                                a.norm4_dual,
                                norm4_pair_eligible[r1][r2]);
                      }

                      if (roots <= 2 || norm4 >= 0) {
                        LOG(INFO)
                            << "TWO_ADD TARGET DET HIT"
                            << " rows=" << r1 << "," << r2
                            << " signs=" << s1 << "," << s2
                            << " h=" << h
                            << " u=" << ui
                            << " v=" << vi
                            << " raw_det=" << raw_det
                            << " index=" << predicted_index
                            << " roots=" << roots
                            << " rootmask="
                          << std::hex << new_root_mask << std::dec
                          << " norm4=" << norm4;
                      }

                      const bool promising =
                          roots <= FLAGS_save_max_roots &&
                          norm4 >= 0 &&
                          std::abs(
                              norm4 -
                              FLAGS_target_norm4) <=
                              FLAGS_save_n4_band;

                      if (promising) {
                        bool save_this = true;

                        if (FLAGS_two_add_rootwalk) {
                          ++rootwalk_considered;
                          save_this = false;

                          /*
                           * Main plateau:
                           *   roots=1, norm4=1311
                           *
                           * Keep a handful of distinct embeddings for
                           * each resulting surviving-root mask.
                           */
                          if (roots == 1 &&
                              norm4 == FLAGS_target_norm4) {
                            int& count =
                                rootwalk_main_saved[new_root_mask];

                            if (count <
                                FLAGS_two_add_rootwalk_per_mask) {
                              ++count;
                              save_this = true;
                            }
                          }

                          /*
                           * Small sideways excursions:
                           *   roots=1, norm4=1310 or 1312.
                           *
                           * These let the graph walk briefly leave the
                           * exact 1311 plateau.
                           */
                          else if (
                              roots == 1 &&
                              std::abs(
                                  norm4 -
                                  FLAGS_target_norm4) == 1) {
                            int& count =
                                rootwalk_side_saved[new_root_mask];

                            if (count <
                                FLAGS_two_add_rootwalk_side_per_mask) {
                              ++count;
                              save_this = true;
                            }
                          }

                          /*
                           * Two-root bridge states exactly at 1311.
                           */
                          else if (
                              roots == 2 &&
                              norm4 == FLAGS_target_norm4) {
                            int& count =
                                rootwalk_two_saved[new_root_mask];

                            if (count <
                                FLAGS_two_add_rootwalk_two_per_mask) {
                              ++count;
                              save_this = true;
                            }
                          }
                        }

                        if (save_this) {
                          SaveTwoRowHit(
                              candidate,
                              roots,
                              norm4,
                              r1,
                              r2,
                              ui,
                              vi,
                              h,
                              raw_det,
                              predicted_index,
                              det_hits,
                              new_root_mask);

                          if (FLAGS_two_add_rootwalk) {
                            ++rootwalk_saved;
                          }
                        }
                      }

                      if (roots == 0 &&
                          norm4 ==
                              FLAGS_target_norm4) {

                        LOG(INFO)
                            << "**********************************************";
                        LOG(INFO)
                            << "*** JACKPOT TWO-ADD det=948 roots=0 norm4=1311 ***";
                        LOG(INFO)
                            << "**********************************************";

                        SaveTwoRowHit(
                            candidate,
                            roots,
                            norm4,
                            r1,
                            r2,
                            ui,
                            vi,
                            h,
                            raw_det,
                            predicted_index,
                            det_hits,
                            new_root_mask);

                        return 0;
                      }
                    }
                  }
                }
              }
            }
          }
        }

        LOG(INFO)
            << "TWO_ADD_FILTER rows="
            << r1 << "," << r2
            << " left_signatures="
            << left_buckets.size()
            << " right_signatures="
            << right_buckets.size()
            << " survivors="
            << row_survivors
            << " pair_upper_bound="
            << row_pair_upper;

        for (int h = -3; h <= 3; ++h) {
          const int hi = h + 3;
          if (h_survivors[hi] == 0)
            continue;

          LOG(INFO)
              << "TWO_ADD_H rows="
              << r1 << "," << r2
              << " h=" << h
              << " survivors="
              << h_survivors[hi]
              << " pair_upper_bound="
              << h_pair_upper[hi];
        }
      }
    }

    out.close();

    
    if (FLAGS_two_add_rootwalk) {
      LOG(INFO)
          << "TWO_ADD ROOTWALK SUMMARY"
          << " considered=" << rootwalk_considered
          << " saved=" << rootwalk_saved
          << " main_masks=" << rootwalk_main_saved.size()
          << " side_masks=" << rootwalk_side_saved.size()
          << " two_masks=" << rootwalk_two_saved.size();
    }

LOG(INFO)
        << "TWO_ADD_FILTER SUMMARY"
        << " signature_tests="
        << total_sig_tests
        << " survivors="
        << total_survivors
        << " pair_upper_bound="
        << total_pair_upper
        << " actual_h_pairs="
        << actual_h_pairs
        << " old_root_killing_pairs="
        << old_root_killing_pairs
        << " pre_rootless_pairs="
        << pre_rootless_pairs
        << " concrete_pairs="
        << concrete_pairs
        << " det948_hits="
        << det_hits
        << " file="
        << survivor_path;
  }

LOG(INFO)
      << "================ FINAL ================";

  LOG(INFO)
      << "tested="
      << tested;

  LOG(INFO)
      << "square-filter survivors="
      << square_survivors;

  LOG(INFO)
      << "det948 hits="
      << det948_hits
      << " distinct Grams="
      << distinct_grams.size()
      << " fingerprints="
      << fingerprints.size();

  for (const auto& [roots, norm4] :
       fingerprints) {

    LOG(INFO)
        << "fingerprint roots="
        << roots
        << " norm4="
        << norm4;
  }

  return 0;
}

} // namespace
} // namespace antipode


int main(
    int argc,
    char** argv) {

  gflags::ParseCommandLineFlags(
      &argc,
      &argv,
      true);

  const antipode::ScopedLogging
      logging(argv[0]);

  CHECK_EQ(argc, 1);

  return antipode::Run();
}
