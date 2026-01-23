#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

namespace {

enum Mutation : std::size_t { S = 0, D = 1, I = 2, N_MUTS = 3 };

inline std::size_t score_index(std::size_t state,
                               std::size_t j,
                               std::size_t k,
                               std::size_t lseq_plus1,
                               std::size_t lref_plus1) {
    return ((state * lseq_plus1) + j) * lref_plus1 + k;
}

inline std::size_t backtrack_index(std::size_t from_state,
                                   std::size_t to_state,
                                   std::size_t j,
                                   std::size_t k,
                                   std::size_t lseq_plus1,
                                   std::size_t lref_plus1) {
    return ((((from_state * Mutation::N_MUTS) + to_state) * lseq_plus1) + j) * lref_plus1 + k;
}

inline std::array<double, 3> get_state_candidates(const std::vector<double>& score,
                                                  std::size_t j,
                                                  std::size_t k,
                                                  std::size_t lseq_plus1,
                                                  std::size_t lref_plus1,
                                                  std::size_t state,
                                                  double subst,
                                                  const double* open_penalty,
                                                  const double* close_penalty) {
    const auto idx = [&](std::size_t st, std::size_t jj, std::size_t kk) {
        return score_index(st, jj, kk, lseq_plus1, lref_plus1);
    };

    if (state == Mutation::S) {
        return {
            score[idx(Mutation::S, j - 1, k - 1)] + subst,
            score[idx(Mutation::D, j - 1, k - 1)] - close_penalty[k - 1] + subst,
            score[idx(Mutation::I, j - 1, k - 1)] - close_penalty[k] + subst
        };
    }

    if (state == Mutation::D) {
        return {
            score[idx(Mutation::S, j, k - 1)] - open_penalty[k],
            score[idx(Mutation::D, j, k - 1)],
            score[idx(Mutation::I, j, k - 1)]
        };
    }

    // state == Mutation::I
    return {
        score[idx(Mutation::S, j - 1, k)] - open_penalty[k],
        score[idx(Mutation::D, j - 1, k)],
        score[idx(Mutation::I, j - 1, k)]
    };
}

inline std::array<bool, 3> argmax3(const std::array<double, 3>& values, double& maxval) {
    maxval = std::max({values[0], values[1], values[2]});
    return {
        values[0] == maxval,
        values[1] == maxval,
        values[2] == maxval
    };
}

}  // namespace

py::tuple cas9_align_cpp(py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast> seq,
                         py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast> ref,
                         py::array_t<double, py::array::c_style | py::array::forcecast> open_penalty,
                         py::array_t<double, py::array::c_style | py::array::forcecast> close_penalty,
                         py::array_t<double, py::array::c_style | py::array::forcecast> sub_score) {
    auto seq_buf = seq.request();
    auto ref_buf = ref.request();
    auto open_buf = open_penalty.request();
    auto close_buf = close_penalty.request();
    auto sub_buf = sub_score.request();

    if (seq_buf.ndim != 1 || ref_buf.ndim != 1) {
        throw std::invalid_argument("seq and ref must be 1-D arrays");
    }

    const std::size_t Lseq = static_cast<std::size_t>(seq_buf.shape[0]);
    const std::size_t Lref = static_cast<std::size_t>(ref_buf.shape[0]);

    if (open_buf.ndim != 1 || close_buf.ndim != 1) {
        throw std::invalid_argument("open_penalty and close_penalty must be 1-D arrays");
    }

    if (open_buf.shape[0] != static_cast<py::ssize_t>(Lref + 1) ||
        close_buf.shape[0] != static_cast<py::ssize_t>(Lref + 1)) {
        throw std::invalid_argument("penalty arrays must have length Lref + 1");
    }

    if (sub_buf.ndim != 1 || sub_buf.shape[0] != 25) {
        throw std::invalid_argument("sub_score must be a 1-D array of length 25");
    }

    const auto* seq_ptr = static_cast<std::uint8_t*>(seq_buf.ptr);
    const auto* ref_ptr = static_cast<std::uint8_t*>(ref_buf.ptr);
    const auto* open_ptr = static_cast<double*>(open_buf.ptr);
    const auto* close_ptr = static_cast<double*>(close_buf.ptr);
    const auto* sub_ptr = static_cast<double*>(sub_buf.ptr);

    const std::size_t lseq_plus1 = Lseq + 1;
    const std::size_t lref_plus1 = Lref + 1;

    const std::size_t score_size = Mutation::N_MUTS * lseq_plus1 * lref_plus1;
    const std::size_t backtrack_size = Mutation::N_MUTS * Mutation::N_MUTS * lseq_plus1 * lref_plus1;

    const double neg_inf = -std::numeric_limits<double>::infinity();

    std::vector<double> score(score_size, neg_inf);
    std::vector<std::uint8_t> backtrack(backtrack_size, 0);

    // base case
    score[score_index(Mutation::S, 0, 0, lseq_plus1, lref_plus1)] = 0.0;

    // initialize first row (insertions)
    for (std::size_t j = 1; j <= Lseq; ++j) {
        score[score_index(Mutation::I, j, 0, lseq_plus1, lref_plus1)] = -open_ptr[0];
        backtrack[backtrack_index(Mutation::I, Mutation::I, j, 0, lseq_plus1, lref_plus1)] = 1;
    }

    // initialize first column (deletions)
    for (std::size_t k = 1; k <= Lref; ++k) {
        score[score_index(Mutation::D, 0, k, lseq_plus1, lref_plus1)] = -open_ptr[0];
        backtrack[backtrack_index(Mutation::D, Mutation::D, 0, k, lseq_plus1, lref_plus1)] = 1;
    }

    // DP filling
    for (std::size_t j = 1; j <= Lseq; ++j) {
        for (std::size_t k = 1; k <= Lref; ++k) {
            std::uint8_t seq_nt = seq_ptr[j - 1];
            std::uint8_t ref_nt = ref_ptr[k - 1];
            double subst = sub_ptr[seq_nt * 5 + ref_nt];

            for (std::size_t state = 0; state < Mutation::N_MUTS; ++state) {
                auto candidates = get_state_candidates(score, j, k, lseq_plus1, lref_plus1,
                                                       state, subst, open_ptr, close_ptr);
                double best_val;
                auto argmax = argmax3(candidates, best_val);
                score[score_index(state, j, k, lseq_plus1, lref_plus1)] = best_val;

                backtrack[backtrack_index(Mutation::S, state, j, k, lseq_plus1, lref_plus1)] = argmax[Mutation::S];
                backtrack[backtrack_index(Mutation::D, state, j, k, lseq_plus1, lref_plus1)] = argmax[Mutation::D];
                backtrack[backtrack_index(Mutation::I, state, j, k, lseq_plus1, lref_plus1)] = argmax[Mutation::I];
            }
        }
    }

    // Determine final state
    std::array<double, 3> final_scores = {
        score[score_index(Mutation::S, Lseq, Lref, lseq_plus1, lref_plus1)],
        score[score_index(Mutation::D, Lseq, Lref, lseq_plus1, lref_plus1)],
        score[score_index(Mutation::I, Lseq, Lref, lseq_plus1, lref_plus1)]
    };

    double best_score;
    auto final_mask = argmax3(final_scores, best_score);

    std::size_t cur_state = Mutation::S;
    if (final_mask[Mutation::S]) {
        cur_state = Mutation::S;
    } else if (final_mask[Mutation::D]) {
        cur_state = Mutation::D;
    } else if (final_mask[Mutation::I]) {
        cur_state = Mutation::I;
    } else {
        throw std::runtime_error("Failed to determine final alignment state");
    }

    std::vector<std::uint8_t> al_seq;
    std::vector<std::uint8_t> al_ref;
    al_seq.reserve(Lseq + Lref);
    al_ref.reserve(Lseq + Lref);

    std::size_t cur_j = Lseq;
    std::size_t cur_k = Lref;

    while (cur_j > 0 || cur_k > 0) {
        if (cur_state >= Mutation::N_MUTS) {
            throw std::runtime_error("Invalid backtrack state encountered");
        }

        std::array<bool, 3> bt = {
            backtrack[backtrack_index(Mutation::S, cur_state, cur_j, cur_k, lseq_plus1, lref_plus1)] != 0,
            backtrack[backtrack_index(Mutation::D, cur_state, cur_j, cur_k, lseq_plus1, lref_plus1)] != 0,
            backtrack[backtrack_index(Mutation::I, cur_state, cur_j, cur_k, lseq_plus1, lref_plus1)] != 0
        };

        if (cur_state == Mutation::S) {
            al_seq.push_back(cur_j > 0 ? seq_ptr[cur_j - 1] : 0);
            al_ref.push_back(cur_k > 0 ? ref_ptr[cur_k - 1] : 0);
            if (cur_j > 0) { --cur_j; }
            if (cur_k > 0) { --cur_k; }

            if (bt[Mutation::S]) {
                cur_state = Mutation::S;
            } else if (bt[Mutation::D]) {
                cur_state = Mutation::D;
            } else if (bt[Mutation::I]) {
                cur_state = Mutation::I;
            } else {
                throw std::runtime_error("Backtrack terminated unexpectedly in substitution state");
            }
        } else if (cur_state == Mutation::D) {
            al_seq.push_back(0);
            al_ref.push_back(cur_k > 0 ? ref_ptr[cur_k - 1] : 0);
            if (cur_k > 0) { --cur_k; }

            if (bt[Mutation::D]) {
                cur_state = Mutation::D;
            } else if (bt[Mutation::S]) {
                cur_state = Mutation::S;
            } else if (bt[Mutation::I]) {
                cur_state = Mutation::I;
            } else {
                throw std::runtime_error("Backtrack terminated unexpectedly in deletion state");
            }
        } else { // Mutation::I
            al_seq.push_back(cur_j > 0 ? seq_ptr[cur_j - 1] : 0);
            al_ref.push_back(0);
            if (cur_j > 0) { --cur_j; }

            if (bt[Mutation::I]) {
                cur_state = Mutation::I;
            } else if (bt[Mutation::S]) {
                cur_state = Mutation::S;
            } else if (bt[Mutation::D]) {
                cur_state = Mutation::D;
            } else {
                throw std::runtime_error("Backtrack terminated unexpectedly in insertion state");
            }
        }
    }

    std::reverse(al_seq.begin(), al_seq.end());
    std::reverse(al_ref.begin(), al_ref.end());

    // Pre-allocate vectors for better performance
    std::vector<std::uint8_t> seq_vec;
    std::vector<std::uint8_t> ref_vec;
    seq_vec.reserve(al_seq.size());
    ref_vec.reserve(al_ref.size());
    
    for (std::uint8_t v : al_seq) {
        seq_vec.push_back(v);
    }
    for (std::uint8_t v : al_ref) {
        ref_vec.push_back(v);
    }
    
    // Convert to pybind11 lists
    py::list seq_list = py::cast(seq_vec);
    py::list ref_list = py::cast(ref_vec);

    return py::make_tuple(best_score, seq_list, ref_list);
}

PYBIND11_MODULE(_cas9_align, m) {
    m.doc() = "C++ acceleration for cas9_align";
    m.def("cas9_align", &cas9_align_cpp,
          py::arg("seq"),
          py::arg("ref"),
          py::arg("open_penalty"),
          py::arg("close_penalty"),
          py::arg("sub_score"),
          "Compute Cas9 alignment using optimized C++ implementation");
}

