
#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// auxiliary function: construct the mapping
std::unordered_map<int, std::vector<int>> build_label_map(const IntegerVector& labels) {
  std::unordered_map<int, std::vector<int>> label_map;
  for (int i = 0; i < labels.size(); ++i) {
    label_map[labels[i]].push_back(i);
  }
  return label_map;
}

// [[Rcpp::export]]
IntegerMatrix generate_sg_pairs_int(const IntegerVector& c_labels,
                                    const IntegerVector& r_labels) {
  auto label_map = build_label_map(c_labels);
  std::vector<int> i_vec, j_vec;

  for (const auto& kv : label_map) {
    const std::vector<int>& indices = kv.second;
    for (size_t i = 0; i < indices.size(); ++i) {
      int idx_i = indices[i];
      for (size_t j = i + 1; j < indices.size(); ++j) {
        int idx_j = indices[j];
        if (r_labels[idx_i] != r_labels[idx_j]) {
          i_vec.push_back(idx_i);
          j_vec.push_back(idx_j);
        }
      }
    }
  }

  int n = i_vec.size();
  IntegerMatrix result(n, 2);
  for (int k = 0; k < n; ++k) {
    result(k, 0) = i_vec[k];
    result(k, 1) = j_vec[k];
  }
  return result;
}

// [[Rcpp::export]]
IntegerMatrix generate_gs_pairs_int(const IntegerVector& c_labels,
                                    const IntegerVector& r_labels) {
  auto label_map = build_label_map(r_labels);
  std::vector<int> i_vec, j_vec;

  for (const auto& kv : label_map) {
    const std::vector<int>& indices = kv.second;
    for (size_t i = 0; i < indices.size(); ++i) {
      int idx_i = indices[i];
      for (size_t j = i + 1; j < indices.size(); ++j) {
        int idx_j = indices[j];
        if (c_labels[idx_i] != c_labels[idx_j]) {
          i_vec.push_back(idx_i);
          j_vec.push_back(idx_j);
        }
      }
    }
  }

  int n = i_vec.size();
  IntegerMatrix result(n, 2);
  for (int k = 0; k < n; ++k) {
    result(k, 0) = i_vec[k];
    result(k, 1) = j_vec[k];
  }
  return result;
}
