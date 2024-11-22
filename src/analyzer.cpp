#include "analyzer.hpp"

#include <iostream>
#include <random>
#include <set>
#include <vector>

namespace mc {

Analyzer::Analyzer() {
  std::random_device rd;
  random_generator = std::mt19937(rd());
  dist = std::uniform_int_distribution<unsigned>(settings.min_size,
                                                 settings.max_size);
}

Analyzer::Analyzer(const unsigned min_size, const unsigned max_size)
    : settings{min_size, max_size} {
  std::random_device rd;
  random_generator = std::mt19937(rd());
  dist = std::uniform_int_distribution<unsigned>(min_size, max_size);
}

Instance Analyzer::randomInstance(const unsigned n) {
  Instance instance(n + 1);
  for (unsigned i = 0U; i < n + 1; i++) instance[i] = dist(random_generator);
  return instance;
}

std::vector<Instance> Analyzer::randomInstances(const unsigned n,
                                                const unsigned n_instances) {
  std::vector<Instance> instances;
  instances.reserve(n_instances);
  for (unsigned i = 0U; i < n_instances; i++)
    instances.emplace_back(randomInstance(n));
  return instances;
}

std::map<Permutation, unsigned> getMapPerm2Index(
    const std::vector<Algorithm>& A) {
  std::map<Permutation, unsigned> perm2index;
  for (unsigned i = 0; i < A.size(); i++) perm2index[A[i].getPermutation()] = i;
  return perm2index;
}

std::vector<double> FLOPsOnInstances(std::vector<Algorithm>& A,
                                     const std::vector<Instance>& S) {
  const unsigned M = A.size();
  const unsigned N = S.size();
  std::vector<double> cost_matrix(M * N);

  for (unsigned j = 0U; j < N; j++) {
    for (unsigned i = 0U; i < M; i++) {
      cost_matrix[j * M + i] = A[i].computeFlops(S[j]);
    }
  }
  return cost_matrix;
}

std::vector<double> getMinA(const unsigned M, const unsigned N,
                            const std::vector<double>& cost_matrix) {
  std::vector<double> min_A(N, std::numeric_limits<double>::max());
  for (unsigned j = 0; j < N; j++) {
    for (unsigned i = 0; i < M; i++) {
      if (cost_matrix[j * M + i] < min_A[j]) min_A[j] = cost_matrix[j * M + i];
    }
  }
  return min_A;
}

std::vector<double> getMinZ(const unsigned M, const unsigned N,
                            const std::vector<double>& cost_matrix,
                            const std::set<unsigned>& Z) {
  std::vector<double> min_Z(N, std::numeric_limits<double>::max());
  for (unsigned j = 0U; j < N; j++) {
    for (const auto& id : Z) {
      if (cost_matrix[j * M + id] < min_Z[j])
        min_Z[j] = cost_matrix[j * M + id];
    }
  }
  return min_Z;
}

double penalty(const double min_A, const double min_Z) {
  return (min_Z / min_A) - 1.0;
}

std::vector<double> getPenaltyZ(const unsigned N,
                                const std::vector<double>& min_A,
                                const std::vector<double>& min_Z) {
  std::vector<double> penalty_Z(N, 0.0);
  for (unsigned j = 0; j < N; j++) {
    penalty_Z[j] = penalty(min_A[j], min_Z[j]);
  }
  return penalty_Z;
}

void printMetrics(const std::vector<double>& penalty_Z,
                  const std::string name_exp) {
  double N = static_cast<double>(penalty_Z.size());
  double max_penalty = -1.0;
  double nnz = 0.0;
  double avg_penalty = 0.0;

  for (unsigned i = 0; i < penalty_Z.size(); i++) {
    if (penalty_Z[i] > 0.0) {
      nnz += 1.0;
      avg_penalty += penalty_Z[i];
      if (penalty_Z[i] > max_penalty) max_penalty = penalty_Z[i];
    }
  }

  std::cout << "================================================\n"
            << name_exp << '\n'
            << "================================================\n"
            << "max_penalty: " << max_penalty << "\n"
            << "freq_penalty: " << nnz / N << "\n"
            << "avg_penalty: " << avg_penalty / N << "\n"
            << "avg_penalty (only non-zero): " << avg_penalty / nnz << "\n"
            << "================================================\n\n";
}

Permutation getPermFromApprx(
    const Instance& instance,
    std::function<Permutation(const Instance&)> apprx) {
  return apprx(instance);
}

std::vector<Permutation> getPermFromApprx(
    const std::vector<Instance>& S,
    std::function<Permutation(const Instance&)> apprx) {
  std::vector<Permutation> perms(S.size());
  for (unsigned i = 0; i < S.size(); i++) perms[i] = apprx(S[i]);

  return perms;
}

std::vector<double> getCostFromApprx(
    std::vector<Algorithm>& A, const std::vector<Instance>& S,
    const std::vector<double>& cost_matrix,
    const std::map<Permutation, unsigned>& perm2index,
    std::function<Permutation(const Instance&)> apprx) {
  const unsigned M = A.size();
  const unsigned N = S.size();
  std::vector<double> cost_apprx(N);
  unsigned idx;
  for (unsigned j = 0; j < N; j++) {
    idx = perm2index.at(apprx(S[j]));
    cost_apprx[j] = cost_matrix[j * M + idx];
  }
  return cost_apprx;
}

}  // namespace mc
