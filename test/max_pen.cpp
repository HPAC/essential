#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include "../src/algorithm.hpp"
#include "../src/analyzer.hpp"
#include "../src/apprx_algorithms.hpp"
#include "../src/definitions.hpp"
#include "../src/generator.hpp"

std::ostream& operator<<(std::ostream& os, const mc::Permutation& perm) {
  os << "[";
  for (const auto& p : perm) {
    os << p << ' ';
  }
  os << "]";
  return os;
}

int main(int argc, char** argv) {
  unsigned n, n_samples;
  if (argc < 3) {
    std::cerr << "Usage: ./max_pen n n_samples\n";
    exit(-1);
  } else {
    n = std::stoi(argv[1]);
    n_samples = std::stoi(argv[2]);
  }

  auto A = mc::generateAlgorithms(n);
  mc::Analyzer analyzer(1U, 1000U);
  std::vector<mc::Instance> S = analyzer.randomInstances(n, n_samples);

  const unsigned M = A.size();
  const unsigned N = S.size();

  auto cost_matrix = mc::FLOPsOnInstances(A, S);
  auto perm2index = mc::getMapPerm2Index(A);
  auto min_A = mc::getMinA(M, N, cost_matrix);

  using AlgMCP = std::function<mc::Permutation(const mc::Instance&)>;

  // Chin's
  AlgMCP chin = mc::chin;
  auto cost_chin = mc::getCostFromApprx(A, S, cost_matrix, perm2index, chin);
  auto penalty_chin = mc::getPenaltyZ(N, min_A, cost_chin);
  mc::printMetrics(penalty_chin, "Chin's:");

  auto it_max = std::max_element(penalty_chin.begin(), penalty_chin.end());
  auto idx_max = std::distance(penalty_chin.begin(), it_max);
  auto instance_max = S[idx_max];
  std::cout << "Chin's max penalty on: " << instance_max << '\n';

  // // Reduce and minimize
  AlgMCP rnm = mc::reduceMin;
  auto cost_rnm = mc::getCostFromApprx(A, S, cost_matrix, perm2index, rnm);
  auto penalty_rnm = mc::getPenaltyZ(N, min_A, cost_rnm);
  mc::printMetrics(penalty_rnm, "Algorithm 3:");

  it_max = std::max_element(penalty_rnm.begin(), penalty_rnm.end());
  idx_max = std::distance(penalty_rnm.begin(), it_max);
  instance_max = S[idx_max];
  std::cout << "Algorithm 3's max penalty on: " << instance_max << '\n';
}