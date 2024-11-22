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
    std::cerr << "Usage: ./experiment n n_samples\n";
    exit(-1);
  } else {
    n = std::stoi(argv[1]);
    n_samples = std::stoi(argv[2]);
  }

  auto A = mc::generateAlgorithms(n);
  mc::Analyzer analyzer(1U, 1000U);

  std::vector<mc::Instance> S = analyzer.randomInstances(n, n_samples);

  auto essential_perms = mc::getEssentialPerms(n);
  std::set<unsigned> E;  // set of indices of essential parenthesisations.
  for (const auto& perm : essential_perms) E.insert(mc::getID(A, perm));

  const unsigned M = A.size();
  const unsigned N = S.size();
  std::cout << "M: " << M << "\n";
  std::cout << "N: " << N << "\n";

  auto cost_matrix = mc::FLOPsOnInstances(A, S);
  auto perm2index = mc::getMapPerm2Index(A);
  auto min_A = mc::getMinA(M, N, cost_matrix);

  // Essentials
  auto min_E = mc::getMinZ(M, N, cost_matrix, E);
  auto penalty_E = mc::getPenaltyZ(N, min_A, min_E);
  mc::printMetrics(penalty_E, "Essentials:");

  using AlgMCP = std::function<mc::Permutation(const mc::Instance&)>;

  // Chandra's
  AlgMCP chandra = mc::chandra;
  auto cost_chandra =
      mc::getCostFromApprx(A, S, cost_matrix, perm2index, chandra);
  auto penalty_chandra = mc::getPenaltyZ(N, min_A, cost_chandra);
  mc::printMetrics(penalty_chandra, "Chandra's:");

  // Chin's
  AlgMCP chin = mc::chin;
  auto cost_chin = mc::getCostFromApprx(A, S, cost_matrix, perm2index, chin);
  auto penalty_chin = mc::getPenaltyZ(N, min_A, cost_chin);
  mc::printMetrics(penalty_chin, "Chin's:");

  // // Reduce and minimize
  AlgMCP rnm = mc::reduceMin;
  auto cost_rnm = mc::getCostFromApprx(A, S, cost_matrix, perm2index, rnm);
  auto penalty_rnm = mc::getPenaltyZ(N, min_A, cost_rnm);
  mc::printMetrics(penalty_rnm, "Algorithm 3:");
}