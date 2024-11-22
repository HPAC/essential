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
  unsigned n;
  mc::Instance instance;

  if (argc < 4) {
    std::cerr << "Usage: ./single_instance n k0 k1 ...\n";
    exit(-1);
  } else {
    n = std::stoi(argv[1]);
    for (int i = 2; i < argc; i++) {
      instance.push_back(std::stoi(argv[i]));
    }
    if (n != instance.size() - 1) {
      std::cerr << "Mismatch between n and instance size\n";
      exit(-1);
    }
  }

  auto A = mc::generateAlgorithms(n);

  std::vector<mc::Instance> S = {instance};

  auto essential_perms = mc::getEssentialPerms(n);
  std::set<unsigned> E;  // set of indices of essential parenthesisations.
  for (const auto& perm : essential_perms) E.insert(mc::getID(A, perm));

  const unsigned M = A.size();
  const unsigned N = S.size();

  auto cost_matrix = mc::FLOPsOnInstances(A, S);
  auto perm2index = mc::getMapPerm2Index(A);
  auto min_A = mc::getMinA(M, N, cost_matrix);

  // Essentials - Algorithm 1.
  auto min_E = mc::getMinZ(M, N, cost_matrix, E);
  auto penalty_E = mc::getPenaltyZ(N, min_A, min_E);
  std::cout << "Penalty Algorithm 1: " << penalty_E[0] << "\n";

  using AlgMCP = std::function<mc::Permutation(const mc::Instance&)>;
  // Chandra's.
  AlgMCP chandra = mc::chandra;
  auto cost_chandra =
      mc::getCostFromApprx(A, S, cost_matrix, perm2index, chandra);
  auto penalty_chandra = mc::getPenaltyZ(N, min_A, cost_chandra);
  std::cout << "Penalty Chandra: " << penalty_chandra[0] << "\n";

  // Chin's.
  AlgMCP chin = mc::chin;
  auto cost_chin = mc::getCostFromApprx(A, S, cost_matrix, perm2index, chin);
  auto penalty_chin = mc::getPenaltyZ(N, min_A, cost_chin);
  std::cout << "Penalty Chin: " << penalty_chin[0] << "\n";

  // // Reduce and minimize - Algorithm 3.
  AlgMCP rnm = mc::reduceMin;
  auto cost_rnm = mc::getCostFromApprx(A, S, cost_matrix, perm2index, rnm);
  auto penalty_rnm = mc::getPenaltyZ(N, min_A, cost_rnm);
  std::cout << "Penalty Algorithm 3: " << penalty_rnm[0] << "\n";
}