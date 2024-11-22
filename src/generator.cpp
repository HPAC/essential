#include "generator.hpp"

#include <algorithm>
#include <numeric>
#include <vector>

#include "algorithm.hpp"
#include "definitions.hpp"
#include "permutation.hpp"

namespace mc {

std::vector<Algorithm> generateAlgorithms(const unsigned n) {
  unsigned n_algs = catalanNumber(n - 1U);
  std::vector<Algorithm> algorithms;
  algorithms.reserve(n_algs);

  Permutation permutation(n - 1U);
  std::iota(permutation.begin(), permutation.end(), 1U);

  do {
    if (isCanonical(permutation)) {
      algorithms.emplace_back(permutation);
    }
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  return algorithms;
}

std::vector<Permutation> getEssentialPerms(const unsigned n) {
  std::vector<Permutation> essentials;
  unsigned n_dims = n + 1U;
  essentials.reserve(n_dims);
  for (unsigned i = 0U; i < n_dims; i++)
    essentials.push_back(getEssentialPerm(n, i));

  return essentials;
}

Permutation getEssentialPerm(const unsigned n, const unsigned h) {
  Permutation perm;
  perm.reserve(n - 1);
  for (unsigned i = h - 1; i > 0 and i < n; i--) perm.push_back(i);
  for (unsigned i = h + 1; i < n; i++) perm.push_back(i);

  if (perm.size() < n - 1) perm.push_back(h);
  return perm;
}

unsigned getID(const std::vector<Algorithm>& algs, const Permutation& perm) {
  auto iter = std::find_if(
      algs.begin(), algs.end(),
      [perm](const Algorithm& alg) { return alg.getPermutation() == perm; });
  return std::distance(algs.begin(), iter);
}

std::vector<unsigned> getIDs(const std::vector<Algorithm>& algs,
                             const std::vector<Permutation>& perms) {
  std::vector<unsigned> IDs;
  IDs.reserve(perms.size());
  for (const auto& perm : perms) {
    IDs.push_back(getID(algs, perm));
  }
  return IDs;
}

unsigned factorial(const unsigned n) {
  unsigned fact = 1U;
  for (unsigned i = 2U; i <= n; i++) {
    fact *= i;
  }
  return fact;
}

unsigned catalanNumber(const unsigned n) {
  return factorial(2 * n) / (factorial(n + 1) * factorial(n));
}

bool isCanonical(const Permutation& perm) {
  return PermutationTransformer::canonicalize(perm) == perm;
}

}  // namespace mc
