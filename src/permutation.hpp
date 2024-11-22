#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <vector>

#include "definitions.hpp"

namespace mc {

struct PermutationTransformer {
  struct InfoEntry {
    int left{-1}, right{-1};
    bool computed{false};
  };

  static std::vector<InfoEntry> table_info;

  /**
   * @brief Returns the canonical form of the input permutation.
   *
   * @param perm          Permutation of which to obtain the canonical form.
   * @return Permutation  Canonical form of the input permutation.
   */
  static Permutation canonicalize(const Permutation& perm);

 private:
  static void clearTable();

  static void buildRepresentation(const Permutation& perm);

  static void addDependencies(const unsigned& p);

  static Permutation buildPermutation(const Permutation& perm);

  static void buildRecursive(const unsigned p, Permutation& canonical_perm);
};

}  // namespace mc

#endif