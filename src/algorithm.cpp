#include "algorithm.hpp"

#include <cstdint>
#include <vector>

#include "definitions.hpp"

namespace mc {

Algorithm::Algorithm(const Permutation& permutation)
    : permutation{permutation} {
  buildTree();
}

double Algorithm::computeFlops(const Instance& instance) {
  assignSizes(instance);

  double flops = 0.0;
  for (unsigned i = permutation.size() + 1; i < tree.size(); i++) {
    propagateSizes(i);
    flops += costMult(tree[tree[i].left]._rows, tree[tree[i].left]._cols,
                      tree[tree[i].right]._cols);
  }

  return flops;
}

void Algorithm::buildTree() {
  tree.reserve(2 * permutation.size() + 1);
  createInputNodes();

  int8_t left, right;
  for (const auto& p : permutation) {
    left = getRoot(p - 1);
    right = getRoot(p);

    tree.emplace_back(left, right);

    tree[left].parent = tree.size() - 1;
    tree[right].parent = tree.size() - 1;
  }
}

void Algorithm::createInputNodes() {
  for (unsigned i = 0U; i < permutation.size() + 1; i++) tree.emplace_back();
}

int8_t Algorithm::getRoot(const int8_t id) const {
  return (tree[id].parent == -1) ? id : getRoot(tree[id].parent);
}

void Algorithm::assignSizes(const Instance& instance) {
  for (unsigned i = 0U; i <= permutation.size(); i++) {
    tree[i]._rows = instance[i];
    tree[i]._cols = instance[i + 1];
  }
}

void Algorithm::propagateSizes(const int8_t id) {
  tree[id]._rows = tree[tree[id].left]._rows;
  tree[id]._cols = tree[tree[id].right]._cols;
}

double Algorithm::costMult(const unsigned m, const unsigned k,
                           const unsigned n) const {
  return static_cast<double>(m) * static_cast<double>(k) *
         static_cast<double>(n);
}

}  // namespace mc
