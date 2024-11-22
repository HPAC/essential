#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <cstdint>
#include <vector>

#include "definitions.hpp"

namespace mc {

class Node {
 public:
  unsigned _rows{}, _cols{};
  int8_t left{-1}, right{-1}, parent{-1};

  // Default Constructor. Constructor for input nodes.
  Node() = default;

  Node(const int8_t left, const int8_t right) : left{left}, right{right} {}
};

class Algorithm {  // An algorithm is a parenthesisation.
 private:
  Permutation permutation;  // Order of computation for the algorithm.
  std::vector<Node> tree;   // Proper in-memory representation of the algorithm.

 public:
  Algorithm() = delete;

  /**
   * @brief Parametrised constructor taking in a permutation.
   *
   * @param permutation vector<unsigned>.
   */
  Algorithm(const Permutation& permutation);

  ~Algorithm() = default;

  // Getter for the permutation.
  inline Permutation getPermutation() const noexcept { return permutation; }

  /**
   * @brief Returns the number of FLOPs for the algorithm on the given instance.
   *
   * @param instance  vector<unsigned>.
   * @return double   number of FLOPs.
   */
  double computeFlops(const Instance& instance);

 private:
  /**
   * @brief Does the actual work when constructing the Algorithm.
   */
  void buildTree();

  /**
   * @brief Creates N (chain's length) nodes - one per matrix in the input
   * chain.
   */
  void createInputNodes();

  /**
   * @brief Returns the topmost node's ID of the node whose ID is passed.
   *
   * @param id      int8_t - index of the node in the tree vector.
   * @return int8_t topmost node's ID.
   */
  int8_t getRoot(const int8_t id) const;

  /**
   * @brief Assigns sizes in the instance to the input nodes.
   *
   * @param instance vector<unsigned>.
   */
  void assignSizes(const Instance& instance);

  /**
   * @brief Propagates the sizes for the node with the passed ID.
   *
   * @param id  int8_t - ID of the node to which sizes are propagated.
   */
  void propagateSizes(const int8_t id);

  /**
   * @brief Returns the cost of a multiplication of two matrices with sizes
   * (m,k,n). For the purpose of this work: costMult(m,k,n) = m * k * n.
   *
   * @param m       unsigned - number of rows of the left matrix.
   * @param k       unsigned - number of columns of the left matrix.
   * @param n       unsigned - number of columns of the right matrix.
   * @return double cost of the multiplication.
   */
  double costMult(const unsigned m, const unsigned k, const unsigned n) const;
};

}  // namespace mc

#endif