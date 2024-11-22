#ifndef GENERATOR_H
#define GENERATOR_H

#include "algorithm.hpp"
#include "definitions.hpp"

namespace mc {

/**
 * @brief Generates all parenthesisations for a given length of the chain.
 *
 * @param n length of the chain.
 * @return std::vector<Algorithm>
 */
std::vector<Algorithm> generateAlgorithms(const unsigned n);

/**
 * @brief Generates the permutation of the essential parenthesisations.
 *
 * @param n length of the chain.
 * @return std::vector<Permutation>
 */
std::vector<Permutation> getEssentialPerms(const unsigned n);

/**
 * @brief Generates the permutation of the essential parenthesisation that fans
 * out from the dimension with id h.
 *
 * @param n   length of the chain
 * @param h   idx of the dimension to fan out from.
 * @return Permutation
 */
Permutation getEssentialPerm(const unsigned n, const unsigned h);

/**
 * @brief Returns the index of the algorithm with the passed permutation.
 *
 * @param algs
 * @param perm
 * @return unsigned
 */
unsigned getID(const std::vector<Algorithm>& algs, const Permutation& perm);

/**
 * @brief Returns the indices of the algorithms with the passed permutations.
 *
 * @param algs
 * @param perms
 * @return std::vector<unsigned>
 */
std::vector<unsigned> getIDs(const std::vector<Algorithm>& algs,
                             const std::vector<Permutation>& perms);

/**
 * @brief Computes the factorial of n.
 *
 * @param n
 * @return unsigned
 */
unsigned factorial(const unsigned n);

/**
 * @brief Computes the nth Calatan number.
 *
 * @param n
 * @return unsigned
 */
unsigned catalanNumber(const unsigned n);

/**
 * @brief Checks whether a permutation is in canonical form.
 *
 * @param perm
 * @return true
 * @return false
 */
bool isCanonical(const Permutation& perm);

}  // namespace mc

#endif