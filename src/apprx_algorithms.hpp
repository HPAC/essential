#ifndef APPRX_ALGORITHMS_H
#define APPRX_ALGORITHMS_H

#include <vector>

#include "definitions.hpp"

namespace mc {

/**
 * @brief Executes Chandra's approximation algorithm on the passed instance.
 *
 * Reference: A.K. Chandra. Computing matrix chain products in near-optimal
 * time. Tech. Rep. RC-5625. IBM Thomas J. Watson Research Center, P. O. Box
 * 218, Yorktown Heights, New York, USA (1975).
 *
 * @param k             Instance.
 * @return Permutation  Yielded order of computation.
 */
Permutation chandra(const Instance& k);

/**
 * @brief Finds the index of the dimension of the essential parenthesisation
 * with minimal cost for the given instance.
 *
 * @param k           Instance.
 * @return unsigned   Index of the dimension of the essential parenthesisation
 * with minimal cost.
 */
unsigned minEssential(const Instance& k);

/**
 * @brief Executes Chin's algorithm on the passed instance.
 *
 * Reference: F.Y. Chin. An O(n) algorithm for determining a near-optimal
 * computation order of matrix chain products. Communications of the ACM. 1978.
 *
 * @param k             Instance.
 * @return Permutation  Yielded order of computation (in our convention).
 */
Permutation chin(const Instance& k);

/**
 * @brief Executes Algorithm 3 in the paper.
 *
 * First finds multiplications that must be in the optimal order, given the
 * condition by Lemma 1 in (Chin 1978). Then, finds the essential
 * parenthesisation on the reduced chain with minimal cost.
 *
 * @param k             Instance.
 * @return Permutation  Yielded order of computation (in our convention).
 */
Permutation reduceMin(const Instance& k);

/**
 * @brief Converts Chin's notation of the order of computation to our notation.
 *
 * Used to be able to find Chin's yielded parenthesisation in our convention.
 *
 * @param chins_perm    vector<unsigned>.
 * @return Permutation  vector<unsigned> - order of computation in our
 * convention.
 */
Permutation chin2Canonical(const Permutation& chins_perm);

}  // namespace mc

#endif