#ifndef ANALYZER_H
#define ANALYZER_H

#include <functional>
#include <map>
#include <random>
#include <set>
#include <vector>

#include "algorithm.hpp"
#include "apprx_algorithms.hpp"
#include "definitions.hpp"

namespace mc {

struct ConfigAnalyzer {
  unsigned min_size{1U};
  unsigned max_size{1'000U};

  ConfigAnalyzer() = default;

  ConfigAnalyzer(const unsigned min_size, const unsigned max_size)
      : min_size{min_size}, max_size{max_size} {}
};

class Analyzer {
 public:
  ConfigAnalyzer settings{};
  std::uniform_int_distribution<unsigned> dist{};
  std::mt19937 random_generator{};

  Analyzer();

  /**
   * @brief Parametrised constructor.
   *
   * @param min_size minimum size that can be generated.
   * @param max_size maximum size that can be generated.
   */
  Analyzer(const unsigned min_size, const unsigned max_size);

  /**
   * @brief Returns a random instance for a chain of length n.
   *
   * @param n           length of the chain.
   * @return Instance   produced instance.
   */
  Instance randomInstance(const unsigned n);

  /**
   * @brief Returns a vector containing n_instances.
   *
   * @param n                       length of the chain.
   * @param n_instances             number of instances.
   * @return std::vector<Instance>  generated instances.
   */
  std::vector<Instance> randomInstances(const unsigned n,
                                        const unsigned n_instances);
};

/**
 * @brief Generates a map from permutations to indices for faster retrieval.
 *
 * @param A  vector containing all parenthesisations (= algorithms).
 * @return std::map<Permutation, unsigned>
 */
std::map<Permutation, unsigned> getMapPerm2Index(
    const std::vector<Algorithm>& A);

/**
 * @brief Computes the cost of the parenthesisations in A on the instances in S.
 *
 * Returns a matrix of dimensions (A.size() * S.size()) where each row
 * corresponds to the same parenthesisation and each column corresponds to the
 * same instance. Let v = FLOPsOnInstances(A, S); then v(i,j) holds the cost of
 * the i-th parenthesisation on the j-th instance.
 *
 * @param A     vector containing all the parenthesisations.
 * @param S     vector containing the instances.
 * @return std::vector<double> matrix with the cost of all parenthesisations on
 * every instance.
 */
std::vector<double> FLOPsOnInstances(std::vector<Algorithm>& A,
                                     const std::vector<Instance>& S);

/**
 * @brief Returns a vector holding the minimum cost for each instance.
 *
 * @param M                    number of parenthesisations in A.
 * @param N                    number of instances in S.
 * @param cost_matrix          MxN matrix in a vector.
 * @return std::vector<double> vector of size N holding the minimum cost for
 * each instance.
 */
std::vector<double> getMinA(const unsigned M, const unsigned N,
                            const std::vector<double>& cost_matrix);

/**
 * @brief Returns a vector holding the minimum cost across the parenthesisations
 * in Z for all instances.
 *
 * This function is primarily used to determine the minimum cost amongst the
 * essential parenthesisations.
 *
 * @param M                     number of parenthesisations in A.
 * @param N                     number of instances in S.
 * @param cost_matrix           MxN matrix holding the cost of every
 * parenthesisation in A on every instance in S.
 * @param Z                     set of parenthesisations' IDs.
 * @return std::vector<double>  minimum cost across Z.
 */
std::vector<double> getMinZ(const unsigned M, const unsigned N,
                            const std::vector<double>& cost_matrix,
                            const std::set<unsigned>& Z);

/**
 * @brief Computes the penalty of one instance given the overall cheapest cost
 * (min_A) and the cost of the parenthesisation of interest (min_Z).
 *
 * @param min_A   overall cheapest cost.
 * @param min_Z   cost of the parenthesisation of interest.
 * @return double ~ ratio between both normalised to 0.
 */
double penalty(const double min_A, const double min_Z);

/**
 * @brief Computes the penalty on a per instance basis.
 *
 * @param N       number of instances (size of min_A and min_Z).
 * @param min_A   vector with the overall cheapest for every instance.
 * @param min_Z   vector with the cheapest across the parenths in Z.
 * @return std::vector<double>
 */
std::vector<double> getPenaltyZ(const unsigned N,
                                const std::vector<double>& min_A,
                                const std::vector<double>& min_Z);

/**
 * @brief Prints metrics on a vector of penalties.
 *
 * @param penalty_Z   vector of penalties for every instance.
 * @param name_exp    string - name of the experiment (approximation alg being
 * used).
 */
void printMetrics(const std::vector<double>& penalty_Z,
                  const std::string name_exp);

/**
 * @brief Executes the passed approximation algorithm for the given instance.
 *
 * @param instance     vector<unsigned>.
 * @param apprx        approximation algorithm to use.
 * @return Permutation the order of execution the approximation algorithm
 * yields.
 */
Permutation getPermFromApprx(const Instance& instance,
                             std::function<Permutation(const Instance&)> apprx);

/**
 * @brief Executes the passed approximation algorithm on a number of instances.
 *
 * @param instances     vector<Instance>.
 * @param apprx         approximation algorithm to use.
 * @return std::vector<Permutation> orders of execution the approximation
 * algorithm yields for every instance.
 */
std::vector<Permutation> getPermFromApprx(
    const std::vector<Instance>& instances,
    std::function<Permutation(const Instance&)> apprx);

/**
 * @brief Produces a vector with the costs of the parenthesisations yielded by
 * the passed approximation algorithm for all instances in S.
 *
 * @param A           all parenthesisations.
 * @param S           vector<instance>.
 * @param cost_matrix matrix holding the cost of every parenthesisation on every
 * instance.
 * @param perm2index  map from permutation to index. Used for fast retrieval of
 * parenthesisation's cost in cost_matrix.
 * @param apprx       approximation algorithm to use.
 * @return std::vector<double>
 */
std::vector<double> getCostFromApprx(
    std::vector<Algorithm>& A, const std::vector<Instance>& S,
    const std::vector<double>& cost_matrix,
    const std::map<Permutation, unsigned>& perm2index,
    std::function<Permutation(const Instance&)> apprx);

}  // namespace mc

#endif