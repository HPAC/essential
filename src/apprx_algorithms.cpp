#include "apprx_algorithms.hpp"

#include <algorithm>
#include <deque>
#include <iostream>

#include "definitions.hpp"
#include "generator.hpp"
#include "permutation.hpp"

namespace mc {

Permutation chandra(const Instance& k) {
  auto it_min = std::min_element(k.begin(), k.end());
  unsigned idx = static_cast<unsigned>(std::distance(k.begin(), it_min));
  return getEssentialPerm(k.size() - 1U, idx);
}

unsigned minEssential(const Instance& k) {
  const unsigned n = k.size() - 1U;
  std::vector<double> z(k.size());
  z[0] = static_cast<double>(k[0]) * static_cast<double>(k[n]);

  for (unsigned i = 1U; i <= n; i++)
    z[i] = static_cast<double>(k[i - 1U]) * static_cast<double>(k[i]);

  double x{0.0};
  for (unsigned i = 0U; i <= n; i++) x += z[i];

  std::vector<double> t(k.size());
  for (unsigned i = 0U; i < n; i++)
    t[i] = static_cast<double>(k[i]) * (x - z[i] - z[i + 1]);
  t[n] = static_cast<double>(k[n]) * (x - z[n] - z[0]);

  unsigned h = 0U;
  double t_h = t[h];
  for (unsigned i = 1; i <= n; i++) {
    if (t[i] < t_h) {
      t_h = t[i];
      h = i;
    }
  }
  return h;
}

Permutation chin(const Instance& k) {
  const int n = k.size() - 1U;

  // Initialise.
  auto it_min = std::min_element(k.begin(), k.end());
  int m = static_cast<int>(std::distance(k.begin(), it_min));

  std::vector<double> r(k.size());
  for (int i = 0; i < static_cast<int>(k.size()); i++)
    r[i] = 1.0 / static_cast<double>(k[i]);

  Permutation v(n - 1U);  // Chin's order of computation.

  std::deque<int> Q;
  // Scan forward.
  Q.push_back(0U);
  int a = 1U;

  for (int i = 1U; i < n; i++) {
    Q.push_back(i);
    while (Q.size() >= 2U and
           (r[Q.back()] + r[m] < r[Q[Q.size() - 2U]] + r[i + 1U])) {
      v[Q.back() - 1] = a;
      a++;
      Q.pop_back();
    }
  }
  Q.push_back(n);
  // Nibble at both ends.
  int b = n - 1U;

  while (Q.size() >= 3U) {
    if (r[Q.back()] + r[m] < r[Q[Q.size() - 2U]] + r[Q.front()]) {
      Q.pop_back();
      v[Q.back() - 1] = b;
      b--;
    } else if (r[Q.front()] + r[m] < r[Q[1]] + r[Q.back()]) {
      Q.pop_front();
      v[Q.front() - 1] = b;
      b--;
    } else {
      break;
    }
  }

  // Associate out from smallest dimension.
  for (int i = m - 1; i > Q.front(); i--) {
    if (v[i - 1] == 0U) {
      v[i - 1] = a;
      a++;
    }
  }
  for (int i = m + 1; i < Q.back(); i++) {
    if (v[i - 1] == 0U) {
      v[i - 1] = a;
      a++;
    }
  }

  // Final multiplication, if any.
  if (m != 0U and m != n and v[m - 1] == 0U) {
    v[m - 1] = a;
  }

  // If Chin's format is to be used. Simply return.
  // Otherwise, convert from Chin's format to ours.
  return chin2Canonical(v);
}

Permutation reduceMin(const Instance& k) {
  const int n = k.size() - 1U;

  // Initialise.
  auto it_min = std::min_element(k.begin(), k.end());
  int m = static_cast<int>(std::distance(k.begin(), it_min));

  std::vector<double> r(k.size());
  for (int i = 0U; i < static_cast<int>(k.size()); i++)
    r[i] = 1.0 / static_cast<double>(k[i]);

  Permutation v(n - 1U);  // Chin's order of computation.

  std::deque<int> Q;
  // Scan forward.
  Q.push_back(0U);
  int a = 1U;

  for (int i = 1U; i < n; i++) {
    Q.push_back(i);
    while (Q.size() >= 2U and
           (r[Q.back()] + r[m] < r[Q[Q.size() - 2U]] + r[i + 1U])) {
      v[Q.back() - 1] = a;
      a++;
      Q.pop_back();
    }
  }
  Q.push_back(n);

  // Nibble at both ends.
  int b = n - 1U;

  while (Q.size() >= 3U) {
    if (r[Q.back()] + r[m] < r[Q[Q.size() - 2U]] + r[Q.front()]) {
      Q.pop_back();
      v[Q.back() - 1] = b;
      b--;
    } else if (r[Q.front()] + r[m] < r[Q[1]] + r[Q.back()]) {
      Q.pop_front();
      v[Q.front() - 1] = b;
      b--;
    } else {
      break;
    }
  }

  // minimise over the essential parenth of the remaining chain.
  if (Q.size() >= 3) {
    Instance q(Q.size());  // reduced instance;
    for (int i = 0; i < static_cast<int>(Q.size()); i++) q[i] = k[Q[i]];
    int h_q = minEssential(q);
    int h = Q[h_q];

    // Apply the selected essential parenthesisation
    for (int i = h - 1; i > Q.front(); i--) {
      if (v[i - 1] == 0U) {
        v[i - 1] = a;
        a++;
      }
    }
    for (int i = h + 1; i < Q.back(); i++) {
      if (v[i - 1] == 0U) {
        v[i - 1] = a;
        a++;
      }
    }

    // Final multiplication, if any.
    if (h != 0U and h != n and v[h - 1] == 0U) {
      v[h - 1] = a;
    }
  }

  // If Chin's format is to be used. Simply return.
  // Otherwise, convert from Chin's format to ours.
  return chin2Canonical(v);
}

Permutation chin2Canonical(const Permutation& chins_perm) {
  Permutation canonical_perm(chins_perm.size());
  for (unsigned i = 0; i < chins_perm.size(); i++) {
    canonical_perm[chins_perm[i] - 1] = i + 1;
  }

  return PermutationTransformer::canonicalize(canonical_perm);
}

}  // namespace mc
