// MANDEL NOTATION HELPERS
//
// Symmetric second- and fourth-order tensors are represented as 6-vectors and
// 6x6 matrices in MANDEL notation: off-diagonal entries carry sqrt(2) for
// vectors, and the outer product of the sqrt(2) factors for matrices.
//
// The point of the convention is that the plain linear-algebra operations
// below coincide with the tensor operations they represent:
//
//     dotProduct6(mandel(A), mandel(B))  ==  A : B
//     vectorNorm6(mandel(A))             ==  A.L2norm()
//     matVecMult6(mandel(C), mandel(A))  ==  mandel(C : A)
//     invert(mandel(C))                  ==  mandel(C^-1)
//
// This is also the convention the reference UMAT uses internally; it converts
// to it on entry (UMAT lines 1359-1365) and back to Abaqus convention on exit
// (lines 1950-1967).
//
// DO NOT add a "plain copy" variant that drops the sqrt(2) factors. A previous
// version of this file did exactly that, which made every shear compliance
// component wrong by a factor of 4 and the yield gradient wrong by a factor of
// 2 in shear. mandelRoundTripError() at the bottom is the regression test.
//
// NOTE: _F_matrix and _f_lin_vector in OrthotropicPlasticityStressUpdate are
// NOT in this convention. They act on plain stress components (F[5][5] =
// 1/tau_xy^2 multiplying sigma_12), and the 6-vector algebra around them in
// computeYieldFunction / computeYieldGradient / computeYieldHessian is written
// to match. Do not "tidy" those loops into the helpers below: the yield
// surface would silently change.

#pragma once

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include <vector>
#include <cmath>
#include <limits>

// index order: xx yy zz yz xz xy
static const int MANDEL_IJ[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};
static const Real MANDEL_C[6] = {1.0, 1.0, 1.0,
                                 1.41421356237309504880,
                                 1.41421356237309504880,
                                 1.41421356237309504880};

// ============================================================================
// CONVERSIONS
// ============================================================================

inline void tensorToMandel(const RankTwoTensor & t, std::vector<Real> & v)
{
  v.resize(6);
  v[0] = t(0, 0);
  v[1] = t(1, 1);
  v[2] = t(2, 2);
  v[3] = MANDEL_C[3] * t(1, 2);
  v[4] = MANDEL_C[4] * t(0, 2);
  v[5] = MANDEL_C[5] * t(0, 1);
}

// Reads the first 6 entries only, so a 7-vector (stress + dkappa) can be
// passed directly.
inline void mandelToTensor(const std::vector<Real> & v, RankTwoTensor & t)
{
  t.zero();
  t(0, 0) = v[0];
  t(1, 1) = v[1];
  t(2, 2) = v[2];
  t(1, 2) = t(2, 1) = v[3] / MANDEL_C[3];
  t(0, 2) = t(2, 0) = v[4] / MANDEL_C[4];
  t(0, 1) = t(1, 0) = v[5] / MANDEL_C[5];
}

// Requires minor symmetries only; major symmetry is not assumed, so this is
// safe for DRRDS and DNPDS, which lack it.
inline void rankFourToMandel(const RankFourTensor & C, std::vector<std::vector<Real>> & M)
{
  M.assign(6, std::vector<Real>(6, 0.0));
  for (int p = 0; p < 6; p++)
    for (int q = 0; q < 6; q++)
      M[p][q] = MANDEL_C[p] * MANDEL_C[q] *
                C(MANDEL_IJ[p][0], MANDEL_IJ[p][1], MANDEL_IJ[q][0], MANDEL_IJ[q][1]);
}

inline void mandelToRankFour(const std::vector<std::vector<Real>> & M, RankFourTensor & C)
{
  C.zero();
  for (int p = 0; p < 6; p++)
  {
    const int i = MANDEL_IJ[p][0], j = MANDEL_IJ[p][1];
    for (int q = 0; q < 6; q++)
    {
      const int k = MANDEL_IJ[q][0], l = MANDEL_IJ[q][1];
      const Real val = M[p][q] / (MANDEL_C[p] * MANDEL_C[q]);
      C(i, j, k, l) = val;
      if (i != j) C(j, i, k, l) = val;
      if (k != l) C(i, j, l, k) = val;
      if (i != j && k != l) C(j, i, l, k) = val;
    }
  }
}

// ============================================================================
// LINEAR ALGEBRA (convention-agnostic)
// ============================================================================

inline void matVecMult6(const std::vector<std::vector<Real>> & A,
                        const std::vector<Real> & v,
                        std::vector<Real> & result)
{
  result.assign(6, 0.0);
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      result[i] += A[i][j] * v[j];
}

inline Real dotProduct6(const std::vector<Real> & a, const std::vector<Real> & b)
{
  Real result = 0.0;
  for (int i = 0; i < 6; i++)
    result += a[i] * b[i];
  return result;
}

inline Real vectorNorm6(const std::vector<Real> & v)
{
  return std::sqrt(dotProduct6(v, v));
}

inline void matVecMult7(const std::vector<std::vector<Real>> & A,
                        const std::vector<Real> & v,
                        std::vector<Real> & result)
{
  result.assign(7, 0.0);
  for (int i = 0; i < 7; i++)
    for (int j = 0; j < 7; j++)
      result[i] += A[i][j] * v[j];
}

inline Real dotProduct7(const std::vector<Real> & a, const std::vector<Real> & b)
{
  Real result = 0.0;
  for (int i = 0; i < 7; i++)
    result += a[i] * b[i];
  return result;
}

// result_ij = C_ijkl * A_kl
inline RankTwoTensor contractRankFourTwo(const RankFourTensor & C, const RankTwoTensor & A)
{
  RankTwoTensor result;
  result.zero();
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          result(i, j) += C(i, j, k, l) * A(k, l);
  return result;
}

// Gauss-Jordan with partial pivoting. Returns false on a singular matrix
// rather than aborting, so callers can fall back instead of killing the run.
inline bool invertMatrix6x6(const std::vector<std::vector<Real>> & A,
                            std::vector<std::vector<Real>> & Ainv)
{
  const int n = 6;
  std::vector<std::vector<Real>> m = A;
  std::vector<std::vector<Real>> inv(n, std::vector<Real>(n, 0.0));
  for (int i = 0; i < n; i++)
    inv[i][i] = 1.0;

  for (int i = 0; i < n; i++)
  {
    int imax = i;
    Real amax = std::abs(m[i][i]);
    for (int k = i + 1; k < n; k++)
      if (std::abs(m[k][i]) > amax) { amax = std::abs(m[k][i]); imax = k; }
    if (amax < 1e-14) return false;
    std::swap(m[i], m[imax]);
    std::swap(inv[i], inv[imax]);

    const Real piv = m[i][i];
    for (int j = 0; j < n; j++) { m[i][j] /= piv; inv[i][j] /= piv; }

    for (int k = 0; k < n; k++)
      if (k != i)
      {
        const Real f = m[k][i];
        for (int j = 0; j < n; j++) { m[k][j] -= f * m[i][j]; inv[k][j] -= f * inv[i][j]; }
      }
  }
  Ainv = inv;
  return true;
}

// ============================================================================
// REGRESSION TEST
// ============================================================================
//
// rankFourToMandel -> invert -> mandelToRankFour must reproduce C.invSymm().
// Returns the max component difference normalised by the largest compliance
// component. A correct implementation returns < 1e-12; the plain-copy
// convention returns O(1). Call it once on the elasticity tensor after any
// change to the conversions above.
inline Real mandelRoundTripError(const RankFourTensor & C)
{
  std::vector<std::vector<Real>> M, Minv;
  rankFourToMandel(C, M);
  if (!invertMatrix6x6(M, Minv))
    return std::numeric_limits<Real>::max();

  RankFourTensor S_got;
  mandelToRankFour(Minv, S_got);
  const RankFourTensor S_ref = C.invSymm();

  Real diff = 0.0, scale = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
        {
          diff = std::max(diff, std::abs(S_got(i, j, k, l) - S_ref(i, j, k, l)));
          scale = std::max(scale, std::abs(S_ref(i, j, k, l)));
        }
  return (scale > 0.0) ? diff / scale : diff;
}
