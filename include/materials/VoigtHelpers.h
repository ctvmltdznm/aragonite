// VOIGT NOTATION HELPERS
// Proper conversions between tensor and Voigt (vector) notation

#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include <vector>

using namespace libMesh;

// Convert RankTwoTensor (3x3 symmetric) to Voigt 6-vector
// Voigt: [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
void tensorToVoigt(const RankTwoTensor & tensor, std::vector<Real> & voigt)
{
  voigt.resize(6);
  voigt[0] = tensor(0, 0);  // xx
  voigt[1] = tensor(1, 1);  // yy
  voigt[2] = tensor(2, 2);  // zz
  voigt[3] = tensor(1, 2);  // yz (or zy, symmetric)
  voigt[4] = tensor(0, 2);  // xz (or zx, symmetric)
  voigt[5] = tensor(0, 1);  // xy (or yx, symmetric)
}

// Convert Voigt 6-vector to RankTwoTensor (3x3 symmetric)
void voigtToTensor(const std::vector<Real> & voigt, RankTwoTensor & tensor)
{
  tensor.zero();
  tensor(0, 0) = voigt[0];  // xx
  tensor(1, 1) = voigt[1];  // yy
  tensor(2, 2) = voigt[2];  // zz
  tensor(1, 2) = tensor(2, 1) = voigt[3];  // yz, zy
  tensor(0, 2) = tensor(2, 0) = voigt[4];  // xz, zx
  tensor(0, 1) = tensor(1, 0) = voigt[5];  // xy, yx
}

// Convert RankFourTensor (3x3x3x3 with minor symmetries) to Voigt 6x6 matrix
// Uses Mandel notation internally for proper contraction
void rankFourToVoigt(const RankFourTensor & C, std::vector<std::vector<Real>> & voigt)
{
  voigt.resize(6, std::vector<Real>(6, 0.0));
  
  // Voigt index mapping: (i,j) -> k
  auto getVoigtIndex = [](int i, int j) -> int {
    const int map[3][3] = {{0, 5, 4},   // 00->0, 01->5, 02->4
                           {5, 1, 3},   // 10->5, 11->1, 12->3
                           {4, 3, 2}};  // 20->4, 21->3, 22->2
    return map[i][j];
  };
  
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      int p = getVoigtIndex(i, j);
      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          int q = getVoigtIndex(k, l);
          voigt[p][q] = C(i, j, k, l);
        }
      }
    }
  }
}

// Convert Voigt 6x6 matrix to RankFourTensor
void voigtToRankFour(const std::vector<std::vector<Real>> & voigt, RankFourTensor & C)
{
  C.zero();
  
  // Reverse mapping: Voigt index -> (i,j) pair
  const int voigt_to_ij[6][2] = {{0,0}, {1,1}, {2,2}, {1,2}, {0,2}, {0,1}};
  
  for (int p = 0; p < 6; p++) {
    int i = voigt_to_ij[p][0];
    int j = voigt_to_ij[p][1];
    
    for (int q = 0; q < 6; q++) {
      int k = voigt_to_ij[q][0];
      int l = voigt_to_ij[q][1];
      
      C(i, j, k, l) = voigt[p][q];
      
      // Enforce minor symmetries
      if (i != j) C(j, i, k, l) = voigt[p][q];
      if (k != l) C(i, j, l, k) = voigt[p][q];
      if (i != j && k != l) C(j, i, l, k) = voigt[p][q];
    }
  }
}

// Matrix-vector multiply: result = A * v (6x6 * 6x1 = 6x1)
void matVecMult6(const std::vector<std::vector<Real>> & A,
                 const std::vector<Real> & v,
                 std::vector<Real> & result)
{
  result.resize(6, 0.0);
  for (int i = 0; i < 6; i++) {
    result[i] = 0.0;
    for (int j = 0; j < 6; j++)
      result[i] += A[i][j] * v[j];
  }
}

// Dot product: result = a · b (6x1 · 6x1 = scalar)
Real dotProduct6(const std::vector<Real> & a, const std::vector<Real> & b)
{
  Real result = 0.0;
  for (int i = 0; i < 6; i++)
    result += a[i] * b[i];
  return result;
}

// Vector norm: ||v|| (6x1)
Real vectorNorm6(const std::vector<Real> & v)
{
  return std::sqrt(dotProduct6(v, v));
}

// Matrix-vector multiply: result = A * v (7x7 * 7x1 = 7x1)
void matVecMult7(const std::vector<std::vector<Real>> & A,
                 const std::vector<Real> & v,
                 std::vector<Real> & result)
{
  result.resize(7, 0.0);
  for (int i = 0; i < 7; i++) {
    result[i] = 0.0;
    for (int j = 0; j < 7; j++)
      result[i] += A[i][j] * v[j];
  }
}

// Dot product: result = a · b (7x1 · 7x1 = scalar)
Real dotProduct7(const std::vector<Real> & a, const std::vector<Real> & b)
{
  Real result = 0.0;
  for (int i = 0; i < 7; i++)
    result += a[i] * b[i];
  return result;
}

// CRITICAL: Contract RankFourTensor with RankTwoTensor
// result_ij = C_ijkl * A_kl
RankTwoTensor contractRankFourTwo(const RankFourTensor & C, const RankTwoTensor & A)
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

// ============================================================================
// HELPER: 7x7 Matrix Inversion using libMesh DenseMatrix
// ============================================================================
void invertMatrix7x7(const std::vector<std::vector<Real>> & A,
                     std::vector<std::vector<Real>> & Ainv)
{
  Ainv.resize(7, std::vector<Real>(7));
  
  // Solve A * X = I for each column of identity
  DenseVector<Real> rhs(7), sol(7);
  
  for (unsigned int col = 0; col < 7; col++) {
    // Create fresh copy (lu_solve modifies matrix in place)
    DenseMatrix<Real> mat(7, 7);
    for (unsigned int i = 0; i < 7; i++)
      for (unsigned int j = 0; j < 7; j++)
        mat(i, j) = A[i][j];
    
    // RHS = column of identity matrix
    for (unsigned int i = 0; i < 7; i++)
      rhs(i) = (i == col) ? 1.0 : 0.0;
    
    // Solve: mat * x = rhs
    mat.lu_solve(rhs, sol);
    
    // Store result (column of inverse)
    for (unsigned int i = 0; i < 7; i++)
      Ainv[i][col] = sol(i);
  }
}