#include "NTL/ZZX.h"
#include "NTL/RR.h"
#include <random>
#include <utility>
#include <complex.h>
#include <math.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>

using namespace std;
using namespace NTL;
using namespace Eigen;

class CryptoContext_CKKS {
  public: 
    int poly_mod_degree;
    ZZX poly_mod;
    ZZ coeff_mod_degree;

    CryptoContext_CKKS(int a, ZZ b) {
      poly_mod_degree = a;
      coeff_mod_degree = b;

      // Define a polynomial to keep x^poly_mod_degree + 1 because all polynomials are reduced by this
      SetCoeff(poly_mod, poly_mod_degree, 1);
      SetCoeff(poly_mod, 0, 1);
    }
};

/* 
  USED FOR ENCODING
  Input: 
  **  z = a vector of complex values
  **  N = poly_mod_degree
  
  Output:
    A polynomial with complex coefficients reduces modulo X^N + 1
*/
vector<dcomplex> canonical_embedding_inverse(vector<dcomplex> z, long N) {
  if (z.size() != N) {
    throw "Wrong vector size in canonical embedding";
  }

  Matrix<dcomplex, Dynamic, Dynamic> vandermonde(N, N);
  Matrix<dcomplex, Dynamic, Dynamic> z_as_column_vec(N, 1);

  dcomplex root_cyclotomic;
  root_cyclotomic.real(cos(M_PI / N));
  root_cyclotomic.imag(sin(M_PI / N));

  for (long i = 0; i < N; i++) {
    for (long j = 0; j < N; j++) {
      vandermonde(i, j) = pow(root_cyclotomic, (2 * i + 1) * j);
    }
  }

  for (long i = 0; i < N; i++) {
    z_as_column_vec(i, 0) = z.at(i);
  }

  auto inverse_vandermonde = vandermonde.inverse();
  auto poly = inverse_vandermonde * z_as_column_vec;

  vector<dcomplex> poly_as_vec;
  for (long i = 0; i < N; i++) {
    poly_as_vec.push_back(poly(i, 0));
  }

  return poly_as_vec;
}


/* 
  Used for DECODING
  Input: 
  **  m = a polynomial in C[X] / (X^N + 1)
  **  N = poly_mod_degree
  
  Output:
  **  A vector with N complex values
*/
vector<dcomplex> canonical_embedding(vector<dcomplex> m, long N) {
  if (m.size() != N) {
    throw "Wrong polynomial degree in canonical embedding";
  }

  Matrix<dcomplex, Dynamic, Dynamic> vandermonde(N, N);
  Matrix<dcomplex, Dynamic, Dynamic> m_as_column_vec(N, 1);

  dcomplex root_cyclotomic;
  root_cyclotomic.real(cos(M_PI / N));
  root_cyclotomic.imag(sin(M_PI / N));

  for (long i = 0; i < N; i++) {
    for (long j = 0; j < N; j++) {
      vandermonde(i, j) = pow(root_cyclotomic, (2 * i + 1) * j);
    }
  }

  for (long i = 0; i < N; i++) {
    m_as_column_vec(i, 0) = m.at(i);
  }

  auto z_vec = vandermonde * m_as_column_vec;

  vector<dcomplex> z_as_vec;
  for (long i = 0; i < N; i++) {
    z_as_vec.push_back(z_vec(i, 0));
  }

  return z_as_vec;
}

// Expand vector from C^N/2 into H by taking conjugates in reverse order
vector<dcomplex> pi_inv(vector<dcomplex> z, long poly_mod_degree) {
  if (z.size() != poly_mod_degree / 2) {
    throw "Vector must have size " + poly_mod_degree / 2;
  }

  vector<dcomplex> vec_expanded;
  for (int i = 0; i < z.size(); i++) {
    vec_expanded.push_back(z.at(i));
  }

  for (int i = 0; i < z.size(); i++) {
    vec_expanded.push_back(conj(z.at(poly_mod_degree / 2 - 1 - i)));
  }

  return vec_expanded;
}

// If vector is in H, apply pi by only taking its first poly_mod_degree / 2 coordinates (the others are just conjugates)
vector<dcomplex> pi(vector<dcomplex> z, long poly_mod_degree) {
  if (z.size() != poly_mod_degree) {
    throw "Vector in H must have size " + poly_mod_degree;
  }

  vector<dcomplex> vec_shrinked;
  for (int i = 0; i < poly_mod_degree / 2; i++) {
    vec_shrinked.push_back(z.at(i));
  }

  return vec_shrinked;
}

double hermitian_prod(vector<dcomplex> z, vector<dcomplex> b_i) {
  if (z.size() != b_i.size()) {
    throw "Vectors in hermitian product must have the same size";
  }

  dcomplex result;
  result.real(0.0);
  result.imag(0.0);

  for (int i = 0; i < z.size(); i++) {
    result = result + z.at(i) * conj(b_i.at(i));
  }

  // Result should be a real number
  return result.real();
}

vector<long> project_to_ZX(vector<dcomplex> z, long N) {
  if (z.size() != N) {
    throw "Wrong vector size in projection of H to Z[X]";
  }

  vector<long> z_coeff_proj;
  
  dcomplex root_cyclotomic;
  root_cyclotomic.real(cos(M_PI / N));
  root_cyclotomic.imag(sin(M_PI / N));

  for (int k = 0; k < z.size(); k++) {
    dcomplex eps_k = pow(root_cyclotomic, k);
    
    vector<dcomplex> b_k;
    for (int i = 1; i <= 2 * N - 1; i+=2) {
      b_k.push_back(pow(eps_k, i));
    }

    double herm_prod_z_bk = hermitian_prod(z, b_k);
    double bk_norm = hermitian_prod(b_k, b_k);  // == poly_mod_degree

    double coeff = herm_prod_z_bk / bk_norm;
    double prob_weight = coeff - floor(coeff);

    random_device rd;
    mt19937 gen(rd());
    binomial_distribution<> round_prob(1, prob_weight);

    long coeff_k = floor(coeff) + round_prob(gen);
    z_coeff_proj.push_back(coeff_k);
  }

  return z_coeff_proj;
}

void printVector(vector<dcomplex> vec) {
  cout << "Vector: " << endl;
  for (dcomplex d : vec) {
    cout << d.real() << " + i * " << d.imag() << endl;
  }

  cout << endl << endl;
} 





int main() {
  long poly_mod_degree = 16;
  double delta_precision = 64;

  vector<dcomplex> z;
  for (int i = 0; i < poly_mod_degree / 2; i++) {
    dcomplex val;
    val.real(100.1 + i);
    val.imag(20.1 + i);

    // Multiply by delta for precision
    val.real(val.real() * delta_precision);
    val.imag(val.imag() * delta_precision);
    
    z.push_back(val);
  }
  
  // Test example given in paper
  // dcomplex val1;
  // val1.real(3 * delta_precision);
  // val1.imag(4 * delta_precision);

  // dcomplex val2;
  // val2.real(2 * delta_precision);
  // val2.imag(-1 * delta_precision);

  // z in C^N/2
  auto z_expanded = pi_inv(z, poly_mod_degree);
  auto z_projected_on_sigma_r = project_to_ZX(z_expanded, poly_mod_degree);

  vector<dcomplex> z_projected_as_complex;
  for (long val : z_projected_on_sigma_r) {
    dcomplex val_as_complex;
    val_as_complex.real(val);
    val_as_complex.imag(0);

    z_projected_as_complex.push_back(val_as_complex);
  }

  // auto z_encoded_as_complex = canonical_embedding_inverse(z_projected_as_complex, poly_mod_degree);
  auto z_encoded_as_complex = z_projected_as_complex;
  vector<long> z_encoded_final;

  for (dcomplex val : z_encoded_as_complex) {
    z_encoded_final.push_back(floor(val.real()));
  }

  // Decode
  vector<dcomplex> rescaled_z;
  for (long val : z_encoded_final) {
    dcomplex val_rescaled_as_complex;
    val_rescaled_as_complex.real(val * 1.0 / delta_precision);
    val_rescaled_as_complex.imag(0);

    rescaled_z.push_back(val_rescaled_as_complex);
  }

  auto apply_canonical_embedding = canonical_embedding(rescaled_z, poly_mod_degree);
  auto apply_pi = pi(apply_canonical_embedding, poly_mod_degree);

  printVector(apply_pi);

  return 0;
}