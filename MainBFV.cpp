#include "NTL/ZZX.h"
#include "NTL/RR.h"
#include <random>
#include <utility>

using namespace std;
using namespace NTL;

class CryptoContext_BFV {
  public: 
    int poly_mod_degree;
    ZZX poly_mod;
    ZZ coeff_mod_degree;
    ZZ plain_mod_degree;

    CryptoContext_BFV(int a, ZZ b, ZZ c) {
      poly_mod_degree = a;
      coeff_mod_degree = b;
      plain_mod_degree = c;

      // Define a polynomial to keep x^poly_mod_degree + 1 because all polynomials are reduced by this
      SetCoeff(poly_mod, poly_mod_degree, 1);
      SetCoeff(poly_mod, 0, 1);
    }
};

long n_f = 100;
long n_i = 1000;  

void reduceCoeffMod(ZZX& poly, ZZ coeffMod) {
  for (long i = 0; i <= deg(poly); i++) {
    ZZ currentCoeff = coeff(poly, i) % coeffMod;
    if (currentCoeff > coeffMod / 2) {
      currentCoeff = currentCoeff - coeffMod;
    }

    SetCoeff(poly, i, currentCoeff);
  }
}

// @TODO: Check this
vector<ZZX> polyToBase(ZZX& poly, ZZ base, ZZ limit) {
  vector<ZZX> rep;

  for (ZZ i(0); i <= limit; i++) {
    ZZX poly_i;

    // find representation of coefficient j in poly
    for (long j = 0; j <= deg(poly); j++) {
      ZZ coeff_j = coeff(poly, j);

      // @TODO: not sure about this
      for (long k = 0; k < j; k++) {
        coeff_j = coeff_j / base;
      }

      coeff_j = coeff_j % base;
      if (coeff_j > base / 2) {
        coeff_j = coeff_j - base;
      }

      SetCoeff(poly_i, j, coeff_j);
    }

    rep.push_back(poly_i);
  }

  return rep;
}

ZZX genSK(CryptoContext_BFV& cryptoContext) {
  int security_param = cryptoContext.poly_mod_degree;
  ZZX poly_mod = cryptoContext.poly_mod;

  /*
  ** FV.SH.SecretKeyGen:
  **  sk <- Z2[X] / (x^n + 1)
  **  Output: sk
  */
  ZZX secretKey;
  for (int i = 0; i < security_param; i++) {
    ZZ coeff(rand() % 2);
    SetCoeff(secretKey, i, coeff);
  }
  secretKey = secretKey % poly_mod; // deg(secretKey) = 8191

  return secretKey;
}

pair<ZZX, ZZX> genPK(CryptoContext_BFV& cryptoContext, ZZX& secretKey) {
  int poly_mod_degree = cryptoContext.poly_mod_degree;
  ZZX poly_mod = cryptoContext.poly_mod;
  ZZ coeff_mod_degree = cryptoContext.coeff_mod_degree;

  // Gaussian distribution
  default_random_engine engine;
  normal_distribution<double> distribution(0.0, 2.0);

  /*
  ** FV.SH.PublicKeyGen:
  **  a <- Zq[X] / (x^n + 1)
  **  e <- Discrete Gaussian distribution over Z[X] / (x^n + 1)
  **  Output: ([-as + e] % q, a) 
  */
  ZZX a;
  for (int i = 0; i < poly_mod_degree; i++) {
    ZZ coeff = ZZ(rand());
    SetCoeff(a, i, coeff);
  }
  reduceCoeffMod(a, coeff_mod_degree);

  a = a % poly_mod; 

  ZZX e;  
  for (int i = 0; i < poly_mod_degree; i++) {
    long coeff = round(distribution(engine));
    SetCoeff(e, i, coeff);
  }
  e = e % poly_mod;

  ZZX pk0 = - a * secretKey + e;
  ZZX pk1 = a;

  pair<ZZX, ZZX> publicKey;
  publicKey.first = pk0;
  publicKey.second = pk1;

  return publicKey;
}

vector<pair<ZZX, ZZX>> genRLK(CryptoContext_BFV& cryptoContext, ZZX& secretKey, ZZ w) {
  int poly_mod_degree = cryptoContext.poly_mod_degree;
  ZZX poly_mod = cryptoContext.poly_mod;
  ZZ coeff_mod_degree = cryptoContext.coeff_mod_degree;
  ZZ plain_mod_degree = cryptoContext.plain_mod_degree;

  vector<pair<ZZX, ZZX>> rlk;
  
  // l = floor(log_w(q))
  RR q_log = log(MakeRR(coeff_mod_degree, 0));
  RR w_log = log(MakeRR(w, 0));
  ZZ l = FloorToZZ(q_log / w_log);

  ZZ w_power(1);
  for (ZZ i(0); i <= l; i++) {
    //!! w_poer should be w^i

    // Sample a_i from R_coeff_mod_degree
    ZZX a_i;
    for (int i = 0; i < poly_mod_degree; i++) {
      ZZ coeff = ZZ(rand());
      SetCoeff(a_i, i, coeff);
    }
    reduceCoeffMod(a_i, coeff_mod_degree);

    a_i = a_i % poly_mod;

    // Sample e_i from Gaussian distribution
    default_random_engine engine;
    normal_distribution<double> distribution(0.0, 2.0);

    ZZX e_i;  
    for (int i = 0; i < poly_mod_degree; i++) {
      long coeff = round(distribution(engine));
      SetCoeff(e_i, i, coeff);
    }
    e_i = e_i % poly_mod;
    
    // @TODO: Correct this to use some predefined function for power for computing w^i
    ZZX rlk_element_0 = - a_i * secretKey - e_i + w_power * secretKey * secretKey;
    ZZX rlk_element_1 = a_i;

    pair<ZZX, ZZX> rlk_element_i;
    rlk_element_i.first = rlk_element_0;
    rlk_element_i.second = rlk_element_1;

    rlk.push_back(rlk_element_i);

    w_power = w_power * w; 
  }

  return rlk; 
}

ZZX integerEncoder(ZZ value) {
  int sign = value >= 0 ? 1 : -1;
  value = value * sign;

  ZZX value_base_2;

  long counter = 0;
  while (value > 0) {
    long coeff = sign * (value % 2);
    SetCoeff(value_base_2, counter, coeff);

    counter++;
    value = value / 2;
  }

  // No need to reduce coeff mod plain or poly.mod because these polynomials should be small
  return value_base_2;
}

ZZ integerDecoder(ZZX poly) {
  ZZ value(0);

  for (long i = 0; i <= deg(poly); i++) {
    value = value + pow(2, i) * coeff(poly, i);
  }

  return value;
}

// When the rational number can be expressed in a finite binary form (using less than n_f terms for the fractional part)
// There is 0 precision loss
// e.g. 22.8125 * 33.8125
ZZX fractionalEncoder(double rational_number, CryptoContext_BFV& cryptoContext) {
  long poly_mod_degree = cryptoContext.poly_mod_degree;
 
  int sign = rational_number >= 0 ? 1 : -1;
  long integer_part = floor(rational_number);
  double fractional_part = rational_number - integer_part;

  ZZX poly = integerEncoder(ZZ(integer_part));

  for (long i = 1; i <= n_f; i++) {
    int bit = floor(2 * fractional_part);
    fractional_part = 2 * fractional_part - bit;

    SetCoeff(poly, poly_mod_degree - i, (-sign)*bit);
  }

  return poly;
}

double fractionalDecoder(ZZX poly, CryptoContext_BFV& cryptoContext) {
  double value = 0;
  
  for (long i = 0; i < n_i; i++) {
    value = value + conv<long>(pow(2, i) * coeff(poly, i));
  }

  for (long i = 1; i <= deg(poly) - n_i; i++) {
    ZZ coeff_i = coeff(poly, conv<long>(deg(poly) + 1 - i));
    value = value - conv<double>(coeff_i) / pow(2, i);  // Subtract because sign was changed when coefficient were assigned for the fractional part
  }

  return value;
}

pair<ZZX, ZZX> encrypt(ZZX& plaintext, CryptoContext_BFV& cryptoContext, pair<ZZX, ZZX>& publicKey) {
  int poly_mod_degree = cryptoContext.poly_mod_degree;
  ZZX poly_mod = cryptoContext.poly_mod;
  ZZ coeff_mod_degree = cryptoContext.coeff_mod_degree;
  ZZ plain_mod_degree = cryptoContext.plain_mod_degree;

  ZZX pk0 = publicKey.first;
  ZZX pk1 = publicKey.second;

  // Gaussian distribution
  default_random_engine engine;
  normal_distribution<double> distribution(0.0, 2.0);

  /* Encrypt:
  **  Plaintext in Zt[X] / (x^n + 1)
  **  Output: pair of polynomials (ct0, ct1) in Zq[X] / (x^n + 1) * Zq[X] / (x^n + 1)
  */

  ZZ delta = coeff_mod_degree / plain_mod_degree;
  ZZ r_q = coeff_mod_degree % plain_mod_degree;

  // Sample u with binary coefficients
  ZZX u;
  for (int i = 0; i < poly_mod_degree; i++) {
    int coeff = rand() % 2;
    SetCoeff(u, i, coeff);
  }
  u = u % poly_mod;

  // Sample e1 and e2 from Gaussian distribution
  ZZX e1;
  ZZX e2;
  for (int i = 0; i < poly_mod_degree; i++) {
    long coeff1 = round(distribution(engine)); 
    long coeff2 = round(distribution(engine));

    SetCoeff(e1, i, coeff1);
    SetCoeff(e2, i, coeff2);
  }
  e1 = e1 % poly_mod;
  e2 = e2 % poly_mod;

  // Encrypt
  ZZX ct0 = pk0 * u + e1 + delta * plaintext;
  ZZX ct1 = pk1 * u + e2;
  ct0 = ct0 % poly_mod;
  ct1 = ct1 % poly_mod;
  reduceCoeffMod(ct0, coeff_mod_degree);
  reduceCoeffMod(ct1, coeff_mod_degree);

  pair<ZZX, ZZX> ciphertext;
  ciphertext.first = ct0;
  ciphertext.second = ct1;

  return ciphertext;
}

ZZX decrypt(pair<ZZX, ZZX> ciphertext, CryptoContext_BFV& cryptoContext, ZZX secretKey) {
  int poly_mod_degree = cryptoContext.poly_mod_degree;
  ZZX poly_mod = cryptoContext.poly_mod;
  ZZ coeff_mod_degree = cryptoContext.coeff_mod_degree;
  ZZ plain_mod_degree = cryptoContext.plain_mod_degree;

  ZZX ct0 = ciphertext.first;
  ZZX ct1 = ciphertext.second;

  // Decryption
  ZZX mid = ct0 + ct1 * secretKey;
  mid = mid % poly_mod;
  reduceCoeffMod(mid, coeff_mod_degree);

  // Need to multiply each coefficient by t / q then round to the nearest integer and reduce mod t
  ZZX decrypted;

  for (int i = 0; i < poly_mod_degree; i++) {
    RR dec_coeff = MakeRR(plain_mod_degree, 0) * MakeRR(coeff(mid, i), 0) / MakeRR(coeff_mod_degree, 0);
    
    ZZ round_coeff = RoundToZZ(dec_coeff);
    round_coeff = round_coeff % plain_mod_degree;
    if (round_coeff > plain_mod_degree / 2) {
      round_coeff = round_coeff - plain_mod_degree;
    } 

    SetCoeff(decrypted, i, round_coeff);
  }

  return decrypted;
}

pair<ZZX, ZZX> addCiphertexts(pair<ZZX, ZZX>& ciphertext1, pair<ZZX, ZZX>& ciphertext2) {
  pair<ZZX, ZZX> sum;
  sum.first = ciphertext1.first + ciphertext2.first;
  sum.second = ciphertext1.second + ciphertext2.second;

  return sum;
}

pair<ZZX, ZZX> multCiphertexts(pair<ZZX, ZZX>& ciphertext0, pair<ZZX, ZZX>& ciphertext1, CryptoContext_BFV& cryptoContext, vector<pair<ZZX, ZZX>> rlk) {
  int poly_mod_degree = cryptoContext.poly_mod_degree;
  ZZX poly_mod = cryptoContext.poly_mod;
  ZZ coeff_mod_degree = cryptoContext.coeff_mod_degree;
  ZZ plain_mod_degree = cryptoContext.plain_mod_degree;

  ZZX c0 = ciphertext0.first * ciphertext1.first;
  ZZX c1 = ciphertext0.first * ciphertext1.second + ciphertext0.second * ciphertext1.first;
  ZZX c2 = ciphertext0.second * ciphertext1.second;

  c0 = c0 % poly_mod;
  c1 = c1 % poly_mod;
  c2 = c2 % poly_mod;

  // Need to multiply each coefficient by t / q then round to the nearest integer and reduce mod q
  for (int i = 0; i < poly_mod_degree; i++) {
    RR c0_coeff = MakeRR(plain_mod_degree, 0) * MakeRR(coeff(c0, i), 0) / MakeRR(coeff_mod_degree, 0);
    RR c1_coeff = MakeRR(plain_mod_degree, 0) * MakeRR(coeff(c1, i), 0) / MakeRR(coeff_mod_degree, 0);
    RR c2_coeff = MakeRR(plain_mod_degree, 0) * MakeRR(coeff(c2, i), 0) / MakeRR(coeff_mod_degree, 0);
    
    ZZ c0_round_coeff = RoundToZZ(c0_coeff) % coeff_mod_degree;
    ZZ c1_round_coeff = RoundToZZ(c1_coeff) % coeff_mod_degree;
    ZZ c2_round_coeff = RoundToZZ(c2_coeff) % coeff_mod_degree;

    if (c0_round_coeff > coeff_mod_degree / 2) {
      c0_round_coeff = c0_round_coeff - coeff_mod_degree;
    } 

    if (c1_round_coeff > coeff_mod_degree / 2) {
      c1_round_coeff = c1_round_coeff - coeff_mod_degree;
    } 

    if (c2_round_coeff > coeff_mod_degree / 2) {
      c2_round_coeff = c2_round_coeff - coeff_mod_degree;
    } 

    SetCoeff(c0, i, c0_round_coeff);
    SetCoeff(c1, i, c1_round_coeff);
    SetCoeff(c2, i, c2_round_coeff);
  }

  // Express c2 in base w and recompute c0' and c1'
  ZZ base = FloorToZZ(sqrt(MakeRR(coeff_mod_degree, 0)));
  RR q_log = log(MakeRR(coeff_mod_degree, 0));
  RR w_log = log(MakeRR(base, 0));
  ZZ limit = FloorToZZ(q_log / w_log);

  vector<ZZX> c2_base = polyToBase(c2, base, limit);

  ZZX c0_prime = c0;
  ZZX c1_prime = c1;

  for (ZZ i(0); i <= limit; i++) {
    pair<ZZX, ZZX> rlk_i = rlk.at(to_long(i));
    c0_prime = c0_prime + rlk_i.first * c2_base.at(to_long(i));
    c1_prime = c1_prime + rlk_i.second * c2_base.at(to_long(i));

    c0_prime = c0_prime % poly_mod;
    c1_prime = c1_prime % poly_mod;
  }

  c0_prime = c0_prime % poly_mod;
  c1_prime = c1_prime % poly_mod;

  pair<ZZX, ZZX> multEncrypted;
  multEncrypted.first = c0_prime;
  multEncrypted.second = c1_prime;

  return multEncrypted;
}

void testMultiplication(CryptoContext_BFV& cryptoContext, pair<ZZX, ZZX> publicKey, ZZX secretKey, vector<pair<ZZX, ZZX>> rlk) {
  int poly_mod_degree = cryptoContext.poly_mod_degree;
  ZZX poly_mod = cryptoContext.poly_mod;
  ZZ coeff_mod_degree = cryptoContext.coeff_mod_degree;
  ZZ plain_mod_degree = cryptoContext.plain_mod_degree;

  for (int test_nr = 0; test_nr < 50; test_nr++) {
    ZZX plain1;
    for (int i = 0; i < poly_mod_degree / 8; i++) {
      ZZ coeff = ZZ(rand()) % plain_mod_degree;
      if (coeff > plain_mod_degree / 2) {
        coeff = coeff - plain_mod_degree;
      }

      SetCoeff(plain1, i, coeff);
    }
    plain1 = plain1 % cryptoContext.poly_mod;

    ZZX plain2;
    for (int i = 0; i < poly_mod_degree / 8; i++) {
      ZZ coeff = ZZ(rand()) % plain_mod_degree;
      if (coeff > plain_mod_degree / 2) {
        coeff = coeff - plain_mod_degree;
      }

      SetCoeff(plain2, i, coeff);
    }
    plain2 = plain2 % cryptoContext.poly_mod;

    ZZX plain3;
    for (int i = 0; i < poly_mod_degree / 8; i++) {
      ZZ coeff = ZZ(rand()) % plain_mod_degree;
      if (coeff > plain_mod_degree / 2) {
        coeff = coeff - plain_mod_degree;
      }

      SetCoeff(plain3, i, coeff);
    }
    plain3 = plain3 % cryptoContext.poly_mod;

    auto ciphertext1 = encrypt(plain1, cryptoContext, publicKey);
    auto ciphertext2 = encrypt(plain2, cryptoContext, publicKey);
    auto ciphertext3 = encrypt(plain3, cryptoContext, publicKey);

    auto ciphertextMult1 = multCiphertexts(ciphertext1, ciphertext2, cryptoContext, rlk);
    auto ciphertextMult2 = multCiphertexts(ciphertext3, ciphertextMult1, cryptoContext, rlk);
    auto decryptedMult = decrypt(ciphertextMult2, cryptoContext, secretKey);

    auto plainMult = plain1 * plain2 * plain3;
    reduceCoeffMod(plainMult, plain_mod_degree);

    if (plainMult == decryptedMult) {
      cout << "Test " << test_nr << " succeeded" << endl;
    } else {
      cout << "Test " << test_nr << " failed" << endl;
    }
  }
}


int main() {
  srand((unsigned) time(0));

  int poly_mod_degree = 8192;
  ZZ coeff_mod_degree(1);
  coeff_mod_degree = coeff_mod_degree * GenPrime_ZZ(60);
  coeff_mod_degree = coeff_mod_degree * GenPrime_ZZ(40);
  coeff_mod_degree = coeff_mod_degree * GenPrime_ZZ(40);
  coeff_mod_degree = coeff_mod_degree * GenPrime_ZZ(60);
  ZZ plain_mod_degree(1000000); 

  CryptoContext_BFV cryptoContext(poly_mod_degree, coeff_mod_degree, plain_mod_degree);

  auto secretKey = genSK(cryptoContext);
  auto publicKey = genPK(cryptoContext, secretKey);

  // Generate rlk and specify base !!
  ZZ base = FloorToZZ(sqrt(MakeRR(coeff_mod_degree, 0)));
  auto rlk = genRLK(cryptoContext, secretKey, base);

  // cout << endl << endl << endl << "***Test***" << endl;
  // testMultiplication(cryptoContext, publicKey, secretKey, rlk);

  double val1;
  cout << "Val1: ";
  cin >> val1;
  ZZX plain_1 = fractionalEncoder(val1, cryptoContext);

  double val2;
  cout << "Val2: ";
  cin >> val2;
  ZZX plain_2 = fractionalEncoder(val2, cryptoContext);

  auto ciphertext1 = encrypt(plain_1, cryptoContext, publicKey);
  auto ciphertext2 = encrypt(plain_2, cryptoContext, publicKey);

  auto ciphertext_sum = multCiphertexts(ciphertext1, ciphertext2, cryptoContext, rlk);
  ZZX plain_sum = decrypt(ciphertext_sum, cryptoContext, secretKey);
  double val_sum = fractionalDecoder(plain_sum, cryptoContext);

  cout << "Prod computed by library: " << val_sum << endl;
  cout << "Actual prod: " << (val1 * val2) << endl;

  return 0;
}