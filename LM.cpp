#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include "utils.cpp"
#pragma once

using namespace std;

#define DELTA_DRV 0.01
#define DELTA_CONVERGE 0.00001
#define CONST_LM 1.1

double* initial_parameters;




void model_func(int, double*, double*, double*, void*);
void GN(int, int, double*, double*, void*, void (*)(int, double*, double*, double*, void*));
void LM(int, int, double*, double*, void*, void (*)(int, double*, double*, double*, void*));
void multiply_const_matrix(int, int, int, double*, double*);
void add_matrix(int, int, double*, double*, double*);
void subtract_matrix(int, int, double*, double*, double*);
void multiply_matrix(int, int, int, int, double*, double*, double*);
void diagonal_matrix(int, double*, double*);
void inverse_matrix(int, double*, double*);
void LU_decompose(int, double*, double*, double*);
double error_i(int, double*, double*, double*, void*, void(*)(int, double*, double*, double*, void*));
double get_S(int, double*, double*, double*, void*, void(*)(int, double*, double*, double*, void*));
void de_vec_dx_j(int, int, int, double*, double*, double*, double*, void*, void(*)(int, double*, double*, double*, void*));
void Hessian(int, int, double*, double*, double*, double*, void*, void(*)(int, double*, double*, double*, void*));
void dS_dpara(int, int, double*, double*, double*, void*, double*, void(*)(int, double*, double*, double*, void*));
bool check_if_converged(int, double*, double*);
void display_matrix(int, int, double*);

/*
int main() {
  string filename = "data";
  int datanum = get_datanum(filename, 1, 2);
  double X[datanum];
  double Y[datanum];
  void* others;
  void (*model)(int, double*, double*, double*, void*) = model_func;
  read_file(filename, 1, 2, X, Y);

  LM(datanum, 3, X, Y, others, model);
}
*/



void model_func(int datalen, double* X, double* parameters, double* ret, void* others) {
  for (int i=0; i<datalen; i++) {
    ret[i] = parameters[0]*X[i]*X[i] + parameters[1]*X[i] + parameters[2];
  }
}

void GN(int datanum, int num_of_parameters, double* X, double* Y, void* others, void (*model)(int, double*, double*, double*, void*)) {
  
  double parameters[num_of_parameters];
  copy_arr(num_of_parameters, initial_parameters, parameters);
  double H[num_of_parameters][num_of_parameters];
  double H_inverse[num_of_parameters][num_of_parameters];
  double dX[num_of_parameters];
  double dSdp[num_of_parameters];
  
  while (true) {
    Hessian(datanum, num_of_parameters, X, Y, parameters, &H[0][0], others, model);    
    inverse_matrix(num_of_parameters, &H[0][0], &H_inverse[0][0]);
    dS_dpara(datanum, num_of_parameters, X, Y, parameters, others, dSdp,  model);  
    multiply_matrix(num_of_parameters, num_of_parameters, num_of_parameters, 1, &H_inverse[0][0], dSdp, dX);
    display_matrix(1, 3, parameters);
    subtract_matrix(1, num_of_parameters, parameters, dX, parameters);
  }
}

void LM(int datanum, int num_of_parameters, double* X, double* Y, void* others, void (*model)(int, double*, double*, double*, void*)) {

  double parameters[num_of_parameters];
  double parameters_buf[num_of_parameters];
  copy_arr(num_of_parameters, initial_parameters, parameters);
  double H[num_of_parameters][num_of_parameters];
  double H_diagonal[num_of_parameters][num_of_parameters];
  double H_inverse[num_of_parameters][num_of_parameters];
  double dX[num_of_parameters];
  double dSdp[num_of_parameters];
  double C = 0.5;
  int iter_count = 1;
  double S, S_buf;
  
  while (true) {
    Hessian(datanum, num_of_parameters, X, Y, parameters, &H[0][0], others, model);
    diagonal_matrix(num_of_parameters, &H[0][0], &H_diagonal[0][0]);
    multiply_const_matrix(C, num_of_parameters, num_of_parameters, &H_diagonal[0][0], &H_diagonal[0][0]);
    add_matrix(num_of_parameters, num_of_parameters, &H[0][0], &H_diagonal[0][0], &H[0][0]);
    inverse_matrix(num_of_parameters, &H[0][0], &H_inverse[0][0]);
    if (!check_nan_matrix(num_of_parameters, num_of_parameters, &H_inverse[0][0])) {
      cerr << "Error: Cannot calculate inverse matrix of Hessian" << endl;
      exit(0);
    }
    dS_dpara(datanum, num_of_parameters, X, Y, parameters, others, dSdp, model);
    multiply_matrix(num_of_parameters, num_of_parameters, num_of_parameters, 1, &H_inverse[0][0], dSdp, dX);
    subtract_matrix(1, num_of_parameters, parameters, dX, parameters_buf);
    
    S = get_S(datanum, X, Y, parameters, others, model);
    S_buf = get_S(datanum, X, Y, parameters_buf, others, model);
    
    cout << "---------------------------------------------------" << endl;
    cout << "iteration: " << iter_count++ << endl;
    cout << "C value: " << C << endl;
    cout << "current parameters:" << endl;
    display_matrix(1, num_of_parameters, parameters);
    cout << "Next parameters:" << endl;
    display_matrix(1, num_of_parameters, parameters_buf);
    cout << "S(current) S(next)" << endl;
    cout << S << " " << S_buf << endl;

    if (check_if_converged(num_of_parameters, parameters, parameters_buf)) {
      cout << "===========================" << endl;
      cout << "final set of parameters:" << endl;
      display_matrix(1, num_of_parameters, parameters);
      return;
    }
    
    if (S_buf < S && check_nan_matrix(num_of_parameters, num_of_parameters, &H_inverse[0][0])) {
      subtract_matrix(1, num_of_parameters, parameters, dX, parameters);
      C /= CONST_LM;
    } else {
      C *= CONST_LM;
    }

  }
}

/* =========================functions==========================*/


void multiply_const_matrix(int cst, int l, int c, double* X, double* ret) {
  for (int i=0; i<l; i++) {
    for (int j=0; j<c; j++) {
      *(ret + i*c + j) = cst * (*(X + i*c + j));
    }
  }
}

void add_matrix(int l, int c, double* A, double* B, double* ret) {
  for (int i=0; i<l; i++) {
    for (int j=0; j<c; j++) {
      *(ret + i*c + j) = *(A + i*c + j) + *(B + i*c + j);
    }
  }
}

void subtract_matrix(int l, int c, double* A, double* B, double* ret) {
  for (int i=0; i<l; i++) {
    for (int j=0; j<c; j++) {
      *(ret + i*c + j) = *(A + i*c + j) - *(B + i*c + j);
    }
  }
}

/* accepts square matrix only */
/* extracts diagonal components of given matrix */
/* fills 0 for other components */
void diagonal_matrix(int n, double* X, double* ret) {
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (i == j) {
	*(ret + i*n + j) = *(X + i*n + j);
      } else {
	*(ret + i*n + j) = 0;
      }
    }
  }
}

void multiply_matrix(int X_l, int X_c, int Y_l, int Y_c, double* X, double* Y, double* ret) {
  
  /* l is for line, c is for column */
  /* X[line][column]                */
  /* returns ret[X_l][Y_c]          */

  if (X_c != Y_l) {
    cerr << "Cannot multiple matrix due to dimension mismatch in function multiple_matrix" << endl;
    exit(0);
  }

  for (int i=0; i<X_l; i++) {
    for (int j=0; j<Y_c; j++) {
      *(ret + i*Y_c + j) = 0;
      for (int k=0; k<X_c; k++) {
	*(ret + i*Y_c + j) += *(X + i*X_c + k) * (*(Y + k*Y_c + j));
      }
    }
  }
}

/* accepts square matrix only */
/* use LU decomposition */
void inverse_matrix(int n, double* X, double* ret) {

  if (n == 1) {
    *ret = 1/(*X);
    return;
  }
    
  double L[n][n];
  double U[n][n];

  /* AA^(-1)  = I                     */
  /* LUA^(-1) = I                     */
  /* B = UA^(-1)                      */
  /* solve LB = I to find B           */
  /* solve B = UA^(-1) to find A^(-1) */
  
  LU_decompose(n, X, &L[0][0], &U[0][0]);

  /* matrix B should look like this  */
  /* 1/l[0][0]     0         0       */
  /* B[1][0]   1/l[1][1]     0       */
  /* B[2][0]    B[2][1]  1/l[2][2]   */

  double B[n][n];
  double sum = 0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (i < j) {B[i][j] = 0;}
      else if (i == j) {B[i][j] = 1/L[i][j];}
      else {
	sum = 0;
	for (int k=0; k<i; k++) {
	  sum += B[k][j] * L[i][k];
	}
	B[i][j] = -sum / L[i][i];
      }
    }
  }

  for (int i=0; i<n; i++) {
    *(ret + n*(n-1) + i) = B[n-1][i];
  }


  // n=4
  for (int i=n-2; i>=0; i--) { // i = 2->0
    for (int j=0; j<n; j++) {  // j = 0->3
      sum = 0;
      for (int k=i+1; k<n; k++) { // k = 0->2 (i=0)
	sum += *(ret + k*n + j) * U[i][k];
      }
      *(ret + i*n + j) = B[i][j] - sum;
    }
  }
}

void LU_decompose(int n, double* A, double* L, double* U) {
  
  /* A[i][j] = *(A + i*n + j) (assumes A is an n dimensional matrix) */

  /* substitute 0s in matrix L */
  /* like below (L is 4 dimensional for ex.) */

  /* L[0][0]    0       0       0    */
  /* L[1][0] L[1][1]    0       0    */
  /* L[2][0] L[2][1] L[2][2]    0    */
  /* L[3][0] L[3][1] L[3][2] L[3][3] */
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (i < j) {*(L + i*n + j) = 0;}
    }
  }

  /* substitute 0s and 1s in matrix U */
  /* like below (L is 4 dimensional for ex.) */

  /* 1 U[0][1] U[0][2] U[0][3] */
  /* 0    1    U[1][2] U[1][3] */
  /* 0    0       1    U[2][3] */
  /* 0    0       0       1    */

  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (i > j) {*(U + i*n + j) = 0;}
      if (i == j) {*(U + i*n + j) = 1;}
    }
  }
  
  /* find L[0][0] L[1][0] L[2][0] ... */
  for (int i=0; i<n; i++) {
    *(L + i*n) = *(A + i*n);
  }

  /* find u[0][1] u[0][2] u[0][3] ... */
  for (int i=1; i<n; i++) {
    *(U + i) = *(A + i) / *(L);
  }

  double A_tmp[n-1][n-1];
  double L_tmp[n-1][n-1];
  double U_tmp[n-1][n-1];
  double tmp[n-1][n-1];

  double l_0_vector[n-1][1];
  double u_0_vector[1][n-1];

  /* copy matrix A into A_tmp */
  for (int i=0; i<n-1; i++) {
    for (int j=0; j<n-1; j++) {
      *(&A_tmp[0][0] + i*(n-1) + j) = *(A + (i+1)*n + (j+1));
    }
  }

  /* copy l0 and u0 vector */
  for (int i=0; i<n-1; i++) {
    *(&l_0_vector[0][0] + i) = *(L + (i+1)*n);
    *(&u_0_vector[0][0] + i) = *(U + (i+1));
  }
  
  multiply_matrix(n-1, 1, 1, n-1, &l_0_vector[0][0], &u_0_vector[0][0], &tmp[0][0]);
  subtract_matrix(n-1, n-1, &A_tmp[0][0], &tmp[0][0], &A_tmp[0][0]);
  if (n == 2) {
    *(L + 3) = *(A + 3) - tmp[0][0];
    return;
  } else {
    LU_decompose(n-1, &A_tmp[0][0], &L_tmp[0][0], &U_tmp[0][0]);
  }

  /* copy L_tmp and U_tmp into L and U */
  for (int i=0; i<n-1; i++) {
    for (int j=0; j<n-1; j++) {
       *(L + (i+1)*n + (j+1)) = *(&L_tmp[0][0] + i*(n-1) + j);
       *(U + (i+1)*n + (j+1)) = *(&U_tmp[0][0] + i*(n-1) + j);
    }
  }
}

void get_errors(int datalen, double* X, double* Y, double* parameters, double* ret, void* others, void (*model)(int, double*, double*, double*, void*)) {
  double buf[datalen];
  model(datalen, X, parameters, buf, others);
  for (int i=0; i<datalen; i++) {
    ret[i] = Y[i] - buf[i];
  }
}

double get_S(int datalen, double* X, double* Y, double* parameters, void* others, void (*model)(int, double*, double*, double*, void*)) {
  double sum = 0;
  double buf[datalen];
  get_errors(datalen, X, Y, parameters, buf, others, model);
  for (int i=0; i<datalen; i++) {
    sum += (buf[i] * buf[i]);
  }
  return sum;
}

void de_vec_dx_j(int datalen, int j, int num_of_parameters, double* X, double* Y, double* parameters, double* ret, void* others, void (*model)(int, double*, double*, double*, void*)) {
  double tmp[num_of_parameters];
  copy_arr(num_of_parameters, parameters, tmp);
  tmp[j] += DELTA_DRV;
  double buf[datalen];
  double buf2[datalen];
  
  get_errors(datalen, X, Y, tmp, buf, others, model);
  get_errors(datalen, X, Y, parameters, buf2, others, model);
  for (int i=0; i<datalen; i++) {
    ret[i] = (buf[i] - buf2[i]) / DELTA_DRV;
  }
}

void Hessian(int datalen, int num_of_parameters, double* X, double* Y, double* parameters, double* ret, void* others, void (*model)(int, double*, double*, double*, void*)) {
  double buf[num_of_parameters][datalen];

  for (int i=0; i<num_of_parameters; i++) {
    de_vec_dx_j(datalen, i, num_of_parameters, X, Y, parameters, &buf[i][0], others, model);
  }
  
  for (int j=0; j<num_of_parameters; j++) {
    for (int k=0; k<num_of_parameters; k++) {
      ret[j*num_of_parameters + k] = 0;
      for (int i=0; i<datalen; i++) {
	ret[j*num_of_parameters + k] += (buf[j][i] * buf[k][i]);
      }
    }
  }
}

/* This is equivalent to nabla(S) */
void dS_dpara(int datalen, int num_of_parameters, double* X, double* Y, double* parameters, void* others, double* ret, void (*model)(int, double*, double*, double*, void*)) {
  double sum;
  double buf[datalen];
  double buf2[datalen];
  
  for (int i=0; i<num_of_parameters; i++) {
    sum = 0;
    get_errors(datalen, X, Y, parameters, buf, others, model);
    de_vec_dx_j(datalen, i, num_of_parameters, X, Y, parameters, buf2, others, model);
    for (int j=0; j<datalen; j++) {
      sum += (buf[j] * buf2[j]);
    }
    ret[i] = sum;
  }
}

bool check_if_converged(int len, double* para1, double* para2) {
  double buf[len];
  subtract_matrix(1, len, para1, para2, buf);
  return (DELTA_CONVERGE > norm_vector(len, buf));
}

void display_matrix(int l, int c, double* A) {
  for (int i=0; i<l; i++) {
    for (int j=0; j<c; j++) {
      cout << setw(7) << *(A + i*c + j) << " ";
    }
    cout << endl;
  }
}
