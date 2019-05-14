#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#pragma once

using namespace std;

double NaN();
bool isNaN(double);
bool isNum(string);
bool check_nan_matrix(int, int, double*);
void copy_arr(int, double*, double*);
int get_datanum(string, int, int);
int read_file(string, int, int, double*, double*);
vector<string> split(string);
void concat_data(int, int*, double**, double*);
double average_arr(int, double*);
void add_const_array(double, int, double*, double*);
double abs(double);
double norm_vector(int, double*);
void waitEnterKey();


double NaN() {
  return numeric_limits<double>::quiet_NaN();
}

bool isNaN(double num1) {
  double num2 = num1;
  return (num1 != num2);
}

bool isNum(string str) {
  try {
    double num = stod(str);
    return 1;
  } catch (...) {
    return 0;
  }
}

bool check_nan_matrix(int l, int c, double* X) {
  for (int i=0; i<l; i++) {
    for (int j=0; j<l; j++) {
      if (isNaN(*(X + i*c + j))) {return false;}
    }
  }
  return true;
}

void copy_arr(int len, double* source, double* target) {
  while (--len >= 0) target[len] = source[len];
}

int get_datanum(string filename, int c1, int c2) {
  ifstream ifs;
  ifs.open(filename);
  string str;
  int count = 0;
  
  if (ifs.fail()) {
    cerr << "Cannot open file. filename=" << filename << endl;
    exit(0);
  }

  vector<string> strs;
  double x, y;
  while (getline(ifs, str)) {
    strs = split(str);
    try {
      x = stod(strs[c1-1]);
      y = stod(strs[c2-1]);
      count++;
    } catch(...) {
      ;
    }
  }

  return count;
}

int read_file(string filename, int c1, int c2, double* X, double* Y) {
  ifstream ifs;
  ifs.open(filename);
  string str;
  int count = 0;
  
  if (ifs.fail()) {
    cerr << "Cannot open file. filename=" << filename << endl;
    exit(0);
  }

  vector<string> strs;
  double x, y;
  while (getline(ifs, str)) {
    strs = split(str);
    try {
      x = stod(strs[c1-1]);
      y = stod(strs[c2-1]);
      X[count] = x;
      Y[count++] = y;
    } catch(...) {
      ;
    }
  }

  return count;
}

vector<string> split(string str) {
  int count = 0;
  vector<string> strings;
  stringstream ss(str);
  string tmp;

  while (!ss.eof()) {
    ss >> tmp;
    strings.push_back(tmp);
  }
  return strings;
}

void concat_data(int num_of_data, int* datalens, double** data, double* ret) {

  double* ptr;
  double* current = ret;
  for (int i=0; i<num_of_data; i++) {
    ptr = data[i];
    copy_arr(datalens[i], ptr, current);
    current += datalens[i];
  }
}

double average_arr(int len, double* arr) {
  double sum = 0;
  for (int i=0; i<len; i++) sum += arr[i];
  return sum / len;
}

void add_const_array(double C, int len, double* arr, double* ret) {
  while (--len >= 0) ret[len] = arr[len] + C;
}

double abs(double x) {
  return (x < 0 ? -x : x);
}

/* Euclidean norm */
double norm_vector(int len, double* vec) {
  double sum = 0;
  while (--len >= 0) sum += (vec[len] * vec[len]);
  return sum;
}

void waitEnterKey() {
  if (cin.get()) return;
}
