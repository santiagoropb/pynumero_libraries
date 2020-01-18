#ifndef __MA27INTERFACE_HPP__
#define __MA27INTERFACE_HPP__

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cassert>


#define _ASSERT_ assert
#define _ASSERT_MSG_ assert_msg
#define _ASSERT_EXIT_ assert_exit
#define _ASSERTION_FAILURE_ assertion_failure

inline void assert_msg(bool cond, const std::string& msg) {
//void assert_msg(bool cond, const std::string& msg) {
  if (!cond) {
    std::cout << "Assertion Failed: " << msg.c_str() << std::endl;
  }
  assert(msg.c_str() && cond);
}

inline void assert_exit(bool cond, const std::string& msg, int exitcode=1) {
//void assert_exit(bool cond, const std::string& msg, int exitcode=1) {
  if (!(cond)) {
    std::cout << msg << std::endl;
    exit(exitcode);
  }
}

inline void assertion_failure(const std::string& msg) {
  std::cout << "Assertion Failed: " << msg.c_str() << std::endl;
  assert(msg.c_str() && false);
}

class MA27_LinearSolver
{
public:

  MA27_LinearSolver(double pivtol = 0.01,
                    double n_a_factor = 5.0,
                    double n_iw_factor = 5.0,
                    double memory_increase_factor = 2.0);

  ~MA27_LinearSolver();

  enum MA27_FACT_STATUS
    {
      MA27_SUCCESS=0,
      MA27_MATRIX_SINGULAR=1,
      MA27_INCORRECT_INERTIA=2,
      MA27_WARNING=3,
      MA27_ERROR=4
    };

  MA27_FACT_STATUS DoSymbolicFactorization(int nrowcols,
                                           int* irow,
                                           int* jcol,
                                           int nnonzeros);

  MA27_FACT_STATUS DoNumericFactorization(int nrowcols,
                                          int nnonzeros,
                                          double* values,
                                          int desired_num_neg_evals=-1);

  void DoBacksolve(double* rhs, int nrhs, double* sol, int nsol);

  int get_nnz(){return nnz_;}
  int get_dim(){return dim_;}
  int get_n_iw(){return n_iw_;}
  int get_n_a(){return n_a_;}
  int get_num_neg_evals(){return num_neg_evals_;}

private:
  // some members for error checking between calls
  int nnz_;
  int dim_;

  // integer values for options
  int icntl_[30];
  // double values for options
  double cntl_[5];

  // integer workspace
  int n_iw_;
  int* iw_;

  // extracted matrix structure
  int* irows_;
  int* jcols_;

  // factors
  int n_a_;
  double* a_;

  double n_a_factor_;
  double n_iw_factor_;
  double memory_increase_factor_;

  // MA27 parameters
  int* ikeep_;
  int nsteps_;
  int maxfrt_;

  int num_neg_evals_;

};

#endif
