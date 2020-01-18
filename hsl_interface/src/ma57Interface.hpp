#ifndef __MA57INTERFACE_HPP__
#define __MA57INTERFACE_HPP__

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

class MA57_LinearSolver
{
public:

  MA57_LinearSolver(double pivtol = 0.01,
                    double prealocate_factor = 1.05);

  ~MA57_LinearSolver();

  enum MA57_FACT_STATUS
    {
      MA57_SUCCESS=0,
      MA57_MATRIX_SINGULAR=1,
      MA57_INCORRECT_INERTIA=2,
      MA57_WARNING=3,
      MA57_ERROR=4
    };

  MA57_FACT_STATUS DoSymbolicFactorization(int nrowcols, int* irow, int* jcol, int nnonzeros);

  MA57_FACT_STATUS DoNumericFactorization(int nrowcols, int nnonzeros, double* values, int desired_num_neg_evals=-1);

  void DoBacksolve(double* rhs, int nrhs, double* sol, int nsol);

  int get_nnz(){return nnz_;}
  int get_dim(){return dim_;}
  int get_lifact(){return lifact_;}
  int get_lfact(){return lfact_;}
  int get_n_a(){return n_a_;}
  int get_num_neg_evals(){return num_neg_evals_;}

private:
  // some members for error checking between calls
  int nnz_;
  int dim_;

  // double values for options
  double cntl_[5];
  // integer values for options
  int icntl_[20];

  // Factor for estimating initial size of work arrays
  double ma57_pre_alloc_;

  int info_[40];
  double rinfo_[20];

  // MA57 parameters
  int lkeep_;      /* LKEEP >= 5*N + NE + max(N,NE) + 42. */
  int  *keep_;

  // integer workspace
  int  *iwork_;      /* 5 * N. */

  int   lfact_;
  double  *fact_;
  int   lifact_;
  int  *ifact_;

  // extracted matrix structure
  int* irows_;
  int* jcols_;

  // factors
  int n_a_;
  double* a_;

  int num_neg_evals_;
};

#endif
