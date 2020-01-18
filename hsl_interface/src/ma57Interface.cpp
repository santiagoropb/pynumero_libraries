#include "ma57Interface.hpp"

#include <cstdio>
#include <limits>

#define F77_FUNC(name,NAME) name ## _

extern "C"
{
  void F77_FUNC(ma57id,MA57ID)(double* CNTL, int* ICNTL);
  void F77_FUNC(ma57ad,MA57AD)(int *N, int *NZ, const int *IRN, const int* ICN,
			       int *LKEEP, int* KEEP, int* IWORK, int *ICNTL,
			       int* INFO, double* RINFO);
  void F77_FUNC(ma57bd,MA57BD)(int *N, int *NZ, double* A, double* FACT, int* LFACT,
			       int* IFACT, int* LIFACT, int* LKEEP, int* KEEP, int* IWORK,
			       int* ICNTL, double* CNTL, int* INFO, double* RINFO);
  void F77_FUNC(ma57cd,MA57CD)(int* JOB, int *N, double* FACT, int* LFACT,
			       int* IFACT, int* LIFACT, int* NRHS, double* RHS,
			       int* LRHS, double* WORK, int* LWORK, int* IWORK,
			       int* ICNTL, int* INFO);
  void F77_FUNC(ma57ed,MA57ED)(int *N, int* IC, int* KEEP, double* FACT, int* LFACT,
			       double* NEWFAC, int* LNEW, int* IFACT, int* LIFACT,
			       int* NEWIFC, int* LINEW, int* INFO);
}

MA57_LinearSolver::MA57_LinearSolver(double pivtol, double prealocate_factor)
  :
  nnz_(0),
  dim_(0),
  ma57_pre_alloc_(prealocate_factor),
  lkeep_(0),
  keep_(NULL),
  iwork_(NULL),
  lfact_(0),
  fact_(NULL),
  lifact_(0),
  ifact_(NULL),
  irows_(NULL),
  jcols_(NULL),
  n_a_(0),
  a_(NULL),
  num_neg_evals_(-1)
{
  // initialize the options to their defaults
  ma57id_(cntl_, icntl_);
  cntl_[0] = pivtol;

  // The following are set based on ipopt options. Added by Jia Kang on Nov. 9 2011
  icntl_[1-1] = 0; /* Error stream */
  icntl_[2-1] = 0; /* Warning stream. */

  icntl_[4-1] = 0; /* Print statistics.  NOT Used. */
  icntl_[5-1] = 0; /* Print error. */

  icntl_[6-1] = 5; /* Pivoting order. */
  icntl_[7-1] = 1; /* Pivoting strategy. */

  icntl_[11-1] = 16; /* Block size used by Level 3 BLAS in MA57BD - should be a multiple of 8.  Default is 16. */
  icntl_[12-1] = 16; /* Two nodes of the assembly tree are merged only if both involve less than ICNTL(12) eliminations.  Default is 16. */

  icntl_[15-1] = 1; // ma57_automatic_scaling
  icntl_[16-1] = 0; /* If set to 1, small entries are removed and corresponding pivots are placed at the end of factorization.  May be useful for highly rank deficient matrices.  Default is 0. */


}

MA57_LinearSolver::~MA57_LinearSolver()
{
  delete [] iwork_;
  delete [] irows_;
  delete [] jcols_;
  delete [] a_;
  delete [] keep_;
  delete [] fact_;
  delete [] ifact_;
}


MA57_LinearSolver::MA57_FACT_STATUS MA57_LinearSolver::DoSymbolicFactorization(int nrowcols, int* irow, int* jcol, int nnonzeros)
{

  MA57_LinearSolver::MA57_FACT_STATUS status = MA57_SUCCESS;

  delete [] iwork_;
  iwork_ = NULL;

  delete [] irows_;
  irows_ = NULL;

  delete [] jcols_;
  jcols_ = NULL;
  nnz_ = 0;

  delete [] a_;
  a_ = NULL;
  n_a_ = 0;

  delete [] keep_;
  keep_ = NULL;
  lkeep_ = 0;

  delete [] fact_;
  fact_ = NULL;
  lfact_ = 0;

  delete [] ifact_;
  ifact_ = NULL;
  lifact_ = 0;

  dim_ = nrowcols;
  nnz_ = nnonzeros;

  irows_ = new int[nnz_];
  jcols_ = new int[nnz_];

  // copy data
  std::copy(irow, irow + nnonzeros, irows_);
  std::copy(jcol, jcol + nnonzeros, jcols_);

  int MAX_N_NZ = (dim_ < nnz_) ? nnz_ : dim_;
  lkeep_ = 5*dim_ + nnz_ + MAX_N_NZ + 42;

  iwork_ = new int[5*dim_];
  keep_  = new int[lkeep_];

  // Initialize to 0 as otherwise MA57ED can sometimes fail
  for (int k=0; k<lkeep_; k++) {
    keep_[k] = 0;
  }

  int N = dim_;
  int NZ = nnz_;

  ma57ad_(&N, &NZ, irows_, jcols_, &lkeep_, keep_, iwork_, icntl_, info_, rinfo_);

  int retflag = info_[0];
  if (retflag < 0) {
    std::cerr << "Error in ma57ad_ with INFO[0]="
              <<retflag << std::endl;
    status = MA57_ERROR;
  }

  if (retflag > 0) {
    std::cerr << "Warning in ma57ad_ with INFO[0]="
              <<retflag << std::endl;
    status = MA57_WARNING;
  }

  lfact_ = info_[8]*ma57_pre_alloc_;
  lifact_ = info_[9]*ma57_pre_alloc_;

  delete [] fact_;
  fact_ = NULL;
  delete [] ifact_;
  ifact_ = NULL;

  fact_  = new double[lfact_];
  ifact_ = new int[lifact_];

  n_a_ = nnz_;
  delete [] a_;
  a_ = NULL;
  a_ = new double[n_a_];

  return status;
}

MA57_LinearSolver::MA57_FACT_STATUS MA57_LinearSolver::DoNumericFactorization(int nrowcols, int nnonzeros, double* values, int desired_num_neg_evals)
{
  _ASSERT_(nrowcols == dim_);
  _ASSERT_(nnonzeros == nnz_);

  int fact_error = 1;

  MA57_LinearSolver::MA57_FACT_STATUS status = MA57_SUCCESS;
  int N = dim_;
  int NZ = nnz_;

  // reset number of negative eigenvalues
  num_neg_evals_ = - 1;

  // copy data
  std::copy(values, values + nnonzeros, a_);

  int count_iters = 0;
  int max_iters = 20;
  while (fact_error > 0 && count_iters<max_iters) {
    ma57bd_(&N, &NZ, a_, fact_, &lfact_, ifact_, &lifact_,
	    &lkeep_, keep_, iwork_, icntl_, cntl_, info_, rinfo_);


    int retflag = info_[0];
    //  std::cout << "retflag" << retflag << std::endl;

    if (retflag == 0) {
      fact_error = 0;
    }
    else if (retflag == -3) {
      // Failure due to insufficient REAL space on a call to MA57B/BD.
      // INFO(17) is set to a value that may suffice.  INFO(2) is set
      // to value of LFACT.  The user can allocate a larger array and
      // copy the contents of FACT into it using MA57E/ED, and recall
      // MA57B/BD.

      double *temp;
      int ic = 0;

      lfact_ = info_[16]*ma57_pre_alloc_;

      if( lfact_ > std::numeric_limits<int>::max() || lfact_ < 0){
          std::cerr << "lfact array needs too much memory. Quitting" << '\n';
          return MA57_ERROR;
      }

      temp = new double[lfact_];

      printf("Reallocating memory for MA57: lfact (%d)\n", lfact_);

      int idmy;
      ma57ed_(&N, &ic, keep_, fact_, &info_[1], temp, &lfact_, ifact_, &info_[1], &idmy, &lfact_, info_);

      delete [] fact_;
      fact_ = temp;
    }
    else if (retflag == -4) {
      // Failure due to insufficient INTEGER space on a call to
      // MA57B/BD.  INFO(18) is set to a value that may suffice.
      // INFO(2) is set to value of LIFACT.  The user can allocate a
      // larger array and copy the contents of IFACT into it using
      // MA57E/ED, and recall MA57B/BD.

      int *temp;
      int ic = 1;

      lifact_ = info_[17]*ma57_pre_alloc_;

      if( lifact_ > std::numeric_limits<int>::max() || lifact_ < 0){
          std::cerr << "lifact array needs too much memory. Quitting" << '\n';
          return MA57_ERROR;
      }

      temp = new int[lifact_];

      double ddmy;
      ma57ed_(&N, &ic, keep_, fact_, &info_[1], &ddmy, &lifact_, ifact_, &info_[1], temp, &lifact_, info_);

      delete [] ifact_;
      ifact_ = temp;

    }
    else if (retflag < 0) {
      std::cerr << "Error in ma57ad_ with INFO[0]="
                <<retflag << std::endl;
      return MA57_ERROR;
    }
    else if (retflag > 0) {
      printf("Warning in MA57BD:  %d\n", info_[0]);
      return MA57_WARNING;
    }
    else if (retflag == 4) {
      printf("System singular, rank = %d\n", info_[25-1]);
      return MA57_MATRIX_SINGULAR;
    }
    ++count_iters;
  }

  if(count_iters >= max_iters)
  {
    std::cerr << "Reallocated memory "<<count_iters<< "times. Quitting"<< '\n';
    return MA57_ERROR;
  }

  int num_neg_evals = info_[24-1];
  int rank = info_[25-1];
  if (desired_num_neg_evals != -1 && desired_num_neg_evals != num_neg_evals) {

     if (status != MA57_MATRIX_SINGULAR)
     {
         status = MA57_INCORRECT_INERTIA;
         num_neg_evals_ = num_neg_evals;
     }
  }
  return status;
}

void MA57_LinearSolver::DoBacksolve(double* rhs, int nrhs, double* sol, int nsol)
{
  _ASSERT_(nrhs == dim_);
  _ASSERT_(nsol == dim_);

  int N = dim_;
  int JOB = 1;

  int  NRHS = 1; // Ax=b, x is a vector
  int LRHS = N;

  int  lwork;
  double* work;

  lwork = N * NRHS;
  work = new double[lwork];

  double* soln_vals = new double[N];
  std::copy(rhs, rhs + N, soln_vals);

  ma57cd_(&JOB, &N, fact_, &lfact_, ifact_, &lifact_, &NRHS, soln_vals, &LRHS, work, &lwork, iwork_, icntl_, info_);

  if (info_[0] != 0){
    printf("Error in MA57CD:  %d.\n", info_[0]);
  }
  std::copy(soln_vals, soln_vals + N, sol);

  delete [] soln_vals;
  delete [] work;
}

extern "C"
{

MA57_LinearSolver* EXTERNAL_MA57Interface_new(double pivot, double prealocate_factor)
{ return new MA57_LinearSolver(pivot, prealocate_factor);}

int EXTERNAL_MA57Interface_get_nnz(MA57_LinearSolver* p_hi)
{ return p_hi->get_nnz();}

int EXTERNAL_MA57Interface_get_dim(MA57_LinearSolver* p_hi)
{ return p_hi->get_dim();}

int EXTERNAL_MA57Interface_get_lfact(MA57_LinearSolver* p_hi)
{ return p_hi->get_lfact();}

int EXTERNAL_MA57Interface_get_lifact(MA57_LinearSolver* p_hi)
{ return p_hi->get_lifact();}

int EXTERNAL_MA57Interface_get_n_a(MA57_LinearSolver* p_hi)
{ return p_hi->get_n_a();}

int EXTERNAL_MA57Interface_get_num_neg_evals(MA57_LinearSolver* p_hi)
{ return p_hi->get_num_neg_evals();}

int EXTERNAL_MA57Interface_do_symbolic_factorization(MA57_LinearSolver* p_hi,
                                                      int nrowcols,
                                                      int* irow,
                                                      int* jcol,
                                                      int nnonzeros)
{
  int status = p_hi->DoSymbolicFactorization(nrowcols, irow, jcol, nnonzeros);
  return status;
}

int EXTERNAL_MA57Interface_do_numeric_factorization(MA57_LinearSolver* p_hi,
                                                     int nrowcols,
                                                     int nnonzeros,
                                                     double* values,
                                                     int desired_num_neg_evals) {
  int status = p_hi->DoNumericFactorization(nrowcols, nnonzeros, values, desired_num_neg_evals);
  return status;
}

void EXTERNAL_MA57Interface_do_backsolve(MA57_LinearSolver* p_hi, double* rhs, int nrhs, double* sol, int nsol)
{p_hi->DoBacksolve(rhs, nrhs, sol, nsol);}

void EXTERNAL_MA57Interface_free_memory(MA57_LinearSolver* p_hi)
{p_hi->~MA57_LinearSolver();}

}
