#include "ma27Interface.hpp"
#include <algorithm>
#include <limits>

#define F77_FUNC(name,NAME) name ## _

extern "C"
{
  void F77_FUNC(ma27id,MA27ID)(int* ICNTL, double* CNTL);
  void F77_FUNC(ma27ad,MA27AD)(int *N, int *NZ, const int *IRN, const int* ICN,
                               int *IW, int* LIW, int* IKEEP, int *IW1,
                               int* NSTEPS, int* IFLAG, int* ICNTL,
                               double* CNTL, int *INFO, double* OPS);
  void F77_FUNC(ma27bd,MA27BD)(int *N, int *NZ, const int *IRN, const int* ICN,
                               double* A, int* LA, int* IW, int* LIW,
                               int* IKEEP, int* NSTEPS, int* MAXFRT,
                               int* IW1, int* ICNTL, double* CNTL,
                               int* INFO);
  void F77_FUNC(ma27cd,MA27CD)(int *N, double* A, int* LA, int* IW,
                               int* LIW, double* W, int* MAXFRT,
                               double* RHS, int* IW1, int* NSTEPS,
                               int* ICNTL, double* CNTL);
}

MA27_LinearSolver::MA27_LinearSolver(double pivtol,
                                     double n_a_factor,
                                     double n_iw_factor,
                                     double memory_increase_factor)
  :
  nnz_(0),
  dim_(0),
  n_iw_(0),
  iw_(NULL),
  irows_(NULL),
  jcols_(NULL),
  n_a_(0),
  a_(NULL),
  ikeep_(NULL),
  nsteps_(0),
  maxfrt_(0),
  n_a_factor_(n_a_factor),
  n_iw_factor_(n_iw_factor),
  memory_increase_factor_(memory_increase_factor),
  num_neg_evals_(-1)
{
  // initialize the options to their defaults
  F77_FUNC(ma27id,MA27ID)(icntl_, cntl_);
  cntl_[0] = pivtol;
}

MA27_LinearSolver::~MA27_LinearSolver()
{
  delete [] iw_;
  delete [] irows_;
  delete [] jcols_;
  delete [] a_;
  delete [] ikeep_;

  //  numeric_fac_timer_.PrintTotal(std::cout, "numeric factor time");
  //  backsolve_timer_.PrintTotal(std::cout, "backsolve time");
}

MA27_LinearSolver::MA27_FACT_STATUS MA27_LinearSolver::DoSymbolicFactorization(int nrowcols, int* irow, int* jcol, int nnonzeros)
{

  MA27_LinearSolver::MA27_FACT_STATUS status = MA27_SUCCESS;

  delete [] iw_;
  iw_ = NULL;
  n_iw_ = 0;

  delete [] irows_;
  irows_ = NULL;

  delete [] jcols_;
  jcols_ = NULL;
  nnz_ = 0;

  delete [] a_;
  a_ = NULL;
  n_a_ = 0;

  delete [] ikeep_;
  ikeep_ = NULL;

  dim_ = nrowcols;
  nnz_ = nnonzeros;

  irows_ = new int[nnz_];
  jcols_ = new int[nnz_];

  // copy data
  std::copy(irow, irow + nnonzeros, irows_);
  std::copy(jcol, jcol + nnonzeros, jcols_);

  // allocate work space memory
  int overestimate_factor_n_iw = 2; // overestimate size by 2
  n_iw_ = overestimate_factor_n_iw * (2*nnz_+3*dim_+1);
  iw_ = new int[n_iw_];

  ikeep_ = new int[3*dim_];

  int N = dim_;
  int NZ = nnz_;
  int IFLAG = 0;
  double OPS;
  int INFO[20];
  int* IW1 = new int[2*dim_];

  // find sequence of pivots
  F77_FUNC(ma27ad,MA27AD)(&N,
                          &NZ,
                          irows_,
                          jcols_,
                          iw_,
                          &n_iw_,
                          ikeep_,
                          IW1,
                          &nsteps_,
	                        &IFLAG,
                          icntl_,
                          cntl_,
                          INFO,
                          &OPS);
  delete [] IW1;

  // check info pointer
  int retflag = INFO[0];
  int errorflag = INFO[1];
  if (retflag != 0) {
    std::cerr << "Error in ma27ad_ with INFO[0]="
              <<retflag <<" and INFO[1]=" << errorflag << std::endl;
    status = MA27_ERROR;
  }

  // allocate memory for numerical factorization
  int recommended_n_iw = INFO[5];

  int max_n_iw = static_cast<int>(std::numeric_limits<int>::max() / n_iw_factor_);
  if (recommended_n_iw <= max_n_iw){
    n_iw_ = static_cast<int>(n_iw_factor_ * recommended_n_iw);
  }
  else{
    n_iw_ = static_cast<int>(recommended_n_iw);
  }

  // check numeric_limits
  if (n_iw_ < 0 || n_iw_ > std::numeric_limits<int>::max())
  {
    std::cerr << "IW array needs too much memory. Quitting" << '\n';
    status = MA27_ERROR;
  }

  delete [] iw_;
  iw_ = NULL;
  iw_ = new int[n_iw_];

  int recommended_n_a = INFO[4];
  int max_n_a = static_cast<int>(std::numeric_limits<int>::max() / n_a_factor_);
  double factor = 1.0;
  if (recommended_n_a <= max_n_a){
    factor = n_a_factor_;
  }
  n_a_ = std::max(nnz_, static_cast<int>(factor * recommended_n_a));

  // check numeric_limits
  if (n_a_ < 0 || n_a_ > std::numeric_limits<int>::max())
  {
    std::cerr << "A array needs too much memory. Quitting" << '\n';
    status = MA27_ERROR;
  }

  delete [] a_;
  a_ = NULL;
  a_ = new double[n_a_];

  return status;
}

MA27_LinearSolver::MA27_FACT_STATUS MA27_LinearSolver::DoNumericFactorization(int nrowcols, int nnonzeros, double* values, int desired_num_neg_evals)
{

  _ASSERT_(nrowcols == dim_);
  _ASSERT_(nnonzeros == nnz_);

  MA27_LinearSolver::MA27_FACT_STATUS status = MA27_SUCCESS;
  int N = dim_;
  int NZ = nnz_;
  int INFO[20];

  // reset number of negative eigenvalues
  num_neg_evals_ = - 1;

  // copy data
  std::copy(values, values + nnonzeros, a_);

  int count_iters = 0;
  int max_iters = 20;
  int fact_error = 1;
  while(fact_error > 0 && count_iters<max_iters){
    int* IW1 = new int[2*dim_];
    F77_FUNC(ma27bd,MA27BD)(&N,
                            &NZ,
                            irows_,
                            jcols_,
                            a_,
                            &n_a_,
                            iw_,
                            &n_iw_,
                            ikeep_,
                            &nsteps_,
  	                        &maxfrt_,
                            IW1,
                            icntl_,
                            cntl_,
                            INFO);
    delete [] IW1;

    int retflag = INFO[0];
    int errorflag = INFO[1];

    if (retflag == 0) {
      fact_error = 0;
    }

    double recomended_value = errorflag;
    if (retflag == -3) {
      delete [] iw_;
      iw_ = NULL;
      delete [] a_;
      a_ = NULL;
      n_iw_ = static_cast<int>(memory_increase_factor_ * recomended_value);
      n_a_ = static_cast<int>(memory_increase_factor_ * n_a_);

      // check numeric_limits
      if (n_iw_ < 0 || n_iw_ > std::numeric_limits<int>::max())
      {
        std::cerr << "IW array needs too much memory. Quitting" << '\n';
        return MA27_ERROR;
      }
      if (n_a_ < 0 || n_a_ > std::numeric_limits<int>::max())
      {
        std::cerr << "A array needs too much memory. Quitting" << '\n';
        return MA27_ERROR;
      }

      std::cout << "Reallocating memory for MA27: liw " << n_iw_ <<"\n";
      std::cout << "Reallocating memory for MA27: la " << n_a_ <<"\n";
      iw_ = new int[n_iw_];
      a_ = new double[n_a_];

    }
    else if (retflag == -4) {
      delete [] iw_;
      iw_ = NULL;
      delete [] a_;
      a_ = NULL;
      n_iw_ = static_cast<int>(memory_increase_factor_ * n_iw_);
      n_a_ = static_cast<int>(memory_increase_factor_ * recomended_value);

      // check numeric_limits
      if (n_iw_ < 0 || n_iw_ > std::numeric_limits<int>::max())
      {
        std::cerr << "IW array needs too much memory. Quitting" << '\n';
        return MA27_ERROR;
      }
      if (n_a_ < 0 || n_a_ > std::numeric_limits<int>::max())
      {
        std::cerr << "A array needs too much memory. Quitting" << '\n';
        return MA27_ERROR;
      }

      std::cout << "Reallocating memory for MA27: liw " << n_iw_ <<"\n";
      std::cout << "Reallocating memory for MA27: la " << n_a_ <<"\n";
      iw_ = new int[n_iw_];
      a_ = new double[n_a_];
    }
    else if (retflag == -5) {
      return MA27_MATRIX_SINGULAR;
    }
    else if (retflag == 3) {
      return MA27_MATRIX_SINGULAR;
    }
    else if (retflag < 0) {
      std::cerr << "Error in ma27ad_ with INFO[0]="
                <<retflag <<" and INFO[1]=" << errorflag << std::endl;
      return MA27_ERROR;
    }
    else if (retflag > 0) {
      std::cerr << "Error in ma27ad_ with INFO[0]="
                <<retflag <<" and INFO[1]=" << errorflag << std::endl;
      return MA27_WARNING;
    }
    ++count_iters;
  }

  if(count_iters >= max_iters)
  {
    std::cerr << "Reallocated memory "<<count_iters<< "times. Quitting"<< '\n';
    return MA27_ERROR;
  }

  int num_int_compressions = INFO[12];
  int num_double_compressions = INFO[11];
  if (num_int_compressions >= 10) {
    std::cerr << "MA27: Number of integer compressions is high - increase n_iw_" << std::endl;
    return MA27_ERROR;
  }

  if (num_double_compressions >= 10) {
    std::cerr << "MA27: Number of double compressions is high - increase n_a_" << std::endl;
    return MA27_ERROR;
  }

  int num_neg_evals = INFO[14];
  if (desired_num_neg_evals != -1 && desired_num_neg_evals != num_neg_evals) {
    if (status != MA27_MATRIX_SINGULAR)
    {
        status = MA27_INCORRECT_INERTIA;
        num_neg_evals_ = num_neg_evals;
    }
  }
  return status;
}

void MA27_LinearSolver::DoBacksolve(double* rhs, int nrhs, double* sol, int nsol)
{

  _ASSERT_(nrhs == dim_);
  _ASSERT_(nsol == dim_);

  int N = dim_;
  double* W = new double[maxfrt_];
  int* IW1 = new int[nsteps_];

  double* soln_vals = new double[N];
  std::copy(rhs, rhs + N, soln_vals);

  F77_FUNC(ma27cd,MA27CD)(&N, a_, &n_a_, iw_, &n_iw_, W, &maxfrt_, soln_vals, IW1, &nsteps_, icntl_, cntl_);

  std::copy(soln_vals, soln_vals + N, sol);

  delete [] soln_vals;
  delete [] W;
  delete [] IW1;

}

extern "C"
{

MA27_LinearSolver* EXTERNAL_MA27Interface_new(double pivot,
                                              double n_a_factor,
                                              double n_iw_factor,
                                              double memory_increase_factor)
{
  return new MA27_LinearSolver(pivot,
                               n_a_factor,
                               n_iw_factor,
                               memory_increase_factor);
}

int EXTERNAL_MA27Interface_get_nnz(MA27_LinearSolver* p_hi)
{ return p_hi->get_nnz();}

int EXTERNAL_MA27Interface_get_n_iw(MA27_LinearSolver* p_hi)
{ return p_hi->get_n_iw();}

int EXTERNAL_MA27Interface_get_n_a(MA27_LinearSolver* p_hi)
{ return p_hi->get_n_a();}

int EXTERNAL_MA27Interface_get_dim(MA27_LinearSolver* p_hi)
{ return p_hi->get_dim();}

int EXTERNAL_MA27Interface_get_num_neg_evals(MA27_LinearSolver* p_hi)
{ return p_hi->get_num_neg_evals();}

int EXTERNAL_MA27Interface_do_symbolic_factorization(MA27_LinearSolver* p_hi,
                                                      int nrowcols,
                                                      int* irow,
                                                      int* jcol,
                                                      int nnonzeros)
{
  int status = p_hi->DoSymbolicFactorization(nrowcols, irow, jcol, nnonzeros);
  return status;
}

int EXTERNAL_MA27Interface_do_numeric_factorization(MA27_LinearSolver* p_hi,
                                                     int nrowcols,
                                                     int nnonzeros,
                                                     double* values,
                                                     int desired_num_neg_evals) {
  int status = p_hi->DoNumericFactorization(nrowcols, nnonzeros, values, desired_num_neg_evals);
  return status;
}

void EXTERNAL_MA27Interface_do_backsolve(MA27_LinearSolver* p_hi, double* rhs, int nrhs, double* sol, int nsol)
{p_hi->DoBacksolve(rhs, nrhs, sol, nsol);}

void EXTERNAL_MA27Interface_free_memory(MA27_LinearSolver* p_hi)
{p_hi->~MA27_LinearSolver();}

}
