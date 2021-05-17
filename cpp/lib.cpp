
#include <iostream>
#include <valarray> // No nice indexing.
// We load Eigen just to have nice indexing & printing!
// Download Eigen & place Eigen/ -> ./Eigen/
//                        unsupported/ -> Eigen_unsupported/
#include <Eigen/Dense>
#include <Eigen_unsupported/Eigen/CXX11/Tensor>

using namespace std;
using namespace Eigen;

typedef TensorMap<Tensor<int64_t, 3>> TensorMap3d;
typedef Matrix<int64_t, 3, 1> Vector3i64;
typedef Matrix<int64_t, Dynamic, 1> VectorXi64;
typedef Matrix<int64_t, Dynamic, Dynamic> MatrixXi64;

extern "C" {
  int density_gradient(
                       double* point_ptr,
                       const int64_t nmo,
                       const int64_t natm,
                       const int64_t nprims,
                       const int64_t mgrp,
                       const int64_t ngto_h,
                       int64_t* ngroup_ptr,
                       int64_t* ityp_ptr,
                       int64_t* nzexp_ptr,
                       int64_t* nlm_ptr,
                       int64_t* nuexp_ptr,
                       double* occ_ptr,
                       double* oexp_ptr,
                       double* xyz_ptr,
                       double* rcutte_ptr,
                       double* coef_ptr,
                       double* grad_ptr
                      ) ;
}

int density_gradient(
                     double* point_ptr,
                     const int64_t nmo,
                     const int64_t natm,
                     const int64_t nprims,
                     const int64_t mgrp,
                     const int64_t ngto_h,
                     int64_t* ngroup_ptr,
                     int64_t* ityp_ptr,
                     int64_t* nzexp_ptr,
                     int64_t* nlm_ptr,
                     int64_t* nuexp_ptr,
                     double* occ_ptr,
                     double* oexp_ptr,
                     double* xyz_ptr,
                     double* rcutte_ptr,
                     double* coef_ptr,
                     double* grad_ptr
										) {

  // input & output variables
  Map<Vector3d>   point(point_ptr, 3);
  Map<VectorXi64> ngroup(ngroup_ptr, natm);
  Map<VectorXi64> ityp(ityp_ptr, nprims);
  Map<MatrixXi64> nzexp(nzexp_ptr, natm, mgrp);
  Map<MatrixXi64> nlm(nlm_ptr, 56, 3); // could be static ... => () -> []
  TensorMap3d     nuexp(nuexp_ptr, natm, mgrp, ngto_h);
  Map<VectorXd>   occ(occ_ptr, nmo);
  Map<VectorXd>   oexp(oexp_ptr, nprims);
  Map<MatrixXd>   xyz(xyz_ptr, natm, 3);
  Map<MatrixXd>   rcutte(rcutte_ptr, natm, mgrp);
  Map<MatrixXd>   coef(coef_ptr, 2*nmo, nprims);
  Map<Vector3d>   grad(grad_ptr, 3);
  // local variables
  Vector3d fun   = Vector3d::Zero(3);
  Vector3d fun1  = Vector3d::Zero(3);
  Vector3d xcoor = Vector3d::Zero(3);
  VectorXd gun   = VectorXd::Zero(nmo);
  MatrixXd gun1  = MatrixXd::Zero(nmo,3);

  grad.setZero();

  // Run over centers
  for ( int ic=0; ic<natm; ic++) {
    // Atomic coordinates of this center
    xcoor(0) = point(0) - xyz(ic,0);
    xcoor(1) = point(1) - xyz(ic,1);
    xcoor(2) = point(2) - xyz(ic,2);
    double dis2 = xcoor.cwiseAbs2().sum();
    // Loop over different shell in this atom
    for ( int m=0; m<ngroup(ic); m++) {
      int k = nuexp(ic, m, 0)-1;
      // Skip to compute this primitive if distance is too big.
      if (dis2 > rcutte(ic, m) * rcutte(ic, m) ) {
        continue;
      }
      double ori = -oexp(k);
      double dp2 = 2.0*ori;
      // All primitives in a shell share the same exponent.
      double aexp = exp( ori * dis2 );
      // Loop over the different primitives in this shell.
      for ( int jj=0; jj<nzexp(ic,m); jj++) {
        // "i" is the original index of the primitive in the WFN.
        int i = nuexp(ic, m, jj)-1;
        int itip = ityp(i)-1;
        // Integer coefficients.
        Vector3i64 it = Vector3i64::Zero();
        it(0) = nlm(itip,0);
        it(1) = nlm(itip,1);
        it(2) = nlm(itip,2);
        double x2;

        for ( int j=0; j<3; j++) {
          int n = it(j);
          double x = xcoor(j);
          if (n == 0) {
                  fun1(j) = dp2 * x;
                  fun(j) = 1.0;
          } else if (n == 1) {
                  fun1(j) = 1.0 + dp2 * x * x;
                  fun(j) = x;
          } else if (n == 2) {
                  x2 = x * x;
                  fun1(j) = x * ( 2.0 + dp2 * x2 );
                  fun(j) = x2;
          } else if (n == 3) {
                  x2 = x * x;
                  fun1(j) = x2 * ( 3.0 + dp2 * x2 );
                  fun(j) = x * x2;
          } else if (n == 4) {
                  x2 = x * x;
                  fun1(j) = x2 * x * ( 4.0 + dp2 * x2 );
                  fun(j) = x2 * x2;
          } else if (n == 5) {
                  x2 = x * x;
                  fun1(j) = x2 * x2 * ( 5.0 + dp2 * x2 );
                  fun(j) = x2 * x2 * x;
          }
        } // endfor j
        double f12 = fun(0) * fun(1) * aexp;
        double f123 = f12 * fun(2);
        double fa = fun1(0) * fun(1) * fun(2) * aexp;
        double fb = fun1(1) * fun(0) * fun(2) * aexp;
        double fc = fun1(2) * f12;

        // run over orbitals
        for ( int j=0; j<nmo; j++) {
          double cfj = coef(j,i);
          gun(j) = gun(j) + cfj * f123;
          gun1(j,0) = gun1(j, 0) + cfj * fa;
          gun1(j,1) = gun1(j, 1) + cfj * fb;
          gun1(j,2) = gun1(j, 2) + cfj * fc;
        }

      } // endfor jj
    } // endfor m
  } // endfor ic

  // Run again over orbitals
  for ( int i=0; i<nmo; i++) {
    double fac = occ(i);
    double facgun = fac * gun(i);
    for ( int j=0; j<3; j++) {
      grad(j) += facgun * gun1(i,j);
    }
  }
  grad *= 2.0;

  return EXIT_SUCCESS;
}
