
#include "lib.cpp"

int main()
{ // Dumb test
	double point[3] = {1.0, 1.0, 1.0 };
	int64_t nmo    = 5;
	int64_t natm   = 5;
	int64_t nprims = 93;
	int64_t mgrp   = 19;
	int64_t ngto_h = 21;
	std::valarray<int64_t> ngroup(natm);
	std::valarray<int64_t> ityp(nprims);
	std::valarray<int64_t> nzexp(natm * mgrp);
	int64_t nlm[56*3];
	std::valarray<int64_t> nuexp(natm * mgrp * ngto_h);
	std::valarray<double> occ(nmo);
	std::valarray<double> oexp(nprims);
	std::valarray<double> xyz(natm * 3);
	std::valarray<double> rcutte(natm * mgrp);
	std::valarray<double> coef(2*nmo * nprims);
	std::valarray<double> grad(3);

  //ngroup = zero<int64_t>(ngroup);
  //cout << ngroup[1] << endl;
  ngroup = ngroup.apply( [](int64_t x){ return x = 0; });
  ityp   = ityp.apply(   [](int64_t x){ return x = 0; });
  nzexp  = nzexp.apply(  [](int64_t x){ return x = 0; });
  nuexp  = nuexp.apply(  [](int64_t x){ return x = 0; });
  occ    = occ.apply(    [](double x){ return x  = 0.0; });
  oexp   = oexp.apply(   [](double x){ return x  = 0.0; });
  xyz    = xyz.apply(    [](double x){ return x  = 0.0; });
  rcutte = rcutte.apply( [](double x){ return x  = 0.0; });
  coef   = coef.apply(   [](double x){ return x  = 0.0; });
  grad   = grad.apply(   [](double x){ return x  = 0.0; });

	density_gradient(&point[0],
									 nmo,
									 natm,
									 nprims,
									 mgrp,
									 ngto_h,
									 &ngroup[0],
									 &ityp[0],
									 &nzexp[0],
									 &nlm[0],
									 &nuexp[0],
									 &occ[0],
									 &oexp[0],
									 &xyz[0],
									 &rcutte[0],
									 &coef[0],
									 &grad[0]
									 );

  cout << grad[0] << endl;
  cout << grad[1] << endl;
  cout << grad[2] << endl;
}
