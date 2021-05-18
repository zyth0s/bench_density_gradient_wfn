
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
//static GLOBAL: jemallocator::Jemalloc = jemallocator::Jemalloc;

extern crate ndarray;
use ndarray::prelude::*;

#[no_mangle]
pub extern fn density_gradient(
                        point_ptr: *const f64,
                        nmo    : i64,
                        natm   : i64,
                        nprims : i64,
                        mgrp   : i64,
                        ngto_h : i64,
                        ngroup_ptr : *const i64,
                        ityp_ptr   : *const i64,
                        nzexp_ptr  : *const i64,
                        nlm_ptr    : *const i64,
                        nuexp_ptr  : *const i64,
                        occ_ptr    : *const f64,
                        oexp_ptr   : *const f64,
                        xyz_ptr    : *const f64,
                        rcutte_ptr : *const f64,
                        coef_ptr   : *const f64,
                        grad_ptr   : *mut   f64
                        ) {
    // input & output variables
    let point    = unsafe { ArrayView::from_shape_ptr(3, point_ptr) };
    let ngroup   = unsafe { ArrayView::from_shape_ptr(natm as usize, ngroup_ptr) };
    let ityp     = unsafe { ArrayView::from_shape_ptr(nprims as usize, ityp_ptr) };
    let nzexp    = unsafe { ArrayView::from_shape_ptr((natm as usize,mgrp as usize).f(), nzexp_ptr) };
    let nlm      = unsafe { ArrayView::from_shape_ptr((56,3).f(), nlm_ptr) };
    let nuexp    = unsafe { ArrayView::from_shape_ptr((natm as usize, mgrp as usize, ngto_h as usize).f(), nuexp_ptr) };
    let occ      = unsafe { ArrayView::from_shape_ptr(nmo as usize, occ_ptr) };
    let oexp     = unsafe { ArrayView::from_shape_ptr(nprims as usize, oexp_ptr) };
    let xyz      = unsafe { ArrayView::from_shape_ptr((natm as usize, 3).f(), xyz_ptr) };
    let rcutte   = unsafe { ArrayView::from_shape_ptr((natm as usize, mgrp as usize).f(), rcutte_ptr) };
    let coef     = unsafe { ArrayView::from_shape_ptr((2*nmo as usize, nprims as usize).f(), coef_ptr) };
    let mut grad = unsafe { ArrayViewMut::from_shape_ptr(3, grad_ptr) };
    // local variables
    let mut fun   = Array1::<f64>::zeros(3);
    let mut fun1  = Array1::<f64>::zeros(3);
    let mut xcoor = Array1::<f64>::zeros(3);
    let mut gun   = Array1::<f64>::zeros(nmo as usize);
    let mut gun1  = Array2::<f64>::zeros((nmo as usize, 3).f());
    let mut it    = Array1::<usize>::zeros(3);
    let mut k      : usize;
    let mut i      : usize;
    let mut itip   : usize;
    let mut n      : usize;
    let mut dis2   : f64;
    let mut ori    : f64;
    let mut dp2    : f64;
    let mut aexp   : f64;
    let mut x      : f64;
    let mut f12    : f64;
    let mut f123   : f64;
    let mut fa     : f64;
    let mut fb     : f64;
    let mut fc     : f64;
    let mut cfj    : f64;
    let mut fac    : f64;
    let mut facgun : f64;
    // output
    grad.fill(0.0);

    for ic in (0 as usize)..(natm as usize) {
        // Atomic coordinates of this center
        xcoor[0] = point[0] - xyz[[ic,0]];
        xcoor[1] = point[1] - xyz[[ic,1]];
        xcoor[2] = point[2] - xyz[[ic,2]];
        dis2 = xcoor.mapv(|a| a.powi(2)).sum();
        // Loop over different shells in this atom
        for m in (0 as usize)..(ngroup[ic] as usize) {
            k = (nuexp[[ic,m,0]]-1) as usize;
            // Skip to compute this primitive if distance is too big.
            if dis2 > (rcutte[[ic,m]]).powi(2) { continue };
            ori = -oexp[k];
            dp2 = 2.0*ori;
            // All primitives in a shell share the same exponent.
            aexp = (ori*dis2).exp();
            // Loop over the different primitives in this shell.
            for jj in (0 as usize)..(nzexp[[ic,m]] as usize) {
                // "i" is the original index of the primitive in the WFN.
                i = (nuexp[[ic,m,jj]]-1) as usize;
                itip = (ityp[i]-1) as usize;
                // Integer coeficients
                it[0] = nlm[[itip,0]] as usize;
                it[1] = nlm[[itip,1]] as usize;
                it[2] = nlm[[itip,2]] as usize;

                for j in (0 as usize)..(3 as usize) {
                    n = it[j] as usize;
                    x = xcoor[j];
                    if n == 0 {
                        fun1[j] = dp2 * x;
                        fun[j] = 1.0;
                    } else if n == 1 {
                        fun1[j] = 1.0 + dp2 * x * x;
                        fun[j] = x;
                    } else if n == 2 {
                        let x2 = x * x;
                        fun1[j] = x * (2.0 + dp2 * x2);
                        fun[j] = x2;
                    } else if n == 3 {
                        let x2 = x * x;
                        fun1[j] = x2 * (3.0 + dp2 * x2);
                        fun[j] = x * x2;
                    } else if n == 4 {
                        let x2 = x * x;
                        fun1[j] = x2 * x * (4.0 + dp2 * x2);
                        fun[j] = x2 * x2;
                    } else if n == 5 {
                        let x2 = x * x;
                        fun1[j] = x2 * x2 * (5.0 + dp2 * x2);
                        fun[j] = x2 * x2 * x;
                    }
                }
                f12 = fun[0] * fun[1] * aexp;
                f123 = f12 * fun[2];
                fa = fun1[0] * fun[1] * fun[2] * aexp;
                fb = fun1[1] * fun[0] * fun[2] * aexp;
                fc = fun1[2] * f12;

                // Run over orbitals
                for j in (0 as usize)..(nmo as usize) {
                    cfj = coef[[j,i]];
                    gun[j] += cfj * f123;
                    gun1[[j,0]] += cfj * fa;
                    gun1[[j,1]] += cfj * fb;
                    gun1[[j,2]] += cfj * fc;
                }
            }
        }
    }
    // Run again over orbitals
    for i in (0 as usize)..(nmo as usize) {
        fac = occ[i];
        facgun = fac * gun[i];
        for j in (0 as usize)..(3 as usize) {
            grad[j] += facgun * gun1[[i,j]];
        }
    }
    grad[0] *= 2.0;
    grad[1] *= 2.0;
    grad[2] *= 2.0;

}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_density_gradient() { // Dumb test
        let point = array![0.0,0.0,0.0];
        let nmo: i64    = 5;
        let natm: i64   = 5;
        let nprims: i64 = 93;
        let mgrp: i64   = 19;
        let ngto_h: i64 = 21;
        let ngroup = Array1::<i64>::zeros(natm as usize);
        let ityp   = Array1::<i64>::zeros(nprims as usize);
        let nzexp  = Array2::<i64>::zeros((natm as usize, mgrp as usize).f());
        let nlm    = Array2::<i64>::zeros((56, 3).f());
        let nuexp  = Array3::<i64>::zeros((natm as usize, mgrp as usize, ngto_h as usize).f());
        let occ    = Array1::<f64>::zeros(nmo as usize);
        let oexp   = Array1::<f64>::zeros(nprims as usize);
        let xyz    = Array2::<f64>::zeros((natm as usize,3).f());
        let rcutte = Array2::<f64>::zeros((natm as usize,mgrp as usize).f());
        let coef   = Array2::<f64>::zeros((2*nmo as usize,nprims as usize).f());
        let mut grad = Array1::<f64>::zeros(3);

        density_gradient(point.as_ptr(),
                         nmo,
                         natm,
                         nprims,
                         mgrp,
                         ngto_h,
                         ngroup.as_ptr(),
                         ityp.as_ptr(),
                         nzexp.as_ptr(),
                         nlm.as_ptr(),
                         nuexp.as_ptr(),
                         occ.as_ptr(),
                         oexp.as_ptr(),
                         xyz.as_ptr(),
                         rcutte.as_ptr(),
                         coef.as_ptr(),
                         grad.as_mut_ptr()
                         );
    }
}
