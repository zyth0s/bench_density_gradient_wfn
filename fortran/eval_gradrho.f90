
module eval_density_gradient

contains

  subroutine density_gradient(        &
                              point,  &
                              nmo,    &
                              natm,   &
                              nprims, &
                              mgrp,   &
                              ngtoH,  &
                              ngroup, &
                              ityp,   &
                              nzexp,  &
                              nlm,    &
                              nuexp,  &
                              occ,    &
                              oexp,   &
                              xyz,    &
                              rcutte, &
                              coef,   &
                              grad_rho) &
                              bind(C, name='density_gradient')

     use, intrinsic :: iso_c_binding, only: c_long, c_double
     implicit none
     ! input variables:
     real(c_double),  dimension(3),               intent(in) :: point
     integer(c_long),                             intent(in) :: nmo, natm, nprims, mgrp, ngtoH
     integer(c_long), dimension(natm),            intent(in) :: ngroup
     integer(c_long), dimension(nprims),          intent(in) :: ityp
     integer(c_long), dimension(natm,mgrp),       intent(in) :: nzexp
     integer(c_long), dimension(56,3),            intent(in) :: nlm ! magic number
     integer(c_long), dimension(natm,mgrp,ngtoH), intent(in) :: nuexp
     real(c_double),  dimension(nmo),             intent(in) :: occ
     real(c_double),  dimension(nprims),          intent(in) :: oexp
     real(c_double),  dimension(natm,3),          intent(in) :: xyz
     real(c_double),  dimension(natm,mgrp),       intent(in) :: rcutte
     real(c_double),  dimension(2*nmo,nprims),    intent(in) :: coef
     ! output variables:
     real(c_double),  dimension(3),               intent(inout) :: grad_rho
     ! local variables:
     integer(c_long)                   :: n, j, k, jj, i, itip, ic, m
     integer(c_long), dimension(3)     :: it
     real(c_double)                    :: ori, aexp, dp2, x, dis2, f12, f123, fa, fb, fc, fac, facgun, x2, cfj
     real(c_double),  dimension(3)     :: fun, fun1, xcoor
     real(c_double),  dimension(nmo)   :: gun
     real(c_double),  dimension(nmo,3) :: gun1

     grad_rho(:) = 0.0
     fun(:)      = 0.0
     fun1(:)     = 0.0
     gun(:)      = 0.0
     gun1(:,:)   = 0.0
     
     do ic = 1,natm
        ! Atomic coordinates of this center
        xcoor = point - xyz(ic,:)
        dis2 = dot_product(xcoor, xcoor)
        ! Loop over different shells in this atom
        do m = 1,ngroup(ic)
           k = nuexp(ic,m,1)
           ! Skip to compute this primitive if distance is too big.
           if (dis2 > rcutte(ic,m)**2) continue

           ori = -oexp(k)
           dp2 = 2*ori
           ! All primitives in a shell share the same exponent.
           aexp = exp(ori*dis2)
           ! Loop over the different primitives in this shell.
           do jj = 1,nzexp(ic,m)
              ! "i" is the original index of the primitive in the WFN.
              i = nuexp(ic,m,jj)
              itip = ityp(i)
              ! Integer coefficients.
              it = nlm(itip,:)
              do j = 1,3
                 n = it(j)
                 x = xcoor(j)
                 if (n == 0) then ! TODO select case (n); case (0); ... end select
                    fun1(j) = dp2 * x
                    fun(j) = 1.0
                 else if (n == 1) then
                    fun1(j) = 1 + dp2 * x * x
                    fun(j) = x
                 else if (n == 2) then
                    x2 = x * x
                    fun1(j) = x * ( 2.0 + dp2 * x2 )
                    fun(j) = x2
                 else if (n == 3) then
                    x2 = x * x
                    fun1(j) = x2 * ( 3.0 + dp2 * x2 )
                    fun(j) = x * x2
                 else if (n == 4) then
                    x2 = x * x
                    fun1(j) = x2 * x * ( 4.0 + dp2  * x2 )
                    fun(j) = x2 * x2
                 else if (n == 5) then
                    x2 = x * x
                    fun1(j) = x2 * x2 * ( 5.0 + dp2 * x2 )
                    fun(j) = x2 * x2 * x
                 end if
              enddo

              f12  = fun(1) * fun(2) * aexp
              f123 = f12 * fun(3)
              fa   = fun1(1) * fun(2) * fun(3) * aexp
              fb   = fun1(2) * fun(1) * fun(3) * aexp
              fc   = fun1(3) * f12

              ! run over orbitals
              do j = 1,nmo
                 cfj = coef(j,i)
                 gun(j)    = gun(j)    + cfj * f123
                 gun1(j,1) = gun1(j,1) + cfj * fa
                 gun1(j,2) = gun1(j,2) + cfj * fb
                 gun1(j,3) = gun1(j,3) + cfj * fc
              enddo
           enddo
        enddo
     enddo
     ! Run again over orbitals
     do i = 1,nmo
        fac = occ(i)
        facgun = fac * gun(i)
        do j = 1,3
           grad_rho(j) = grad_rho(j) + facgun * gun1(i,j)
        enddo
     enddo
     grad_rho(:) = 2*grad_rho(:)
  end subroutine

end module eval_density_gradient
