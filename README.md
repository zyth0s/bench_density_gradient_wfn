# bench_density_gradient_wfn

Analysis of the performance of evaluating the density gradient
employing various language implementations (C++/Fortran/Rust/Julia)

Requires a private package.


```
Calculating the density gradient of CH₄ at a point (x10) [f]:
  17.677 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [cpp]:
  19.745 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [rs]:
  31.248 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [jl]:
  64.193 μs (862 allocations: 95.38 KiB)


Calculating the density gradient of C₂H₄ at a point (x10) [f]:
  26.512 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [cpp]:
  45.543 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [rs]:
  79.007 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [jl]:
  141.856 μs (1662 allocations: 183.66 KiB)


Calculating the density gradient of imidazol at a point (x10) [f]:
  66.627 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [cpp]:
  123.802 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [rs]:
  184.032 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [jl]:
  315.007 μs (2872 allocations: 319.44 KiB)

```


Impressions
===========

1. **Fortran** is the fastest. Perfect for this task. The routine is simple enough to not
          experience any of its drawbacks. No dependencies. No extra effort is needed to
          handle arrays. Declaration of variables is archaic. Picked this.

3. **Rust** is 2x slower than Fortran. Suffers with large array sizes. Would need to link to BLAS/LAPACK.
       It required also some additional work [Fortran < Rust < C++]
       because multidimensional arrays are not first
       order citizens. But it is trivial to add ndarray.
       More verbose for the same reason, also due
       to forced type conversions. Inference helps anyway.
       Too, nalgebra library has some very nice things (but <3D).
       Wish they join together. Might be still a risky alternative!?
       Provided some pre-training, and having to develop a large library...
       Even tried another allocator, why not? Saw no effect.
       Fortran compilers are more widespread and are faster after all.
       Will stay with Fortran for this reason (I do not want more issues)
       but it is worth keeping an eye on Rust, definetely.

2. **C++** is a bit worse than Rust [due to location of declarations?].
      Suffers with larger matrix sizes. It is tedious.
      Does C++ have multidimensional dynamic arrays? Yes,
      valarrays (or vectors if you wish but not for numerics).
      The former are poorly documented and have terrible slicing.
      Ok, are there 3rd party libs? Yes, Eigen[try this; header-only], Armadillo[known], ...
      Eigen, despite being more common, IMHO has worse documentation.
      Oh, Eigen only has 1D/2D. If you want >2D => unsupported but avail?!.
      Armadillo has up to 3D. Also, Eigen interface is not very uniform.
      How to set them up? Manually. Download and place them somewhere,
      set up the includes, adapt the build script. Buff... fine. Once it is
      done it is like Rust. Well, you cannot use inference (auto)
      with those libraries [messy]. And the compiler is more silent [bug land].
      Would require storing a copy of Eigen (I do not want more 3rd party code here).

4. **Julia** is a joy. Interactive, expressive, but
        not fast enough [although loable]. That is crucial
        for this routine. Looking for ways to speed it up.
        If speed equals this stays.


INIT: declare variables at the beginning in all implementations
      to see if they catch up with Fortran.
      Already made Rust faster than C++. Differences between 
      compiled languages should be due to implementation differences.

TODO: try with a large system

TODO: link to LAPACK with ndalgebra-linalg in Rust library.

TODO: link to LAPACK in the Fortran library.

