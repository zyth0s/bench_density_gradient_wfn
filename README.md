# bench_density_gradient_wfn

Analyze the performance of evaluating the density gradient
with various language implementations (C++/Fortran/Rust/Julia).

```
Calculating the density gradient of CH₄ at a point (x10) [f]:
  17.897 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [rs]:
  18.343 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [cpp]:
  19.656 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [jl]:
  54.734 μs (862 allocations: 95.38 KiB)


Calculating the density gradient of C₂H₄ at a point (x10) [f]:
  26.353 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [rs]:
  41.708 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [cpp]:
  45.770 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [jl]:
  117.139 μs (1662 allocations: 183.66 KiB)


Calculating the density gradient of imidazol at a point (x10) [f]:
  66.975 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [rs]:
  115.989 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [cpp]:
  123.619 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [jl]:
  240.226 μs (2872 allocations: 319.44 KiB)

```


Impressions
===========

1. **Fortran** is the fastest. Perfect for this task. The routine is simple
   enough to not experience any of its drawbacks. Good for numerical kernels.
   No dependencies. No extra effort is needed to handle arrays. Compiles in the
   blink of an eye. Declaration of variables is archaic. Picked this.

3. **Rust** is 2x slower than Fortran. Suffers with large array sizes.
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

2. **C++** is a bit worse than Rust.
    Also suffers with larger matrix sizes. It is tedious.
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
    Would require storing a copy of Eigen (I do not want more 3rd party
    code in my repository).

4. **Julia** is a joy. Interactive, expressive, but
    not fast enough [although loable]. That is crucial
    for this routine. Looking for ways to speed it up.
    If speed equals this stays.


DONE: declare variables at the beginning in all implementations
      to see if they catch up with Fortran.
      Already made Rust faster than C++. Differences between 
      compiled languages should be due to implementation differences.

TODO: try with a large system
