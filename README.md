# bench_density_gradient_wfn

Analyze the performance of evaluating the density gradient
with various language implementations (C++/Fortran/Rust/Julia).

```
Calculating the density gradient of CH₄ at a point (x10) [f]:
  17.897 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [rs]:
  17.283 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [cpp]:
  19.656 μs (2 allocations: 224 bytes)

Calculating the density gradient of CH₄ at a point (x10) [jl]:
  20.021 μs (122 allocations: 14.44 KiB)


Calculating the density gradient of C₂H₄ at a point (x10) [f]:
  26.353 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [rs]:
  39.857 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [cpp]:
  45.770 μs (2 allocations: 224 bytes)

Calculating the density gradient of C₂H₄ at a point (x10) [jl]:
  41.737 μs (142 allocations: 17.41 KiB)


Calculating the density gradient of imidazol at a point (x10) [f]:
  66.975 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [rs]:
  109.139 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [cpp]:
  123.619 μs (2 allocations: 224 bytes)

Calculating the density gradient of imidazol at a point (x10) [jl]:
  112.843 μs (202 allocations: 27.41 KiB)

```


Impressions
===========

1. **Fortran** is the fastest. Perfect for this task. The routine is simple
   enough to not experience any of its drawbacks. Good for numerical kernels.
   No dependencies. No extra effort is needed to handle arrays. Compiles in the
   blink of an eye. Declaration of variables is archaic. Picked this.
   Optimization steps:
   1. Predeclare variables (forced to do so).
   2. Flags(same as for others): `-fPIC -O3 -march=native -funroll-loops`

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
   Optimization steps:
   1. Predeclare variables.
   2. Flags: `RUSTFLAGS="-C target-cpu=native" cargo b --release`.
      Release: full LTO, 1 codegen unit, and abort on panic.
   3. Unchecked indexing with `uget` and `uget_mut`.

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
    Optimization steps:
    1. Predeclare variables.
    2. Flags(same as for others): `-O3 -march=native -fPIC -funroll-loops`

4. **Julia** is a joy. Interactive, expressive, and quite fast with some modifications.
    Scales better than Rust & C++. There is still a (small but crucial) difference with Fortran.
    If speed equals this stays.
    Optimization steps:
    1. StaticArrays.
    2. Predeclare variables to avoid type instability accessing struct elements.
    3. Flags(same as for others): `-O3 -C skylake --check-bounds=no`
    4. Avoid array copies with `@views`.
    5. Guarantee `@inbounds` access.


DONE: declare variables at the beginning in all implementations
      to see if they catch up with Fortran.
      Already made Rust faster than C++. Differences between 
      compiled languages should be due to implementation differences.

Compilers:
```
GNU Fortran (Homebrew GCC 10.2.0) 10.2.0

rustc 1.52.1 (9bc8c42bb 2021-05-09)

g++
Apple clang version 11.0.3 (clang-1103.0.32.62)
Target: x86_64-apple-darwin19.4.0
Thread model: posix


Julia Version 1.6.0
Commit f9720dc2eb (2021-03-24 12:55 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin19.6.0)
  CPU: Intel(R) Core(TM) i5-7360U CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, skylake)
Environment:
  JULIA_NUM_THREADS = 1

```
