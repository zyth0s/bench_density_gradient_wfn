# bench_density_gradient_wfn

Analysis of the performance of evaluating the density gradient
employing various language implementations (C++/Fortran/Rust/Julia)

Requires a private package.


```
Calculating the density gradient of CH₄ at a point x10 [cpp]:
  61.601 μs (102 allocations: 1.78 KiB)

Calculating the density gradient of CH₄ at a point x10 [f]:
  73.489 μs (202 allocations: 3.34 KiB)

Calculating the density gradient of CH₄ at a point x10 [rs]:
  75.102 μs (102 allocations: 1.78 KiB)

Calculating the density gradient of CH₄ at a point x10 [jl]:
  215.357 μs (6159 allocations: 113.61 KiB)
```


Feelings
========

* **C++** is the fastest [small margin] but also more tedious.
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
      Would require storing a copy of Eigen (I do not like adding 3rd party code).

* **Fortran** is fast. Perfect for this task. The routine is simple enough to not
          experience any of its drawbacks. No dependencies. Best choice.

* **Rust** is equally fast. It requires a bit more of work [Fortran < Rust < C++]
       because multidimensional arrays are not first
       order citizens. But it is trivial to add ndarray.
       More verbose for the same reason, also due
       to forced type conversions. Inference helps anyway.
       FWIW, nalgebra library has some very nice things (but <3D).
       Wish they join together. Rust might be still a risky alternative!?
       Provided some pre-training, and having to develop a large library...
       Even tried another allocator, why not? Saw no effect.
       Fortran compilers are more widespread and faster after all.
       Will stay with Fortran for this reason (no more issues)
       but it is worth keeping an eye on Rust, definetely.


* **Julia** is a joy. Interactive, expressive, but
        slightly slower (2-3x) [loable]. That is crucial
        for this routine. Looking for ways to speed it up.
        If speed equals... this stays.
