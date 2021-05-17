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

