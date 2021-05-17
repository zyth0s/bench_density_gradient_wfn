
# Benchmark

const RUSTLIB = joinpath(@__DIR__, "rust/target/release/libexp_chemint_rust.dylib")
const CPPLIB = joinpath(@__DIR__, "cpp/lib.so")
const FORTRANLIB = joinpath(@__DIR__, "fortran/libeval_gradrho.so")

function density_gradient_rust(wfn; nsamples=10)

   point = ones(3)
   ngtoH = size(wfn.nuexp, 3)
   ∇ρ = zeros(3)

   for i in 1:nsamples
      @ccall RUSTLIB.density_gradient(
                 point      :: Ptr{Cdouble},
                 wfn.nmo    :: Clong,
                 wfn.natm   :: Clong,
                 wfn.nprims :: Clong,  
                 wfn.maxgrp :: Clong,
                 ngtoH      :: Clong,
                 wfn.ngroup :: Ptr{Clong},
                 wfn.ityp   :: Ptr{Clong},
                 wfn.nzexp  :: Ptr{Clong},
                 wfn.nlm    :: Ptr{Clong},
                 wfn.nuexp  :: Ptr{Clong},
                 wfn.occ    :: Ptr{Cdouble},
                 wfn.oexp   :: Ptr{Cdouble},
                 wfn.coords :: Ptr{Cdouble},
                 wfn.rcutte :: Ptr{Cdouble},
                 wfn.coef   :: Ptr{Cdouble},
                 #ρ          :: Ref{Cdouble}, # inout
                 ∇ρ         :: Ptr{Cdouble}, # inout
             )::Cvoid
   end
end

function density_gradient_cpp(wfn; nsamples=10)

   point = ones(3)
   ngtoH = size(wfn.nuexp, 3)
   ∇ρ = zeros(3)

   for i in 1:nsamples
      @ccall CPPLIB.density_gradient(
                 point      :: Ptr{Cdouble},
                 wfn.nmo    :: Clong,
                 wfn.natm   :: Clong,
                 wfn.nprims :: Clong,  
                 wfn.maxgrp :: Clong,
                 ngtoH      :: Clong,
                 wfn.ngroup :: Ptr{Clong},
                 wfn.ityp   :: Ptr{Clong},
                 wfn.nzexp  :: Ptr{Clong},
                 wfn.nlm    :: Ptr{Clong},
                 wfn.nuexp  :: Ptr{Clong},
                 wfn.occ    :: Ptr{Cdouble},
                 wfn.oexp   :: Ptr{Cdouble},
                 wfn.coords :: Ptr{Cdouble},
                 wfn.rcutte :: Ptr{Cdouble},
                 wfn.coef   :: Ptr{Cdouble},
                 #ρ          :: Ref{Cdouble}, # inout
                 ∇ρ         :: Ptr{Cdouble}, # inout
             )::Cvoid
   end
end

function density_gradient_fortran(wfn; nsamples=10)

   point = ones(3)
   ngtoH = size(wfn.nuexp, 3)
   ∇ρ = zeros(3)

   for i in 1:nsamples
      @ccall FORTRANLIB.density_gradient(
                 point      :: Ptr{Cdouble},
                 wfn.nmo    :: Ref{Clong},
                 wfn.natm   :: Ref{Clong},
                 wfn.nprims :: Ref{Clong},  
                 wfn.maxgrp :: Ref{Clong},
                 ngtoH      :: Ref{Clong},
                 wfn.ngroup :: Ptr{Clong},
                 wfn.ityp   :: Ptr{Clong},
                 wfn.nzexp  :: Ptr{Clong},
                 wfn.nlm    :: Ptr{Clong},
                 wfn.nuexp  :: Ptr{Clong},
                 wfn.occ    :: Ptr{Cdouble},
                 wfn.oexp   :: Ptr{Cdouble},
                 wfn.coords :: Ptr{Cdouble},
                 wfn.rcutte :: Ptr{Cdouble},
                 wfn.coef   :: Ptr{Cdouble},
                 #ρ          :: Ref{Cdouble}, # inout
                 ∇ρ         :: Ptr{Cdouble}, # inout
             )::Cvoid
   end
end

# Pure Julia
function density_gradient_julia(wfn; nsamples=10)

   point = ones(3)
   ngtoH = size(wfn.nuexp, 3)
   ∇ρ = zeros(3)

   for i in 1:nsamples
      ChemInt.jl_density_gradient!(∇ρ, point, wfn, 0)
   end
   nothing
end

# ------------------------------------------------
#
# Prepare

cd("fortran")
run(`make`)
cd("../cpp")
run(`make`)
cd("../rust")
run(`make`)
cd("..")


using ChemInt
using BenchmarkTools

@doc """
    bench_library(lang=:f)

Benchmark C++ (:cpp), Rust (:rs), Fortran(:f), or Julia (:jl)
"""
function bench_library(lang=:f,system=:ch4)

   wfn = if system == :ch4
            parse_wavefunction(joinpath(@__DIR__, "ch4.wfn"))
         elseif system == :c2h4
            parse_wavefunction(joinpath(@__DIR__, "c2h4.wfn"))
         elseif system == :imidazol
            parse_wavefunction(joinpath(@__DIR__, "imidazol.wfn"))
         end

   formula = if system == :ch4
            "CH₄"
         elseif system == :c2h4
            "C₂H₄"
         elseif system == :imidazol
            "imidazol"
         end

   println()
   println("Calculating the density gradient of $formula at a point (x10) [$(string(lang))]:")
   if lang == :cpp
      return @btime density_gradient_cpp($wfn)
   elseif lang == :rs
      return @btime density_gradient_rust($wfn)
   elseif lang == :f
      return @btime density_gradient_fortran($wfn)
   elseif lang == :jl
      return @btime density_gradient_julia($wfn)
   else
      @warn "lang must be one of :cpp, :rs, :f, or :jl"
   end
end

for system in [:ch4, :c2h4, :imidazol]
   for lang in [:f,:cpp,:rs,:jl]
      bench_library(lang, system)
   end
   println()
end

