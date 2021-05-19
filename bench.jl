
# Benchmark

using Libdl: dlopen, dlsym, dlext
find_grad_fn(path_basename,extension) = dlsym(
               dlopen(joinpath(@__DIR__,
                      "$path_basename.$extension")),
               :density_gradient)

function density_gradient_rust(wfn; nsamples=10)

   point = ones(3)
   ngtoH = size(wfn.nuexp, 3)
   ∇ρ = zeros(3)

   for i in 1:nsamples
      @ccall $RUST_gradient(
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
                 ∇ρ         :: Ptr{Cdouble}, # inout
             )::Cvoid
   end
end

function density_gradient_cpp(wfn; nsamples=10)

   point = ones(3)
   ngtoH = size(wfn.nuexp, 3)
   ∇ρ = zeros(3)

   for i in 1:nsamples
      @ccall $CPP_gradient(
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
                 ∇ρ         :: Ptr{Cdouble}, # inout
             )::Cvoid
   end
end

function density_gradient_fortran(wfn; nsamples=10)

   point = ones(3)
   ngtoH = size(wfn.nuexp, 3)
   ∇ρ = zeros(3)

   for i in 1:nsamples
      @ccall $FORTRAN_gradient(
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

const RUST_gradient    = find_grad_fn("rust/target/release/libexp_chemint_rust", dlext)
const CPP_gradient     = find_grad_fn("cpp/lib", "so")
const FORTRAN_gradient = find_grad_fn("fortran/libeval_gradrho", "so")

struct AIMWFN_Mol_Gamess_HF
   nmo       :: Int64
   nocc      :: Int64
   nprims    :: Int64
   natm      :: Int64
   atnam     :: Vector{String}
   coords    :: Matrix{Float64}
   charges   :: Vector{Float64}
   icen      :: Vector{Int64}
   ityp      :: Vector{Int64}
   oexp      :: Vector{Float64}
   occ       :: Vector{Float64}
   eorb      :: Vector{Float64}
   coef      :: Matrix{Float64}
   tote      :: Float64
   virial    :: Float64
   npc       :: Vector{Int64}
   icenat    :: Matrix{Int64}
   nlm       :: Matrix{Int64}
   ngroup    :: Vector{Int64}
   nzexp     :: Matrix{Int64}
   nuexp     :: Array{Int64,3}
   npcant    :: Int64
   maxgrp    :: Int64
   numshells :: Int64
   cuttz     :: Float64
   rcutte    :: Matrix{Float64}
   chkfile   :: String
end



using Serialization
using BenchmarkTools
#using ChemInt

@doc """
    bench_library(lang=:f)

Benchmark C++ (:cpp), Rust (:rs), Fortran(:f), or Julia (:jl)
"""
function bench_library(lang=:f,system=:ch4)

   wfnpath = joinpath(@__DIR__, "$(string(system)).wfn")
   wfn = if :ChemInt in names(Main, imported = true)
            parse_wavefunction(wfnpath)
         else
            deserialize(wfnpath * ".ser")
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
      if :ChemInt in names(Main, imported = true)
         return @btime density_gradient_julia($wfn)
      else
         println("  ChemInt not present: skipping Julia benchmark")
      end
   else
      @warn "lang must be one of :cpp, :rs, :f, or :jl"
   end
end

for system in [:ch4, :c2h4, :imidazol]
   for lang in [:f,:rs,:cpp,:jl]
      bench_library(lang, system)
   end
   println()
end

