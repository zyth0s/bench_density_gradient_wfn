
function density_gradient!(∇ρ, point, wfn)


    ρ    = 0.0
   ∇ρ   .= 0.0
   ϕ, ∇ϕ = eval_mo1(wfn, point)

   occ = wfn.occ

   # Run again over orbitals
   @inbounds for imo in eachindex(occ)

      @views begin
         n    = occ[imo]
         nϕ   = n   *  ϕ[imo]
          ρ  += nϕ  *  ϕ[imo]
         ∇ρ .+= nϕ .* ∇ϕ[imo,:]
      end
   end
   ∇ρ .*= 2 # + c.c.
   ∇ρ, ρ
end


function eval_mo1(wfn, point) #derivative=1)

   nmo    = wfn.nmo
   natm   = wfn.natm
   coords = wfn.coords
   ngroup = wfn.ngroup
   nuexp  = wfn.nuexp
   rcutte = wfn.rcutte
   oexp   = wfn.oexp
   nzexp  = wfn.nzexp
   ityp   = wfn.ityp
   nlm    = wfn.nlm
   coef   = wfn.coef

    ϕ   = zeros(nmo)
   ∇ϕ   = zeros(nmo,3)
    G1D = @MVector zeros(3) # Gaussian 1D functions x, y, z
   ∇G1D = @MVector zeros(3)

   @views begin
      @inbounds for ic in 1:natm

         # Atomic coordinates of this center
         xcoor = point - coords[ic,:]
         r² = sum( abs2.(xcoor) ) # ≡ dot(xcoor, xcoor)
         # Loop over different shell in this atom
         @inbounds for m in 1:ngroup[ic]

            k = nuexp[ic,m,1]
            # Skip to compute this primitive if distance is too big.
            r² > rcutte[ic,m]^2 && continue

            α = -oexp[k]
            𝟐α = 2α
            # All primitives in a shell share the same exponent.
            aexp = exp(α * r²)
            # Loop over the different primitives in this shell.
            @inbounds for jj in 1:nzexp[ic,m]

               # "i" is the original index of the primitive in the WFN.
               i = nuexp[ic,m,jj]
               itip = ityp[i]
               # Integer coefficients.
               it = nlm[itip,:]
               @inbounds for j in 1:3

                  n = it[j]
                  # r⃗ = ∑ⱼ³ rⱼ e⃗ⱼ # cartesian components
                  rⱼ = xcoor[j]
                  # G1D = r^n exp(α * r²)
                  if n == 0
                     # ∇G1D = 2αr * G1D
                     ∇G1D[j] = 𝟐α * rⱼ
                      G1D[j] = 1.0
                  elseif n == 1
                     ∇G1D[j] = 1 + 𝟐α * rⱼ * rⱼ
                      G1D[j] = rⱼ
                  elseif n == 2
                     rⱼ²   = rⱼ * rⱼ
                     ∇G1D[j] = rⱼ * ( 2.0 + 𝟐α * rⱼ² )
                      G1D[j] = rⱼ²
                  elseif n == 3
                     rⱼ²   = rⱼ  * rⱼ
                     ∇G1D[j] = rⱼ² * ( 3.0 + 𝟐α * rⱼ² )
                      G1D[j] = rⱼ  * rⱼ²
                  elseif n == 4
                     rⱼ²   = rⱼ  * rⱼ
                     ∇G1D[j] = rⱼ² * rⱼ * ( 4.0 + 𝟐α  * rⱼ² )
                      G1D[j] = rⱼ² * rⱼ²
                  elseif n == 5
                     rⱼ²   = rⱼ  * rⱼ
                     ∇G1D[j] = rⱼ² * rⱼ² * ( 5.0 + 𝟐α * rⱼ² )
                      G1D[j] = rⱼ² * rⱼ² * rⱼ
                  end
               end

               # exp(-αr₁²) * exp(-αr₂²) * exp(-αr₃²) = exp(-αr²) ≡ aexp
                 G3D =  prod(G1D) * aexp # Gxyz = Gx * Gy * Gz
               ∇xG3D = ∇G1D[1] *  G1D[2] *  G1D[3] * aexp
               ∇yG3D =  G1D[1] * ∇G1D[2] *  G1D[3] * aexp
               ∇zG3D =  G1D[1] *  G1D[2] * ∇G1D[3] * aexp

               # AO -> MO basis; run over orbitals
               @inbounds for imo in 1:nmo

                  c = coef[imo,i]
                   ϕ[imo]   += c *   G3D
                  ∇ϕ[imo,1] += c * ∇xG3D
                  ∇ϕ[imo,2] += c * ∇yG3D
                  ∇ϕ[imo,3] += c * ∇zG3D
               end
            end
         end
      end
   end
   ϕ, ∇ϕ
end
