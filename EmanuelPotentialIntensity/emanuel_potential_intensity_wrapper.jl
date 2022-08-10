function get_pcmin(SST :: Float32, PSL :: Float32, P::Vector{Float32}, T::Vector{Float32}, R::Vector{Float32})
    @assert length(P) == length(T) == length(R)
    N = Ref{Int32}(length(P))
    NA = Ref{Int32}(length(P))
    PMIN = Ref{Float32}(0.0)
    VMAX = Ref{Float32}(0.0)
    IFLAG = Ref{Int32}(1)
    ccall((:__emanuelpotentialintensity_MOD_pcmin, joinpath(@__DIR__,"EmanuelPotentialIntensity.so")),
                    Cvoid,
                    (Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Int32}, Ref{Int32}, Ref{Float32}, Ref{Float32}, Ref{Int32}),
          Ref(SST), Ref(PSL), P, T, R, NA, N, PMIN, VMAX, IFLAG)
    if IFLAG[] == 0
        error("PI routine did not converge (hypercane)")
    elseif IFLAG[] == 2
        error("Cape routine did not converge")
    end
    return PMIN[], VMAX[]
end


function get_cape(tparcel :: Float32, rparcel :: Float32, pparcel:: Float32, T::Vector{Float32}, R::Vector{Float32}, P::Vector{Float32})
    @assert length(P) == length(T) == length(R)
    N = Ref{Int32}(length(P))
    ND = Ref{Int32}(length(P))
    CAPED = Ref{Float32}(0.0)
    SIG = Ref{Float32}(0.0)
    TOB = Ref{Float32}(0.0)
    IFLAG = Ref{Int32}(1)
    ccall((:__emanuelpotentialintensity_MOD_cape, joinpath(@__DIR__,"EmanuelPotentialIntensity.so")),
                    Cvoid,
                    (Ref{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}, Ptr{Float32}, Ref{Int32}, Ref{Int32},Ref{Float32}, Ref{Float32}, Ref{Float32}, Ref{Int32}),
          Ref(tparcel), Ref(rparcel),Ref(pparcel), T, R,P, ND, N,SIG, CAPED, TOB, IFLAG)
    if IFLAG[] == 0
        error("Cape routine did not run due to improper sounding (e.g. no water vapor at parcel level)")
    elseif IFLAG[] == 2
        error("Cape routine did not converge")
    end
    return CAPED[], TOB[]
end
