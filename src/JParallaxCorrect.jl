module JParallaxCorrect

const lib = joinpath(@__DIR__, "..", "lib", "libparallax.so")

export ParallaxCorrect, c_parallax_correct, c_tri_interpolation_2d

function ParallaxCorrect()
end

function J_parallax_correct(lon_cld::Matrix{Float64},
                            lat_cld::Matrix{Float64}, 
                            h_cld::Matrix{Float64},
                            re::Float64,
                            rp::Float64,
                            hsat::Float64,
                            psat::Float64,
                            lsat::Float64,
                            missing_value::Float64)

    n, m = size(lon_cld)

    @assert stride(lon_cld,1) == 1   # column-major保証

    lon_cor = Matrix{Float64}(undef,n,m)
    lat_cor = Matrix{Float64}(undef,n,m)

    ccall((:c_parallax_correct, lib),
          Cvoid,
          (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Float64, Float64, Float64, 
           Float64, Float64, Float64),
          n, m, lon_cld, lat_cld, h_cld, lon_cor, lat_cor,
          re, rp, hsat, psat, lsat, missing_value)

    return lon_cor, lat_cor

end

!-- 上までチェックした
function J_tri_interpolation_2d(x_in::Matrix{Float64},
                                y_in::Matrix{Float64}, 
                                h_cld::Matrix{Float64},
				iv::Matrix{Float64},
				ivad::Matrix{Float64},
				x_out::Vector{Float64},
				y_out::Vector{Float64},
				ov::Matrix{Float64},
				ovad::Matrix{Float64},
                                missing_value::Float64)

    n, m = size(lon_cld)

    @assert stride(lon_cld,1) == 1   # column-major保証

    lon_cor = Matrix{Float64}(undef,n,m)
    lat_cor = Matrix{Float64}(undef,n,m)

    ccall((:c_tri_interpolation_2d, lib),
          Cvoid,
          (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Float64, Float64, Float64, 
           Float64, Float64, Float64),
          n, m, lon_cld, lat_cld, h_cld, lon_cor, lat_cor,
          re, rp, hsat, psat, lsat, missing_value)

    return lon_cor, lat_cor

end

function mat_scale!(A::Matrix{Float64}, alpha::Float64)
    n, m = size(A)

    @assert stride(A,1) == 1   # column-major保証

    ccall((:mat_scale, lib),
          Cvoid,
          (Cint, Cint, Ptr{Float64}, Float64),
          n, m, A, alpha)

    return A
end

end
