"""
Author: Satoki Tsujino
Date: 2026/03/04
"""
module JParallaxCorrect

const lib = joinpath(@__DIR__, "..", "lib", "libparallax.so")

export ParallaxCorrect, parallax_correct_core, tri_interpolation_2d, convert_Tbb2Zph

"""
"""
function ParallaxCorrect(lon_pix::Matrix{Float64},
                         lat_pix::Matrix{Float64}, 
                         h_pix::Matrix{Float64},
                         val_pix::Matrix{Float64},
                         lon_grid::Vector{Float64},
                         lat_grid::Vector{Float64},
                         re::Float64,
                         rp::Float64,
                         hsat::Float64,
                         psat::Float64,
                         lsat::Float64,
                         missing_value::Float64)

    n = size(lon_grid)[1]
    m = size(lat_grid)[1]

    h_grid = Matrix{Float64}(undef,n,m)
    val_grid = Matrix{Float64}(undef,n,m)

    #-- Performing the parallax correction from prescribed lon/lat/height at each pixel 
    #--   lon_cor, lat_cor: lon/lat after the parallax correction at each pixel
    lon_cor, lat_cor = parallax_correct_core(lon_pix, lat_pix, h_pix,
                          re, rp, hsat, psat, lsat, missing_value)

    #-- Assigning the height/brightness temperature/other variables to prescribed lon/lat gridpoints, 
    #    based on triangle interpolation from the parallax-corrected lon/lat at each pixel
    #--   lon_cor, lat_cor: lon/lat after the parallax correction at each pixel
    h_grid, val_grid = tri_interpolation_2d(lon_cor, lat_cor, h_pix, val_pix,
                           lon_grid, lat_grid, missing_value)

    return h_grid, val_grid

end

#-- 以下は c_interface からのラッパー関数

function parallax_correct_core(lon_pix::Matrix{Float64},
                               lat_pix::Matrix{Float64}, 
                               h_pix::Matrix{Float64},
                               re::Float64,
                               rp::Float64,
                               hsat::Float64,
                               psat::Float64,
                               lsat::Float64,
                               missing_value::Float64)

    n, m = size(lon_pix)

    @assert stride(lon_pix,1) == 1   # column-major保証

    lon_cor = Matrix{Float64}(undef,n,m)
    lat_cor = Matrix{Float64}(undef,n,m)

    ccall((:c_parallax_correct, lib),
          Cvoid,
          (Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Float64, Float64, Float64, 
           Float64, Float64, Float64),
          n, m, lon_pix, lat_pix, h_pix, lon_cor, lat_cor,
          re, rp, hsat, psat, lsat, missing_value)

    return lon_cor, lat_cor

end

function tri_interpolation_2d(x_in::Matrix{Float64},
                              y_in::Matrix{Float64}, 
                              iv::Matrix{Float64},
                              ivad::Matrix{Float64},
                              x_out::Vector{Float64},
                              y_out::Vector{Float64},
                              missing_value::Float64)

    n, m = size(x_in)
    l = size(x_out)[1]
    k = size(y_out)[1]

    @assert stride(x_in,1) == 1   # column-major保証

    ov = Matrix{Float64}(undef,l,k)
    ovad = Matrix{Float64}(undef,l,k)

    ccall((:c_tri_interpolation_2d, lib),
          Cvoid,
          (Cint, Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Float64),
          n, m, l, k, x_in, y_in, iv, ivad, x_out, y_out, ov, ovad, missing_value)

    return ov, ovad

end

function convert_Tbb2Zph(tval::Matrix{Float64},
                         t1d::Vector{Float64},
                         z1d::Vector{Float64},
                         missing_value::Float64)

    n, m = size(tval)
    l = size(t1d)[1]

    @assert stride(tval,1) == 1   # column-major保証

    zval = Matrix{Float64}(undef,n,m)

    ccall((:c_convert_Tbb2Zph, lib),
          Cvoid,
          (Cint, Cint, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
           Float64),
          n, m, l, tval, zval, t1d, z1d, missing_value)

    return zval

end

end
