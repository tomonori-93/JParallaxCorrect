# How to use JParallaxCorrect
This part will be moved to a separately repository for the wrapper package of JuliaLang.

You can use the subroutines and functions in the library (`libparallax.so`) by using JuliaLang scripts which are wrapper functions for `libparallax.so`.

* You can find [example programs](./) in `sample/`
  * `sample.jl`: Read infrared brightness temperature from NetCDF data captured by the Himawari-8/9 satellite, convert the brightness temperature to geopotential height by `../data/bb2Zgph.dat`, and perform the parallax correction. The corrected results are output in a new NetCDF file named `[input_filename].para.nc`. 
    * Running: `$ julia sample.jl`

<!-- # Example programs -->
