# JParallaxCorrect
Fortran library to perform the parallax correction for images captured by geostationary meteorological satellites.

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17490992.svg)](https://doi.org/10.5281/zenodo.17490992) -->

# Features
* Calculation of the parallax based on mathematically exact formulations
  * Only approximation: The earth on which original images of satellites are projected is a rotational ellipsoid.
* The core library for [JParallaxCorrect.jl](https://github.com/tomonori-93/JParallaxCorrect.jl) which is a JuliaLang package to easily perform the parallax correction
<!--  * Supporting NetCDF: It is commonly used in meteorological and oceanographical communities. -->

# References
* Tsujino, S., Horinouchi, T., Tsukada, T., Kuo, H.-C., Yamada, H., & Tsuboki, K. (2021). Inner-core wind field in a concentric eyewall replacement of Typhoon Trami (2018): A quantitative analysis based on the Himawari-8 satellite. _Journal of Geophysical Research: Atmospheres_, **126**, e2020JD034434. https://doi.org/10.1029/2020JD034434
  * The original paper using the method.
* Tsujino, S., Tsuboki, K., Kanada, S., Hirano, S., Takahashi, T., Yamaguchi, M., Kato, M., Shimada, U., & Wada, A. (2026): Title. _SOLA_, **Volume**, pp-pp, [https://CHECK](https://CHECK)
  * The forumation was specified in Supporting Information [SI1](https://CHECKKKK)

* [Method descriptions](https://tomonori-93.github.io/juparallax/CHECKKK)

* [日本語による手法の定式化とその導出: Formulation and derivation (Japanese document)](https://www.gfd-dennou.org/library/davis/stpk/manual.pdf)

# Images
![Test Image 1](image/image1.png)



# USAGE
### [Install](install.md)
### [How to use and examples](tools/README.md)
### [Demo](demo/sample.md)

# Future works
* Validation with dual-Doppler analysis
* Development of objective methods to find the vortex center
* Application to various vortices and cases
* Dynamical analyses by using this method

<!--# Cite as
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17490992.svg)](https://doi.org/10.5281/zenodo.17490992)

The above DOI corresponds to the latest versioned release as published to Zenodo, where you will find all earlier releases. To cite ro-crate-py independent of version, use https://doi.org/10.5281/zenodo.17490992, which will always redirect to the latest release.
-->
