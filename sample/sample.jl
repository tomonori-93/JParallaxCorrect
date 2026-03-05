include("../src/JParallaxCorrect.jl")
include("../src/wrap_netcdf.jl")

# For Himawari-8/9 satellite, a sample script of parallax correction
using NCDatasets  # Reading NetCDF
using DelimitedFiles
using DataStructures: OrderedDict  # Creating NetCDF attributes
using .JParallaxCorrect  # Parallax correction
using .wrap_netcdf

#-- setting section
infile_list = "input_files.txt"  # 1st column: NetCDF, 2nd column: Temperature vs height profile data (ASCII)
hskip_zprof = 4  # array number to skip reading in the temp/height profile data
tord = 2  # Column order of temperature in infile_list[2]
zord = 1  # Column order of height in infile_list[2]
read_vname = "tbb"
height_vname = "gph"
lon_vname = "longitude"
lat_vname = "latitude"
_fillval_name = "_FillValue"
rename = "equator_radius"
rpname = "polar_radius"
hsatname = "satellite_earth_center_distance"
lsatname = "satellite_longitude"
psatname = "satellite_latitude"
dlon = 0.02  # Horizontal grid spacing [degree] for longitude
dlat = 0.02  # Horizontal grid spacing [degree] for latitude
#-----------------------

d2r = π/180.0
r2d = 180.0/π

# Read infile_list
infiles = readlines(infile_list)
nl = size(infiles)[1]

# Start loop to read and perform parallax correction
for i in 1:nl
    local infile_nc, infile_zp = split(infiles[i])

    # input temperature/height profiles
    temperature_1d = Float64.(readdlm(infile_zp)[hskip_zprof+1:end,tord])
    height_1d = Float64.(readdlm(infile_zp)[hskip_zprof+1:end,zord])

    # Input NetCDF satellite data
    local al = NCDataset(infile_nc,"r")
    println("Read : $infile_nc")

    # Set arrays and parameters for reading Himawari data
    R_e = al[read_vname].attrib[rename] * 1.0e3  # km -> m
    R_p = al[read_vname].attrib[rpname] * 1.0e3  # km -> m
    h_sat = al[read_vname].attrib[hsatname] * 1.0e3  # km -> m
    l_sat = al[read_vname].attrib[lsatname] * d2r  # degree -> radian
    p_sat = al[read_vname].attrib[psatname] * d2r  # degree -> radian
    # In NCDatasets, _FillValue is forced to missing type in JuliaLang, so if set Float, you need it as follows.
    missing_value = al[read_vname].attrib[_fillval_name]
    ncol, nlin = size(al[read_vname].var)
    ds_val = fill(missing_value,ncol,nlin)  # NOTE: x,y
    ds_gph = fill(missing_value,ncol,nlin)  # NOTE: x,y
    ds_lon = fill(convert(Float64,missing_value),ncol,nlin)
    ds_lat = fill(convert(Float64,missing_value),ncol,nlin)

    # Substitute NetCDF data into ds_{val,lon,lat}
    ds_val = map(x -> convert(Float64,x), al[read_vname].var)
    ds_lon = map(x -> convert(Float64,x), al[lon_vname].var)
    ds_lat = map(x -> convert(Float64,x), al[lat_vname].var)

    # Allocate horizontal grid for output file (lon_grid/lat_grid)
    reshape_lon = reshape(ds_lon,(ncol * nlin,))  # temporary
    reshape_lat = reshape(ds_lat,(ncol * nlin,))  # temporary
    lonmin = minimum(reshape_lon)
    lonmax = maximum(reshape_lon)
    latmin = minimum(reshape_lat)
    latmax = maximum(reshape_lat)
    nlon = floor((lonmax - lonmin) / dlon) + 1
    nlat = floor((latmax - latmin) / dlat) + 1
    lon_grid = lonmin .+ dlon .* (Vector(1:nlon) .- 1)
    lat_grid = latmin .+ dlat .* (Vector(1:nlat) .- 1)

    # Convert Tbb to Geopotential height (Zph)
    ds_gph = convert_Tbb2Zph(ds_val, temperature_1d, height_1d, missing_value)  # assuming ds_val = tbb

    # Perform parallax correction and assign the corrected values to lon_grid/lat_grid
    # lon/lat: [radian], gph: [m]
    gph_grid, val_grid = ParallaxCorrect(ds_lon.*d2r, ds_lat.*d2r, ds_gph, ds_val,
                            lon_grid.*d2r, lat_grid.*d2r, R_e, R_p, h_sat, p_sat, l_sat, missing_value)

    lon_grid = Array{Real}(lon_grid)
    lat_grid = Array{Real}(lat_grid)
    gph_grid = Array{Real}(gph_grid)
    val_grid = Array{Real}(val_grid)

    outfile = chopsuffix(infile_nc,".nc") * ".para.nc"
    ncdump_parallax(outfile,String(infile_nc),lon_grid,lat_grid,gph_grid,val_grid,read_vname)

    println("Output file: $outfile ...")

end

