"""
Author: Satoki Tsujino
Date: 2026/03/04
"""
module wrap_netcdf  # sub module for handling NetCDF
# For Himawari-8/9 satellite, wrapper library
using NCDatasets  # Reading NetCDF
using DataStructures: OrderedDict  # Creating NetCDF attributes

export ncdump_parallax

#-- Define function to create NetCDF for VTTrac
"""
   ncdump_vttrac(fname_nc,ref_fname_nc,n_grid,lon_grid,lat_grid,gph_grid,u_grid,v_grid,score_grid,_fillvalue)

Create a NetCDF file and store the VTTrac results

# Arguments
- `fname_nc::String` : Created NetCDF file name
- `ref_fname_nc::String` : Input NetCDF file name (used to copy global attributes and obs_time)
- `lon_grid::Array{Real,1}` : [nx_grid] longitude [degree]
- `lat_grid::Array{Real,1}` : [ny_grid] latitude [degree]
- `gph_grid::Array{Real,2}` : [nx_grid,ny_grid] height [m] at each lon_grid/lat_grid
- `val_grid::Array{Real,2}` : [nx_grid,ny_grid] variable [unit] at each lon_grid/lat_grid
- `read_vname::String` : varname to be output
"""
function ncdump_parallax(fname_nc::String,ref_fname_nc::String,lon_grid::Array{Real,1},lat_grid::Array{Real,1},gph_grid::Array{Real,2},val_grid::Array{Real,2},read_vname::String)

    ds_new = NCDataset(fname_nc,"c")
    ds_old = NCDataset(ref_fname_nc,"r")

    dlon_grid = lon_grid[2] - lon_grid[1]
    dlat_grid = lat_grid[2] - lat_grid[1]

    # Define the dimensions.
    defDim(ds_new,"longitude",size(lon_grid[:])[1])
    defDim(ds_new,"latitude",size(lat_grid[:])[1])
    defDim(ds_new,"obs_time",1)

    # Define global attributes
    ds_new.attrib["title"] = "Parallax-corrected satellite image"
    ds_new.attrib["source"] = "test"
    ds_new.attrib["institution"] = "Meteorological Research Institute"
    ds_new.attrib["satellite_name"] = ds_old[read_vname].attrib["satellite_name"]
    ds_new.attrib["used_band"] = ds_old[read_vname].attrib["band_number"]
    ds_new.attrib["equator_radius"] = ds_old[read_vname].attrib["equator_radius"]
    ds_new.attrib["polar_radius"] = ds_old[read_vname].attrib["polar_radius"]
    ds_new.attrib["satellite_longitude"] = ds_old[read_vname].attrib["satellite_longitude"]
    ds_new.attrib["satellite_latitude"] = ds_old[read_vname].attrib["satellite_latitude"]
    ds_new.attrib["satellite_earth_center_distance"] = ds_old[read_vname].attrib["satellite_earth_center_distance"]
    ds_new.attrib["nadir_longitude"] = ds_old[read_vname].attrib["nadir_longitude"]
    ds_new.attrib["nadir_latitude"] = ds_old[read_vname].attrib["nadir_latitude"]
    ds_new.attrib["wavelength"] = ds_old[read_vname].attrib["wavelength"]
    _fillvalue = ds_old[read_vname].attrib["_FillValue"]

    # Define the variables with the attribute units
    #println(typeof(ds_old["obs_time"].attrib["units"]))
    #println(typeof(ds_old["obs_time"].attrib["long_name"]))
    obs_time = 0.5 .* (ds_old["start_time"].var .+ ds_old["end_time"].var)
    tim = defVar(ds_new,"obs_time",Float64,("obs_time",), attrib = OrderedDict(
            "units" => ds_old["start_time"].attrib["units"], 
            "long_name" => ds_old["start_time"].attrib["long_name"]
            ))
    lon = defVar(ds_new,"longitude",Float32,("longitude",), attrib = OrderedDict(
            "units" => "degree_east",
            "long_name" => "longitude",
            ))
    lat = defVar(ds_new,"latitude",Float32,("latitude",), attrib = OrderedDict(
            "units" => "degree_north",
            "long_name" => "latitude",
            ))
    gph = defVar(ds_new,"gph",Float32,("longitude","latitude"), attrib = OrderedDict(
            "units" => "m",
            "long_name" => "Cloud top height",
            "_FillValue" => Float32(_fillvalue),
            ))
    val = defVar(ds_new,read_vname,Float32,("longitude","latitude"), attrib = OrderedDict(
            "units" => ds_old[read_vname].attrib["units"],
            "long_name" => ds_old[read_vname].attrib["long_name"],
            "_FillValue" => Float32(_fillvalue),
            ))

    # add additional attributes
    #v.attrib["comments"] = "this is a string attribute with Unicode Ω ∈ ∑ ∫ f(x) dx"

    # Generate some example data
    tim[:] = obs_time[:]
    lon[:] = lon_grid
    lat[:] = lat_grid
    gph[:,:] = gph_grid
    val[:,:] = val_grid

    close(ds_new)
end

end
