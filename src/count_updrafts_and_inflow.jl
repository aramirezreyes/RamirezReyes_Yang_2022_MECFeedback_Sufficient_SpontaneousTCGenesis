function days_to_indices(day_begin,day_end; total_days = 100, total_length=1200)
    indices_per_day = total_length ÷ total_days
    index_begin = (day_begin - 1)*indices_per_day + 1
    index_end = (day_end)*indices_per_day
    index_begin:index_end
end

function count_updrafts_and_downdrafts(radiusbins,vertical_windspeed :: Array{T,2},center;gridspacing = 1) where T
    sx,sy = size(vertical_windspeed)
    binspacing = radiusbins[2] - radiusbins[1]
    nupdrafts_per_distance = zeros(Int64,length(radiusbins) - 1)
    ndowndrafts_per_distance = zeros(Int64,length(radiusbins) - 1)
    npoints_in_distance = zeros(Int64,length(radiusbins) - 1)
    for index in CartesianIndices(vertical_windspeed)
        distance_from_center = AvailablePotentialEnergyFramework.distance(index,center,gridspacing)
        if distance_from_center >= 1000
            in_which_bin = ceil(Int, (distance_from_center - 1000) / binspacing) 
            isupdraft = vertical_windspeed[index] >= 2.0
            isdowndraft = vertical_windspeed[index] <= -2.0
            npoints_in_distance[in_which_bin] += 1
            if isupdraft 
                nupdrafts_per_distance[in_which_bin] +=1
            end
            if isdowndraft 
                ndowndrafts_per_distance[in_which_bin] +=1
            end
        end
    end
    return nupdrafts_per_distance, ndowndrafts_per_distance, npoints_in_distance
end

function add_updraft_and_downdrafts_all_cyclones_in_timestep(cyclones, radiusbins,vertical_windspeed :: Array{T,2},center;gridspacing = 1) where T
    count = 0
    nupdrafts_per_distance = zeros(Int64,length(radiusbins) - 1)
    ndowndrafts_per_distance = zeros(Int64,length(radiusbins) - 1)
    npoints_in_distance = zeros(Int64,length(radiusbins) - 1)
    for cyclone in cyclones[1]
        if cyclone[2] != 1000
            cyclone_center = (cyclone[1][1], cyclone[1][2])
            centered_speed = circshift(vertical_windspeed,(center[1] - cyclone_center[1],center[2] - cyclone_center[2]));
            count += 1
            nup, ndown, npoints = count_updrafts_and_downdrafts(radiusbins,centered_speed,center;gridspacing = 2000)
            nupdrafts_per_distance .= nupdrafts_per_distance .+ nup
            ndowndrafts_per_distance .= ndowndrafts_per_distance .+ ndown 
            npoints_in_distance .= npoints
        end
    end
    return count, nupdrafts_per_distance, ndowndrafts_per_distance, npoints_in_distance
end 
function add_updraft_and_downdrafts_all_cyclones_in_all_timesteps_and_do_composite(radiusbins,surface_pressure :: Array{T,3},vertical_windspeed :: Array{T,3},u_speed :: Array{T,4},v_speed :: Array{T,4},center;gridspacing = 1) where T
    cyclonecount = 0
    sx,sy,st = size(vertical_windspeed)
    sxu,syu,szu,stu = size(u_speed)
    u_composite = zeros(Float64,sxu,syu,szu)
    v_composite = zeros(Float64,sxu,syu,szu)
    buf1 = similar(u_composite)
    buf2 = similar(u_composite)
    nupdrafts_per_distance = zeros(Int64,length(radiusbins) - 1)
    ndowndrafts_per_distance = zeros(Int64,length(radiusbins) - 1)
    npoints_in_distance = zeros(Int64,length(radiusbins) - 1)
    for timestep in 1:st
        surface_pressure_anomaly = surface_pressure[:,:,timestep] .- mean(surface_pressure[:,:,timestep],dims = (1,2));
        cyclones = detect_cyclones(AvailablePotentialEnergyFramework.PressureMinima(),surface_pressure_anomaly,-9,2000)
        if cyclones == (nothing,nothing)
            continue
        end
        cyclone_count, nup, ndown, npoints = add_updraft_and_downdrafts_all_cyclones_in_timestep(cyclones, radiusbins, vertical_windspeed[:,:,timestep], (256,256);gridspacing = 2000)
        count = add_allcyclones!(u_composite, buf1, buf2, u_speed[:,:,:,timestep], cyclones[2], cyclones[1]; maskcyclones = false)
        count = add_allcyclones!(v_composite, buf1, buf2, v_speed[:,:,:,timestep], cyclones[2], cyclones[1]; maskcyclones = false)
        nupdrafts_per_distance .= nupdrafts_per_distance .+ nup
        ndowndrafts_per_distance .= ndowndrafts_per_distance .+ ndown
        cyclonecount = cyclonecount + cyclone_count
        npoints_in_distance .= npoints
    end
    return cyclonecount, nupdrafts_per_distance, ndowndrafts_per_distance, npoints_in_distance, u_composite./cyclonecount, v_composite./cyclonecount
end


function get_snapshot_from_3d(file_path, var,level,timestep)
    profile = NetCDF.open(file_path) do f
        wholevar = f[var][:,:,level,timestep]
    end
end
function get_snapshot_from_2d(file_path, var,timestep)
    profile = NetCDF.open(file_path) do f
        wholevar = f[var][:,:,timestep]
    end
end

function readfile_and_count_updrafts_and_downdrafts(radiusbins,data_dir,file_2d,file_3d,z_level,pick_timesteps_3d,center;gridspacing = 2000)
    pick_timesteps_2d = pick_timesteps_3d .* 2 .- 1
    w = get_snapshot_from_3d(joinpath(data_dir,file_3d),"W",z_level,pick_timesteps_3d)
    u = get_snapshot_from_3d(joinpath(data_dir,file_3d),"U",1:50,pick_timesteps_3d)
    v = get_snapshot_from_3d(joinpath(data_dir,file_3d),"V",1:50,pick_timesteps_3d)
    p = get_snapshot_from_2d(joinpath(data_dir,file_2d),"PSFC",pick_timesteps_2d)
    cyclone_count, nup, ndown, npoints = add_updraft_and_downdrafts_all_cyclones_in_all_timesteps_and_do_composite(radiusbins,p,w,u,v,(256,256);gridspacing = 2000)
end


function get_tangential_and_radial_speed(composite)
    tangential = similar(composite[6])
    radial = similar(composite[6]);
    for index in CartesianIndices(tangential)
        center = (256,256)
        index_of_point = (index[1],index[2])
        tangential[index],radial[index] = AvailablePotentialEnergyFramework.velocity_cartesian_to_polar(composite[5][index],composite[6][index],index_of_point,center)
    end
    return tangential,radial
end

function get_azimuthal_average(array :: Array{T,3},radiusbins) where T
    azimuthalaverage = zeros(eltype(array),length(radiusbins) - 1,size(array,3));
    for rindex in 1:(length(radiusbins) - 1)
        azimuthalaverage[rindex,:] .= AvailablePotentialEnergyFramework.averageallindistance((radiusbins[rindex],radiusbins[rindex+1]),array,(256,256),2000.0)   
    end
    return azimuthalaverage
end


function get_inflow_and_rmax(composite, radiusbins)
    tangential_profiles = get_tangential_and_radial_speed(composite)
    radial = get_azimuthal_average(tangential_profiles[2], radiusbins)
    tangential = get_azimuthal_average(tangential_profiles[1], radiusbins)
    max_wind,index_rmax = findmax(tangential[:,1])
    rmax = radiusbins[index_rmax]
    index_final = findmin(abs,radiusbins .- 6rmax)[2]
    index_middle = ceil(Int,(index_final + index_rmax) ÷ 2)
    inflow = mean(radial[index_rmax:index_final,1:32])
    return inflow,index_rmax,max_wind
end

function normalize_updrafts_and_downdrafts(composite)
    cyclonecount, nupdrafts_per_distance, ndowndrafts_per_distance, npoints_in_distance, u_composite, v_composite = composite
    normalized_updraft = nupdrafts_per_distance ./ cyclonecount ./ npoints_in_distance
    normalized_downdraft = ndowndrafts_per_distance ./ cyclonecount ./ npoints_in_distance
    return normalized_updraft, normalized_downdraft
end


function compute_updraft_and_downdraft_ratio(updrafts,downdrafts,index_rmax, radiusbins)
    cyclone_size = radiusbins[index_rmax]
    index_final = findmin(abs,radiusbins .- 6cyclone_size)[2]
    index_middle = findmin(abs,radiusbins .- 3cyclone_size)[2]
    
    updrafts_first = sum(updrafts[1:index_middle]) 
    updrafts_second =  sum(updrafts[index_middle:index_final]) 
    updrafts_total =  sum(updrafts[1:index_final]) 
    
    return updraft_ratio = updrafts_second / updrafts_first, updrafts_second / updrafts_total
end

function compute_updraft_ratio_and_mean_inflow(composite,radiusbins)
    updrafts, downdrafts  = normalize_updrafts_and_downdrafts(composite)
    inflow,index_rmax,max_wind = get_inflow_and_rmax(composite,radiusbins)
    updraft_ratio = compute_updraft_and_downdraft_ratio(updrafts,downdrafts,index_rmax, radiusbins)
    return updraft_ratio, inflow,max_wind
end


function read_file_and_compute_updraft_ratio_and_mean_inflow(radiusbins,data_dir,file2d,file3d,z_index,days; center = (256,256), gridspacing = 2000)
    composite = readfile_and_count_updrafts_and_downdrafts(radiusbins,data_dir,file2d,file3d,z_index,days_to_indices(days...),center;gridspacing)
    return updraft_ratio, inflow, max_wind = compute_updraft_ratio_and_mean_inflow(composite,radiusbins)
end
