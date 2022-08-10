

function detect_and_plot_cyclones(file2d,exp_name,output_dir; pressure_threshold = -9, offset_in_days = 0, total_days = 100)
    output_timestep = 3600 #seconds
    indices_in_part = 50
    indices_in_day = 86400 ÷ output_timestep
    number_of_parts = total_days*indices_in_day÷indices_in_part
    initial_timestep = offset_in_days*indices_in_day + indices_in_part*(part_number-1) + 1
    final_timestep = offset_in_days*indices_in_day + indices_in_part*(part_number)

    iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
    iterator_time = pre_iterator_time_3d   = initial_timeindex:final_timeindex

    t, surface_pressure = Dataset(file2d) do dataset
        t                   = variable(dataset,"time")[terator_time_2d]     :: Array{Float32,1}
        surface_pressure    = variable(dataset,"PSFC")[:,:,iterator_time_2d] :: Array{Float32,3}
        (t, surface_pressure)
    end
    
    pressure_anomaly    = surface_pressure .- mean(surface_pressure,dims=(1,2))
    
    @info "The surface pressure array is of size" size(surface_pressure) 
    sx,sy,st = size(pressure)
    buf2d_1 = Array{eltype(surface_pressure),2}(undef,sx,sy);
    buf2d_2 = Array{eltype(surface_pressure),2}(undef,sx,sy);
    centers_labels_and_cyclones = @views [
        detect_cyclones!(buf_2d_1, buf_2d_2, pressure_anomaly[:,:,timeindex],pressure_threshold,2000) for
        timeindex in 1:st ];
    centers_labels_and_cyclones,surface_pressure
end
