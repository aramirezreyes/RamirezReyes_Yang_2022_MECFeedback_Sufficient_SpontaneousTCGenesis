"""

        composite_and_radialaverage_withmask(azimuthal_average,composite_addition,buf,radiusbins,surface_pressure,variable_of_interest)
    """


function composite_and_radialaverage_in_chunks!(azimuthal_average,composite_addition,buf1,buf2,radiusbins,variable_name,pre_iterator_time,nchunks,file_data,diagnostics_data_dir,exp_name,centers_labels_and_cyclones;maskcyclones = true)
    totalcyclonecount = 0
    variables_diag = ["convec_heating_anomaly","rad_heating_anomaly","buoyancy_anomaly"]
    variables_orig_3d = ["U","V","QV","TABS","QRAD","PP","W"]
    variables_orig_2d = ["PSFC","USFC","VSFC","PW","Prec","LHF","SHF"]
    chunk_size_2d = 2*length(pre_iterator_time)÷nchunks 
    chunk_size_3d = length(pre_iterator_time)÷nchunks
    
    nd = ndims = if (variable_name ∈ variables_orig_3d)
            nd = ndims = 4
        elseif (variable_name ∈ variables_diag) 
            nd = ndims = 4
        elseif (variable_name ∈ variables_orig_2d) 
            nd = ndims = 3
    end
    for chunk in 1:nchunks
        iterator_time_2d = (1202 + (chunk - 1)*chunk_size_2d): 2 : (1202 + chunk*chunk_size_2d - 1)
        iterator_time_data =  ((chunk - 1)*chunk_size_3d + 1 ) : (chunk*chunk_size_3d)
        iterator_time_3d = 600 .+ iterator_time_data
        if (variable_name ∈ variables_orig_3d)
            variable_of_interest = Dataset(file_data) do ds
                variable(ds,variable_name)[:,:,:,iterator_time_3d]     :: Array{Float32,4}
            end
        elseif (variable_name ∈ variables_diag) 
            variable_of_interest = gather_diagnostic_3d_in_chunks(diagnostics_data_dir,exp_name,variable_name,chunk,chunk_size_3d,nchunks)
        elseif (variable_name ∈ variables_orig_2d) 
            variable_of_interest = Dataset(file_data) do ds
                variable(ds,variable_name)[:,:,iterator_time_2d]     :: Array{Float32,3}
            end
        end
        

        
        for timeindex_in in 1:size(variable_of_interest,nd) #TODO this time index must correspond to an "inside" index
            timeindex_out = iterator_time_data[timeindex_in]#corresponds to the index in the centers_labels_and_cyclones array
            centers_and_labels,cyclones = centers_labels_and_cyclones[timeindex_out]
            if !isnothing(centers_and_labels)
                count = add_allcyclones!(composite_addition,buf1,buf2,selectdim(variable_of_interest,nd,timeindex_in),cyclones,centers_and_labels;maskcyclones)
                totalcyclonecount += count
            end
        end
    end
    
    if !iszero(totalcyclonecount)
        composite_addition ./= totalcyclonecount
        if nd == 4
            for idx in 1:(length(radiusbins) - 1)            
                azimuthal_average[idx,:] .= averageallindistance((radiusbins[idx],radiusbins[idx+1]),composite_addition,(256,256), 2000 )
            end
        elseif nd == 3
            for idx in 1:(length(radiusbins) - 1)            
                azimuthal_average[idx] = averageallindistance((radiusbins[idx],radiusbins[idx+1]),composite_addition,(256,256), 2000 )
            end
        end
        
        
    end
    return totalcyclonecount
end



function gather_diagnostic_3d_in_chunks(data_dir,exp_name,desired_var,chunk,chunksize,nchunks)
    diagnostic = zeros(Float32,512,512,80,chunksize)
    howmanyfiles = 12÷nchunks
    @info "how many files" howmanyfiles
    for iter in 1:howmanyfiles
        filenumber = ((chunk - 1)*howmanyfiles + iter)
        
        @info "this is iter" iter
        @info "file number and chunk:" filenumber chunk
        file_name = string(exp_name,"_",filenumber+12,".jld")
        @info "I am opening the file $file_name and reading the variable $desired_var"
        file_path = joinpath(data_dir,file_name)
        file = jldopen(file_path, "r")
        index_initial = (iter - 1)*(50)+1
        index_end = iter*50
        @info "indices to write:" index_initial,index_end
        var_from_file = read(file, desired_var)
        diagnostic[:,:,:,index_initial:index_end] .=  var_from_file[:,:,:,:]
        close(file)
    end
    return diagnostic
end



"""
        getcomposites(file2d,file3d,diagnostics_dir,exp_,initial_timeindex,final_timeindex,partnumber,outputInterval)


"""
function get_composites_in_chunks(file2d,file3d,exp_name)

    @info "Starting composite routine"
    diagnostics_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/MECFeedback/ApeBudgetOutputs_nosmoothing/"
    composites_data_dir =  "/global/cscratch1/sd/aramreye/for_postprocessing/MECFeedback/CompositeOutputs_50d_9hpa/"
    azimuthal_averages_data_dir = "/global/cscratch1/sd/aramreye/for_postprocessing/MECFeedback/AzimuthualAverageOutputs_50d_9hpa/"
    centers_data_file = string(composites_data_dir,"centers_detection_",exp_name,".jld")
    
    variables_diag = ["convec_heating_anomaly","rad_heating_anomaly","buoyancy_anomaly"]
    variables_orig_3d = ["U","V","QV","TABS","QRAD","PP","W"]
    variables_orig_2d = ["PSFC","USFC","VSFC","PW","Prec","LHF","SHF"]

    #variables_diag = []
    #variables_orig_3d = ["U","V"]
    #variables_orig_2d = ["PSFC","USFC","VSFC"]
    
    initial_timeindex = 601
    #    initial_timeindex = 1195
    final_timeindex = 1200
    #final_timeindex = 650
    
    
    pre_iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
    pre_iterator_time = pre_iterator_time_3d   = initial_timeindex:final_timeindex


    
    radiusbins = 1000:2000:512000; 
    ds3d                = Dataset(file3d)
    ds2d                = Dataset(file2d)
    z                   = variable(ds3d,"z")[:]                       :: Array{Float32,1}
    t                   = variable(ds3d,"time")[pre_iterator_time_3d]     :: Array{Float32,1}
    surface_pressure    = variable(ds2d,"PSFC")[:,:,pre_iterator_time_2d] :: Array{Float32,3}
    close(ds2d)
    close(ds3d)
    pressure_anomaly    = surface_pressure .- mean(surface_pressure,dims=(1,2))

    totalcyclonecount = 0
    
    @info size(surface_pressure) 


    ##Create buffers to reuse on allocating functions

    size_var_of_interest_2d = (512,512)
    size_var_of_interest_3d = (size_var_of_interest_2d...,length(z))
    ## Create 3D buffers
    composite_addition_3d   = zeros(Float32,size_var_of_interest_3d)
    buf_3d_2                = similar(composite_addition_3d)
    buf_3d_1                = similar(composite_addition_3d)
    azimuthal_average_3d    = zeros(Float32,length(radiusbins) - 1,length(z))
    ## Create 2D buffers
    composite_addition_2d   = zeros(Float32,size_var_of_interest_2d)
    buf_2d_1                = similar(composite_addition_2d)
    buf_2d_2                = similar(composite_addition_2d)
    azimuthal_average_2d    = zeros(Float32,length(radiusbins) - 1)
    @info size_var_of_interest_2d size_var_of_interest_3d

    centers_labels_and_cyclones = @views [
        detect_cyclones!(AvailablePotentialEnergyFramework.PressureMinima(),buf_2d_1, buf_2d_2, pressure_anomaly[:,:,timeindex],-5,2000) for
        timeindex in 1:size(surface_pressure,3)]; pressure_anomaly = []; surface_pressure = [];

    composite_filename_nomask   = string(composites_data_dir,exp_name        ,"_nomask.jld")
    azimuthal_filename_nomask   = string(azimuthal_averages_data_dir,exp_name,"_nomask.jld")

    ### *** add allcyclones ***
    @info "Creating composite for surface pressure "
    flush(stdout)
    composite_addition_2d .= 0.0
    buf_2d_1              .= 0.0
    buf_2d_2              .= 0.0
    azimuthal_average_2d  .= 0.0

    nchunks=6   ## this is for 50 days!
     
    #create files:
    ds =  jldopen(composite_filename_nomask, "w") do file end
    ds =  jldopen(azimuthal_filename_nomask, "w") do file end

    for current_variable in variables_orig_3d
           
     #     ## Create products without mask
         @info "Creating composite without mask for: " current_variable
         flush(stdout)
         composite_addition_3d .= 0.0
         buf_3d_1              .= 0.0
         buf_3d_2              .= 0.0
         azimuthal_average_3d  .= 0.0
         totalcyclonecount = composite_and_radialaverage_in_chunks!(azimuthal_average_3d,composite_addition_3d,buf_3d_1,buf_3d_2,radiusbins,current_variable,pre_iterator_time_3d,nchunks,file3d,diagnostics_data_dir,exp_name,centers_labels_and_cyclones;maskcyclones = false)
         @info "Writing"
         flush(stdout)
         jldopen(composite_filename_nomask, "r+",mmaparrays=true) do file
             write(file,current_variable,composite_addition_3d)
         end        
         jldopen(azimuthal_filename_nomask, "r+",mmaparrays=true) do file
             write(file,current_variable,azimuthal_average_3d)
         end

        if current_variable == first(variables_orig_3d)
            jldopen(composite_filename_nomask, "r+",mmaparrays=true) do file
             write(file,"n_cyclones_averaged",totalcyclonecount)
         end
        end
         GC.gc()
    end
    

    for current_variable in variables_orig_2d
         composite_addition_2d .= 0.0
         buf_2d_1              .= 0.0
         buf_2d_2              .= 0.0
         azimuthal_average_2d  .= 0.0

         ## Create products without mask
         @info "Creating composite without mask for: " current_variable
         flush(stdout)
         composite_addition_2d .= 0.0
         buf_2d_1              .= 0.0
         buf_2d_2              .= 0.0
         azimuthal_average_2d  .= 0.0
    
         composite_and_radialaverage_in_chunks!(azimuthal_average_2d,composite_addition_2d,buf_2d_1,buf_2d_2,radiusbins,current_variable,pre_iterator_time_3d,nchunks,file2d,diagnostics_data_dir,exp_name,centers_labels_and_cyclones;maskcyclones = false)
         @info "Writing"
         flush(stdout)
         jldopen(composite_filename_nomask, "r+",mmaparrays=true) do file
             write(file,current_variable,composite_addition_2d)
         end        
         jldopen(azimuthal_filename_nomask, "r+",mmaparrays=true) do file
             write(file,current_variable,azimuthal_average_2d)
         end
         GC.gc()
    end


    for current_variable in variables_diag
    
         ## Create products without mask
         @info "Creating composite without mask for: " current_variable
         flush(stdout)
         composite_addition_3d .= 0.0
         buf_3d_1              .= 0.0
         buf_3d_2              .= 0.0
        azimuthal_average_3d  .= 0.0
        composite_and_radialaverage_in_chunks!(azimuthal_average_3d,composite_addition_3d,buf_3d_1,buf_3d_2,radiusbins,current_variable,pre_iterator_time_3d,nchunks,file3d,diagnostics_data_dir,exp_name,centers_labels_and_cyclones;maskcyclones = false)
    
         @info "Writing"
         flush(stdout)
         jldopen(composite_filename_nomask, "r+",mmaparrays=true) do file
             write(file,current_variable,composite_addition_3d)
         end        
         jldopen(azimuthal_filename_nomask, "r+",mmaparrays=true) do file
             write(file,current_variable,azimuthal_average_3d)
         end
         GC.gc()
         GC.gc()
    end

    

  
    
    return nothing

end
