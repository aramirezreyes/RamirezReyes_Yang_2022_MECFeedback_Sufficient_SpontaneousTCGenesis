"""
    getapeanalysis(file2d,file3d,outputfile,initial_timeindex,final_timeindex,partnumber,outputInterval,position,FloatType::Type=Float64)
"""
function getapeanalysis(file2d,file3d,outputfile,initial_timeindex,final_timeindex,partnumber,outputInterval,position,FloatType::Type=Float32;tropopause_height = 15000,dosmoothing=true,smooth_x = 15, smooth_time_in_days = 5)

    @info "Starting APE analysis routine"
    day           = 86400
    sst           = 300
    dt            = outputInterval
    
    dayLength     = 60*60*24÷outputInterval; #How many data points make one day
    #final_time    = floor(Int,final_time_days*dayLength) 
    #initial_time  = floor(Int,initial_time_days*dayLength)+1
    
    Pref          = 1000*1f2                 #Pa
    Ts            = sst                      #Sea surface temperature
    qs            = 25.7*1f-3 
    c1            = (ustrip(Dryair.R)/ustrip(Dryair.cp))
    

    #iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
    #iterator_time_3d    = initial_timeindex:1:final_timeindex
    if dosmoothing
        #smooth_x      = smooth_y = 15 #it was 11
        smooth_time   = floor(Int,dayLength*smooth_time_in_days)+1 

        if position == 1
            iterator_time_2d    = initial_timeindex*2-1:2:(final_timeindex+div(smooth_time-1,2))*2  
            iterator_time_3d    = initial_timeindex:1:final_timeindex+div(smooth_time-1,2)
        elseif position == 2
            iterator_time_2d    = (initial_timeindex-div(smooth_time-1,2))*2-1:2:(final_timeindex+div(smooth_time-1,2))*2  
            iterator_time_3d    = initial_timeindex-div(smooth_time-1,2):1:final_timeindex+div(smooth_time-1,2)
        elseif position == 3
            iterator_time_2d    = (initial_timeindex-div(smooth_time-1,2))*2-1:2:final_timeindex*2  
            iterator_time_3d    = initial_timeindex-div(smooth_time-1,2):1:final_timeindex
        elseif position == 4
        end
        @info " Reading files: " file2d file3d
        
        ds3d                = Dataset(file3d)
        ds2d                = Dataset(file2d)
        x                   = variable(ds3d,"x")[:]    :: Array{Float32,1}
        y                   = variable(ds3d,"y")[:]    :: Array{Float32,1}
        z                   = variable(ds3d,"z")[:]    :: Array{Float32,1}
        t                   = variable(ds3d,"time")[iterator_time_3d] :: Array{Float32,1}
        P0                  = variable(ds3d,"p")[:]    :: Array{Float32,1}
        U                   = variable(ds3d,"U")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 10%"    
        V                   = variable(ds3d,"V")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 20%"    
        W                   = variable(ds3d,"W")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 30%"    
        RAD                 = variable(ds3d,"QRAD")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 40%"    
        T                   = variable(ds3d,"TABS")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 50%"    
        Tv                  = variable(ds3d,"QV")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 60%"    
        PP                  = variable(ds3d,"PP")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 70%"                  
        SHF                 = variable(ds2d,"SHF")[:,:,iterator_time_2d] :: Array{Float32,3}
        @info "Read 80%"                 
        LHF                 = variable(ds2d,"LHF")[:,:,iterator_time_2d] :: Array{Float32,3}
        @info "Read 90%"

        close(ds3d)   
        close(ds2d)   
        


        #if 3 nodes:
        # node 1: workers 1, 2, 5,8
        # node 2: workers 3, 6, 9
        # node 3: workers 4, 7, 10
        @info "Smoothing batch 1"
        fU   = @spawnat 3 cutborders!(filter_array(U,smooth_x,smooth_time,1),smooth_time,position)
        fV   = @spawnat 4 cutborders!(filter_array(V,smooth_x,smooth_time,1),smooth_time,position)
        fW   = @spawnat 6 cutborders!(filter_array(W,smooth_x,smooth_time,1),smooth_time,position)
        fTv  = @spawnat 7 cutborders!(filter_array(Tv,smooth_x,smooth_time,1),smooth_time,position)
        U = fetch(fU)
        V = fetch(fV)
        W = fetch(fW)
        Tv = fetch(fTv)
        
        fetch(@spawnat 3 GC.gc())
        fetch(@spawnat 4 GC.gc())
        fetch(@spawnat 6 GC.gc())
        fetch(@spawnat 7 GC.gc())
        @info "smoothing batch 2"
        fT   = @spawnat 3 cutborders!(filter_array(T,smooth_x,smooth_time,1),smooth_time,position)
        fPP  = @spawnat 4 cutborders!(filter_array(PP,smooth_x,smooth_time,1),smooth_time,position)
        fRAD = @spawnat 6 cutborders!(filter_array(RAD,smooth_x,smooth_time,1),smooth_time,position)
        SHF = cutborders!(filter_array(SHF,smooth_x,smooth_time,1),smooth_time,position)
        LHF = cutborders!(filter_array(LHF,smooth_x,smooth_time,1),smooth_time,position)
        t = cutborders!(t,smooth_time,position)
        
        T = fetch(fT)
        PP = fetch(fPP)
        RAD = fetch(fRAD)

        @everywhere GC.gc()
        @everywhere GC.gc()

    elseif !dosmoothing


        iterator_time_2d    = initial_timeindex*2-1:2:final_timeindex*2
        iterator_time_3d    = initial_timeindex:1:final_timeindex

        ### Array allocation

                @info " Reading files: " file2d file3d
        
        ds3d                = Dataset(file3d)
        ds2d                = Dataset(file2d)
        x                   = variable(ds3d,"x")[:]    :: Array{Float32,1}
        y                   = variable(ds3d,"y")[:]    :: Array{Float32,1}
        z                   = variable(ds3d,"z")[:]    :: Array{Float32,1}
        t                   = variable(ds3d,"time")[iterator_time_3d] :: Array{Float32,1}
        P0                  = variable(ds3d,"p")[:]    :: Array{Float32,1}
        U                   = variable(ds3d,"U")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 10%"    
        V                   = variable(ds3d,"V")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 20%"    
        W                   = variable(ds3d,"W")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 30%"    
        RAD                 = variable(ds3d,"QRAD")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 40%"    
        T                   = variable(ds3d,"TABS")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 50%"    
        Tv                  = variable(ds3d,"QV")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 60%"    
        PP                  = variable(ds3d,"PP")[:,:,:,iterator_time_3d] :: Array{Float32,4}
        @info "Read 70%"                  
        SHF                 = variable(ds2d,"SHF")[:,:,iterator_time_2d] :: Array{Float32,3}
        @info "Read 80%"                 
        LHF                 = variable(ds2d,"LHF")[:,:,iterator_time_2d] :: Array{Float32,3}
        @info "Read 90%"

        close(ds3d)   
        close(ds2d)   
        
    end

    ThetaV       = similar(T)
    B            = similar(T)
    var_Tv       = similar(T)
    var_ThetaV   = similar(T)
    xBar_Pt      = Array{eltype(T),4}(undef,1,1,size(PP,3),size(PP,4))
    rho0         = Array{eltype(T),4}(undef,1,1,size(PP,3),size(PP,4))
    xBar_Tv      = Array{eltype(T),4}(undef,1,1,size(Tv,3),size(Tv,4))
    xBar_ThetaV  = Array{eltype(T),4}(undef,1,1,size(ThetaV,3),size(ThetaV,4))
    #******
    dx           = x[2]-x[1]
    dy           = y[2]-y[1]                            # Grid size
    kz           = length(z)                            # vertical levels
    kx           = length(x)                            # # of horizonal grid points
    ky           = length(y)                            # # of horizonal grid points


    ###### Processing
    
    @. SHF    = surface_sensible_heat_flux_to_buoyancy(Ts, SHF)                            
    @. LHF    = surface_latent_heat_flux_to_buoyancy(Ts, LHF)
    @. SHF   += LHF                                 # Now it is transformed
    @. RAD    = RAD/day                              # K/s Heating rate per second
    @. Tv     = get_virtual_temperature(T,Tv)        # Virtual temperature
    PP       .= PP .+ 1f2reshape(P0,(1,1,kz,1))
    @. ThetaV = get_potential_temperature(Tv,PP,P0[1]) # Virtual potential temp

    mean!(xBar_Pt,PP)                                     
    mean!(xBar_Tv,Tv)              
    mean!(xBar_ThetaV,ThetaV)   
    @. var_Tv       = Tv     - xBar_Tv
    @. var_ThetaV   = ThetaV - xBar_ThetaV
    @. rho0         = xBar_Pt ./xBar_Tv ./ Float32(ustrip(Dryair.R))
    @. B            = get_buoyancy(var_ThetaV,xBar_ThetaV) # g .* var_ThetaV./xBar_ThetaV
    @. RAD          = radiative_heating_rate_to_buoyancy(xBar_Tv,RAD)
    N2              = compute_N2(xBar_Tv,z)
    PP              = []
    LHF             = []
    T               = []
    ThetaV          = []
    Tv              = []
    var_ThetaV      = []
    var_Tv          = []
    @info "Computing buoyancy budget"
    ### NOTE that SHF is now the sum, saving memory
    # Buoyancy budget
    dz          = 50
    @info size(B), size(RAD), size(SHF), size(U),size(V) ,size(W), size(N2), size(dx),size(dy), size(dz), size(dt), size(x),size(y), size(z), size(t)
    Diabatic_other = fetch(@spawnat 3  get_diabatic_as_residual_buoyancy(B, RAD, SHF, U,V ,W, N2, dx,dy, dz, dt))
    @everywhere GC.gc()
    @everywhere GC.gc()
    @info "Computing ape budget for the whole troposphere"

    # APE budget
    z_up        = tropopause_height
    z_BL        = 2000

     (int_mass,
     int_KE,
     int_APE,
     int_APE_rate,
     int_APE_Ub2,
     int_APE_Vb2,
     int_APE_WN2,
     int_APE_RAD,
     int_APE_DIA,
     xBar_APE_Fs,
     residual) = fetch( @spawnat 3 getapebudget(B, U,V, W, N2, RAD, SHF, Diabatic_other, rho0, x,y, z, t, dx,dy, dz, dt, z_up))
    @everywhere GC.gc()
    @everywhere GC.gc()

    Diabatic_other .= Diabatic_other .- mean(Diabatic_other,dims=(1,2)) #They are now perturbations
    RAD            .= RAD .- mean(RAD,dims=(1,2))
    #dia_ape         = Diabatic_other.*(B .- mean(B,dims=(1,2))) 
    #rad_ape         = RAD .* (B .- mean(B,dims=(1,2))) 

    @info "WritingData"


    jldopen(string(outputfile,"_",partnumber,".jld"), "w", compress=true) do file
        write(file,"rho0",rho0)
        write(file,"N2",N2)
        write(file,"int_APE",int_APE)
        write(file,"int_KE",int_KE)
        write(file,"int_RAD",int_APE_RAD)
        write(file,"int_DIA",int_APE_DIA)
        write(file,"int_WN2",int_APE_WN2)
        write(file,"int_Ub2",int_APE_Ub2)
        write(file,"int_Vb2",int_APE_Vb2)
        write(file,"int_APE_rate",int_APE_rate)
        write(file,"APE_Fs",xBar_APE_Fs)
        write(file,"convec_heating_anomaly",Diabatic_other)
        write(file,"rad_heating_anomaly",RAD)
        write(file,"buoyancy_anomaly",B)
    end
    GC.gc()
    @everywhere GC.gc()
    return nothing

end

function computebudgets(exp_name; total_days = 100, offset_in_days = 0, dosmoothing = true, smooth_x = 15, smooth_time_in_days = 5, outputsuffix = "")
    outputfilename = string(exp_name)
    outputfolder = if dosmoothing
        if smooth_x != 15
            outputfilename = string(outputfilename,"_smoothx_$(smooth_x)")
        end
        if smooth_time_in_days != 5
            outputfilename = string(outputfilename,"_smootht_$(smooth_time_in_days)")
        end
        "/global/cscratch1/sd/aramreye/for_postprocessing/MECFeedback/ApeBudgetOutputs/"        
    elseif !dosmoothing
        "/global/cscratch1/sd/aramreye/for_postprocessing/MECFeedback/ApeBudgetOutputs_nosmoothing/"
    end
    if total_days != 100 || offset_in_days != 0
        outputfilename = string(outputfilename,"_$(total_days)days_startingfromday_$(offset_in_days)")
    end

    outputfile = joinpath(outputfolder,outputfilename)
    input_folder = "/global/cscratch1/sd/aramreye/for_postprocessing/largencfiles/"
    file2d = string(input_folder,exp_name,"_2d.nc")
    file3d = string(input_folder,exp_name,"_3d.nc")
    #total_days = 100
    output_timestep = 7200 #seconds
    indices_in_part = 50
    indices_in_day = 86400 ÷ output_timestep
    number_of_parts = total_days*indices_in_day÷indices_in_part

    for part_number in 1:number_of_parts
        if part_number == 1
            position = 1
        elseif part_number == number_of_parts
            position = 3
        else
            position = 2
        end
        initial_timestep = offset_in_days*indices_in_day + indices_in_part*(part_number-1) + 1
        final_timestep = offset_in_days*indices_in_day + indices_in_part*(part_number)

        @info "Starting part $part_number with indices from" initial_timestep final_timestep
        flush(stdout)
        flush(stderr)
        getapeanalysis(file2d,file3d,outputfile,initial_timestep,final_timestep,part_number,output_timestep,position,Float32; dosmoothing, smooth_x, smooth_time_in_days)

    end

end

function computebudgets_nosmoothing(exp_name)

    #file2d = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_2d.nc"
    #file3d = "/global/cscratch1/sd/aramreye/sam3d/subsetsForLin/f5e-4_2km_1000km_control_3d.nc"
    outputfile = string("/global/cscratch1/sd/aramreye/for_postprocessing/MECFeedback/ApeBudgetOutputs_nosmoothing/",exp_name)
    input_folder = "/global/cscratch1/sd/aramreye/for_postprocessing/largencfiles/"
    file2d = string(input_folder,exp_name,"_2d.nc")
    file3d = string(input_folder,exp_name,"_3d.nc")
    total_days = 100
    output_timestep = 7200 #seconds
    indices_in_part = 50
    indices_in_day = 86400 ÷ output_timestep
    number_of_parts = total_days*indices_in_day÷indices_in_part

    for part_number in 1:number_of_parts
        if part_number == 1
            position = 1
        elseif part_number == number_of_parts
            position = 3
        else
            position = 2
        end
        initial_timestep = indices_in_part*(part_number-1) + 1
        final_timestep = indices_in_part*(part_number)
        @info "Starting part $part_number with indices from" initial_timestep final_timestep
        flush(stdout)
        flush(stderr)
        getapeanalysis(file2d,file3d,outputfile,initial_timestep,final_timestep,part_number,output_timestep,position,Float32; dosmoothing=false)

    end

end

