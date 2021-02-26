## run_c.jl --- routines to run the coupled photochemistry model
## to calculate the evolution of atmosphere saving data every time point

folder_directory_run = "/Users/shungokoyama/programming/photochem_evolution/"
including_file = folder_directory_run * "photochemistry_c.jl"
include(including_file)
update!(n_current,0.)

# we can run with logarithmic time steps
timepts = 10 .^ range( log10(0.1), log10(3.15e12), length = 1300 )
timepts_added = 10 .^ range( log10(3.15e12), log10(4.2*3.15e14), length = 2700 )
append!(timepts,timepts_added)
timediff = timepts[2:end]-timepts[1:end-1]
#append!(timediff,1.e12*ones(Float64,500)) #add

core_filename = "CO2_" * string(CO2_pressure) * "mbar_Tsurf_" * string(surface_temperature) *
                "_Tian_BD.h5"
output_folder = "output"
mkdir(output_folder)
run_filename = folder_directory_run * output_folder * "/"* core_filename
hfile = folder_directory_run * output_folder * "/" * "Hesc_depfluxes_" * core_filename

#stock C-related BD fluxes
CO2outgassing_list = Float64[]
Cesc_list = Float64[]
Oreturn_list = Float64[]
Hesc_list = Float64[]
Othermal_list = Float64[]
Cthermal_list = Float64[]
Cnonthermal_list = Float64[]
Ononthermal_list = Float64[]


function runprofile(n_current, dtlist, filename)

    #このn_internalをこの関数の中で更新していく。
    #n_internal = deepcopy(n_current)
    elapsed_time = 0.0

    # create a matrix to contain the data in n_current, which is a Dict
    #それぞれのspeciesに対する高度プロファイル
    n_internal_mat = Array{Float64}(undef,length(alt)-2,length(collect(keys(n_current))));
    for ispecies in 1:length(collect(keys(n_current)))
        for ialt in 1:length(alt)-2
            n_internal_mat[ialt,ispecies] = n_current[collect(keys(n_current))[ispecies]][ialt]
        end
    end

    #write initial input of the density profile
    h5write(filename,"n_current/init",n_internal_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string,collect(keys(n_current))))
    h5write(filename,"n_current/timelist",cumsum(dtlist)) #cumulative sum along this vector

    #stock t-dependent fluxes, 2020/6/18, koyama added
    append!(CO2outgassing_list, eval(speciesbclist[:CO2][1,2]))
    append!(Cesc_list, eval(speciesbclist[:CO][2,2]) * n_current[:CO][end])
    append!(Oreturn_list, eval(speciesbclist[:CO][2,2]) * -1* n_current[:CO][end] )
    append!(Hesc_list, eval(speciesbclist[:H][2,2])*n_current[:H][end] + 2*eval(speciesbclist[:H2][2,2])*n_current[:H2][end])
    append!(Othermal_list, calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:O][end],16.)*O_effusion_velocity*modified_factor(n_current,Texo_tian,Zexo_tian,16.0))
    append!(Cthermal_list, calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:CO][end]*C_CO2_ratio,12.)*C_CO_ratio*flex_effusion_velocity(Texo_tian,12.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,12.0))
    append!(Cnonthermal_list, C_nonthermal)
    append!(Ononthermal_list, O_nonthermal)

    # TIME LOOP - simulate over time
    # This section was changed so that it only prints statements as it starts
    # and finishes new tasks, rather than printing a statement at each dt.
    # the old lines for printing at each dt have been commented out and can be
    # restored if desired.
    thisi=0
    for dt in dtlist
        #println(filename*": iteration = "* string(thisi+=1)*" "*Libc.strftime(time()))
        thisi += 1
        #println("dt = "* string(dt::Float64))
        elapsed_time+=dt
        #println("elapsed_time = "*string(elapsed_time))

        #stock t-dependent fluxes, 2020/6/18, koyama added
        append!(CO2outgassing_list, eval(speciesbclist[:CO2][1,2]))
        append!(Cesc_list, eval(speciesbclist[:CO][2,2]) * n_current[:CO][end])
        append!(Oreturn_list, eval(speciesbclist[:CO][2,2]) * -1* n_current[:CO][end] )
        append!(Hesc_list, eval(speciesbclist[:H][2,2])*n_current[:H][end] + 2*eval(speciesbclist[:H2][2,2])*n_current[:H2][end])
        append!(Othermal_list, calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:O][end],16.)*O_effusion_velocity*modified_factor(n_current,Texo_tian,Zexo_tian,16.0))
        append!(Cthermal_list, calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:CO][end]*C_CO2_ratio,12.)*C_CO_ratio*flex_effusion_velocity(Texo_tian,12.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,12.0))
        append!(Cnonthermal_list, C_nonthermal)
        append!(Ononthermal_list, O_nonthermal)
        # where the action happens - update n_current for new timesteps
        update!(n_current,dt)
        #n_current = n_current
        # save the concentrations to history
        # write n_current into n_current_mat
        n_internal_mat = Array{Float64}(undef,length(alt)-2,length(collect(keys(n_current))));
        for ispecies in 1:length(collect(keys(n_current)))
            for ialt in 1:length(alt)-2
                n_internal_mat[ialt,ispecies] = n_current[collect(keys(n_current))[ispecies]][ialt]
            end
        end

        if thisi%10 == 0
            #println("ela_time: ", elapsed_time)
            println("totaltime: ", totaltime)
            println("CO2 outgassing: ",CO2outgassing_list[end])
            println("C esc total: ",Cesc_list[end])
            println("C thermal in run: ",Cthermal_list[end])
            println("C nonthermal in run: ",Cnonthermal_list[end])
            #println("Oreturn: ",Oreturn_list[end])
            println("O thermal esc: ", Othermal_list[end])
            println("O nonthermal in run: ",Ononthermal_list[end])
            println("H esc in run: ",Hesc_list[end])
            #println("C effusion: ",eval(speciesbclist[:CO][2,2]))
            plotatm(n_current)
        end

        # write n_current_mat to file
        h5write(filename, string("n_current/iter_",thisi), n_internal_mat)


    end
    return n_current
end


function read_ncurrent_from_file(readfile,tag)
    thisalt = h5read(readfile,"n_current/alt")
    if thisalt != alt
        throw("altitudes in file do not match altitudes in memory!")
    end
    n_current_tag_list = map(Symbol,h5read(readfile,"n_current/species"))
    n_current_mat = h5read(readfile,tag);
    n_current = Dict{Symbol,Array{Float64,1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]]=reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

function get_H_fluxes(readfile)
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)

    Hfluxes=fill(0.,timelength)
    n_current=read_ncurrent_from_file(readfile,string("n_current/init"))
    Hfluxes[1]=(n_current[:H][end]*speciesbcs(:H)[2,2]
                  +2*n_current[:H2][end]*speciesbcs(:H2)[2,2])

    for i in 1:(timelength-1)
        n_current=read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        Hfluxes[i+1]=(n_current[:H][end]*speciesbcs(:H)[2,2]
                      +2*n_current[:H2][end]*speciesbcs(:H2)[2,2])
    end
    Hfluxes
end

#to plot deposition fluxes of O3,H2O2,HO2
function get_all_outfluxes(readfile)
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)

    #create dominant deposition specise flux arrays
    Hfluxes=fill(0.,timelength)
    H2O2_dep=fill(0.,timelength)
    O3_dep=fill(0.,timelength)
    HO2_dep=fill(0.,timelength)
    H_dep=fill(0., timelength)
    totalH_outflux = fill(0.,timelength)
    totalO_outflux = fill(0.,timelength)

    n_current=read_ncurrent_from_file(readfile,string("n_current/init"))
    Hfluxes[1]=(n_current[:H][end]*eval(speciesbcs(:H)[2,2])
                  +2*n_current[:H2][end]*eval(speciesbcs(:H2)[2,2]))
    H2O2_dep[1] = n_current[:H2O2][1]*speciesbclist[:H2O2][1,2]
    O3_dep[1] = n_current[:O3][1]*speciesbclist[:O3][1,2]
    HO2_dep[1] = n_current[:HO2][1]*speciesbclist[:HO2][1,2]
    H_dep[1] = n_current[:H][1]*speciesbclist[:H][1,2]
    totalH_outflux[1] = Hesc_list[1]+2*H2O2_dep[1]+HO2_dep[1]+H_dep[1] + n_current[:OH][1]*speciesbclist[:OH][1,2]
    totalO_outflux[1] = 2*H2O2_dep[1] + 3*O3_dep[1] + 2*HO2_dep[1]  + n_current[:OH][1]*speciesbclist[:OH][1,2] +n_current[:O][1]*speciesbclist[:O][1,2] +
                        n_current[:O1D][1]*speciesbclist[:O1D][1,2] - 2*CO2outgassing_list[1] - Oreturn_list[1] + 2*Cesc_list[1] #koyama 2020/06/18 changed

    for i in 1:(timelength-1)
        n_current=read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        Hfluxes[i+1]=(n_current[:H][end]*eval(speciesbcs(:H)[2,2])
                      +2*n_current[:H2][end]*eval(speciesbcs(:H2)[2,2]))
        H2O2_dep[1+i] = n_current[:H2O2][1]*speciesbclist[:H2O2][1,2]
        O3_dep[1+i] = n_current[:O3][1]*speciesbclist[:O3][1,2]
        HO2_dep[1+i] = n_current[:HO2][1]*speciesbclist[:HO2][1,2]
        H_dep[1+i] = n_current[:H][1]*speciesbclist[:H][1,2]
        totalH_outflux[i+1] = Hesc_list[i+1]+2*H2O2_dep[i+1]+HO2_dep[i+1] + H_dep[i+1]+n_current[:OH][1]*speciesbclist[:OH][1,2]
        totalO_outflux[i+1] = 2*H2O2_dep[i+1] + 3*O3_dep[i+1] + 2*HO2_dep[i+1] + n_current[:OH][1]*speciesbclist[:OH][1,2] +n_current[:O][1]*speciesbclist[:O][1,2] +
                            n_current[:O1D][1]*speciesbclist[:O1D][1,2] - 2*CO2outgassing_list[i+1] - Oreturn_list[i+1] + 2*Cesc_list[i+1] #koyama 2020/06/18 changed
    end
    (Hfluxes, H2O2_dep, O3_dep, HO2_dep, H_dep, totalH_outflux, totalO_outflux)
end

#ファイルを読み込んで, reactionratesとfluxes関数を使って入れて行くだけ
function get_rates_and_fluxes(readfile)
    mydset = h5open(readfile,"r")
    mydata = read(mydset)
    timelength = length(mydata["n_current"]["timelist"])+1
    close(mydset)
    reactionrateshist = fill(convert(Float64,NaN),timelength,length(intaltgrid),length(reactionnet))
    fluxhist = fill(convert(Float64,NaN),timelength,length(intaltgrid),length(specieslist))
    n_current = read_ncurrent_from_file(readfile,string("n_current/init"))
    reactionrateshist[1,:,:] = reactionrates(n_current)
    fluxhist[1,:,:] = fluxes(n_current,dz)
    for i in 1:(timelength-1)
        n_current = read_ncurrent_from_file(readfile,string("n_current/iter_",i))
        reactionrateshist[i+1,:,:] = reactionrates(n_current)
        fluxhist[i+1,:,:] = fluxes(n_current,dz)
    end
    (reactionrateshist,fluxhist)
end

function get_all_rates_and_fluxes(readfile)
    (reactionrateshist,fluxhist)=get_rates_and_fluxes(readfile)
    h5write(readfile,"fluxes/flux_history",fluxhist)
    h5write(readfile,"rates/reaction_rates_history",reactionrateshist)
    return
end

#pmap(x->println(string("parmsvec[i][1]=",x[1],", parmsvec[i][2]=",x[2],", filename=",x[3])),[[p,f;] for (p,f) in zip(parmsvec,filenamevec)])

##This runs the simulation for all added ppms and altitudes
#pmap( x->runwaterprofile(n_current,x[1],x[2],timediff,x[3]),[[p,f;] for (p,f) in zip(parmsvec,filenamevec)])
println("start to run!")
result = runprofile(n_current, timediff, run_filename)
println("finished running then start to make files for plot")
get_all_rates_and_fluxes(run_filename)


#This gets deposition get_H_fluxes
(Hfluxes,H2O2dep,O3dep,HO2dep,Hdep,totalH_outflux,totalO_outflux)=get_all_outfluxes(run_filename)
HOratio = totalH_outflux./totalO_outflux


# write out the H fluxes and deposition fluxes =======================================================
println("Writing escape and depostion fluxes file")
h5open(hfile, isfile(hfile) ? "r+" : "w") do file
   #write(file,"fluxes/fluxvals",Hfluxes)
   write(file,"fluxes/fluxvals",Hesc_list)
   write(file,"fluxes/times",h5read(run_filename,"n_current/timelist"))
   write(file, "depositions/H2O2",H2O2dep)
   write(file,"depositions/O3",O3dep)
   write(file,"depositions/HO2",HO2dep)
   write(file, "depositions/H", Hdep)
   write(file,"outfluxes/H",totalH_outflux)
   write(file,"outfluxes/O",totalO_outflux)
   write(file,"HOratio",HOratio)
   write(file, "outfluxes/CO2out",CO2outgassing_list)
   write(file, "outfluxes/Oreturn",Oreturn_list)
   write(file, "outfluxes/CO2esc",Cesc_list)
   write(file, "outfluxes/Hesc",Hesc_list)
   write(file, "outfluxes/Othermal",Othermal_list)
   write(file, "outfluxes/Cthermal",Cthermal_list)
   write(file, "outfluxes/Ononthermal",Ononthermal_list)
   write(file, "outfluxes/Cnonthermal",Cnonthermal_list)

   #write(file,"waterprofs/ppm",writewaterprof)
   #write(file,"waterprofs/alt",alt[2:end-1])
end

timelist = cumsum(timediff)
regtimeindex_rev = findfirst(x ->((x<1.99)|(2.01<x)), reverse(HOratio)) -1
regtimeindex = length(HOratio)-regtimeindex_rev + 1
println(run_filename)

write_ncurrent(n_current, output_folder*"/atmos_final.h5")

#redox balance
Closs = 4*Cesc_list[end]
Hloss = Hesc_list[end]
Oloss = 2*(Othermal_list[end]+Ononthermal_list[end]+3*O3dep[end]+H2O2dep[end])
println("C/O: ",Closs/Oloss)
println("H/O: ",Hloss/Oloss)
println("C&H/O: ",(Hloss+Closs)/Oloss)

println("done!")
