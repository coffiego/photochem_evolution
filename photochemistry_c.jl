# Originally
# This julia code contains a photochemical model of the Martian
# atmosphere prepared by Mike Chaffin, originally for the 8th
# International Conference on Mars, held at CalTech in July 2014.

#2020/6/6 Koyama
#Investigate the effect of C escape on CO2 main atmospheres
#Assuming boundary conditions of Tian+2009 from 4.1 Ga
#This is a following work of Yagi's M thesis 2020

####################################################################################



####################################################################################
#here is parameters to investigate regulation regulationtimescale
CO2_pressure = 1#mbar
surface_temperature = 210#K
Oxygen_escape_rate = 1.2e8#cm-2s-1
Oxygen_escape_rate_before = 1.2e8 #cm-2s-1
c_escape_rate = 0.
get_converged = false#to get converged: true, to see timechange from converged: false
depv = 0.02
N2_Ar = false #if N2 and Ar are added or not.
feuv=15.5
#Ttropo = surface_temperature - 84.0
folder_directory = "/Users/shungokoyama/programming/photochem_evolution/"
global totaltime=0.
dt = 0.

###ALL distance/area/volume units must be in CM!!!!

using PyPlot
using HDF5, JLD
using DelimitedFiles
using SparseArrays
using LinearAlgebra

################################################################################
################################# Species LIST #################################
################################################################################

#array of symbols for each species
#const fullspecieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
#                         :O1D, :H,:N2,:Ar,:CO2pl,:HOCO];
const fullspecieslist = [:CO2, :O2, :O3, :H2, :OH, :HO2, :H2O, :H2O2, :O, :CO,
                        :O1D, :H,:N2,:Ar,:CO2pl,:HOCO];
specieslist=fullspecieslist;
const nochemspecies = [:N2, :Ar, :CO2pl,:H2O];
const chemspecies = setdiff(specieslist, nochemspecies);
const notransportspecies = [:CO2pl,:H2O];
const transportspecies = setdiff(specieslist, notransportspecies);
const speciesmolmasslist = Dict(:CO2=>44, :O2=>32, :O3=>48, :H2=>2, :OH=>17,
                                :HO2=>33, :H2O=>18, :H2O2=>34, :O=>16, :CO=>28,
                                :O1D=>16, :H=>1, :N2=>28, :Ar=>40, :CO2pl=>44,
                                :HOCO=>45)

function fluxsymbol(x)
    Symbol(string("f",string(x)))
end
const fluxlist = map(fluxsymbol, fullspecieslist)

# array of species for which photolysis is important. All rates should
# start with J, end with a lowercase letter, and contain a species in
# specieslist above, which is used to compute photolysis. Different
# letters at the end correspond to different products of photolytic
# destruction or to photoionization.
const Jratelist=[:JCO2ion,:JCO2toCOpO,:JCO2toCOpO1D,:JO2toOpO,:JO2toOpO1D,
                 :JO3toO2pO,:JO3toO2pO1D,:JO3toOpOpO,:JH2toHpH,:JOHtoOpH,
                 :JOHtoO1DpH,:JHO2toOHpH,:JH2OtoHpOH,:JH2OtoH2pO1D,:JH2OtoHpHpO,
                 :JH2O2to2OH,:JH2O2toHO2pH,:JH2O2toH2OpO1D];


################################################################################
###################### Load Converged Test Case from File ######################
################################################################################

# the test case was created by hand by Mike Chaffin and saved for automated use.
# Change following line as needed depending on local machine
# readfile is changed according to the condition of "get_converged" koyama
if get_converged
    readfile =  folder_directory * "initial/" *
                "initial_condition_CO2_" * string(CO2_pressure) * "mbar.h5"
    #readfile = "/Users/shungokoyama/Programming/result/koyama/initial_condition_CO2_1bar.h5"
else
    #readfile = folder_directory * "conv_Tsurf"*string(surface_temperature)*"_Oesc"*string(Oxygen_escape_rate_before)*
    #            "_Tmod/converged_Tsurf_" * string(surface_temperature) * "_CO2_" *
    #            string(CO2_pressure) * "mbar_Oesc_" * string(Oxygen_escape_rate_before) *
    #            "_depv_" * string(depv) * "_nofixedCO2_noHOCO_Tmod_1billion.h5"
    #readfile = "/Users/shungokoyama/Programming/result/koyama/converged_Tsurf_250_CO2_1bar_Oesc_1.2e8_depv_0.02_nofixedCO2.h5"
    #readfile = "/Users/shungokoyama/Programming/result/stat/converged_Tsurf_" * string(surface_temperature) * "_CO2_" *
    #    string(CO2_pressure) * "mbar_Oesc_" * string(Oxygen_escape_rate_before) *
    #    "_depv_" * string(depv) * "_nofixedCO2_noHOCO_Tmod_1billion.h5"
    readfile = folder_directory * "converged_Tsurf_" * string(surface_temperature) * "_CO2_" *
        string(CO2_pressure) * "mbar_Oesc_" * string(Oxygen_escape_rate_before) *
        "_depv_" * string(depv) * "_1billion_ribas_4.1Ga.h5"
end
#const alt=h5read(readfile,"n_current/alt")

#here should be changed
#1. start from 1 mbar
readfile = folder_directory * "converged_Tsurf_" * string(surface_temperature) * "_CO2_" *
    string(CO2_pressure) * "mbar_Oesc_" * string(Oxygen_escape_rate_before) *
    "_depv_" * string(depv) * "_Tinf_480.h5"

#2. start from the previous euv calculation result
#readfile = folder_directory*"Tinf480_lowestCO2_euv13.5_1526K_flexo_1mbar/atmos_final.h5"

#readfile = folder_directory*"Tinf480_lowestCO2_euv11/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv12/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv12.5_1254K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv13_1382K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv13.5_1526K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv14_1670K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv14.5_1826K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv15_2026K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowestCO2_euv15.5_2252K/atmos_final.h5"

#readfile = folder_directory*"Tinf480_middleCO2_euv10/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv11/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv12/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv12.5_1254K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv13_1382K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv13.5_1526K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv14_1670K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv14.5_1826K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv15_2026K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_middleCO2_euv15.5_2252K/atmos_final.h5"

#readfile = folder_directory*"Tinf480_lowerCO2_euv13_1382K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowerCO2_euv13.5_1526K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowerCO2_euv14_1670K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowerCO2_euv14.5_1826K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowerCO2_euv15_2026K/atmos_final.h5"
#readfile = folder_directory*"Tinf480_lowerCO2_euv15.5_2252K/atmos_final.h5"


function get_ncurrent(readfile)
    alt=h5read(readfile,"n_current/alt")
    n_current_tag_list=map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat=h5read(readfile,"n_current/n_current_mat");
    n_current=Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]]=reshape(n_current_mat[:,ispecies],length(alt)-2)
    end
    n_current
end

#n_current=get_ncurrent(readfile)

function write_ncurrent(n_current, filename)
    n_current_mat=Array{Float64}(undef,length(alt)-2, length(collect(keys(n_current))));
    for ispecies in [1:length(collect(keys(n_current)));]
        for ialt in [1:length(alt)-2;]
            n_current_mat[ialt, ispecies]=n_current[collect(keys(n_current))[ispecies]][ialt]
        end
    end
    h5write(filename,"n_current/n_current_mat",n_current_mat)
    h5write(filename,"n_current/alt",alt)
    h5write(filename,"n_current/species",map(string, collect(keys(n_current))))
end

################################################################################
# finding first value which is equal to val in arr, then return its index
# if it is not found, returns 0

# inputs:
# val: value you want to find in the array
# arr: the array in which the value is searched
# output:
# the index of first value in the array which is equal to val
# or 0 (if not found)
function findfirst_or_zero(val,arr)
    output = findfirst(isequal(val),arr)
    if output === nothing
        return 0
    else
        return output
    end
end

################################################################################
############################### REACTION NETWORK ###############################
################################################################################

#function to replace three body rates with the recommended expression
threebody(k0, kinf) = :($k0*M/(1+$k0*M/$kinf)*0.6^((1+(log10($k0*M/$kinf))^2)^-1))

threebodyca(k0, kinf) = :($k0/(1+$k0/($kinf/M))*0.6^((1+(log10($k0/($kinf*M)))^2)^-1))

reactionnet=[
             #Photodissociation
             [[:CO2],[:CO,:O],:JCO2toCOpO],
             [[:CO2],[:CO,:O1D],:JCO2toCOpO1D],
             [[:O2],[:O,:O],:JO2toOpO],
             [[:O2],[:O,:O1D],:JO2toOpO1D],
             [[:O3],[:O2,:O],:JO3toO2pO],
             [[:O3],[:O2,:O1D],:JO3toO2pO1D],
             [[:O3],[:O,:O,:O],:JO3toOpOpO],
             [[:H2],[:H,:H],:JH2toHpH],
             [[:OH],[:O,:H],:JOHtoOpH],
             [[:OH],[:O1D,:H],:JOHtoO1DpH],
             [[:HO2],[:OH,:O],:JHO2toOHpH], # other branches should be
                                            # here, but have not been
                                            # measured
             [[:H2O],[:H,:OH],:JH2OtoHpOH],
             [[:H2O],[:H2,:O1D],:JH2OtoH2pO1D],
             [[:H2O],[:H,:H,:O],:JH2OtoHpHpO],
             [[:H2O2],[:OH,:OH],:JH2O2to2OH],
             [[:H2O2],[:HO2,:H],:JH2O2toHO2pH],
             [[:H2O2],[:H2O,:O1D],:JH2O2toH2OpO1D],

             # recombination of O
             [[:O,:O,:M],[:O2,:M],:(1.8*3.0e-33*(300/T)^3.25)],
             [[:O,:O2,:N2],[:O3,:N2],:(5e-35*exp(724. /T))],
             [[:O,:O2,:CO2],[:O3,:CO2],:(2.5*6.0e-34*(300. /T)^2.4)],
             [[:O,:O3],[:O2,:O2],:(8.0e-12*exp(-2060. /T))],
             [[:O,:CO,:M],[:CO2,:M],:(2.2e-33*exp(-1780. /T))],

             # O1D attack
             [[:O1D,:O2],[:O,:O2],:(3.2e-11*exp(70. /T))],
             [[:O1D,:O3],[:O2,:O2],:(1.2e-10)],
             [[:O1D,:O3],[:O,:O,:O2],:(1.2e-10)],
             [[:O1D,:H2],[:H,:OH],:(1.2e-10)],
             [[:O1D,:CO2],[:O,:CO2],:(7.5e-11*exp(115. /T))],
             [[:O1D,:H2O],[:OH,:OH],:(1.63e-10*exp(60. /T))],

             # loss of H2
             [[:H2,:O],[:OH,:H],:(6.34e-12*exp(-4000/T))],
             [[:OH,:H2],[:H2O,:H],:(9.01e-13*exp(-1526/T))],

             # recombination of H
             [[:H,:H,:CO2],[:H2,:CO2],:(1.6e-32*(298. /T)^2.27)],
             [[:H,:OH,:CO2],[:H2O,:CO2],:(1.9*6.8e-31*(300. /T)^2)],
             [[:H,:HO2],[:OH,:OH],:(7.2e-11)],
             [[:H,:HO2],[:H2O,:O1D],:(1.6e-12)],#O1D is theoretically mandated
             [[:H,:HO2],[:H2,:O2],:(0.5*6.9e-12)],
             [[:H,:H2O2],[:HO2,:H2],:(2.8e-12*exp(-1890/T))],
             [[:H,:H2O2],[:H2O,:OH],:(1.7e-11*exp(-1800/T))],

             # Interconversion of odd H
             [[:H,:O2],[:HO2],threebody(:(2.0*4.4e-32*(T/300)^-1.3),
                                        :(7.5e-11*(T/300)^0.2))],
             [[:H,:O3],[:OH,:O2],:(1.4e-10*exp(-470/T))],
             [[:O,:OH],[:O2,:H],:(1.8e-11*exp(180/T))],
             [[:O,:HO2],[:OH,:O2],:(3.0e-11*exp(200. /T))],
             [[:O,:H2O2],[:OH,:HO2],:(1.4e-12*exp(-2000/T))],
             [[:OH,:OH],[:H2O,:O],:(1.8e-12)],
             [[:OH,:OH],[:H2O2],threebody(:(1.3*6.9e-31*(T/300)^-1.0),:(2.6e-11))],
             [[:OH,:O3],[:HO2,:O2],:(1.7e-12*exp(-940/T))],
             [[:OH,:HO2],[:H2O,:O2],:(4.8e-11*exp(250/T))],
             [[:OH,:H2O2],[:H2O,:HO2],:(1.8e-12)],
             [[:HO2,:O3],[:OH,:O2,:O2],:(1.0e-14*exp(-490/T))],
             [[:HO2,:HO2],[:H2O2,:O2],:(3.0e-13*exp(460/T))],
             [[:HO2,:HO2,:M],[:H2O2,:O2,:M],:(2*2.1e-33*exp(920/T))],

             # CO2 recombination due to odd H (with HOCO intermediate)
             [[:CO,:OH],[:CO2,:H],threebodyca(:(1.5e-13*(T/300)^0.6),:(2.1e9*(T/300)^6.1))],
             #[[:OH,:CO],[:HOCO],threebody(:(5.9e-33*(T/300)^-1.4),:(1.1e-12*(T/300)^1.3))], #koyama delete it for early mars case
             [[:OH,:CO],[:HOCO],:(0.0)], #koyama delete it for early mars case
             #[[:HOCO,:O2],[:HO2,:CO2],:(2.0e-12)], #koyama deleted early mars case
             [[:HOCO,:O2],[:HO2,:CO2],:(0.0)],  #koyama deleted early mars case

             # CO2+ attack on molecular hydrogen
             [[:CO2pl,:H2],[:H,:H],:(8.7e-10)] #I delete CO2+ to introduce no CO2 fixed BD
             ]

################################################################################
############################# FUNDAMENTAL CONSTANTS ############################
################################################################################

# fundamental constants
const boltzmannK=1.38e-23;    # J/K
const bigG=6.67e-11;          # N m^2/kg^2
const mH=1.67e-27;            # kg

# mars parameters
const marsM=0.1075*5.972e24;  # kg
const radiusM=3396e5;         # cm

################################################################################
################################ TEMPERATURE ###################################
################################################################################

function Tspl(z::Float64, lapserate=-1.4e-5, Tsurf=211, ztropo=50e5, zexo=200e5, Texo=300) #change koyama
    # DO NOT MODIFY! If you want to change the temperature, define a
    # new function or select different arguments and pass to Temp(z)

    # spline interpolation between adaiabatic lapse rate in lower atm
    # (Zahnle 2008) and isothermal exosphere (300K, typical?)  Too hot
    # at 125km?? should be close to 170 (Krasnopolsky)
    if z < ztropo
        return Tsurf + lapserate*z
    end
    if z > zexo
        return Texo
    end
    Ttropo = Tsurf+lapserate*ztropo
    Tret = ((z - zexo)^2 * (Ttropo*(2*z + zexo - 3*ztropo) +
            lapserate*(z - ztropo)*(zexo - ztropo))
            - Texo*(z - ztropo)^2*(2*z - 3*zexo + ztropo))/(zexo - ztropo)^3
    # I got the above expression from Mathematica; seems to do what I
    # asked it to, matching the temp and lapse rate at the tropopause
    # and the temp with no derivative at the exobase.
end

#it's not used here
function Tpiecewise(z::Float64, Tinf=240.0, Ttropo=Ttropo, lapserate=-1.4e-5, ztropo=90e5, ztropowidth=30e5) #Koyama flexible
    # DO NOT MODIFY! If you want to change the temperature, define a
    # new function or select different arguments and pass to Temp(z)

    # a piecewise function for temperature as a function of altitude,
    # using Krasnopolsky's 2010 temperatures for altitudes
    # >htropo=90km, fixed at Ttropo=125K between htropo and
    # htropo-htropowidth=60km, and rising at a constant lapse rate
    # (1.4K/km) below.
    if z >= ztropo
        return Tinf - (Tinf - Ttropo)*exp(-((z-ztropo)^2)/(8e10*Tinf))
    end
    if ztropo > z >= ztropo - ztropowidth
        return Ttropo
    end
    if ztropo-ztropowidth > z
        return Ttropo-lapserate*(ztropo-ztropowidth-z)
    end
end

############koyama Temperature#######################
const lapse = -1.4e-5 #K/km
CO2mm = speciesmolmasslist[:CO2]
function calcrho(T,P)
    rho = P*1e2/(boltzmannK*T)*1e-6
    return rho
end
function scaleH(z, T::Float64, mm::Real)
    boltzmannK*T/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
    # constants are in MKS. Convert to m and back to cm.
end

###def####
dz_temp = 1.0e5
P_array =Float64[]
T_array =Float64[]
rho_array=Float64[]
z_array=Float64[]
P = CO2_pressure
Ts = surface_temperature*1.0
T = Ts
Ttropo = 125.0 #K lower boundary of the thermosphere
#Tinf = 240.0 #K
Tinf = 480.
z = 0.0
rho = calcrho(T,P)

### EUV temperature ####
if feuv==16.5
    start_time=0.75
    Texo_tian = 2777.
elseif feuv==16
    start_time=0.77
    Texo_tian = 2503.
elseif feuv == 15.5
    start_time = 0.789
    Texo_tian = 2252.
elseif feuv==15.
    start_time=0.806
    Texo_tian = 2026.
elseif feuv==14.5
    start_time=0.824
    Texo_tian = 1826.
elseif feuv==14.
    start_time=0.842
    Texo_tian = 1670.
elseif feuv==13.5
    start_time=0.860
    Texo_tian = 1526.
elseif feuv==13.
    start_time=0.879
    Texo_tian = 1382.
elseif feuv==12.5
    start_time=0.898
    Texo_tian =1254.
elseif feuv==12.
    start_time=0.915
    Texo_tian = 1154.
elseif feuv==11.
    start_time=0.958
    Texo_tian = 959.
elseif feuv==10.
    start_time=1.0
    Texo_tian = 810.
elseif feuv==8.
    start_time=1.18
    Texo_tian = 588.
elseif feuv==5.
    start_time=1.23 #rough estimate
    Texo_tian = 387.
    Tinf = 320.
elseif feuv==3.
    start_time=1.3 #rough estimate
    Texo_tian = 305.
    Tinf = 270.
else
    start_time = 0.5
    Texo_tian = 3000.
end

################################################################################
#Calc rho from the bottom up till rho = rho(exo)=3.5e12 cm^-3
#T decreases by lapse rate of -1.4K/km till T=Ttropo(125K)
while(rho>3.5e12)
    global z,P,Ts,T,Ttropo,Tinf,rho,dz_temp,lapse
    append!(z_array, z)
    append!(P_array,P)
    append!(T_array,T)
    append!(rho_array,rho)
    z = z+dz_temp
    T = Ts + lapse*z
    if T<Ttropo
        T = Ttropo
    end
    P = P*exp(-1/scaleH(z,T,CO2mm)*dz_temp) #using scaleH to calc P
    rho = calcrho(T,P)
end

# alt is composed of even numbers spaced by 2km
if z%2e5==0.
    const z_upper_tropo = z
else
    const z_upper_tropo = z-dz_temp
end
if z_array[findfirst(T_array.==Ttropo)]%2e5 == 0
    const z_lower_tropo = z_array[findfirst(T_array.==Ttropo)]
else
    const z_lower_tropo = z_array[findfirst(T_array.==Ttropo)] - dz_temp
end

z = z-dz_temp
z_CO2pl_index = Int((z_upper_tropo-1e6)/2e5)+1 #lower boundary of CO2+ is below z_upper_tropo by 10km
#z_CO2pl_index = Int(floor(z_array[findfirst(rho_array.<1.56e13)]/2e5))+1#position of CO2+ lower boundary
zexo = z_upper_tropo + 1.1e7


########to calc T above the bottom of thermosphere.#####
for i in 1:((zexo-z)/dz_temp)
    global z
    z=z+dz_temp
    #T = T + (Tinf-Ttropo)/(zexo-z_upper_tropo)*dz_temp
    #T = Tinf - (Tinf - Ttropo)*exp(-((z-z_upper_tropo)^2)/(8e10*Tinf))
    T = Tinf - (Tinf - Ttropo)*exp(-((z-z_upper_tropo)^2)/(4e10*Tinf)) #change 8e10->4e10
    append!(z_array,z)
    append!(T_array,T)
end


function Tflex(z::Float64, T_array)
    return T_array[Int(z/1e5)+1]
end
###########################koyama original temperature#############################
#If changes to the temperature are needed, they should be made here
#Temp(z::Float64) = Tpiecewise(z)#originally Tpiecewise koyama
Temp(z::Float64) = Tflex(z,T_array)

################################################################################
########################### discretization parameters ##########################
################################################################################

#set grid to adapt to different atmospheric pressure
const alt = z_array[1:2:end]
const zmin=alt[1]
const zmax=alt[end];
const dz=alt[2]-alt[1];

################################################################################
############## additional CO2 when you get a converged state###################
if get_converged
    n_current_tag_list=map(Symbol, h5read(readfile,"n_current/species"))
    n_current_mat=h5read(readfile,"n_current/n_current_mat");
    n_current=Dict{Symbol, Array{Float64, 1}}()
    for ispecies in [1:length(n_current_tag_list);]
        n_current[n_current_tag_list[ispecies]]= zeros(length(alt)-2)
    end
    if length(alt)==101 #exobase is at 200km like original
        n_current[:CO2][1:end] = n_current_mat[:,12]
    elseif length(alt)>101 #exobase is higher than 200km
        n_current[:CO2][1:99] = n_current_mat[:,12]
        for i in 100:(length(alt)-2)
            n_current[:CO2][i] = n_current[:CO2][i-1]*exp(-1/scaleH(alt[i-1],Temp(alt[i-1]),CO2mm)*dz)
        end
    else#zexo is lower than 200km
        n_current[:CO2][1:end] = n_current_mat[1:(length(alt)-2),12]
    end

else
    n_current=get_ncurrent(readfile)
end



################################################################################
############################# AUX DENSITY FUNCTIONS ############################
################################################################################

n_alt_index=Dict([z=>clamp((i-1),1, length(alt)-2) for (i, z) in enumerate(alt)])
# used in combination with n_current. Gets the index corresponding to a given altitude

function n_tot(n_current, z)
    # get the total number density at this altitude
    thisaltindex = n_alt_index[z]
    return sum( [n_current[s][thisaltindex] for s in specieslist] )
end

# mean molecular mass at a given altitude
function meanmass(n_current, z)
    #find the mean molecular mass at a given altitude
    thisaltindex = n_alt_index[z]

    c = [n_current[sp][thisaltindex] for sp in specieslist]
    m = [speciesmolmasslist[sp] for sp in specieslist]

    return sum(c.*m)/sum(c)
end


################################################################################
############################# BOUNDARY CONDITIONS ##############################
################################################################################

# water saturation vapor pressure, T in K, Psat in mmHg
# (from Washburn 1924)
Psat(T::Float64)=10^(-2445.5646/T+8.2312*log10(T)-0.01677006*T+1.20514e-5*T^2-6.757169)


#######koyama version H2O profile#################################
#this H2O profile is 9.5 pr micrometer
#z[1:15](<30km) を相対湿度22%にする
#全て飽和水蒸気を混合率に直し, 最小になる高度より上の高度を全て最小の飽和水蒸気の値にする
#30km~60km(飽和水蒸気量が最小になる高度)はlogで繋ぐ。
# convert to H2O satuation vapor pressure:
# 133.3 N/m2 = 1 mmHg, divided by k*T=J=N*m = 1/m3, convert to 1/cm3... 1/cm3 = 10^6/m3
H2Osat = map(x->(1e-6*133.3*Psat(x))/(boltzmannK*x),map(Temp, alt))
H2Osat_lowalt=H2Osat[1:15]*0.22 #15, 0.22
H2Osatfrac_lowalt = H2Osat_lowalt./map(z->n_tot(n_current, z),alt[1:15])
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z),alt)
H2Oinitfrac = H2Osatfrac[1:findfirst_or_zero(minimum(H2Osatfrac),H2Osatfrac)]
##log10(fraction of H2O=y) = aa * z + bb で繋げる(30-60kmのH2O profile)

H2Osatfrac_midalt = Float64[]
# this connects from 30km to z_lower_tropo by "log scale"
aa = log10(H2Osatfrac_lowalt[end]/minimum(H2Osatfrac)) / (length(H2Osatfrac_lowalt)-length(H2Oinitfrac))
bb = log10(H2Osatfrac_lowalt[end]) - aa*length(H2Osatfrac_lowalt)
for i in [length(H2Osatfrac_lowalt)+1:length(H2Oinitfrac)-1;]
    append!(H2Osatfrac_midalt, 10^(aa*i+bb))#y = 10^(aa*z + bb)
end


#[0-30km](H2Osatの22%) + [30km-z_lower_tropo km](y = 10^(aa*z + bb))+ [60-200km](H2Osatの最小値)
H2Ofrac = vcat(H2Osatfrac_lowalt, H2Osatfrac_midalt, minimum(H2Oinitfrac)*ones(length(alt)-length(H2Oinitfrac)-1))
#飽和水蒸気量を超えていないか下でチェックする。
#for i in [1:length(H2Oinitfrac);]
#    H2Ofrac[i] = H2Ofrac[i] < H2Osatfrac[i+1] ? H2Ofrac[i] : H2Osatfrac[i+1]
#end

#clf()
##plot(H2Ofrac,z_array[3:2:end-1])
#xscale("log")

#############koyama version H2O profile end######################

############Chaffin version H2O profile#######################
#=
# convert to H2O satuation vapor pressure:
# 133.3 N/m2(Pa) = 1 mmHg, ここでmmHgからPaに変える
#次に圧力を密度に変えるP=rho*k*Tよりrho=P/(kT)
#最後にk*T=J=N*m = 1/m3, convert to 1/cm3... 1/cm3 = 10^6/m3なので1/cm3に変えるために1e-6をかける
H2Osat = map(x->(1e-6*133.3*Psat(x))/(boltzmannK*x),map(Temp, alt))#これで温度分布に対応した飽和水蒸気密度
H2Osatfrac = H2Osat./map(z->n_tot(n_current, z),alt)#mixing ratioに変換
H2Oinitfrac = H2Osatfrac[1:findfirst(H2Osatfrac, minimum(H2Osatfrac))]#mixing ratioで一番小さくなるところまでを取り出す
H2Oinitfrac = [H2Oinitfrac, fill(minimum(H2Osatfrac),
               length(alt)-2-length(H2Oinitfrac));]#mixing ratioが最小になった高度より上を全て最小の値にして200kmまで追加(cold trapの効果)
H2Oinitfrac[find(x->x<30e5, alt)]=1e-4#30km以下の高度のmixing ratioを1e-4にする, find function is not valid after ver0.7 julia
for i in [1:length(H2Oinitfrac);]
    H2Oinitfrac[i] = H2Oinitfrac[i] < H2Osatfrac[i+1] ? H2Oinitfrac[i] : H2Osatfrac[i+1] #a ? b : cの用法　a=trueだったらb, a=falseでcを返す
    #つまり、一つ上の高度の混合比より小さければそのままで、大きければ一つ上の高度の混合比と同じにする<-なぜ？小山
    #ここは不明だが、高度の高い方から行うと不連続がなくなるので本当はそれがやりたかったのでは？
end
=#
##################################################################


#=
H2Olowfrac = H2Osatfrac[1:findfirst(H2Osatfrac, minimum(H2Osatfrac))]
H2Olowfrac = [H2Olowfrac, fill(minimum(H2Osatfrac),length(alt)-2-length(H2Olowfrac));]
H2Olowfrac[find(x->x<30e5, alt)]=0.5e-4
for i in [1:length(H2Olowfrac);]
    H2Olowfrac[i] = H2Olowfrac[i] < H2Osatfrac[i+1] ? H2Olowfrac[i] : H2Osatfrac[i+1]
end

H2Onohighalt = deepcopy(H2Oinitfrac)
H2Onohighalt[find(x->x>30e5, alt[2:end-1])]=0.0

# add in a detached water vapor layer, which looks kinda like a
# gaussian packet floating at 60km (Maltagliati 2013)
detachedlayer = 1e-6*map(x->80.*exp(-((x-60)/12.5)^2),alt[2:end-1]/1e5)+H2Oinitfrac
=#

# this computes the total water column in precipitable microns
# n_col(molecules/cm2)/(molecules/mol)*(gm/mol)*(1cc/gm) = (cm)*(10^4μm/cm)
# = precipitable μm
#(sum(([1e-4, H2Oinitfrac;]).*map(z->n_tot(n_current, z),alt[1:end-2]))*dz)/6.02e23*18*1e4)
#(sum(([1e-4, detachedlayer;]).*map(z->n_tot(n_current, z),alt[1:end-1]))*dz)/6.02e23*18*1e4
#(sum(([1e-4, H2Ofrac;]).*map(z->n_tot(n_current, z),alt[1:end-2]))*dz)/6.02e23*18*1e4)


# effusion velocities need to be computed from the atmospheric temperature at
# the upper boundary:
function effusion_velocity(Texo::Float64, m::Float64)
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+zmax))
    ## println("lambda = ",lambda)
    vth=sqrt(2*boltzmannK*Texo/(m*mH))
    ## println("vth = ", vth)
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)
    return v
end

#H_effusion_velocity = effusion_velocity(Temp(zmax),1.0)
#H2_effusion_velocity = effusion_velocity(Temp(zmax),2.0)

### C escape ###################################################################
# flexibile c effsion velocity
# exobase can be changed
# to check exobase does not have a large impact
function flex_effusion_velocity(Texo::Float64, m::Float64, Zexobase::Float64)
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+Zexobase))
    ## println("lambda = ",lambda)
    vth=sqrt(2*boltzmannK*Texo/(m*mH))
    ## println("vth = ", vth)
    v = 1e2*exp(-lambda)*vth*(lambda+1)/(2*pi^0.5)
    return v
end


#T profile above the thermosphere (previous ver)
#function T_thermos(z,Texo,Tinf,zexo_pre)
#    return Texo - (Texo - Tinf)*exp(-((z-zexo_pre)^2)/(3e11*Texo))
#end

#it's approximate to Kaori-san's model
function T_thermos(z,Texo,Tinf,zexo_pre)
    return Texo - (Texo - Tinf)*exp(-((z-zexo_pre))/(1e4*Texo))
end


#to estimate exobase
# altitude where mean free path = scale Height
function estimate_exobase(Texo,Tinf,rhoexo,sp_mass)
    zs=[]
    rhos=[]
    mfps=[]
    Hs=[]
    Ts=[]
    z = alt[end]
    zexo_pre = alt[end]
    zexo_max = 1.e9 #should be lower thatn 10000 km
    dz_temp = 1.e5 #1km
    rho = rhoexo
    xs = 3.e-15 #Oxygen atom collision cross section cm^2 from Kasting and Catling,2017, p137

    #to find exobase -> min(|H-mfp|)
    min_distance = 1.e8 #initial large value
    zexo_new = 0.
    zexo_index=1
    rho_exo_new = 0.
    for i in 1:((zexo_max-z)/dz_temp)
        z=z+dz_temp
        T = T_thermos(z,Texo,Tinf,zexo_pre) #TO do change!!!
        H = scaleH(z,T,sp_mass) #cm
        rho = rho*exp(-1/H*dz_temp) #cm-3
        mfp = 1/(rho*xs) #mean free path [cm]
        append!(zs,z)
        append!(rhos,rho)
        append!(Hs,H)
        append!(mfps,mfp)
        append!(Ts,T)
        #find minimun value of |H-mfp| = exobase
        if abs(H-mfp) < min_distance
            min_distance = abs(H-mfp)
            zexo_new = z
            rho_exo_new = rho
            zexo_index=i
        end
    end
    figure()
    plot(Hs,zs/1e5,label="scale Height")
    plot(mfps,zs/1e5,label="mean free path")
    figure()
    Tss = vcat(T_array,Ts[1:Int(zexo_index)])
    zss = vcat(z_array, zs[1:Int(zexo_index)])
    println(zexo_index)
    plot(Tss,zss/1e5,label="Temperature",color="black")
    yscale("log")
    xlim(0,3000)
    #xscale("log")
    xlabel("Temperature [K]",size=20)
    ylabel("Altitude [km]",size=20)
    xticks(size = 15)
    yticks(size = 15)
    #legend()
    println("n(exo): ",rho_exo_new)
    println("exobase: ", zexo_new)

    return zexo_new
end

# calc species density of the exobase with scale Height
# to estimate escape flux at the exobase
function calc_rho_exo(zexo,Texo,Tinf,rhoexo,sp_mass)
    z = alt[end]
    zexo_pre = alt[end]
    dz_temp = 1.e5
    rho = rhoexo
    for i in 1:((zexo-z)/dz_temp)
        z=z+dz_temp
        T = T_thermos(z,Texo,Tinf,zexo_pre) #TO do change!!!
        rho = rho*exp(-1/scaleH(z,T,sp_mass)*dz_temp)
    end
    return rho
end


#4.6Ga = 0 is defined as start time(t=0)
#Approximate flux equation is taken from Fig.4 of Tian+2009
#4.815e10*(t)^(-1.317) = flux(t) Tian
#start_time = 0.479677377 # Fast : 0.869864494,  Tian : 0.479677377 #単位はGyr
#start_time = 0.65 # 20 EUV at 3.95Ga assuming moderate rotator in Amerstorfer+2017
#start_time = 0.75 # 16.5 EUV at Ga
#start_time = 0.77 # 16 EUV at 3.83 Ga
#start_time = 0.789 # 15.5 EUV at 3.815 Ga
#start_time = 0.806 # 15 EUV at 3.8 Ga
#start_time = 0.824 # 14.5 EUV at 3.78 Ga
#start_time = 0.842 # 14 EUV at 3.76 Ga
#start_time = 0.860 #13.5 EUV at 3.740 Ga
#start_time = 0.879 #13 EUV at 3.721 Ga
#start_time = 0.85 # 13.7 EUV at 3.75 Ga

#Texo_tian = 1382. #fEUV=13, Kaorisan model
#Texo_tian = 1526. #fEUV=13.5, Kaorisan model
#Texo_tian = 1670. #fEUV=14, Kaorisan model
#Texo_tian = 1826. #fEUV=14.5, Kaorisan model
#Texo_tian = 2026. #fEUV=15, Kaorisan model
#Texo_tian = 2252. #fEUV=15.5, Kaorisan model
#Texo_tian = 2503. #fEUV=16, Kaorisan model
#Texo_tian = 2777. #fEUV=16.5, Kaorisan model

#start_time=0.
#Texo_tian = 2000. #fEUV=15.5, Kaorisan model

#Non-thermal escape flux of O & C
#Linealy extrapolated from Amerstorfer+2017
O_nonthermal = (1+0.16*feuv)*1.e8 #3<EUV<20
C_nonthermal = 0.
if 10 <= feuv <= 20
    C_nonthermal = (-2.7+0.43*feuv)*1.e8
elseif 3 <= feuv < 10
    C_nonthermal = (-0.12+0.17*feuv)*1.e8
else
    C_nonthermal = 0.
end

#C_nonthermal = (-2.7+0.43*feuv)*1.e8 #10<EUV<20
#C_nonthermal = (-0.12+0.17*feuv)*1.e8 #3<EUV<10
println("non-thermal O escape: ",O_nonthermal)
println("non-thermal C escape: ",C_nonthermal)


println("Texo: ",Texo_tian)


#flux_end_time_C = 0.65 # Fast : 1.116267056,  Tian : 0.65 #単位はGyr

#EUV time evolution
#reuv = 20. #EUV flux change
#from Amerstorfer et al. 2017
#moderate rotater
#20 EUV at elapse_time = 0.65 (3.95 Ga)
#10 EUV at elapse_time = 1.00 (3.6 Ga)
#approximate expression: y = a*log(x) + b
#a_euv*log(0.65) + b_euv = 20
#a_euv*log(0.80) + b_euv = 10
a_euv = 10. /log(0.65/1.)
b_euv = 10 - a_euv*log(1.)
reuv = a_euv*log(start_time) + b_euv
println("EUV: ", reuv)
#plot
#linspace(0.65,0.8,10)
function reuv_amers(a,b,t)
    return a*log(t) + b
end

function plot_reuv(a_euv,b_euv)
    figure()
    reuv_tlist=linspace(0.65,1.0,100)
    reuv_list=map(t->reuv_amers(a_euv,b_euv,t),reuv_tlist)
    reuv_tlist=4.6-reuv_tlist

    plot(reuv_tlist,reuv_list,label="moderate rotator",color="black")
    xlim(4.0,3.55)
    xlabel("Ga",size=20)
    ylabel("Leuv/Leuv,sun",size=20)
    xticks(size = 15)
    yticks(size = 15)
    legend()
end

function plot_cesc()
    figure()
    #CO2 outgassing
    tlist=linspace(0.5,1.,100)
    CO2out_lowest = map(t->1.4e10*exp(-2.624*t),tlist) #low
    CO2out_low = map(t->2.8e10*exp(-2.624*t),tlist) #low
    CO2out_middle = map(t->7.0e10*exp(-2.624*t),tlist) #middle
    CO2out_high = map(t->1.4e11*exp(-2.624*t),tlist) #high
    tlist=4.6-tlist
    plot(tlist, CO2out_lowest,color="black")
    #plot(tlist, CO2out_low,color="black")
    plot(tlist, CO2out_middle,color="black",label="CO2 outgassing")
    #plot(tlist,CO2out_high,color="black", label="CO2 outgassing")

    #Amerstorfer+17 moderate rotator
    #euv [20,16.5,16,15.5,15,14.5,14,13.5,13]
    #Texo [2800,2777,2503,2252,2026,1826,1670,1526,1382]
    #C escape C:CO = 0.1
    #Cesc = [4.31e10,3.69e10,5.62e9,1.35e9,4.49e8,1.69e8,7.24e7,2.91e7,9.72e6]
    #Oesc = [2.74e11, 2.24e11, 1.86e10,2.72e9,6.17e8,1.65e8,5.24e7,1.53e7,3.49e6]

    #C escape C:CO=0.3
    #Cesc = [1.95e11,1.69e11,2.42e10,4.93e9,1.47e9,5.24e8,2.19e8,8.70e7,2.90e7]
    #Oesc = [4.59e11,3.84e11,3.01e10,3.54e9,6.93e8,1.72e8,5.30e7,1.52e7,3.46e6]
    #C escape C:CO=0.5
    #Cesc = [4.18e11,3.66e11,5.23e10,9.54e9,2.62e9,8.94e8,3.67e8,1.45e8,4.81e7]
    #Oesc = [6.25e11,5.3e11,4.27e10,4.33e9,7.58e8,1.78e8,5.34e7,1.52e7,3.44e6]
    Cesc = [3.66e11,5.23e10,9.54e9,2.62e9,8.94e8,3.67e8,1.45e8,4.81e7] #remove 20 EUV
    Oesc = [5.3e11,4.27e10,4.33e9,7.58e8,1.78e8,5.34e7,1.52e7,3.44e6] #remove 20 EUV


    #C escape C:CO = 1.0
    #Cesc = [1.179e12,1.058e12,1.60e11,2.47e10,5.91e9,1.87e9,7.43e8,2.89e8,9.54e7]
    #Oesc = [9.33e11,9.27e11,7.50e10,6.13e9,8.89e8,1.88e8,5.42e7,1.51e7,3.40e6]


    #elapsed=[0.65,0.756,0.772,0.789,0.806,0.824,0.842,0.860,0.879]
    elapsed=[0.756,0.772,0.789,0.806,0.824,0.842,0.860,0.879] #remove 20 EUV

    plot(4.6-elapsed,2*Cesc,marker="o",fillstyle="none",linestyle="--",label="C escape x2",color="orange")
    plot(4.6-elapsed,Oesc,marker="^", fillstyle="none",linestyle="--",label="O escape",color="g")
    yscale("log")
    xlabel("Time [Ga]",size=20)
    ylabel(L"Flux $[cm^{-2}s^{-1}]$",size=20)
    xticks(size = 15)
    yticks(size = 15)
    legend()
end


#Texo_tian = 1382. #fEUV=13, Kaorisan model
#Texo_tian = 1526. #fEUV=13.5, Kaorisan model
#Texo_tian = 1670. #fEUV=14, Kaorisan model
#Texo_tian = 1826. #fEUV=14.5, Kaorisan model
#Texo_tian = 2026. #fEUV=15, Kaorisan model
#Texo_tian = 2252. #fEUV=15.5, Kaorisan model
#Texo_tian = 2503. #fEUV=16, Kaorisan model
#Texo_tian = 2777. #fEUV=16.5, Kaorisan model
#Texo_tian = 2800. #3000K from Tian 2009, figure 4
#Texo_tian = 3000. #3000K from Tian 2009, figure 4

C_CO_ratio = 1.0 #at exobase in the atmosphere ~1000km
C_CO2_ratio = 0.5 #at exobase in the model ~175km estimated from C:CO ratio in Amerstorfer+2017, figure1

#Zexo_tian = 1.0e8 #[cm] = 1000km
Zexo_tian = estimate_exobase(Texo_tian,Tinf,n_current[:O][end],16.)

#C_effusion_velocity = effusion_velocity(Texo_tian,12.0)
C_effusion_velocity = flex_effusion_velocity(Texo_tian,12.0, Zexo_tian) #1000km from Tian+2009
H_effusion_velocity = flex_effusion_velocity(Texo_tian,1.0,Zexo_tian)
H2_effusion_velocity = flex_effusion_velocity(Texo_tian,2.0,Zexo_tian)
O_effusion_velocity = flex_effusion_velocity(Texo_tian,16.0,Zexo_tian)
#print(C_effusion_velocity)

##### modified jeans escape ####################################################
#####
# Tian+2009 calculated escape flux by modified jeans escape
# the effect of bulk velocity is taken into account in modified jeans escape
# Step
# we followed the way of Koskinen+2013
# using equation (9) from Volkov+2011b
# bulk velocity is determined by (6) from Koskinen+2013, Tian+2009
# bulk velocity = sum(m_i * n_i * effusion_v_i)/sum(m_i*n_i), i=C,H,H2,O
#####
#bulk velocity
function bulk_velocity(n_current, Texo, Zexobase)
    #number density -> density[kg/cm3]
    #rho_H = 1*mH*n_current[:H][end]
    #rho_H2 = 1*mH*2*n_current[:H2][end]
    #rho_C = 12*mH*n_current[:CO][end]
    #rho_O = 16*mH*n_current[:O][end]
    rho_H = 1*mH*calc_rho_exo(Zexobase,Texo,Tinf,n_current[:H][end],1.)
    rho_H2 = 1*mH*2*calc_rho_exo(Zexobase,Texo,Tinf,n_current[:H2][end],2.)
    rho_C = 12*mH*calc_rho_exo(Zexobase,Texo,Tinf,n_current[:CO][end]*C_CO2_ratio,12.)
    rho_O = 16*mH*calc_rho_exo(Zexobase,Texo,Tinf,n_current[:O][end],16.)
    sum_rho = rho_H + rho_H2 + rho_C + rho_O #[kg], H,H2,C,O

    #effusion velocity of each species
    H_eff=flex_effusion_velocity(Texo, 1.0, Zexobase) #H
    H2_eff=flex_effusion_velocity(Texo, 2.0, Zexobase) #H2
    C_eff=flex_effusion_velocity(Texo, 12.0, Zexobase) #C
    O_eff=flex_effusion_velocity(Texo, 16.0, Zexobase) #O

    #average
    bulk = (rho_H*H_eff + rho_H2*H2_eff + rho_C*C_eff + rho_O*O_eff)/sum_rho
    return bulk
end

using SpecialFunctions
#(9) from volkov+2009
function modified_factor(n_current,Texo, Zexobase, m)
    S = bulk_velocity(n_current,Texo,Zexobase)/sqrt(2*boltzmannK*Texo/(m*mH))
    lambda = (m*mH*bigG*marsM)/(boltzmannK*Texo*1e-2*(radiusM+Zexobase)) #jeans parameter
    #println(lambda)
    Se = sqrt(lambda)
    nume = 0.5*exp(-S^2-Se^2) + (S*Se+S^2-0.5)*exp(-(Se-S)^2) + sqrt(pi)*S^3*(1-erf(Se-S))
    deno = S^2*(1+Se^2)*exp(-Se^2)
    factor = nume/deno
    return factor
end

function plot_modified_effusion(n_current, Texo_list, Zexobase, m)
    effusion_modified_list = map(x -> flex_effusion_velocity(x,m, Zexobase)*modified_factor(n_current,x,Zexobase,m), Texo_list)
    effusion_list = map(x -> flex_effusion_velocity(x,m, Zexobase), Texo_list)
    figure()
    plot(Texo_list, effusion_modified_list, label="modified jeans", color="black")
    plot(Texo_list, effusion_list, label="jeans", color="black",linestyle="--")
    #xlim(1e5,1e10)
    yscale("log")
    xlabel("Exobase temperature[K]",size=20)
    ylabel("Effusion velocity",size=20)
    xticks(size = 15)
    yticks(size = 15)
    legend()
end

function plot_modified_effusion_compare(n_current, Texo_list, Zexobase, m, name)
    effusion_modified_list = map(x -> flex_effusion_velocity(x,m, Zexobase)*modified_factor(n_current,x,Zexobase,m), Texo_list)
    effusion_list = map(x -> flex_effusion_velocity(x,m, Zexobase), Texo_list)
    #figure()
    #plot(Texo_list, effusion_modified_list, label=name)
    #plot(Texo_list, effusion_list, label=name,linestyle="--")
    plot(Texo_list, effusion_list, label=name)
    #xlim(1e5,1e10)
    yscale("log")
    xlabel("Exobase temperature[K]",size=20)
    ylabel("Effusion velocity",size=20)
    xticks(size = 15)
    yticks(size = 15)
    legend()
end


################################################################################

#CO escape ver
function O_upper_boundary_CO(n_current,elapsed_time_boundary,dt)

    #flux_upper_boundary_O = - C_effusion_velocity*modified_factor(n_current,Texo_tian,1.0e8,12.0) * n_current[:CO][end] + 4.2e8 #CO escape ver, non-thermal is taken from Amerstorfer+2017
    #flux_upper_boundary_O = - C_effusion_velocity*modified_factor(n_current,Texo_tian,1.0e8,12.0) * n_current[:CO][end] + n_current[:O][end]*O_effusion_velocity*modified_factor(n_current,Texo_tian,Zexo_tian,16.0) + 4.2e8
    flux_upper_boundary_O = - calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:CO][end]*C_CO2_ratio,12.)*C_CO_ratio * flex_effusion_velocity(Texo_tian,12.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,12.0) +calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:O][end],16.)*flex_effusion_velocity(Texo_tian,16.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,16.0) +O_nonthermal-C_nonthermal
    return flux_upper_boundary_O
end

function C_upper_boundary(n_current,elapsed_time_boundary,dt)
    return calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:CO][end]*C_CO2_ratio,12.)*C_CO_ratio*flex_effusion_velocity(Texo_tian,12.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,12.0)+C_nonthermal
end

####### CO2 outgassing #####################################################
#outgassing rate is also taken from Tian+2009
function CO2_outgas_linear_flux_lower_boundary(elapsed_time_boundary,dt)
    #f_outgas= 10^(log10(1.0e10)-log10(9)*(elapsed_time_boundary/(2*3.14e7*1.0e9))) # ← 時間の関数 初期値を設定する (20億年で一桁落ち)

    #Tian 2009, CO2 outgassing rate, 手入力
    #### middle outgassing ######################
    t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    f_outgas = 7.0e10*exp(-2.624*t1_Gyr) #koyama calc again from Tian+2009

    ##### high_outgas##############
    #t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    #f_outgas= 1.4e11*exp(-2.624*t1_Gyr)
    #########################

    ###higher outgas##
    #t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    #f_outgas= 2.8e11*exp(-2.624*t1_Gyr)

    ### low CO2 outgassing rate
    #t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    #f_outgas= 2.8e10*exp(-2.624*t1_Gyr)

    # low and lowest
    #t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    #f_outgas= 2.e10*exp(-2.624*t1_Gyr)

    ### lowest CO2 outgassing
    #t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    #f_outgas= 1.4e10*exp(-2.624*t1_Gyr)

    ### try lower than lowest
    #t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    #f_outgas= 1.4e9*exp(-2.624*t1_Gyr)

    ### between low and middle ####
    #t1_Gyr=(elapsed_time_boundary/3.15e16) + start_time #単位はGyr
    #f_outgas= 5.0e10*exp(-2.624*t1_Gyr)

    #println("f_outgas=",f_outgas)
    return f_outgas
end

function CO2_outgas_flux_lower_boundary(elapsed_time_boundary,dt)
    outgas_flux_CO2 = CO2_outgas_linear_flux_lower_boundary(elapsed_time_boundary,dt)
    #outgas_flux_CO = CO_outgas_rectangle_flux_lower_boundary(elapsed_time_boundary,dt)
    #outgas_flux_CO = 1.0e10 #脱ガス一定の場合
    return outgas_flux_CO2
end


#koyama 2020/06/17
#境界条件がタイムステップ毎に変わるようにするため関数ではなく、マクロで保存するように変更した。
speciesbclist=Dict(
              #CO escape ver
              :CO2=>["f" :(CO2_outgas_flux_lower_boundary(totaltime,dt)); "f" 0.], #CO esc ,normal:2.1e17 koyama
              :CO=>["f" 0.; "v" :(C_upper_boundary(n_current,totaltime,dt)/n_current[:CO][end])], #八木くんの境界条件が "f" になっていたが、boundaryconditionsで解決していた。w
              :O=>["v" 0.; "v" :(O_upper_boundary_CO(n_current,totaltime,dt)/n_current[:O][end])],#normal:1.2e8
              #:O=>["v" 0.; "v" :(calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:O][end],16.)/n_current[:O][end]*flex_effusion_velocity(Texo_tian,16.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,16.0))], #H deposition exists?
              #:CO=>["f" 0.; "v" 0.], #八木くんの境界条件が "f" になっていたが、boundaryconditionsで解決していた。w


              #others BD
              :H2=>["f" 0.; "v" :(calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:H2][end],2.)/n_current[:H2][end]*flex_effusion_velocity(Texo_tian,2.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,2.0))],
              :H=>["v" 0.; "v" :(calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:H][end],1.)/n_current[:H][end]*flex_effusion_velocity(Texo_tian,1.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,1.0))], #H deposition exists?
              :H2O2=>["v" depv; "f" 0.],
              :O3=>["v" depv; "f" 0.],
              :OH=>["v" 0.; "f" 0.],
              #:HOCO=>["v" depv; "f" 0.],
              :HO2=>["v" 0.; "f" 0.],
              :O1D=>["v" 0.; "f" 0.],

              #non-related
              #:Ar=>["n" 2.0e-2*2.1e17; "f" 0.],
              #:N2=>["n" 1.9e-2*2.1e17; "f" 0.],
              #:H2O=>["n" H2Osat[1]; "f" 0.], # boundary condition does
                                             # not matter if H2O is
                                             # fixed
              );

function speciesbcs(species)
    get(speciesbclist,
        species,
        ["f" 0.; "f" 0.])
end


################################################################################
############################ DIFFUSION COEFFICIENTS ############################
################################################################################

#Koyama change 60km to z_lower_tropo, to adapt to different atmospheric pressure
function Keddy(n_current, z)
    # eddy diffusion coefficient, stolen from Krasnopolsky (1993).
    # Scales as the inverse sqrt of atmospheric number density
    z <= z_lower_tropo ? 10^6 : 2e13/sqrt(n_tot(n_current, z))
end

function Keddy(z::Real, nt::Real)
    # eddy diffusion coefficient, stolen from Krasnopolsky (1993).
    # Scales as the inverse sqrt of atmospheric number density
    z <= z_lower_tropo ? 10^6 : 2e13/sqrt(nt)
end

# molecular diffusion parameters. molecular diffusion is different only
# for small molecules and atoms (H and H2), otherwise all species share
# the same values (Krasnopolsky 1993 <- Banks and Kockarts Aeronomy)
# THESE ARE IN cm^-2 s^-2!!!
diffparams(species) = get(Dict(:H=>[8.4, 0.597],:H2=>[2.23, 0.75]),species,[1.0, 0.75])
function Dcoef(T, n::Real, species::Symbol)
    dparms = diffparams(species)
    dparms[1]*1e17*T^(dparms[2])/n
end
Dcoef(z, species::Symbol, n_current) = Dcoef(Temp(z),n_tot(n_current, z),species)

#thermal diffusion factors (from Krasnopolsky 2002)
thermaldiff(species) = get(Dict(:H=>-0.25,:H2=>-0.25,:HD=>-0.25,:D2=>-0.25,
                                :He=>-0.25), species, 0)

################################################################################
################################## TRANSPORT ###################################
################################################################################

# scale height at a given altitude
function scaleH(z, T::Float64, mm::Real)
    boltzmannK*T/(mm*mH*marsM*bigG)*(((z+radiusM)*1e-2)^2)*1e2
    # constants are in MKS. Convert to m and back to cm.
end

function scaleH(z, species::Symbol)
    T=Temp(z)
    mm = speciesmolmasslist[species]
    scaleH(z, T, mm)
end

function scaleH(z, T::Float64, species::Symbol)
    mm = speciesmolmasslist[species]
    scaleH(z, T, mm)
end

function scaleH(z, T::Float64, n_current)
    mm = meanmass(n_current, z)
    scaleH(z, T, mm)
end

function scaleH(z, n_current)
    T=Temp(z)
    scaleH(z, T)
end

# at each level of the atmosphere, density can be transferred up or
# down with a specified rate coefficient.
#
#                          | n_i+1
#       ^  tspecies_i_up   v tspecies_i+1_down
#   n_i |
#       v tspecies_i_down  ^ tspecies_i-1_up
#                          | n_i-1
#
#
# the flux at each cell boundary is the sum of the upward flux from
# the cell below and the downward flux of the cell above. These fluxes
# are determined using flux coefficients that come from the diffusion
# equation. Care must be taken at the upper and lower boundary so that
# tspecies_top_up and tspecies_bottom_down properly reflect the
# boundary conditions of the atmosphere.


# This is handled in the code with the population of the appropriate
# reactions, with variable rate coeffecients that are populated
# between each timestep (similar to the way photolysis rates are
# included). We need to make reactions at each interior altitude
# level:
#          n_i -> n_i+1  tspecies_i_up
#          n_i -> n_i-1  tspecies_i_down
#
# At the upper and lower boundary we omit the species on the RHS, so
# that these reactions are potentially non-conservative:
#
#         n_top    -> NULL  tspecies_top_up
#         n_bottom -> NULL  tspecies_bottom_down
#
# These coefficients then describe the diffusion velocity at the top
# and bottom of ther atmosphere.


# function to generate coefficients of the transport network
function fluxcoefs(z, dz, ntv, Kv, Dv, Tv, Hsv, H0v, species)
    Dl = (Dv[1] + Dv[2])/2.0
    Kl = (Kv[1] + Kv[2])/2.0
    Tl = (Tv[1] + Tv[2])/2.0
    dTdzl = (Tv[2] - Tv[1])/dz
    Hsl = (Hsv[1] + Hsv[2])/2.0
    H0l = (H0v[1] + H0v[2])/2.0

    # we have two flux terms to combine:
    sumeddyl = (Dl+Kl)/dz/dz
    gravthermall = (Dl*(1/Hsl + (1+thermaldiff(species))/Tl*dTdzl) +
                    Kl*(1/H0l + 1/Tl*dTdzl))/(2*dz)

    Du = (Dv[2] + Dv[3])/2.0
    Ku = (Kv[2] + Kv[3])/2.0
    Tu = (Tv[2] + Tv[3])/2.0
    dTdzu = (Tv[3] - Tv[2])/dz
    Hsu = (Hsv[2] + Hsv[3])/2.0
    H0u = (H0v[2] + H0v[3])/2.0

    # we have two flux terms to combine:
    sumeddyu = (Du+Ku)/dz/dz
    gravthermalu = (Du*(1/Hsu + (1 + thermaldiff(species))/Tu*dTdzu) +
                  Ku*(1/H0u + 1/Tu*dTdzu))/(2*dz)

    # this results in the following coupling coefficients:
    return [sumeddyl+gravthermall, # down
            sumeddyu-gravthermalu] # up
end

# overload to generate the coefficients if they are not supplied
function fluxcoefs(z, dz, species, n_current)
    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)
    Tp = Temp(z+dz)
    T0 = Temp(z)
    Tm = Temp(z-dz)
    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    # return the coefficients
    return fluxcoefs(z, dz,
              [ntm, nt0, ntp],
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)
end

# define transport coefficients for a given atmospheric layers for transport
# from a layer to the one above
function lower_up(z, dz, species, n_current)
    ntp = n_tot(n_current, z+dz)
    nt0 = n_tot(n_current, z)
    ntm = 1
    Kp = Keddy(z+dz, ntp)
    K0 = Keddy(z,nt0)
    Km = 1
    Tp = Temp(z+dz)
    T0 = Temp(z)
    Tm = 1
    Dp = Dcoef(Tp, ntp, species)
    D0 = Dcoef(T0, nt0, species)
    Dm = 1
    Hsp = scaleH(z+dz, species)
    Hs0 = scaleH(z,species)
    Hsm = 1
    H0p = scaleH(z+dz, Tp, n_current)
    H00 = scaleH(z,T0, n_current)
    H0m = 1

    #return the coefficients
    return fluxcoefs(z, dz,
              [ntm, nt0, ntp],
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[2]
end

# define transport coefficients for a given atmospheric layers for transport
# from a layer to the one below
function upper_down(z, dz, species, n_current)
    ntp = 1
    nt0 = n_tot(n_current, z)
    ntm = n_tot(n_current, z-dz)
    Kp = 1
    K0 = Keddy(z, nt0)
    Km = Keddy(z-dz, ntm)
    Tp = 1
    T0 = Temp(z)
    Tm = Temp(z-dz)
    Dp = 1
    D0 = Dcoef(T0, nt0, species)
    Dm = Dcoef(Tm, ntm, species)
    Hsp = 1
    Hs0 = scaleH(z, species)
    Hsm = scaleH(z-dz, species)
    H0p = 1
    H00 = scaleH(z, T0, n_current)
    H0m = scaleH(z-dz, Tm, n_current)

    #return the coefficients
    return fluxcoefs(z, dz,
              [ntm, nt0, ntp],
              [Km , K0, Kp],
              [Dm , D0, Dp],
              [Tm , T0, Tp],
              [Hsm, Hs0, Hsp],
              [H0m, H00, H0p],
              species)[1]
end

function boundaryconditions(species, dz, n_current)
    # returns the symbolic transport coefficients that encode the
    # boundary conditions for the null-pointing equations

    # n_1->NULL t_lower_bc
    # n_f->NULL t_upper_bc


    # this defines two additional symbols for each species that need
    # to be resolved in the function call macro:
    #            tspecies_lower_up and
    #            tspecies_upper_down
    # these are found by passing the appropriate values to fluxcoefs
    # and selecting the correct output

    bcs = speciesbcs(species)
    if issubset([species],notransportspecies)
        bcs = ["f" 0.; "f" 0.]
    end

    ####koyama added 2020/06/17
    #これをbcs[1,2]=eval(bcs[1,2])にしてしまうと、オリジナルのspeciesbclistも変更されてしまう！！
    #time-depenedent boundary condition!!
    bcs_lower = eval(bcs[1,2])
    bcs_upper = eval(bcs[2,2])
    ######################################

    #first element returned corresponds to lower BC, second to upper
    #BC transport rate. Within each element, the two rates correspond
    #to the two equations
    # n_b  -> NULL (first rate, depends on species concentration)
    # NULL -> n_b  (second rate, independent of species concentration
    bcvec = Float64[0 0;0 0]

    # LOWER
    if bcs[1, 1] == "n"
        bcvec[1,:]=[fluxcoefs(alt[2], dz, species, n_current)[1],
                    lower_up(alt[1], dz, species, n_current)*bcs_lower]
    elseif bcs[1, 1] == "f"
        bcvec[1,:] = [0.0, bcs_lower/dz]
    elseif bcs[1, 1] == "v"
        bcvec[1,:] = [bcs_lower/dz, 0.0]
    else
        throw("Improper lower boundary condition!")
    end

    # UPPER
    if bcs[2, 1] == "n"
        bcvec[2,:] = [fluxcoefs(alt[end-1],dz, species, n_current)[2],
                    upper_down(alt[end],dz, species, n_current)*bcs_upper] #I fixed it from bcs[1,2] to bcs[2,2]
    elseif bcs[2, 1] == "f"
            bcvec[2,:] = [0.0,-bcs_upper/dz]
    elseif bcs[2, 1] == "v"
        bcvec[2,:] = [bcs_upper/dz, 0.0]
    else
        throw("Improper upper boundary condition!")
    end

    #return the bc vec
    return bcvec
end


################################################################################
################################## CHEMISTRY ###################################
################################################################################

# TODO: again, the Any[] syntax might be removable once [[],[]] concatenation
# is diabled rather than depracated
function getpos(array, test::Function, n=Any[])
    # this function searches through an arbitrarily structured array
    # finding elements that match the test function supplied, and returns a
    # one-dimensional array of the indicies of these elements.
    if !isa(array, Array)
        test(array) ? Any[n] : Any[]
    else
        vcat([ getpos(array[i], test, Any[n...,i]) for i=1:size(array)[1] ]...)
    end
end

function getpos(array, value)
    # overloading getpos for the most common use case, finding indicies
    # corresponding to elements in array that match value.
    getpos(array, x->x==value)
end

function deletefirst(A, v)
    # Returns list A with its first element equal to v removed.
    index = findfirst_or_zero(v, A)
    keep = setdiff([1:length(A);],index)
    A[keep]
end

function loss_equations(network, species)
    # given a network of equations in the form of reactionnet above, this
    # function returns the LHS (reactants) and rate coefficient for all
    # reactions where the supplied species is consumed.

    # get list of all chemical reactions species participates in:
    speciespos = getpos(network, species)
    # find pos where species is on LHS but not RHS:
    lhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1],  # we only need the reaction number
               map(x->speciespos[x],  # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(lhspos, rhspos)
        lhspos = deletefirst(lhspos, i)
    end

    # get the products and rate coefficient for the identified
    # reactions.
    losseqns=map(x->vcat(Any[network[x][1]...,network[x][3]]), lhspos)
    # automatically finds a species where it occurs twice on the LHS
end

function loss_rate(network, species)
    # return a symbolic expression for the loss rate of species in the
    # supplied reaction network.
    leqn=loss_equations(network, species) # get the equations
    lval=:(+($( # and add the products together
               map(x->:(*($(x...))) # take the product of the
                                    # concentrations and coefficients
                                    # for each reaction
                   ,leqn)...)))
end

function production_equations(network, species)
    # given a network of equations in the form of reactionnet above, this
    # function returns the LHS (reactants) and rate coefficient for all
    # reactions where the supplied species is a product.
    speciespos = getpos(network, species)#list of all reactions where species is produced
    # find pos where species is on RHS but not LHS
    lhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], #s elect the appropriate reactions
                   findall(x->x[2]==1, speciespos)))
    rhspos = map(x->x[1], # we only need the reaction number
               map(x->speciespos[x], # select the appropriate reactions
                   findall(x->x[2]==2, speciespos)))
    for i in intersect(rhspos, lhspos)
        rhspos=deletefirst(rhspos, i)
    end

    prodeqns = map(x->vcat(Any[network[x][1]...,network[x][3]]), # get the products and rate
                                                            # coefficient for the identified
                                                            # reactions.
                 # automatically finds and counts duplicate
                 # production for each molecule produced
                 rhspos)

    return prodeqns
end

function production_rate(network, species)
    # return a symbolic expression for the loss rate of species in the
    # supplied reaction network.

    # get the reactants and rate coefficients
    peqn = production_equations(network, species)

    # and add them all up
    # take the product of each set of reactants and coeffecient
    pval = :(+($(map(x -> :(*($(x...))), peqn) ...) ))
end

function chemical_jacobian(chemnetwork, transportnetwork, specieslist, dspecieslist)
    # Compute the symbolic chemical jacobian of a supplied reaction network
    # for the specified list of species. Returns three arrays suitable for
    # constructing a sparse matrix: lists of the first and second indices
    # and the symbolic value to place at that index.

    # set up output vectors: indices and values
    ivec = Int64[] # list of first indices (corresponding to the species being produced and lost)
    jvec = Int64[] # list of second indices (corresponding to the derivative being taken)
    tvec = Any[] # list of the symbolic values corresponding to the jacobian

    nspecies = length(specieslist)
    ndspecies = length(dspecieslist)
    for i in 1:nspecies #for each species
        ispecies = specieslist[i]
        #get the production and loss equations
        peqn = []
        leqn = []
        if issubset([ispecies],chemspecies)
            peqn = [peqn; production_equations(chemnetwork, ispecies)]
            leqn = [leqn; loss_equations(chemnetwork, ispecies)]
        end
        if issubset([ispecies],transportspecies)
            peqn = [peqn; production_equations(transportnetwork, ispecies)]
            leqn = [leqn; loss_equations(transportnetwork, ispecies)]
        end
        for j in 1:ndspecies #now take the derivative with resp`ect to the other species
            jspecies = dspecieslist[j]
            #find the places where the production rates depend on
            #jspecies, and return the list rates with the first
            #occurrance of jspecies deleted. (Note: this seamlessly
            #deals with multiple copies of a species on either side of
            #an equation, because it is found twice wherever it lives)
            ppos = map(x->deletefirst(peqn[x[1]],jspecies),getpos(peqn, jspecies))
            lpos = map(x->deletefirst(leqn[x[1]],jspecies),getpos(leqn, jspecies))
            if length(ppos)+length(lpos)>0 #if there is a dependence
                #make note of where this dependency exists
                append!(ivec,[i])
                append!(jvec,[j])
                #smash the production and loss rates together,
                #multiplying for each distinct equation, adding
                #together the production and loss seperately, and
                #subtracting loss from production.
                if length(ppos)==0
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($lval))
                elseif length(lpos)==0
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    tval = :(+($pval))
                else
                    pval = :(+($(map(x->:(*($(x...))),ppos)...)))
                    lval = :(+($(map(x->:(*($(x...))),lpos)...)))
                    tval = :(-($pval,$lval))
                end
                #attach the symbolic expression to the return values
                append!(tvec,[tval])
            end
        end
    end
    return (ivec, jvec, Expr(:vcat, tvec...))
end

################################################################################
####################### COMBINED CHEMISTRY AND TRANSPORT #######################
################################################################################

# We now have objects that return the list of indices and coefficients
# for transport, assuming no other species in the atmosphere
# (transportmat), and for chemistry, assuming no other altitudes
# (chemical_jacobian). We need to perform a kind of outer product on
# these operators, to determine a fully coupled set of equations for
# all species at all altitudes.

# need to get a list of all species at all altitudes to iterate over
const intaltgrid=round.(Int64, alt/1e5)[2:end-1]
const replacespecies=[fullspecieslist, Jratelist,[:T,:M]]

# the rates at each altitude can be computed using the reaction network
# already in place, plus additional equations describing the transport
# to and from the cells above and below:

upeqns = [Any[Any[[s], [Symbol(string(s)*"_above")],Symbol("t"*string(s)*"_up")],
         Any[[Symbol(string(s)*"_above")],[s],Symbol("t"*string(s)*"_above_down")]]
        for s in specieslist]

downeqns = [Any[Any[[s], [Symbol(string(s)*"_below")],Symbol("t"*string(s)*"_down")],
           Any[[Symbol(string(s)*"_below")],[s],Symbol("t"*string(s)*"_below_up")]]
          for s in specieslist]

# NOTE: next line works but is really ugly. It is no longer possible to transpose
# symbol arrays in Julia 0.6, and permutedims doesn't work either.
local_transport_rates = [[[Symbol("t"*string(s)*"_up") for s in specieslist]
                          [Symbol("t"*string(s)*"_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_above_down") for s in specieslist]
                          [Symbol("t"*string(s)*"_below_up") for s in specieslist]]...;]


transportnet = [[upeqns...;];[downeqns...;]]

# define names for all the species active in the coupled rates:
activespecies = union(chemspecies, transportspecies)
active_above = [Symbol(string(s)*"_above") for s in activespecies]
active_below = [Symbol(string(s)*"_below") for s in activespecies]

inactivespecies = intersect(nochemspecies, notransportspecies)

function getrate(chemnet, transportnet, species)
    # Calculates the rate at which a given species is either produced or lost.
    # Production is from chemical reaction yields or entry from other
    # atmospheric layers. Loss is due to consumption in reactions or migration
    # to other layers. Comment added by EMC 8/29/17
    rate = :(0.0)
    if issubset([species],chemspecies)
        rate = :($rate
               + $(production_rate(chemnet, species))
               - $(      loss_rate(chemnet, species)))
    end
    if issubset([species],transportspecies)
        rate = :($rate
               + $(production_rate(transportnet, species))
               - $(      loss_rate(transportnet, species)))
    end

    return rate
end

# obtain the rates and jacobian for each altitude
const rates_local = Expr(:vcat, map(x->getrate(reactionnet, transportnet, x),activespecies)...);
const chemJ_local = chemical_jacobian(reactionnet, transportnet, activespecies, activespecies);
const chemJ_above = chemical_jacobian(reactionnet, transportnet, activespecies, active_above);
const chemJ_below = chemical_jacobian(reactionnet, transportnet, activespecies, active_below);

arglist_local = [activespecies;
                 active_above;
                 active_below;
                 inactivespecies;
                 Jratelist;
                 :T; :M;
                 local_transport_rates;
                 :dt]

arglist_local_typed=[:($s::Float64) for s in arglist_local]

@eval begin
    function ratefn_local($(arglist_local_typed[1:end-1]...))
        $rates_local
    end
end

@eval begin
    function chemJmat_local($(arglist_local_typed...))
        localchemJi = $(chemJ_local[1])
        localchemJj = $(chemJ_local[2])
        localchemJval = -dt*$(chemJ_local[3])

        abovechemJi = $(chemJ_above[1])
        abovechemJj = $(chemJ_above[2])
        abovechemJval = -dt*$(chemJ_above[3])

        belowchemJi = $(chemJ_below[1])
        belowchemJj = $(chemJ_below[2])
        belowchemJval = -dt*$(chemJ_below[3])

        ((localchemJi, localchemJj, localchemJval),
         (abovechemJi, abovechemJj, abovechemJval),
         (belowchemJi, belowchemJj, belowchemJval))
    end
end

# a function to return chemical reaction rates for specified species
# concentrations
@eval begin
    function reactionrates_local($(specieslist...), $(Jratelist...), T,  M)
        $(Expr(:vcat, map(x->Expr(:call,:*,x[1]..., x[3]), reactionnet)...))
    end
end


function reactionrates(n_current)
    theserates = fill(convert(Float64, NaN),(length(intaltgrid),length(reactionnet)))
    for ialt in 1:length(intaltgrid)
        theserates[ialt,:] = reactionrates_local([[n_current[sp][ialt] for sp in specieslist];
                                                [n_current[J][ialt] for J in Jratelist];
                                                Temp(alt[ialt+1]);
                                                n_tot(n_current, alt[ialt+1])]...)
    end
    return theserates
end



function getflux(n_current, dz, species)
    thesecoefs = [fluxcoefs(a, dz, species, n_current) for a in alt[2:end-1]]
    thesebcs = boundaryconditions(species, dz, n_current)

    thesefluxes = fill(convert(Float64, NaN),length(intaltgrid))

    thesefluxes[1] = (-(n_current[species][2]*thesecoefs[2][1]
                      -n_current[species][1]*thesecoefs[1][2])
                    +(-n_current[species][1]*thesebcs[1, 1]
                      +thesebcs[1, 2]))/2.0
    for ialt in 2:length(intaltgrid)-1
        thesefluxes[ialt] = (-(n_current[species][ialt+1]*thesecoefs[ialt+1][1]
                             - n_current[species][ialt]*thesecoefs[ialt][2])
                             + (-n_current[species][ialt]*thesecoefs[ialt][1]
                             + n_current[species][ialt-1]*thesecoefs[ialt-1][2]))/2.0
    end
    thesefluxes[end] = (-(thesebcs[2, 2]
                        - n_current[species][end]*thesebcs[2, 1])
                        + (-n_current[species][end]*thesecoefs[end][1]
                        + n_current[species][end-1]*thesecoefs[end-1][2]))/2.0
    return dz*thesefluxes
end

function fluxes(n_current, dz)
    thesefluxes = fill(convert(Float64, NaN),(length(intaltgrid),length(specieslist)))
    for isp in 1:length(specieslist)
        thesefluxes[:,isp] = getflux(n_current, dz, specieslist[isp])
    end
    thesefluxes
end


function ratefn(nthis, inactive, Jrates, T, M, tup, tdown, tlower, tupper)
    # at each altitude, get the appropriate group of concentrations,
    # coefficients, and rates to pass to ratefn_local
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive,(length(inactivespecies),length(intaltgrid)))

    returnrates = zero(nthismat)

    # fill the first altitude entry with information for all species
    returnrates[:,1] = ratefn_local([nthismat[:,1]; nthismat[:,2];
                                    fill(1.0, length(activespecies));
                                    inactivemat[:,1]; Jrates[:,1]; T[1];M[1];
                                    tup[:,1]; tlower[:,1]; tdown[:,2];
                                    tlower[:,2]]...)

    # iterate through other altitudes except the last level, filling the info in
    for ialt in 2:(length(intaltgrid)-1)
        returnrates[:,ialt] = ratefn_local([nthismat[:,ialt];
                                          nthismat[:,ialt+1];
                                          nthismat[:,ialt-1];
                                          inactivemat[:,ialt];
                                          Jrates[:,ialt];
                                          T[ialt]; M[ialt];
                                          tup[:,ialt];tdown[:,ialt];
                                          tdown[:,ialt+1];tup[:,ialt-1]]...)
    end

    # fill in the last level of altitude (200 km)
    returnrates[:,end] = ratefn_local([nthismat[:,end];
                                       fill(1.0, length(activespecies));
                                       nthismat[:,end-1];
                                       inactivemat[:,end];
                                       Jrates[:,end];
                                       T[end]; M[end];
                                       tupper[:,1]; tdown[:,end];
                                       tupper[:,2]; tup[:,end-1]]...)
    return [returnrates...;]
end


function chemJmat(nthis, inactive, Jrates, T, M, tup, tdown, tlower, tupper, dt)
    nthismat = reshape(nthis, (length(activespecies), length(intaltgrid)))
    inactivemat = reshape(inactive, (length(inactivespecies), length(intaltgrid)))
    chemJi = Int64[]
    chemJj = Int64[]
    chemJval = Float64[]
    (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,1];
                                              nthismat[:,2];
                                              fill(1.0, length(activespecies));
                                              inactivemat[:,1];
                                              Jrates[:,1];
                                              T[1]; M[1];
                                              tup[:,1]; tlower[:,1];
                                              tdown[:,2]; tlower[:,2];dt]...)
    #add the influence of the local densities
    append!(chemJi, tclocal[1])
    append!(chemJj, tclocal[2])
    append!(chemJval, tclocal[3])
    #and the upper densities
    append!(chemJi, tcupper[1])
    append!(chemJj, tcupper[2].+length(activespecies))
    append!(chemJval, tcupper[3])

    for ialt in 2:(length(intaltgrid)-1)
        (tclocal, tcupper, tclower) = chemJmat_local([nthismat[:,ialt];
                                                      nthismat[:,ialt+1];
                                                      nthismat[:,ialt-1];
                                                      inactivemat[:,ialt];
                                                      Jrates[:,ialt];
                                                      T[ialt]; M[ialt];
                                                      tup[:,ialt];tdown[:,ialt];
                                                      tdown[:,ialt+1];
                                                      tup[:,ialt-1];dt]...)
        #add the influence of the local densities
        append!(chemJi, tclocal[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclocal[2].+(ialt-1)*length(activespecies))
        append!(chemJval, tclocal[3])
        #and the upper densities
        append!(chemJi, tcupper[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tcupper[2].+(ialt  )*length(activespecies))
        append!(chemJval, tcupper[3])
        #and the lower densities
        append!(chemJi, tclower[1].+(ialt-1)*length(activespecies))
        append!(chemJj, tclower[2].+(ialt-2)*length(activespecies))
        append!(chemJval, tclower[3])
    end

    (tclocal, tcupper, tclower)=chemJmat_local([nthismat[:,end];
                                              fill(1.0, length(activespecies));
                                              nthismat[:,end-1];
                                              inactivemat[:,end];
                                              Jrates[:,end];
                                              T[end]; M[end];
                                              tupper[:,1]; tdown[:,end];
                                              tupper[:,2]; tup[:,end-1];dt]...)
    #add the influence of the local densities
    append!(chemJi, tclocal[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclocal[2].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJval, tclocal[3])
    #and the lower densities
    append!(chemJi, tclower[1].+(length(intaltgrid)-1)*length(activespecies))
    append!(chemJj, tclower[2].+(length(intaltgrid)-2)*length(activespecies))
    append!(chemJval, tclower[3])

    #make sure to add 1's along the diagonal
    append!(chemJi,[1:length(nthis);])
    append!(chemJj,[1:length(nthis);])
    append!(chemJval, fill(1.0, length(nthis)))

    sparse(chemJi,
           chemJj,
           chemJval,
           length(nthis),
           length(nthis),
           +);

end

function next_timestep(nstart::Array{Float64, 1},
                       nthis::Array{Float64, 1},
                       inactive::Array{Float64, 1},
                       Jrates::Array{Float64, 2},
                       T::Array{Float64, 1},
                       M::Array{Float64, 1},
                       tup::Array{Float64, 2},
                       tdown::Array{Float64, 2},
                       tlower::Array{Float64, 2},
                       tupper::Array{Float64, 2},
                       dt::Float64)
   # moves to the next timestep using Newton's method on the linearized
   # coupled transport and chemical reaction network.

    eps = 1.0#ensure at least one iteration
    iter = 0
    while eps>1e-8 #koyama change
        nold = deepcopy(nthis)

        # I modified some of the code solving differential equation
        # (taking the inverse of Jacobian) to avoid an error.
        # following the way of (Appendix in Catling & Kasting, 2017).
        # stuff concentrations into update function and jacobian
        fval = dt*ratefn(nthis, inactive, Jrates, T, M, tup,tdown, tlower, tupper) - 0.5*(nthis - nstart)
	#fval = nthis - nstart - dt*ratefn(nthis, inactive, Jrates, T, M, tup,tdown, tlower, tupper) #this is the original of Chaffin
        #if eps>1e-2; updatemat=chemJmat([nthis, nochems, phrates, T, M, dt]...); end;

        #updatemat = chemJmat(nthis, inactive, Jrates, T, M, tup, tdown, tlower,tupper, dt)#this is the original of Chaffin
	    updatemat = chemJmat(nthis, inactive, Jrates, T, M, tup, tdown, tlower,tupper, dt)+ 0.5*Matrix(1.0I,length(nthis),length(nthis))
        # update
        #nthis = nthis - (updatemat \ fval) #chaffin original
	nthis = nthis + (updatemat \ fval)
        #check relative size of update
        eps = maximum(abs.(nthis-nold)./nold)
        if iter>5e2; println("eps=",eps); end; #koyama change
        iter += 1
        if iter>1e3; throw("too many iterations in next_timestep!"); end;
    end
    return nthis
end


function update!(n_current::Dict{Symbol, Array{Float64, 1}},dt)
    # update n_current using the coupled reaction network, moving to
    # the next timestep
    #set auxiliary (not solved for in chemistry) species values
    inactive = deepcopy(Float64[[n_current[sp][ialt] for sp in inactivespecies, ialt in 1:length(intaltgrid)]...])

    #set photolysis rates
    Jrates = deepcopy(Float64[n_current[sp][ialt] for sp in Jratelist, ialt in 1:length(intaltgrid)])

    #extract concentrations
    nstart = deepcopy([[n_current[sp][ialt] for sp in activespecies, ialt in 1:length(intaltgrid)]...])
    M = sum([n_current[sp] for sp in fullspecieslist])

    # TODO: Is is really necessary to do this every time update! is run?
    # set temperature and total atmospheric concentration
    T = Float64[Temp(a) for a in alt[2:end-1]]

    # take initial guess
    nthis = deepcopy(nstart)

    # get the transport rates
    tup = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[2] for sp in specieslist, a in alt[2:end-1]]
    tdown = Float64[issubset([sp],notransportspecies) ? 0.0 : fluxcoefs(a, dz, sp, n_current)[1] for sp in specieslist, a in alt[2:end-1]]

    # put the lower bcs and upper bcs in separate arrays; but they are not the
    # right shape!
    tlower_temp = [boundaryconditions(sp, dz, n_current)[1,:] for sp in specieslist]
    tupper_temp = [boundaryconditions(sp, dz, n_current)[2,:] for sp in specieslist]

    # reshape tlower and tupper into 2x2 arrays
    tlower = zeros(Float64, length(tlower_temp), 2)
    tupper = zeros(Float64, length(tupper_temp), 2)

    # tlower_temp & tupper_temp have same length; OK to use lower for the range
    for r in range(1, length(tlower_temp),step=1)
        tlower[r, :] = tlower_temp[r]
        tupper[r, :] = tupper_temp[r]
    end

    # update to next timestep
    nthis = next_timestep(nstart, nthis, inactive, Jrates, T, M, tup, tdown,
                          tlower, tupper, dt)
    nthismat = reshape(nthis,(length(activespecies),length(intaltgrid)))

    # write found values out to n_current
    for s in 1:length(activespecies)
        for ia in 1:length(intaltgrid)
            tn = nthismat[s, ia]
            n_current[activespecies[s]][ia] = tn > 0. ? tn : 0.
        end
    end

    update_Jrates!(n_current)
    ######################################koyama#########
    #here to check when it becomes CO atmos
    global totaltime #when you access global variable, you also need to declare inside the function
    totaltime += dt #calc total time run
    #println("dt: ", dt)
    #=
    if totaltime > 3.1536e7
        #print("dt: ", dt)
        #Pco = sum(n_current[:CO])*2e5*28/6.02e23 * 1.0e-3 *3.771*100*1.0e3#micro bar
        #Ph2 = sum(n_current[:H2])*2e5*2/6.02e23 * 1.0e-3 *3.771*100*1.0e3 #miro bar
        #Po2 = sum(n_current[:O2])*2e5*32/6.02e23 * 1.0e-3 *3.771*100*1.0e3#micro bar
        #Pox = 2*Po2 - Pco - Ph2
        #Htot = H_effusion_velocity*n_current[:H][end]+2*H2_effusion_velocity*n_current[:H2][end] + n_current[:H2O2][1]*0.02*2 + n_current[:HO2][1]*0.02
        #Otot = speciesbclist[:O][2,2]+n_current[:H2O2][1]*0.02*2+n_current[:O3][1]*0.02*3+n_current[:HO2][1]*0.02*2
        print(", total time: ", totaltime)
        #print(", CO mixing ratio: ",sum(n_current[:CO])/sum(map(z->n_tot(n_current, z),alt[1:end-2])))
        #print(", Pox: ", Pox)
        #println(", H/O: ", Htot/Otot)
        #print(", Htot flux: ", Htot)
        #println(", Otot flux: ", Otot)
    end
    =#
    ######################################koyama##########



end #update!

speciescolor =Dict(
   :HOCO => "#dead91",
   :H2O2 => "#bebddb",
   :HO2 => "#9e9ac8",
   :O3 => "#83d0c1",
   :OH => "#b593e2",
   :H2O => "#1f78b4",
   :O1D => "#269e56",
   :CO2pl => "#eed2d4",
   :H => "#ef3b2c",
   :H2 => "#fc9aa1",
   :O2 => "#41ae76",
   :O => "#00702d",
   :CO => "#fd8d3c",
   :Ar => "#808080",
   :N2 => "#cccccc",
   :CO2 => "#fdd0a2"
   );


function plotatm(n_current)
    clf()
    for sp in fullspecieslist
        plot(n_current[sp], alt[2:end-1]/1e5, color = speciescolor[sp],
             linewidth = 4, label=sp)
    end
    ylim(0, Int(zmax/1e5))
    xscale("log")
    xlim(1e-15, 1e22)
    xticks(size = 20)
    yticks(size = 20)
    xlabel("Species concentration [/cm3]",size = 25)
    ylabel("Altitude [km]",size=25)
    grid("on")
    legend(bbox_to_anchor=[1.01,1],loc=2,borderaxespad=0)
end

function plotReaction(n_current)
    figure()
    reactionrateshist = fill(convert(Float64,NaN),length(intaltgrid),length(reactionnet))
    reactionrateshist[:,:] = reactionrates(n_current)
    for i in 1:length(reactionnet)
        plot(reactionrateshist[:,i],alt[2:end-1]/1e5)
    end
    ylim(0, 200)
    xscale("log")
    xticks(size = 20)
    yticks(size = 20)
    legend()
end

function plotTemp()
    figure()
    clf()
    plot(map(z->Tflex(z,T_array),alt[1:end-1]),alt[1:end-1]/1e5)
    xlabel("Temperature[K]",size=25)
    ylabel("Altitude[km]",size=25)
    title("Temperature")
end

function plotH2O()
    figure()
    plot(n_current[:H2O], alt[2:end-1]/1e5, color = speciescolor[:H2O],
        linewidth = 4, label="H2O")
    ylim(0, Int(zmax/1e5))
    xscale("log")
    xlim(1e-5, 1e22)
    xticks(size = 20)
    yticks(size = 20)
    xlabel("Species concentration [/cm3]",size = 25)
    ylabel("Altitude [km]",size=25)
    grid("on")
    #legend(bbox_to_anchor=[1.01,1],loc=2,borderaxespad=0)

end

function plotDiff(n_current)
    figure()
    plot(map(z->Keddy(n_current,z),alt[1:end-1]),alt[1:end-1]/1e5,label="Keddy")
    plot(map(z->Dcoef(z,:H, n_current),alt[1:end-1]),alt[1:end-1]/1e5,label="H")
    plot(map(z->Dcoef(z,:H2, n_current),alt[1:end-1]),alt[1:end-1]/1e5,label="H2")
    plot(map(z->Dcoef(z,:O, n_current),alt[1:end-1]),alt[1:end-1]/1e5,label="O")
    xlim(1e5,1e10)
    xticks(size = 20)
    yticks(size = 20)
    xscale("log")
    xlabel("Diffusion coefficients",size=20)
    ylabel("Altitude [km]",size=20)
end



function plotflux(n_current)
    fluxhist = fill(convert(Float64,NaN),length(intaltgrid),length(specieslist))
    fluxhist[:,:] = fluxes(n_current,dz)
    #koyama I'm lazy  to separate down and up fluxes sorry
    #I'll  make it when free.
end

function plotatm()
    plotatm(n_current)
end


################################################################################
######################### PHOTOCHEMICAL CROSS SECTIONS #########################
################################################################################

# Change following line as needed depending on local machine
xsecfolder=folder_directory * "uvxsect/";

function readandskip(a, delim::Char, T::Type; skipstart=0)
    # function to read in data from a file, skipping zero or more lines at
    # the beginning.
    a = open(a,"r")
    if skipstart>0
        for i in [1:skipstart;]
            readline(a)
        end
    end
    a = readdlm(a, delim, T)
end


#CO2, temperature-dependent between 195-295K
co2xdata = readandskip(xsecfolder*"CO2.dat",'\t',Float64, skipstart = 4)
function co2xsect(T::Float64)
    T = clamp(T, 195, 295)
    #clamp(T, 195, 295) #this is not correct koyama, I changed it to like above
    Tfrac = (T-195)/(295-195)

    arr = [co2xdata[:,1]; (1-Tfrac)*co2xdata[:,2]+Tfrac*co2xdata[:,3]]
    reshape(arr, length(co2xdata[:,1]),2)
end

#CO2 photoionization (used to screen high energy sunlight)
co2exdata = readandskip(xsecfolder*"binnedCO2e.csv",',',Float64, skipstart = 4)

#H2O
h2oxdata = readandskip(xsecfolder*"h2oavgtbl.dat",'\t',Float64, skipstart = 4)

#H2O2
#the data in the table cover the range 190-260nm
h2o2xdata = readandskip(xsecfolder*"H2O2.dat",'\t',Float64, skipstart = 3)
#from 260-350 the following analytic calculation fitting the
#temperature dependence is recommended:
function h2o2xsect_l(l::Float64, T::Float64)
    l = clamp(l, 260, 350)
    T = clamp(T, 200, 400)

    A = [64761., -921.70972, 4.535649,
       -0.0044589016, -0.00004035101,
       1.6878206e-7, -2.652014e-10, 1.5534675e-13]
    B = [6812.3, -51.351, 0.11522, -0.000030493, -1.0924e-7]

    lpowA = map(n->l^n,[0:7;])
    lpowB = map(n->l^n,[0:4;])

    expfac = 1.0/(1+exp(-1265/T))

    1e-21*(expfac*reduce(+,map(*,A, lpowA))+(1-expfac)*reduce(+,map(*,B, lpowB)))
end

function h2o2xsect(T::Float64)
    retl = h2o2xdata[:,1]
    retx = 1e4*h2o2xdata[:,2]#factor of 1e4 b/c file is in 1/m2
    addl = [260.5:349.5;]
    retl = [retl; addl]
    retx = [retx; map(x->h2o2xsect_l(x, T),addl)]
    reshape([retl; retx],length(retl),2)
end



#Ozone, including IR bands which must be resampled from wavenumber
o3xdata = readandskip(xsecfolder*"O3.dat",'\t',Float64, skipstart=3)
o3ls = o3xdata[:,1]
o3xs = o3xdata[:,2]
o3chapxdata = readandskip(xsecfolder*"O3Chap.dat",'\t',Float64, skipstart=3)
o3chapxdata[:,1] = map(p->1e7/p, o3chapxdata[:,1])
for i in [round(Int, floor(minimum(o3chapxdata[:,1]))):round(Int, ceil(maximum(o3chapxdata))-1);]
    posss = getpos(o3chapxdata, x->i<x<i+1)
    dl = diff([map(x->o3chapxdata[x[1],1],posss);i])
    x = map(x->o3chapxdata[x[1],2],posss)
    ax = reduce(+,map(*,x, dl))/reduce(+,dl)
    global o3ls = [o3ls; i+0.5]
    global o3xs = [o3xs; ax]
end
o3xdata = reshape([o3ls; o3xs],length(o3ls),2)




#Oxygen, including temperature-dependent Schumann-Runge bands.
o2xdata = readandskip(xsecfolder*"O2.dat",'\t',Float64, skipstart = 3)
function binupO2(list)
    ret = Float64[];
    for i in [176:203;]
        posss = getpos(list[:,1],x->i<x<i+1)
        dl = diff([map(x->list[x[1],1],posss);i])
        x0 = map(x->list[x[1],2],posss)
        x1 = map(x->list[x[1],3],posss)
        x2 = map(x->list[x[1],4],posss)
        ax0 = reduce(+,map(*,x0, dl))/reduce(+,dl)
        ax1 = reduce(+,map(*,x1, dl))/reduce(+,dl)
        ax2 = reduce(+,map(*,x2, dl))/reduce(+,dl)
        append!(ret,[i+0.5, ax0, ax1, ax2])
    end
    return transpose(reshape(ret, 4, 203-176+1))
end
o2schr130K = readandskip(xsecfolder*"130-190.cf4",'\t',Float64, skipstart = 3)
o2schr130K[:,1] = map(p->1e7/p, o2schr130K[:,1])
o2schr130K = binupO2(o2schr130K)
o2schr190K = readandskip(xsecfolder*"190-280.cf4",'\t',Float64, skipstart = 3)
o2schr190K[:,1] = map(p->1e7/p, o2schr190K[:,1])
o2schr190K = binupO2(o2schr190K)
o2schr280K = readandskip(xsecfolder*"280-500.cf4",'\t',Float64, skipstart = 3)
o2schr280K[:,1] = map(p->1e7/p, o2schr280K[:,1])
o2schr280K = binupO2(o2schr280K)

function o2xsect(T::Float64)
    o2x = deepcopy(o2xdata);
    #fill in the schumann-runge bands according to Minschwaner 1992
    T = clamp(T, 130, 500)
    if 130<=T<190
        o2schr = o2schr130K
    elseif 190<=T<280
        o2schr = o2schr190K
    else
        o2schr = o2schr280K
    end

    del = ((T-100)/10)^2

    for i in [176.5:203.5;]
        posO2 = findfirst_or_zero(i,o2x[:,1]) #findfirst(o2x, i)
        posschr = findfirst_or_zero(i,o2schr[:,1])
        o2x[posO2, 2] += 1e-20*(o2schr[posschr, 2]*del^2
                                + o2schr[posschr, 3]*del
                                + o2schr[posschr, 4])
    end

    # add in the herzberg continuum (though tiny)
    # measured by yoshino 1992
    for l in [192.5:244.5;]
        posO2 = findfirst_or_zero(l,o2x[:,1])
        o2x[posO2, 2] += 1e-24*(-2.3837947e4
                            +4.1973085e2*l
                            -2.7640139e0*l^2
                            +8.0723193e-3*l^3
                            -8.8255447e-6*l^4)
    end

    return o2x
end

#HO2
function ho2xsect_l(l::Float64)
    #function to compute HO2 cross-section as a function of wavelength l
    #in nm, as given by Sander 2011 JPL Compilation
    a = 4.91
    b = 30612.0
    sigmamed = 1.64e-18
    vmed = 50260.0
    v = 1e7/l;
    if 190<=l<=250
        return HO2absx = sigmamed / ( 1 - b/v ) * exp( -a * log( (v-b)/(vmed-b) )^2 )
    else
        return 0.0
    end
end
ho2xsect = [190.5:249.5;]
ho2xsect = reshape([ho2xsect; map(ho2xsect_l, ho2xsect)],length(ho2xsect),2)


#H2
h2xdata = readandskip(xsecfolder*"binnedH2.csv",',',Float64, skipstart=4)

#OH
ohxdata = readandskip(xsecfolder*"binnedOH.csv",',',Float64, skipstart=4)
ohO1Dxdata = readandskip(xsecfolder*"binnedOHo1D.csv",',',Float64, skipstart=4)

#SOLAR FLUX
# Change following line as needed depending on local machine

const solarflux=readandskip(folder_directory*"marssolarphotonflux.dat",'\t',Float64,skipstart=4)[1:2000,:]
solarflux[:,2] = solarflux[:,2]/2 #(daytime + night) /2

#### Ribas+2005 correction ######
#The least-squares fit of the power-law exponent beta taken from Claire+2012 but originally from Ribas+2005
#The validated range is limited between 0-120nm but is extrapolated in longer wavelength
#t[Gyr], 4.56Ga=0
function flux_multi(t,l)
    beta = -1.21 + l/200. #least-sqrt fit Claire+2012, Fig1
    return (t/4.56)^beta
end

absorber = Dict(:JCO2ion =>:CO2,
                :JCO2toCOpO =>:CO2,
                :JCO2toCOpO1D =>:CO2,
                :JO2toOpO =>:O2,
                :JO2toOpO1D =>:O2,
                :JO3toO2pO =>:O3,
                :JO3toO2pO1D =>:O3,
                :JO3toOpOpO =>:O3,
                :JH2toHpH =>:H2,
                :JOHtoOpH =>:OH,
                :JOHtoO1DpH =>:OH,
                :JHO2toOHpH =>:HO2,
                :JH2OtoHpOH =>:H2O,
                :JH2OtoH2pO1D =>:H2O,
                :JH2OtoHpHpO =>:H2O,
                :JH2O2to2OH =>:H2O2,
                :JH2O2toHO2pH =>:H2O2,
                :JH2O2toH2OpO1D =>:H2O2);

function padtosolar(crosssection::Array{Float64, 2})
    # a function to take an Nx2 array and pad it with zeroes until it's the
    # same length as the solar flux. Returns the cross sections only, as
    # the wavelengths are shared by solarflux

    positions = map(x->findfirst_or_zero(x,solarflux[:,1]),crosssection[:,1])
    retxsec = fill(0.,length(solarflux[:,1]))
    retxsec[positions] = crosssection[:,2]
    retxsec
end

function quantumyield(xsect::Array, arr)
    #function to assemble cross-sections for a given pathway. Inputs are
    #an Nx2 array xsect with wavelengths and photoabsorption cross
    #sections, and arr, a tuple of tuples with a condition and a quantum
    #yield multiplicative factor, either constant or a function of
    #wavelength in the given regime. Return is an array with all of the
    #matching wavelengths and the scaled cross-sections.
    lambdas = Float64[];
    rxs = Float64[];
    for (cond, qeff) in arr
        places = findall(cond, xsect[:,1])
        append!(lambdas, xsect[places, 1])
        #if we have a number then map to a function
        isa(qeff, Function) ? (qefffn = qeff) : (qefffn = x->qeff)
        append!(rxs, map(*,map(qefffn, xsect[places, 1]),xsect[places, 2]))
    end

    reshape([lambdas; rxs],length(lambdas),2)
end

#build the array of cross-sections
crosssection = Dict{Symbol, Array{Array{Float64}}}()
## #this is a dictionary of the 1-nm photodissociation or photoionization
## #cross-sections important in the atmosphere. keys are symbols found in
## #jratelist. each entry is an array of arrays, yielding the wavelengths
## #and cross-sections for each altitude in the atmosphere.
## #
## #NOTE: jspecies refers to the photodissociation or photoionization
## #cross section for a particular species which produces a UNIQUE SET OF
## #PRODUCTS. In this sense, crosssection has already folded in quantum
## #efficiency considerations.

#now add the cross-sections

#CO2 photoionization
setindex!(crosssection,
          fill(co2exdata, length(alt)),
          :JCO2ion)
#CO2+hv->CO+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->l>167, 1),
                                   (l->95>l, 0.5))), map(t->co2xsect(t),map(Temp, alt))),
          :JCO2toCOpO)
#CO2+hv->CO+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((l->95<l<167, 1),
                                   (l->l<95, 0.5))), map(t->co2xsect(t),map(Temp, alt))),
          :JCO2toCOpO1D)
#O2+hv->O+O
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x>175, 1),)), map(t->o2xsect(t),map(Temp, alt))),
          :JO2toOpO)
#O2+hv->O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<175, 1),)), map(t->o2xsect(t),map(Temp, alt))),
          :JO2toOpO1D)

# The quantum yield of O1D from ozone photolysis is actually
# well-studied! This adds some complications for processing.
function O3O1Dquantumyield(lambda, temp)
    if lambda < 306. || lambda > 328.
        return 0.
    end
    temp=clamp(temp, 200, 320)#expression is only valid in this T range

    X = [304.225, 314.957, 310.737];
    w = [5.576, 6.601, 2.187];
    A = [0.8036, 8.9061, 0.1192];
    v = [0.,825.518];
    c = 0.0765;
    R = 0.695;
    q = exp.(-v/(R*temp))
    qrat = q[1]/(q[1]+q[2])

    (q[1]/sum(q)*A[1]*exp.(-((X[1]-lambda)/w[1])^4.)
     +q[2]/sum(q)*A[2]*(temp/300.)^2. *exp.(-((X[2]-lambda)/w[2])^2.)
     +A[3]*(temp/300.)^1.5*exp.(-((X[3]-lambda)/w[3])^2.)
     +c)
end

# O3+hv->O2+O
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1-(1.37e-2*193-2.16)),
                               (l->193<=l<225, l->(1. -(1.37e-2*l-2.16))),
                               (l->225<=l<306, 0.1),
                               (l->306<=l<328, l->(1. -O3O1Dquantumyield(l, t))),
                               (l->328<=l<340, 0.92),
                               (l->340<=l, 1.0)
                               ))
              ,map(Temp, alt)
              ),
          :JO3toO2pO)

# O3+hv->O2+O1D
setindex!(crosssection,
          map(t->quantumyield(o3xdata,
                              (
                               (l->l<193, 1.37e-2*193-2.16),
                               (l->193<=l<225, l->(1.37e-2*l-2.16)),
                               (l->225<=l<306, 0.9),
                               (l->306<=l<328, l->O3O1Dquantumyield(l, t)),
                               (l->328<=l<340, 0.08),
                               (l->340<=l, 0.0)
                               ))
              ,map(Temp, alt)
              ),
          :JO3toO2pO1D)

# O3+hv->O+O+O
setindex!(crosssection,
          fill(quantumyield(o3xdata,((x->true, 0.),)),length(alt)),
          :JO3toOpOpO)

# H2+hv->H+H
setindex!(crosssection,
          fill(h2xdata, length(alt)),
          :JH2toHpH)

# OH+hv->O+H
setindex!(crosssection,
          fill(ohxdata, length(alt)),
          :JOHtoOpH)

# OH+hv->O1D+H
setindex!(crosssection,
          fill(ohO1Dxdata, length(alt)),
          :JOHtoO1DpH)

# HO2+hv->OH+H
setindex!(crosssection,
          fill(ho2xsect, length(alt)),
          :JHO2toOHpH)

# H2O+hv->H+OH
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.89),(x->x>145, 1))),length(alt)),
          :JH2OtoHpOH)
# H2O+hv->H2+O1D
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->x<145, 0.11),(x->x>145, 0))),length(alt)),
          :JH2OtoH2pO1D)

# H2O+hv->H+H+O
setindex!(crosssection,
          fill(quantumyield(h2oxdata,((x->true, 0),)),length(alt)),
          :JH2OtoHpHpO)

# H2O2+hv->OH+OH
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.85),(x->x>230, 1))), map(t->h2o2xsect(t),map(Temp, alt))),
          :JH2O2to2OH)

# H2O2+hv->HO2+H
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->x<230, 0.15),(x->x>230, 0))), map(t->h2o2xsect(t),map(Temp, alt))),
          :JH2O2toHO2pH)

# H2O2+hv->H2O+O1D
setindex!(crosssection,
          map(xs->quantumyield(xs,((x->true, 0),)), map(t->h2o2xsect(t),map(Temp, alt))),
          :JH2O2toH2OpO1D)

lambdas = Float64[]
for j in Jratelist, ialt in 1:length(alt)
    global lambdas
    lambdas = union(lambdas, crosssection[j][ialt][:,1])
end

if !(setdiff(solarflux[:,1],lambdas)==[])
    throw("Need a broader range of solar flux values!")
end

#pad all cross-sections to solar
for j in Jratelist, ialt in 1:length(alt)
    crosssection[j][ialt] = padtosolar(crosssection[j][ialt])
end

# we need some global objects for the Jrates calculation:
# intensity as a function of wavelength at each altitude

# this is the unitialized array for storing values
# TODO: do we even need this here?
solarabs = fill(fill(0.,size(solarflux, 1)),length(alt)-2);


function update_Jrates!(n_current::Dict{Symbol, Array{Float64, 1}})
    #this function updates the photolysis rates stored in n_current to
    #reflect the altitude distribution of absorbing species
    #    global solarabs::Array{Array{Float64, 1},1}

    solarabs = Array{Array{Float64}}(undef,length(alt)-2)
    for i in range(1, length(alt)-2,step=1)
        solarabs[i] = zeros(Float64, 2000)
    end

    nalt = size(solarabs, 1)
    nlambda = size(solarabs[1],1)

    for jspecies in Jratelist
        species = absorber[jspecies]

        # println(string(jspecies," is absorbed by ",species))
        jcolumn = 0.
        for ialt in [nalt:-1:1;]
            #get the vertical column of the absorbing constituient
            jcolumn += n_current[species][ialt]*dz
            # if jspecies==:JO2toOpO
            #     println(string("At alt = ",alt[ialt+1],
            #                    ", n_",species," = ",n_current[species][ialt],
            #                    ", jcolumn = ",jcolumn))
            #     println("and solarabs[ialt] is $(solarabs[ialt]) before we do axpy")
            #     readline(STDIN)
            # end

            # add the total extinction to solarabs:
            # multiplies air column density at all wavelengths by crosssection
            # to get optical depth
            BLAS.axpy!(nlambda, jcolumn, crosssection[jspecies][ialt+1], 1,
                       solarabs[ialt],1)
        end
    end

    #solarabs now records the total optical depth of the atmosphere at
    #each wavelength and altitude

    # actinic flux at each wavelength is solar flux diminished by total
    # optical depth
    for ialt in [1:nalt;]
        solarabs[ialt] = solarflux[:,2].*exp.(-solarabs[ialt])
    end

    #each species absorbs according to its cross section at each
    #altitude times the actinic flux.
    for j in Jratelist
        for ialt in [1:nalt;]
            n_current[j][ialt] = BLAS.dot(nlambda, solarabs[ialt], 1,
                                          crosssection[j][ialt+1], 1)
        end
    end
    #for ialt in [1:nalt;]#koyama added this part to modify an error that photodissociation rate goes to minus
    #    if n_current[:JCO2toCOpO][ialt] < 0
    #        n_current[:JCO2toCOpO][ialt] = -1*n_current[:JCO2toCOpO][ialt]
    #    end
    #end #koyama until here
end

function timeupdate(mytime)
    for i = 1:50
        #println("dt: ", mytime)
        update!(n_current, mytime)
        #plotatm()
    end
    println("dt: ", mytime)
    println(totaltime)
    print_c()
    plotatm()
    #show()
    ## yield()
end

function timeupdate100(mytime)
    for i = 1:100
        #plotatm()
        #println("dt: ", mytime)
        if i%10 == 0
            println("total: ",totaltime)
            plotatm()
        end
        update!(n_current, mytime)
    end
    println("dt: ", mytime)
    println(totaltime)
    print_c()
    plotatm()
    ## show()
    ## yield()
end

function timeupdate1000(mytime)
    for i = 1:1000
        #plotatm()
        #println("dt: ", mytime)
        if i%100 == 0
            println("total: ",totaltime)
            plotatm()
            print_c()
        end
        update!(n_current, mytime)
    end
    println("dt: ", mytime)
    println(totaltime)
    print_c()
    plotatm()
    ## show()
    ## yield()
end

function timeupdate500(mytime)
    for i = 1:500
        #plotatm()
        #println("dt: ", mytime)
        if i%100 == 0
            println("total: ",totaltime)
            plotatm()
            print_c()
        end
        update!(n_current, mytime)
    end
    println("dt: ", mytime)
    println(totaltime)
    print_c()
    plotatm()
    ## show()
    ## yield()
end

function print_c()
    println("totaltime: ", totaltime)
    println("CO2 outgassing: ",eval(speciesbclist[:CO2][1,2]))
    println("C esc in run: ",eval(speciesbclist[:CO][2,2])*n_current[:CO][end])
    #println("Oreturn: ",Oreturn_list[end])
    println("O thermal esc: ", calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:O][end],16.)*flex_effusion_velocity(Texo_tian,16.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,16.0))
    println("H esc in run: ",eval(speciesbclist[:H][2,2])*n_current[:H][end] + 2*eval(speciesbclist[:H2][2,2])*n_current[:H2][end] )
    println("C effusion: ",eval(speciesbclist[:CO][2,2]))
    println("O2 mixing: ", sum(n_current[:O2])/sum(map(z->n_tot(n_current, z),alt[1:end-2])))
    println("O effusion: ",eval(speciesbclist[:O][2,2]))
    println("O total flux: ",eval(speciesbclist[:O][2,2])*n_current[:O][end])
    println("O2 pressure:", calc_O2_pressure(n_current[:O2]))
    println("CO2 pressure:", calc_CO2_pressure(n_current[:CO2]))

    #plotatm(n_current)
end

#Grott et al., 2011
#R19
#fp=0.01 case
# t(Myr): time from the birth of Mars
function grott(t)
    IW=1.
    A = 224.39
    a=1505.1
    alpha=2.7606
    beta=3.3600
    xi=0.4
    pCO2 = 10^(IW)*xi*A*(tanh((t/a)^alpha))^(1/beta)
    return pCO2
end


function printkey()
    Hesc = eval(speciesbclist[:H][2,2])*n_current[:H][end]+2*eval(speciesbclist[:H2][2,2])*n_current[:H2][end]
    COmixing = sum(n_current[:CO])/sum(map(z->n_tot(n_current, z),alt[1:end-2]))
    Oesc = calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:O][end],16.)*O_effusion_velocity*modified_factor(n_current,Texo_tian,Zexo_tian,16.0)
    Cesc = calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:CO][end]*C_CO2_ratio,12.)*C_CO_ratio*flex_effusion_velocity(Texo_tian,12.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,12.0)
    #println(totaltime)
    H2O2dep = n_current[:H2O2][1]*speciesbclist[:H2O2][1,2]
    O3dep =  n_current[:O3][1]*speciesbclist[:O3][1,2]
    HO2dep = n_current[:HO2][1]*speciesbclist[:HO2][1,2]
    #HOCOdep = n_current[:HOCO][1]*speciesbclist[:HOCO][1,2]
    OHdep = n_current[:OH][1]*speciesbclist[:OH][1,2]
    Hdep = n_current[:H][1]*speciesbclist[:H][1,2]
    O1Ddep = n_current[:O1D][1]*speciesbclist[:O1D][1,2]
    Odep =  n_current[:O][1]*speciesbclist[:O][1,2]
    fO2 = sum(n_current[:O2])/sum(map(z->n_tot(n_current, z),alt[1:end-2]))
    fH2 = sum(n_current[:H2])/sum(map(z->n_tot(n_current, z),alt[1:end-2]))
    fCO = sum(n_current[:CO])/sum(map(z->n_tot(n_current, z),alt[1:end-2]))

    # C-involved boundary
    CO2outgassing = eval(speciesbclist[:CO2][1,2])
    #Cesc = eval(speciesbclist[:CO2][2,2]) * n_current[:CO2][end]
    #Cesc = eval(speciesbclist[:CO][2,2]) * n_current[:CO][end] #CO escape ver
    Oreturn = -eval(speciesbclist[:O][2,2])

    println("CO2 outgassing: ", CO2outgassing)
    println("C thermal escape: ", Cesc)
    println("C nonthermal: ", C_nonthermal)
    println("O thermal: ", Oesc)
    println("O nonthermal: ", O_nonthermal)

    println("O2 mixing ratio: ", fO2)
    println("H2 mixing ratio: ", fH2)
    println("CO mixing ratio: ", fCO)
    println("Hesc: ", Hesc)
    println("CO mixing ratio: ",COmixing)
    println("H2O2 dep flux: ", H2O2dep)
    println("O3 dep flux: ", O3dep)
    println("HO2 dep flux:", HO2dep)
    #println("HOCO dep flux: ", HOCOdep)
    println("H dep flux: ", Hdep)
    Htotalflux = Hesc + H2O2dep*2 + HO2dep + OHdep + Hdep
    Ototalflux = Oesc + O_nonthermal+ H2O2dep*2 + O3dep*3 + HO2dep*2 + OHdep + Odep + O1Ddep
    Ctotalflux = Cesc + C_nonthermal
    #redox H or O flux is net flux considering H2O is neutral species. ex) H2O2 can be interpreted as one oxygen atom.
    redox_Hflux = Hesc + Hdep
    redox_Oflux = eval(speciesbclist[:O][2,2]) + H2O2dep + O3dep*3 + 1.5*HO2dep + 0.5*OHdep + Odep + O1Ddep
    println("H total outflux: ", Htotalflux)
    println("O total outflux: ", Ototalflux)
    println("C total ouflux: ",Ctotalflux)
    #println("HOratio: ", Htotalflux/Ototalflux)
    println("H/2O: ",Htotalflux/(2*Ototalflux))
    println("2C/O: ",2*Ctotalflux/(Ototalflux))
    println("(H+4C)/2O: ",(Htotalflux+4*Ctotalflux)/(2*Ototalflux))
    #println("reodx H outflux: ", redox_Hflux)
    #println("redox O outflux: ", redox_Oflux)
    #println("redox HO ratio: ", redox_Hflux/redox_Oflux)
    return COmixing, Htotalflux/Ototalflux
end

function writeredox(file_name)
    Hesc = eval(speciesbclist[:H][2,2])*n_current[:H][end]+2*eval(speciesbclist[:H2][2,2])*n_current[:H2][end]
    #COmixing = sum(n_current[:CO])/sum(map(z->n_tot(n_current, z),alt[1:end-2]))
    Oesc = calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:O][end],16.)*O_effusion_velocity*modified_factor(n_current,Texo_tian,Zexo_tian,16.0)
    Cesc = calc_rho_exo(Zexo_tian,Texo_tian,Tinf,n_current[:CO][end]*C_CO2_ratio,12.)*C_CO_ratio*flex_effusion_velocity(Texo_tian,12.0, Zexo_tian)*modified_factor(n_current,Texo_tian,Zexo_tian,12.0)

    H2O2dep = n_current[:H2O2][1]*speciesbclist[:H2O2][1,2]
    O3dep =  n_current[:O3][1]*speciesbclist[:O3][1,2]

    # C-involved boundary
    CO2outgassing = eval(speciesbclist[:CO2][1,2])
    pCO2 = calc_CO2_pressure(n_current[:CO2])

    Htotalflux = Hesc + H2O2dep*2
    Ototalflux = Oesc + O_nonthermal+ H2O2dep*2 + O3dep*3
    Ctotalflux = Cesc + C_nonthermal
    H_O = Htotalflux/(2*Ototalflux)
    C_O = 2*Ctotalflux/(Ototalflux)
    H_C_O = (Htotalflux+4*Ctotalflux)/(2*Ototalflux)
    println("H/2O: ",Htotalflux/(2*Ototalflux))
    println("2C/O: ",2*Ctotalflux/(Ototalflux))
    println("(H+4C)/2O: ",(Htotalflux+4*Ctotalflux)/(2*Ototalflux))
    key_parameters=[feuv, pCO2,surface_temperature,Hesc,Oesc,Cesc,O_nonthermal,C_nonthermal,H2O2dep,O3dep,Htotalflux,Ototalflux,Ctotalflux,H_O,C_O,H_C_O]
    key_parameters=reshape(key_parameters,1,length(key_parameters))
    file_name = folder_directory * file_name
    open(file_name, "a") do io
        writedlm(io,key_parameters,',')
    end
end


#convert density into CO2 pressure
function calc_CO2_pressure(CO2_density_array)
    G = 6.674e-11
    M = 0.64171e24
    p = 0
    #alt = [(j+1)*2*1000 for j in 1:length(CO2_density_array)]
    for i in 1:length(CO2_density_array)
        p += CO2_density_array[i]*2e5*44/(6e23)*1e-3*G*M/((3389500+alt[i])^2) #N/cm^2
    end
    return p*100 #mbar
end

#calc O2 pressure
function calc_O2_pressure(O2_density_array)
    G = 6.674e-11
    M = 0.64171e24
    p = 0
    #alt = [(j+1)*2*1000 for j in 1:length(CO2_density_array)]
    for i in 1:length(O2_density_array)
        p += O2_density_array[i]*2e5*32/(6e23)*1e-3*G*M/((3389500+alt[i])^2) #N/cm^2
    end
    return p*100 #mbar
end

#calc any species pressure
function calc_pressure(density_array,sp_mass)
    G = 6.674e-11
    M = 0.64171e24
    p = 0
    #alt = [(j+1)*2*1000 for j in 1:length(CO2_density_array)]
    for i in 1:length(density_array)
        p += density_array[i]*2e5*sp_mass/(6e23)*1e-3*G*M/((3389500+alt[i])^2) #N/cm^2
    end
    return p*100 #mbar
end

#header = ["CO2pressure", "Ts","H2Oppm","Oesc","depv","N2Ar","CO_01G","HOratio_01G","CO_1G","HOratio_1G"]
#header = reshape(header,1,length(header))
#This function saves parameters and results of COmixing and HOratio
function writekeyparameters(filename,COmixing_01G,HOratio_01G,COmixing_1G,HOratio_1G)
    fO2 = sum(n_current[:O2])/sum(map(z->n_tot(n_current, z),alt[1:end-2]))
    fH2 = sum(n_current[:H2])/sum(map(z->n_tot(n_current, z),alt[1:end-2]))
    fOx = 2*fO2 - fH2 - COmixing_1G
    key_parameters=[CO2_pressure,surface_temperature,H2Oppm,Oxygen_escape_rate, depv, Int(N2_Ar),COmixing_01G,HOratio_01G,COmixing_1G,HOratio_1G,fO2,COmixing_1G,fH2,fOx]
    key_parameters=reshape(key_parameters,1,length(key_parameters))
    file_name = folder_directory * filename
    open(file_name, "a") do io
        writedlm(io,key_parameters,',')
    end
end

function write_rates(n_current,filename)
    reactionrateshist = fill(convert(Float64,NaN),length(intaltgrid),length(reactionnet))
    reactionrateshist[:,:] = reactionrates(n_current)
    # reaction rates are saved in the same order of reactionnet
    h5write(filename,"rates/reaction_rates",reactionrateshist)
end

function write_fluxes(n_current,filename)
    fluxeshist = fill(convert(Float64,NaN),length(intaltgrid),length(specieslist))
    fluxeshist[:,:] = fluxes(n_current,dz)
    h5write(filename,"fluxes/flux",fluxeshist)
end

################################################################################
############################### CONVERGENCE CODE ###############################
################################################################################
# Extra code to reach convergence to equilibrium over millions of years
# STANDARD WATER CASE ----------------------------------------------------------
#n_current[:H2O]=H2Oinitfrac.*map(z->n_tot(n_current, z),alt[2:end-1])

n_current[:H2O]=H2Ofrac.*map(z->n_tot(n_current, z),alt[1:end-2]) #koyama H2O profile relative humidity: 25% below 30km

#if N2_Ar is false, remove N2 and Ar.
if N2_Ar == false
    n_current[:N2][:] .= 0.
    n_current[:Ar][:] .= 0.
end
H2Oppm = sum(n_current[:H2O])*2e5/3.34e18
println(sum(n_current[:H2O])*2e5/3.34e18)

####shift CO2+ according to altitude of ztropo######
#CO2+

if get_converged
    n_CO2pl_temp = get_ncurrent(readfile)[:CO2pl]
    n_current[:CO2pl][z_CO2pl_index:(length(alt)-2)] = n_CO2pl_temp[41:(length(alt)-2)-(z_CO2pl_index-41)]#it shifts CO2+ to upper according to each tropopause
end


#n_current[:CO2pl][41:z_CO2pl_index] = zero(0.0)
#####################################################

global totaltime=0.
plotTemp()
plotDiff(n_current)

#this part is to get a converged state
if get_converged
    [timeupdate(t) for t in [10.0^(1.0*i) for i in -3:13]]
    for i in 1:26
        update!(n_current,1.0e14)
        plotatm()
    end
    written_file_name =  folder_directory *
                        "converged_Ts_" * string(surface_temperature) * "_CO2_" *
                        string(CO2_pressure) * "mbar_Oesc_" * string(Oxygen_escape_rate) *
                        "_depv_" * string(depv)
    #write_ncurrent(n_current,written_file_name*"_0.1billion.h5") #here save ncurrent of 0.1G years model run
    #println(written_file_name)
    println("0.1Gyr")
    #COmixing_01G, HOratio_01G = printkey()
    for i in 1:24
        update!(n_current,1.0e14)
        plotatm()
    end
    for i in 1:26
        update!(n_current,1.0e15)
        plotatm()
    end
    println("1Gyr")
    #COmixing_1G, HOratio_1G = printkey()
    write_ncurrent(n_current,written_file_name*"_c.h5") #save 1 Gyr model run
    write_rates(n_current,written_file_name*"_c.h5")
    printkey()
    #plotReaction(n_current)
    #plotTemp()
    #plotDiff(n_current)
    plotatm()

    #writekeyparameters("summary_conv_Tmod.csv", COmixing_01G,HOratio_01G,COmixing_1G,HOratio_1G)
end
##ここで1億年


#timeupdate(1.0e14)
#timeupdate(1.0e15)
#write_ncurrent(n_current,"converged_Tsurf_220_CO2_250mbar_Oesc_1.2e8_depv_0.02_nofixedCO2.h5")
#write_ncurrent(n_current,"converged_Tsurf_270_CO2_1bar_Oesc_1.2e8_depv_0.02_nofixedCO2_noHOCO.h5")

#[timeupdate(t) for t in [10.0^(1.0*i) for i in 11:16]]
#write_ncurrent(n_current,"converged_Tsurf_230_CO2_7bar_Oesc_1e7_depv_0.02_tempdt8.h5")
#println(readfile)
#printkey()
#write_ncurrent(n_current,"converged_Tsurf_250_CO2_1bar_Oesc_1e9_depv_0.02_nofixedCO2_tempdt8_HOCO1e-4.h5")

## n_current[:H2O]=detachedlayer.*map(z->n_tot(n_current, z),alt[2:end-1])
## [timeupdate(t) for t in [10.0^(1.0*i) for i in -3:14]]
#for i in 1:100
#    update!(n_current, 1e14)
#end
## write_ncurrent(n_current,"converged_highwater.h5")
