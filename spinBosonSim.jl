import Pkg
Pkg.activate(".")
Pkg.instantiate()
using ITensors, ITensorMPS
using DelimitedFiles
using LinearAlgebra
using MKL

#Includes
include("Lib/convert_to_Vidal.jl")
include("Lib/tebd.jl")

#System parameters

#System
ϵ = -200. #ϵ/2 σ_z
Δ = 100. #Δ/2 σ_x

#Chain
chain_size = 400 #100 hosc
local_dim = 10 #local dimension

#Times

#Integration step
#The integration step should be set by looking at the 
#asymptotic coupling of the chain modes.
τ= 1 / 1000 / 50
#Simulation time
tmax=0.2
#Measurement graining
mStep = 50

#MPS pars
minBondDim = 10

Threads.nthreads()

#Load chain coefficients
freqs = readdlm("Data/FMO_T77_freqs.dat",'\n' ,Float64); #chain fequencies
coups = readdlm("Data/FMO_T77_coups.dat",'\n',Float64);  #chain couplings; first element:sys-bath coupling
println(length(freqs))
println(length(coups))

#Define indices
#System
sys = siteinds("S=1/2",1);
#Environment
env = [Index(local_dim,"Boson,Site,n=$(i+1)") for i in 1:chain_size]
#System+environment
sysenv = vcat(sys,env);

stateSys = ["Up"]   

#Standard approach: chain always in the vacuum state
stateEnv = ["0" for n=1:chain_size];

stateSE = vcat(stateSys,stateEnv);

psi0 = productMPS(ComplexF64,sysenv,stateSE);
psi = psi0;

Lambdas, Gammas = convert_to_Vidal(psi)

@show length(Lambdas)
@show length(Gammas)

#Ready to go!

SpinBoson_evolution_TEBD(Gammas, Lambdas, sysenv;
ϵ = ϵ,
Δ = Δ , 
sysenvInt = "Z",
ChainLength = chain_size, 
tau = τ, 
ttotal = tmax,
measStep = mStep,
freqs = freqs, 
coups = coups, 
minBondDim = 5, 
cutoff=1E-12)