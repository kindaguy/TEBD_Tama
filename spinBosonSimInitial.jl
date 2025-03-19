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
include("Lib/prepareState.jl")


#System parameters

#System
ϵ = -120. #ϵ/2 σ_z
Δ = -175.4 #Δ/2 σ_x

#Chain
chain_size = 400 #100 hosc
local_dim = 5 #local dimension

#Times

#Integration step
#The integration step should be set by looking at the 
#asymptotic coupling of the chain modes.
τ= 1 / 1000 / 50
#Simulation time
tmax=0.05
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

#Define indices and initial state
sysenv,psi0 = prepareStateMonomer(chain_size,local_dim,filename="Input/specialPairHigh.dat");
#System

psi = psi0;

Lambdas, Gammas = convert_to_Vidal(psi)

@show length(Lambdas)
@show length(Gammas)

#Ready to go!

newGammas,newLambdas = SpinBoson_evolution_TEBD(Gammas, Lambdas, sysenv;
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
cutoff=1E-12,
performProj=true)