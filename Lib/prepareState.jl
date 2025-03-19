# import Pkg
# Pkg.activate(".")
# Pkg.instantiate()
# using ITensors, ITensorMPS
# using DelimitedFiles
# using LinearAlgebra

"""Prepares the state of the monomer; if kwarg filename the state of the monomer,i.e. two complex values, from a file.
The environment starts from the factorized vaccum state.
Parameters:
oscDim: local dimension of the oscillators; 
nchain: number of chain modes; 

kwargs:
filename: file containing initial state
state: standard state "Up","Dn",



an instantiates the corresponding tensor."""

function prepareStateMonomer(
   nchain::Integer,
   oscDim::Integer;
   filename::Union{AbstractString,Nothing}=nothing, 
   state::Union{AbstractString,Nothing}=nothing)
   
   #Environment is always the same
   #NOTE:n-index set consistently with single site system
   env = [Index(oscDim,"Boson,Site,n=$(i+1)") for i in 1:nchain]
   #Set the initial state (vacuum)
   stateenv = ["0" for _ in 1:nchain]

   #Prepare an appo system: initial state according to parameters
   apposys = siteinds("S=1/2",1)
   appostate = [state !==nothing ? state : "Up"]
   res = productMPS(ComplexF64,vcat(apposys,env),vcat(appostate,stateenv))

   if filename !== nothing
      values = readdlm(filename,',',ComplexF64,'\n')
      #Index
      i = Index(2,"S=1/2,Site,n=1")
      #Turn to ITensor
      sys = ITensor(i)
      #Assign values
      sys[i=>1] = values[1]
      sys[i=>2] = values[2]

      #Replace the first tensor
      res[1] = sys
      #WHATCH out: this changes the indices!
      indices = siteinds(res,tags="Site")
      return indices,res

      #Now the tensor we need is ready.
      #The easiest way to build the tensor product between this (system) tensor and the remaining (env) tensor
      #is described here:
      #https://itensor.discourse.group/t/tensor-product-together-mps/824/2
      
      

      #Finally: prepare the state we need.
      #Create MPS for the given initial state
      
   else

      return vcat(apposys,env),res
      
   end
end

