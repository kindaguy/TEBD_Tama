function envZeroBasisUp(inds,nchain)
   envstate0 = ["0" for _ in 1:nchain]
   refUp = productMPS(ComplexF64,inds,vcat(["Up"],envstate0));
   return convert_to_Vidal(refUp);
end


function envZeroBasisDn(inds,nchain)
   envstate0 = ["0" for _ in 1:nchain]
   refDn = productMPS(ComplexF64,inds,vcat(["Dn"],envstate0));
   return convert_to_Vidal(refDn);
end

"""
   twoSiteMod(inds,refGammas,site1,op1,site2,op2)

Prepares the state |ψ'> = op1 op2 |ψ> with op1 and op2 applied to sites
site1 and site2 respectively
"""

"""createSingleEx(inds,refGammas,site): Returns Vidal's Gammas matrices for the site-th element of the single excitation subspace basis
   Note that Lambda matrices are always trivial (factorized state)
"""
function createSingleEx(inds,refGammas, site)
   appo = [a for a in refGammas]
   appo[site+1] = noprime!(op("Adag",inds[site+1])*refGammas[site+1])
   return appo
end

"""createSingleExBasis(inds,refGammas,chain_length):returns the single excitation subspace basis starting from the seed refGammas,
that model an environment in the vacuum state
"""
function createSingleExBasis(inds,refGammas,chain_length)
   return [createSingleEx(inds,refGammas,i) for i in 1:chain_length]
end


"""createDoubleExBasis(inds,refGammas,chain_length):returns the single excitation subspace basis starting from the seed refGammas,
that model an environment in the vacuum state
"""
function createDoubleExBasis(inds,refGammas,chain_length)
   doubleExBasis = Vector{Vector{ITensor}}([])
   for i in 1:chain_length
      appo = createSingleEx(inds,refGammas,i)
      for j in i:chain_length
         pippo = createSingleEx(inds,appo,j)
         push!(doubleExBasis,pippo)
      end
   end
   return doubleExBasis
end

"""indDoubleExBasis(chain_length):returns a k =>(i,j) dictionary
"""
function indDoubleExBasis(chain_length)
   Dict{Int64,Tuple{Int64,Int64}}()
   veci = Vector{Int64}([])
   for i in 1:chain_length
      for j in i:chain_length
         push!(veci,1)
         a=length(veci)
         push!(inds,a=>(i,j))
      end
   end
   return inds
end




"""project(LambdasPsi,GammasPsi,LambdasPhi,GammasPhi) performs the projection <ϕ|ψ>
"""
#Note: this function remains always the same: does not matter what's ϕ
function project(stateLambdas,stateGammas,projLambdas,projGammas)
   ll = size(projLambdas)[1] 
   appo = stateGammas[1]*stateLambdas[1]*dag(projGammas[1]*projLambdas[1])
   for i in 2:ll  
      appo = appo*stateGammas[i]*stateLambdas[i]*dag(projGammas[i]*projLambdas[i])
   end
   appo = appo * stateGammas[ll+1]*dag(projGammas[ll+1])
   return scalar(appo)
end
