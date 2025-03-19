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

"""project(LambdasPsi,GammasPsi,LambdasPhi,GammasPhi) performs the projection <ϕ|ψ>
"""

function project(stateLambdas,stateGammas,projLambdas,projGammas)
   ll = size(projLambdas)[1] 
   appo = stateGammas[1]*stateLambdas[1]*dag(projGammas[1]*projLambdas[1])
   for i in 2:ll  
      appo = appo*stateGammas[i]*stateLambdas[i]*dag(projGammas[i]*projLambdas[i])
   end
   appo = appo * stateGammas[ll+1]*dag(projGammas[ll+1])
   return scalar(appo)
end
