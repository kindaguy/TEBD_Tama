
# We define generic functions to evaluate <ψ| A_i ⊗B_j |ψ>
#with i,j arbitrary sites of the overall system.

"""
   applySingleSiteOp(inds,refGammas,site,what)

   Returns the Vidal's Gammas matrices 
for the site-th element (of the chain) modified by what. Note that Lambda matrices are remain unmodified.
"""
function applySingleSiteOp(inds,refGammas, site::Int64, what::String)
   appo = [a for a in refGammas]
   appo[site] = noprime!(op(what,inds[site])*refGammas[site])
   return appo
end

"""
   applyTwoSiteOp(inds,refGammas,site1,what1,site2,what2)

Returns the Vidal's Gammas matrices for the site-th element (of the chain)
modified by what1/what2 at sites site1, site2. 
Note that Lambda matrices are remain unmodified.
"""
function applyTwoSiteOp(inds,refGammas, site1::Int64, what1::String,site2::Int64, what2::String)
   appo = [a for a in refGammas]
   appo[site1] = noprime!(op(what1,inds[site1])*refGammas[site1])
   appo[site2] = noprime!(op(what2,inds[site2])*appo[site2])
   return appo
end


"""
projectBetween(LambdasPsi,GammasPsi,LambdasPhi,GammasPhi,site1, site2)

performs the projection <ϕ|ψ>, assuming that ϕ and ψ differ only at the sites
site1 and site2
"""
#Note: this function remains always the same: does not matter what's ϕ
function projectBetween(stateLambdas,stateGammas,projLambdas,projGammas,site1,site2)
   
   if site1<=site2
      pL = site1
      pR = site2
   else
      pL = site2
      pR = site1
   end

   #Quantifies the distance between the two sites
   #ll = pR-pL 
   
   #Operators both on first site
   if pL == pR == 1
      #This is a rank zero object
       appo = stateGammas[pL] * stateLambdas[pL] *  dag(projGammas[pL]*projLambdas[pL])
       return scalar(appo)
   #else if last site
   elseif pL == pR == size(stateGammas,1) 
      println("hit last site:",pL," ", pR)
      appo = stateLambdas[pL-1] * stateGammas[pR]   *  dag(projLambdas[pL-1] * projGammas[pR])
      return scalar(appo)
   elseif pL == pR && pL>1 && pR<size(stateGammas,1)
      #both operators on same site
      appo1 = stateLambdas[pR-1]*stateGammas[pR]* stateLambdas[pR]
      appo2 = dag(projLambdas[pR-1]*projGammas[pR]* projLambdas[pR])
      appo = appo1 * appo2
      return scalar(appo)
   else
      #pL and pR have right and left bonds and at least one bond in between
      #the bonds pL and pR must not be primed; the remaining ones does
      
      #This needed as to have tounprime visible outside for loop
      tounprimeR = commonind(projGammas[pL],projLambdas[pL])
      if pL==1
         prevState =  stateGammas[pL]*dag(prime(projGammas[pL],tounprimeR))    
      else
         prevState = stateLambdas[pL-1]* stateGammas[pL]*dag(projLambdas[pL-1]*prime(projGammas[pL],tounprimeR)) 
      end

      if pR < size(stateGammas,1)
         for i in pL+1:pR
            tounprimeL = commonind(projLambdas[i-1],prevState)
            tounprimeR = commonind(projGammas[i],projLambdas[i])
            prevState = prevState * stateLambdas[i-1]* stateGammas[i] * prime(dag(projLambdas[i-1]*projGammas[i]),[tounprimeL,tounprimeR])
         end
         appo = prevState * stateLambdas[pR] * prime(projLambdas[pR],tounprimeR)
         return scalar(appo)
      else
          for i in pL+1:pR-1
            tounprimeL = commonind(projLambdas[i-1],prevState)
            tounprimeR = commonind(projGammas[i],projLambdas[i])
            prevState = prevState * stateLambdas[i-1]* stateGammas[i] * prime(dag(projLambdas[i-1]*projGammas[i]),[tounprimeL,tounprimeR])
         end
         appo = prevState * stateLambdas[pR-1]*stateGammas[pR] * dag(prime(projLambdas[pR-1],tounprimeR)*projGammas[pR])
         return scalar(appo)   
      end
      #Gamma[pR] exits primed
      
      appo = prevState * stateLambdas[pR] * prime(projLambdas[pR],tounprimeR)
      return scalar(appo)
   end
   
end

function indDoubleMeas(chain_length)
   inds = Dict{Int64,Tuple{Int64,Int64}}()
   veci = Vector{Int64}([])
   for i in 1:chain_length
      for j in i:chain_length
         push!(veci,1)
         a=length(veci)
         push!(inds,a=>(i+1,j+1))
      end
   end
   return inds
end