using ITensors, ITensorMPS


abstract type MPSGauge end

"""
    VidalGauge

A finite size matrix product state type in Vidal / Gamma Lambda gauge.

Can be represented as

g1 -l1 -g2 -l2 --
|       |

"""

mutable struct VidalGauge <: MPSGauge
        gammasV :: Vector{ITensor}
        lambdasV :: Vector{ITensor}
end


"""
   VidalGauge()

Construct an empty VidalGauge.

"""

VidalGauge() = VidalGauge(ITensor[],ITensor[])


"""
   convertToVidal!(v :: VidalGauge, psi :: MPS )

Save MPS in mixed canonical gauge in Vidal gauge 

"""

function convertToVidal!(v :: VidalGauge, psi:: MPS)

    #begin from a right canonical form
    orthogonalize!(psi, 1)

    N = length(psi)
    M = psi[1] #initialize the first M as the first site, as it is the orthogonality centre

    for i in 1:N-1

        #DEBUG: println("SVD$i")

        svd_inds=nothing
        
        svd_inds=uniqueinds(M,psi[i+1])
            
        A, Lambda, V = svd(M, svd_inds; lefttags="u=$i", righttags="v=$i")
        
        #DEBUG:
        #println("A: ", inds(A))
        #println("Lambda: ", inds(Lambda))
        #println("V: ", inds(V))
        #println(A * Lambda * dag(V) ≈ M)

        # push the Lambda matrix into the Lambdas vector  
        push!(v.lambdasV, Lambda)

        if i != 1 #do from the second site (NOTE: Lambdas will have N-1 elements) 
            # Assign the same indices (same ids) to the inverse of Lambdas, 
            # so they don't have just the same values but they can be contracted in the right way
            #LambdaInv = randomITensor(inds(Lambdas[i-1]))
            LambdaInv = diag_itensor(inds(v.lambdasV[i-1]))
            #println("debug ", LambdaInv)
            #println("debug2 ", Lambdas[i-1])
            #LambdaInv .= 0
            #println("debug3 pre assegnazione", LambdaInv)
            #LambdaInv .= Lambdas[i-1] .*2
            LambdaInv .= (v.lambdasV[i-1]) .^ (-1)
            #println("debug4 post assegnazione", LambdaInv)
            #CheckLambdaInv = prime(LambdaInv, commonind(Lambdas[i-1],LambdaInv))
            #println("debug check ", CheckLambdaInv)
            #ide = delta(noncommoninds(Lambdas[i-1],CheckLambdaInv))
            #println(" ide ", ide)
            #println(Lambdas[i-1]*CheckLambdaInv ≈ ide)
            # define the Gammas as per Vidal's form
            gamma = LambdaInv * A
        else #The first gamma is just the first A, 
            #note it is left normalized A A^dag =id (giusto?) 
            gamma = A
        end
        #push the gamma tensor into the Gammas vector
        push!(v.gammasV, gamma)

        if i != N - 1 #do unless it is the last site, 
            # in the last site there is no need to define M using the lambda
            M = v.lambdasV[i] * V * psi[i+1]
            # DEBUG: println(psi0)

        else #in the last site I leave out the Lambda, as it needs to stay out
            M = V * psi[i+1]
            # I don't need to make another SVD, so the M is the last gamma
            # Note that as V and psi[i+1] are both right orthogonal(RO), M is also RO
            push!(v.gammasV, M)
        end
    end
        
end


function convertToLeftCanonicalGauge(v :: VidalGauge)
    
    psiLeft=MPS(length(v.gammasV))
    psiLeft[1]=v.gammasV[1]

    for i in 1:length(v.lambdasV)
        psiLeft[i+1]=v.lambdasV[i]*v.gammasV[i+1]
    end
    
    return psiLeft

end


function convertToRightCanonicalGauge(v :: VidalGauge)
    
    psiRight=MPS(length(v.gammasV))

    for k in 1:length(v.lambdasV)
        psiRight[k]=v.gammasV[k]*v.lambdasV[k]
    end

    psiRight[end]=v.gammasV[end]
    
    return psiRight

end


function convertToMixedCanonicalGauge(v :: VidalGauge, i :: Int)
    
    N=length(v.gammasV)
    
    i == 1 && return psi = convertToRightCanonicalGauge(v)
    i == N && return psi = convertToLeftCanonicalGauge(v)
    
    psi=MPS(N)

    #left can until i-1
    psi[1]=v.gammasV[1]

    for k in 1:i-2
        psi[k+1]=v.lambdasV[k]*v.gammasV[k+1]
    end
    
    #oc in site i
    psi[i]=v.lambdasV[i-1]*v.gammasV[i]*v.lambdasV[i]
    
    #right can until
    for k in i+1:length(v.lambdasV)
        psi[k]=v.gammasV[k]*v.lambdasV[k]
    end

    psi[end]=v.gammasV[end]
    
    return psi

end


function checkOrthogonalityRight(v :: VidalGauge)
    
    psi_Right=convertToRightCanonicalGauge(v)
    
    for i in 2:length(psi_Right)
    
        M = dag(psi_Right[i])
        prime!(M; tags="v=$(i-1)")
        #println("M $(inds(M))")
    
        result = psi_Right[i]*M
        indices=inds(result)
        identity= delta(indices) 

        if (result ≈ identity) continue
        else 
            return false
            break
        end

        #La prima non va contratta, non ha indice v
        #println(psi_Left[i]*M)
        #println(psi_Left[i]*M - delta(indk,ind))
    end
    
    return true

end

#=
function convert_to_Vidal_until_OC(psi, OC)

    #begin from a right canonical form
    orthogonalize!(psi, 1)

    Gammas = ITensor[]
    Lambdas = ITensor[]
    N = length(psi)
    M = psi[1] #initialize the first M as the first site, as it is the orthogonality centre
    link_index = nothing

    for i in 1:OC

        #DEBUG: println("SVD$i")

        #find the site indices of the site
        site_index = siteind(psi, i)
        indices = inds(M)

        #DEBUG: 
        #println("Site$i indices:$site_index")
        #println("M = Gamma$(i-1) indices:$indices")
        #println(hastags(indices[1], "u"))

        # do from the second site (the first site does not have a left index) 
        if i != 1

            #find the left index coming from the previous Lambda  
            #to reshape it with the site index for the SVD
            for idx in indices
                if hastags(idx, "u=$(i-1)")

                    #DEBUG: println("in!")

                    link_index = idx
                    break
                end
            end

            # check if the "u" has been found
            if link_index === nothing
                error("The index with tag 'u=$(i-1)' associated to site $i has not been found")
            end

            #DEBUG: println(link_index, ",", site_index)

            #Do the SVD, renaming the left index linked to the Lamda as u=$i
            A, Lambda, V = svd(M, (link_index, site_index); lefttags="u=$i", righttags="v=$i")
        else
            #SVD for the first site
            A, Lambda, V = svd(M, site_index; lefttags="u=$i")
        end

        #DEBUG:
        #println("A: ", inds(A))
        #println("Lambda: ", inds(Lambda))
        #println("V: ", inds(V))
        #println(A * Lambda * dag(V) ≈ M)

        # push the Lambda matrix into the Lambdas vector  
        push!(Lambdas, Lambda)

        if i != 1 #do from the second site (NOTE: Lambdas will have N-1 elements) 
            # Assign the same indices (same ids) to the inverse of Lambdas, 
            # so they don't have just the same values but they can be contracted in the right way
            LambdaInv = ITensor(inds(Lambdas[i-1]))
            LambdaInv .= (Lambdas[i-1]) .^ (-1)
            # define the Gammas as per Vidal's form
            gamma = A * LambdaInv
        else #The first gamma is just the first A, 
            #note it is left normalized A A^dag =id (giusto?) 
            gamma = A
        end
        #push the gamma tensor into the Gammas vector
        push!(Gammas, gamma)

        if i != OC #do unless it is the last site, 
            # in the last site there is no need to define M using the lambda
            M = Lambdas[i] * V * psi0[i+1]
            # DEBUG: println(psi0)

        else #in the last site I leave out the Lambda, as it needs to stay out
            M = V * psi0[i+1]
            # I don't need to make another SVD, so the M is the last gamma
            # Note that as V and psi[i+1] are both right orthogonal(RO), M is also RO
            #push!(Gammas, M) in questo caso NO
        end
    end

    return Lambdas, Gammas, M
end=#