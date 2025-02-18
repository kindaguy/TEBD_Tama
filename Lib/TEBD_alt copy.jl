using DelimitedFiles
using ITensors

#import just for debugging purposes 
include("printing_functs.jl")

function MPS_evolution_TEBD(Gammas, Lambdas, s; sysEnergy, ChainLength, tau, ttotal, cutoff=1E-14, minBondDim=5, maxBondDim=100, freqs, coups)

    #system energy
    eps = sysEnergy

    tot_freqs = vcat([eps], freqs)
    tot_freqs[1] == eps && tot_freqs[2] == freqs[1]
    freqs = tot_freqs

    #total_size
    ntot = ChainLength + 1

    println("ChainLength ", ntot)
    println("length lambda ", length(Lambdas))

    #println(psi0)
    psi0 = convert_to_MPS(Gammas, Lambdas, ntot)

    #println(psi0)

    #TO DO: carry out this computation using directly the Vidal Form
    initmag = expect(psi0, "Z")
    println("first value of exp Z: ", initmag)


    # Define the gates as in the trotter formula e^[(F+G)t]=e^[F*t/2]*e^[G*t]*e^[F*t/2]
    # F odd gate, G even gate

    gates = ITensor[]

    #NOTE: Factor 2 in front of "Sz" ops are due to Sz=1 / 2 * σz
    h_start = 0.5 * freqs[1] * op("Z", s[1]) * op("Id", s[2]) +
              0.5 * 0.5 * freqs[2] * op("Z", s[2]) * op("Id", s[1]) +
              coups[1] * op("S+", s[1]) * op("S-", s[2]) +
              coups[1] * op("S+", s[2]) * op("S-", s[1])
    push!(gates, exp(-im * tau / 2 * h_start)) #first one is always divided by 2

    for j in 2:(ntot-2)

        t = isodd(j) ? tau / 2 : tau


        s1 = s[j]
        s2 = s[j+1]

        hj = coups[j] * op("S+", s1) * op("S-", s2) +
             coups[j] * op("S-", s1) * op("S+", s2) +
             0.5 * 0.5 * freqs[j] * op("Z", s1) * op("Id", s2) +
             0.5 * 0.5 * freqs[j+1] * op("Z", s2) * op("Id", s1)

        Gj = exp(-im * t * hj)
        push!(gates, Gj)
    end

    t = isodd(ntot) ? tau / 2 : tau

    h_end = 0.5 * freqs[ntot] * op("Z", s[ntot]) * op("Id", s[ntot-1]) +
            0.5 * 0.5 * freqs[ntot-1] * op("Z", s[ntot-1]) * op("Id", s[ntot]) +
            coups[ntot-1] * op("S+", s[ntot]) * op("S-", s[ntot-1]) +
            coups[ntot-1] * op("S+", s[ntot-1]) * op("S-", s[ntot])

    push!(gates, exp(-im * t * h_end))

    # Open files to save data 
    ioMagMeas = open("Data/magMeas_TEBD.dat", "w")
    ioTv = open("Data/tv_TEBD.dat", "w")
    ioNormCheck = open("Data/normCheck_TEBD.dat", "w")

    # If needed save them to be stored in arrays
    #magMeas = Vector{Float64}[]
    #tv = Vector{Float64}()
    #normCheck = Vector{Float64}()
    psiV = MPS[]
    LambdasV = Vector{Vector{ITensor}}()
    GammasV = Vector{Vector{ITensor}}()


    for t in 0.0:tau:tmax
        #println("tempo: ", t)
        psi = convert_to_MPS(Gammas, Lambdas, ntot)
        println(ioTv, t)
        writedlm(ioMagMeas, [expect(psi, "Z")], ',')
        println(ioNormCheck, norm(psi))
        #push!(tv, t)
        #push!(magMeas, expect(psi, "Z"))
        #push!(normCheck, norm(psi))
        push!(psiV, psi)
        push!(GammasV, Gammas)
        push!(LambdasV, Lambdas)

        Gammas, Lambdas = apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim=maxBondDim)

    end

    close(ioTv)
    close(ioMagMeas)
    close(ioNormCheck)

    return psiV, GammasV, LambdasV

end

function SpinBoson_evolution_TEBD(Gammas, Lambdas, s; ϵ,Δ, sysenvInt::String, ChainLength, tau, ttotal, cutoff=1E-14, minBondDim=5, maxBondDim=100, freqs, coups)

    #system energy
    #eps = sysEnergy

    # tot_freqs = vcat([eps], freqs)
    # tot_freqs[1] == eps && tot_freqs[2] == freqs[1]
    # freqs = tot_freqs

    #total_size
    ntot = ChainLength + 1

    println("ChainLength ", ntot)
    println("length lambda ", length(Lambdas))

    #println(psi0)
    #psi0 = convert_to_MPS(Gammas, Lambdas, ntot)

    #println(psi0)

    #TO DO: carry out this computation using directly the Vidal Form
    #initmag = expect(psi0, "Z")
    initmag = noprime!(op("Z",s[1])*Gammas[1])*Lambdas[1]*Lambdas[1]*conj(Gammas[1])
    println("first value of exp Z: ", scalar(initmag))


    # # Define the gates as in the trotter formula e^[(F+G)t]=e^[F*t/2]*e^[G*t]*e^[F*t/2]
    # # F odd gate, G even gate

    gates = ITensor[]

    #NOTE: generalize as to include (1+σ_z)/2 required, for example, in dimer sims
    h_start = 0.5 * ϵ * op("Z", s[1]) * op("Id", s[2]) +
              0.5 * Δ * op("X",s[1]) * op("Id",s[2]) +  
              0.5 * 0.5 * freqs[1] * op("N", s[2]) * op("Id", s[1]) +
              coups[1] * op(sysenvInt, s[1]) * op("A", s[2]) +
              coups[1] * op("Adag", s[2]) * op(sysenvInt, s[1])
    push!(gates, exp(-im * tau / 2 * h_start)) #first one is always divided by 2

        #Debug
        #println("hstart: ",h_start)
    for j in 2:(ntot-2)

        t = isodd(j) ? tau / 2 : tau


        s1 = s[j]
        s2 = s[j+1]

        hj = coups[j] * op("Adag", s1) * op("A", s2) +
             coups[j] * op("A", s1) * op("Adag", s2) +
             0.5 * freqs[j] * op("N", s1) * op("Id", s2) +
             0.5 * freqs[j+1] * op("N", s2) * op("Id", s1)
        
        # #Debug
        # if j==2
        #     println("h2:",hj)
        # end

        Gj = exp(-im * t * hj)
        push!(gates, Gj)
    end

    t = isodd(ntot) ? tau / 2 : tau

    h_end = freqs[ntot-1] * op("N", s[ntot]) * op("Id", s[ntot-1]) +
            0.5  * freqs[ntot-2] * op("N", s[ntot-1]) * op("Id", s[ntot]) +
            coups[ntot-1] * op("Adag", s[ntot]) * op("A", s[ntot-1]) +
            coups[ntot-1] * op("Adag", s[ntot-1]) * op("A", s[ntot])

    push!(gates, exp(-im * t * h_end))

    # Open files to save data 
    ioMagMeas = open("Data/magMeas_TEBD.dat", "w")
    ioTv = open("Data/tv_TEBD.dat", "w")
    ioNormCheck = open("Data/normCheck_TEBD.dat", "w")

    # If needed save them to be stored in arrays
    magMeas = Vector{ComplexF64}()
    tv = Vector{Float64}()
    normCheck = Vector{ComplexF64}()
    #psiV = MPS[]
    #LambdasV = Vector{Vector{ITensor}}()
    #GammasV = Vector{Vector{ITensor}}()


    for t in 0.0:tau:tmax
        println("tempo: ", t)
        #psi = convert_to_MPS(Gammas, Lambdas, ntot)
        println(ioTv, t)
        push!(tv, t)
        
        
        #Magnetization measure
        appo =  noprime!(op("Z",s[1])*Gammas[1])*Lambdas[1]*Lambdas[1]*conj(Gammas[1])
        writedlm(ioMagMeas, [scalar(appo)], ',')
        push!(magMeas, scalar(appo))

        #Norm measure
        appoNorm = Gammas[1]*Lambdas[1]*Lambdas[1]*conj(Gammas[1])
        println(ioNormCheck, scalar(appoNorm))
        push!(normCheck, scalar(appoNorm))
        #push!(psiV, psi)
        #push!(GammasV, Gammas)
        #push!(LambdasV, Lambdas)

        Gammas, Lambdas = apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim=maxBondDim)

    end

    close(ioTv)
    close(ioMagMeas)
    close(ioNormCheck)

    return Gammas, Lambdas

end

function apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim)

    #println("first odd")
    Gammas, Lambdas = apply_odd(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim=maxBondDim)
    #println("Gammas[1]",Gammas[1])
    #println("Gammas[2]",Gammas[2])
    #println("After first gate:")
    #println(Lambdas)
    #print_alternated_inds(Gammas, Lambdas)
    #pippo = readline()
    #println("the even step")
    Gammas, Lambdas = apply_even(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim=maxBondDim)

    #println("second odd")
    #println("Gammas[1]: ", Gammas[1])
    Gammas, Lambdas = apply_odd(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim=maxBondDim)    
    #print_alternated_inds(Gammas, Lambdas) 
    return Gammas, Lambdas
end

# Application of ODD GATES
function apply_odd(gates, Gammas, Lambdas; cutoff, mindim, maxBondDim)

    #Number of links:
    #we are updating on links, so this is the relevant Number
    N = length(Lambdas)
    
    #Determine number of threads
    nt = Threads.nthreads()
    
    #For all variables instantiate a number of copies
    #equal to the number of threads

    #List of variables to pool

    #There are things that are already vectors.
    #The "problem" is that:
    #- sequential approach: we push one element at a time.
    #-multi-threading: we need to record things in the right place. 
    #This place is typically determined by the loop index
    
    #Try to define a single code that can be executed both
    #sequentially than in multi-threading

        LambdasInvEven = ITensor[]
        LambdasOdd = ITensor[]
        GammasOdd = ITensor[]
        SVDLeftIndices = Index[]
    
    #print_inds(gates)

    # compute the inverse just for the right Lambdas (even), 
    # the inversion is well defined as in the previous SVD there is a cutoff that
    # cuts the smallest singular values decreasing the bond dimension
    # Q: perdo tanta informazione troncando?
    # Q: mindim non è compatibile, ma non capisco perchè
    for l in 1:length(Lambdas)
        if (l % 2 == 0)
            #LambdaInv = ITensor(inds(Lambdas[l]))
            LambdaInv = diag_itensor(inds(Lambdas[l]))
            LambdaInv .= inv.(Lambdas[l])
            push!(LambdasInvEven, LambdaInv)
        end
    end

    # println(LambdasInvEven)

Threads.@threads  for i in 1:2:N
    println(Threads.threadid())
        
        #println("odd", i)

        empty!(SVDLeftIndices)
    

        if i == 1 #first site
            #println("first site")
            #Debug
            # @show Gammas[1]
            # @show Lambdas[1]
            # @show Gammas[2]
            # @show Lambdas[2]
            double_site_psi = Gammas[i] * Lambdas[i] * Gammas[i+1] * Lambdas[i+1]
            #println("first site completed")
        elseif i == N #last site, if it is an odd one, otherwise it will not be evolved 
            #println("Last site")
            double_site_psi = Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1]
        else
            #println("site",i)
            double_site_psi = Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1] * Lambdas[i+1]
        end

        #  println("double_site_psi ", double_site_psi)
        #  println("gates[$i] ", gates[i])

        Res = gates[i] * double_site_psi
        
        #Debug
        # if i==1
        #     @show Res
        # end

        # println("Res ", Res)

        indices = inds(Res)

        #println(indices)

        #find the indices for the SVD
        SVDLeftIndices = find_inds(indices, i)
        # if i==1
            # println("first_site SVDLeftIdx",SVDLeftIndices)
        # end
# 

        #NOTA: se uso mindim=mindim si rompre tutto!!
        A, S, B = length(SVDLeftIndices) == 1 ?
                  svd(Res, SVDLeftIndices[1]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim) :
                  svd(Res, SVDLeftIndices[1], SVDLeftIndices[2]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim)

        # return to Vidal form
        A = i != 1 ? LambdasInvEven[Int((i - 1) / 2)] * A : A
        B = i != N ? B * LambdasInvEven[Int((i + 1) / 2)] : B

        if (i != 1)
            push!(LambdasOdd, Lambdas[i-1]) #in between the gates, leave the same lambdas
        end
        push!(LambdasOdd, S)
        push!(GammasOdd, A)
        push!(GammasOdd, B)

    end


    if (N % 2 == 0) #N is even, so the number of sites is odd
        push!(LambdasOdd, Lambdas[N])
        push!(GammasOdd, Gammas[N+1])
    end

    #print_alternated_inds(GammasOdd, LambdasOdd)

    # unprime the indices, maybe it can be done in the for loop
    for i in eachindex(Lambdas)
        noprime!(LambdasOdd[i])
        noprime!(GammasOdd[i])
    end

    noprime!(GammasOdd[end])

    return GammasOdd, LambdasOdd

end


# Application of EVEN GATES
function apply_even(gates, Gammas, Lambdas; cutoff, mindim, maxBondDim)

    #println("Even apply")

    LambdasInvOdd = ITensor[]
    LambdasEven = ITensor[]
    GammasEven = ITensor[]
    N = length(Lambdas)
    SVDLeftIndices = Index[]

    # compute the inverse just for the right Lambdas (odd), 
    # the inversion is well defined as in the previous SVD there is a cutoff that
    # cuts the smallest singular values decreasing the bond dimension 
    for l in 1:2:length(Lambdas)
        #LambdaInv = ITensor(inds(Lambdas[l]))
        #LambdaInv .= inv.(Lambdas[l])
        LambdaInv = diag_itensor(inds(Lambdas[l]))
        LambdaInv .= (Lambdas[l]) .^ (-1)
        push!(LambdasInvOdd, LambdaInv)
    end

    #println(LambdasInvOdd)

    #there is no first site to evolve, so the gamma is the same
    push!(GammasEven, Gammas[1])

Threads.@threads for i in 2:2:N

        #println("Even ",i)

        empty!(SVDLeftIndices)

        #there is no first site to evolve
        if i == N
            double_site_psi = Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1]
        else
            double_site_psi = Lambdas[i-1] * Gammas[i] * Lambdas[i] * Gammas[i+1] * Lambdas[i+1]
        end

        Res = gates[i] * double_site_psi

        indices = inds(Res)

        #println(indices)

        #find the indices for the SVD
        SVDLeftIndices = find_inds(indices, i)

        #println(SVDLeftIndices)

        A, S, B = length(SVDLeftIndices) == 1 ?
                  svd(Res, SVDLeftIndices[1]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim) :
                  svd(Res, SVDLeftIndices[1], SVDLeftIndices[2]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim)

        #println(inds(A))
        #println(inds(S))
        #println(inds(B))

        # return to Vidal form
        A = i != 1 ? LambdasInvOdd[Int(i / 2)] * A : A
        B = i != N ? B * LambdasInvOdd[Int((i + 2) / 2)] : B

        #println("gamma", i, ": ", inds(A))
        #println("gamma", (i + 1), ": ", inds(B))

        if (i != 1)
            push!(LambdasEven, Lambdas[i-1]) #in between the gates, leave the same lambdas
        end
        push!(LambdasEven, S)
        push!(GammasEven, A)
        push!(GammasEven, B)

    end


    if (N % 2 != 0) #N is odd, so the number of sites is even
        push!(LambdasEven, Lambdas[N])
        push!(GammasEven, Gammas[N+1])
    end

    #println("Even apply")
    #print_alternated_inds(GammasEven, LambdasEven)

    # unprime the indices
    for i in eachindex(Lambdas)
        noprime!(LambdasEven[i])
        noprime!(GammasEven[i])
    end

    noprime!(GammasEven[end])

    return GammasEven, LambdasEven

end

# Function to find the right indices (left bond one and site one) for the SVD
function find_inds(indices, i)
    SVDLeftIndices = Index[]
    for idx in indices
         if hastags(idx, "u=$(i-1)")

            push!(SVDLeftIndices, idx)
        end
        if hastags(idx, "n=$i")

            push!(SVDLeftIndices, idx)
        end
    end
    return SVDLeftIndices
end

# Defined to compute the expectation value using ITensors methods 
# left-canonical form
function convert_to_MPS(Gammas, Lambdas, ntot)

    psi = MPS(ntot)
    #println(length(psi))
    psi[1] = Gammas[1] #left normalized (A^dag A = id)

    for i in eachindex(Lambdas)
        psi[i+1] = Lambdas[i] * Gammas[i+1]
    end
    #println(psi)

    return psi
end