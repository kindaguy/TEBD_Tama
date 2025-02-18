using DelimitedFiles
using ITensors, ITensorMPS

#import just for debugging purposes 
include("printing_functs.jl")
include("VidalGauge.jl")

"""
    createGatesForTEBD2(sysEnergy, ChainLenght :: Int, tau , ttotal, freqs, coups)

Define the gates as in the trotter formula e^[(F+G)t]=e^[F*t/2]*e^[G*t]*e^[F*t/2]
F odd gate, G even gate
"""

function createGatesForTEBD2(sysEnergy, ntot, tau, freqs, coups)

    #system energy
    eps = sysEnergy

    tot_freqs = vcat([eps], freqs)
    tot_freqs[1] == eps && tot_freqs[2] == freqs[1]
    freqs = tot_freqs


    gates = ITensor[]

    #NOTE: Factor 2 in front of "Z" ops are due to Z=1 / 2 * σz
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

    return gates

end


"""
    TEBD2(v::VidalGauge, s::siteinds; sysEnergy, ChainLength::Int, tau, ttotal, cutoff=1E-14, minBondDim=5, maxBondDim=100, freqs, coups)

evolve a state in Vidal gauge within the TEBD framework

"""

function TEBD2(v::VidalGauge, s::siteinds; sysEnergy, ChainLength::Int, tau, ttotal, cutoff=1E-14, minBondDim=5, maxBondDim=100, freqs, coups)

    #total_size
    ntot = ChainLength + 1

    #println("ChainLength ", ntot)
    #println("length lambda ", length(v.lambdasV))

    #println(psi0)
    psi0 = convertToLeftCanonicalGauge(v)

    #println(psi0)

    #TODO: carry out this computation using directly the Vidal Form
    initmag = expect(psi0, "Z")
    println("first value of exp Z: ", initmag)

    gates = createGatesForTEBD2(sysEnergy, ntot, tau, freqs, coups)

    # Open files to save data 
    ioMagMeas = open("Data/magMeas_TEBD.dat", "w")
    ioTv = open("Data/tv_TEBD.dat", "w")
    ioNormCheck = open("Data/normCheck_TEBD.dat", "w")

    # If needed save them to be stored in arrays
    #magMeas = Vector{Float64}[]
    #tv = Vector{Float64}()
    #normCheck = Vector{Float64}()
    psiV = MPS[]
    vidalV = VidalGauge[]

    for t in 0.0:tau:ttotal
        #println("tempo: ", t)
        psi = convertToLeftCanonicalGauge(v)
        println(ioTv, t)
        writedlm(ioMagMeas, [expect(psi, "Z")], ',')
        println(ioNormCheck, norm(psi))
        #push!(tv, t)
        #push!(magMeas, expect(psi, "Z"))
        #push!(normCheck, norm(psi))
        push!(psiV, psi)
        push!(vidalV, v)

        v = apply_TEBD(gates, v; cutoff, mindim=5, maxBondDim=maxBondDim)

    end

    close(ioTv)
    close(ioMagMeas)
    close(ioNormCheck)

    return psiV, VidalV

end

function apply_TEBD(gates, v; cutoff, mindim=5, maxBondDim)

    v = apply_odd(gates, v; cutoff, mindim=5, maxBondDim=maxBondDim)

    #println("After first gate:")
    #println(v.lambdasV)
    #print_alternated_inds(Gammas, v.lambdasV)

    v = apply_even(gates, v; cutoff, mindim=5, maxBondDim=maxBondDim)

    v = apply_odd(gates, v; cutoff, mindim=5, maxBondDim=maxBondDim)

    #print_alternated_inds(Gammas, v.lambdasV)

    return v
end

# Application of ODD GATES
function apply_odd(gates, v; cutoff, mindim, maxBondDim)

    #println("Odd apply")

    LambdasInvEven = ITensor[]
    vOdd = VidalGauge()
    N = length(v.lambdasV)
    SVDLeftIndices = Index[]

    #print_inds(gates)

    # compute the inverse just for the right v.lambdasV (even), 
    # the inversion is well defined as in the previous SVD there is a cutoff that
    # cuts the smallest singular values decreasing the bond dimension
    # Q: perdo tanta informazione troncando?
    # Q: mindim non è compatibile, ma non capisco perchè
    for l in 1:length(v.lambdasV)
        if (l % 2 == 0)
            #LambdaInv = ITensor(inds(v.lambdasV[l]))
            LambdaInv = diag_itensor(inds(v.lambdasV[l]))
            LambdaInv .= inv.(v.lambdasV[l])
            push!(LambdasInvEven, LambdaInv)
        end
    end

    # println(LambdasInvEven)

    for i in 1:2:N

        empty!(SVDLeftIndices)

        if i == 1 #first site
            double_site_psi = v.gammasV[i] * v.lambdasV[i] * v.gammasV[i+1] * v.lambdasV[i+1]
        elseif i == N #last site, if it is an odd one, otherwise it will not be evolved 
            double_site_psi = v.lambdasV[i-1] * v.gammasV[i] * v.lambdasV[i] * v.gammasV[i+1]
        else
            double_site_psi = v.lambdasV[i-1] * v.gammasV[i] * v.lambdasV[i] * v.gammasV[i+1] * v.lambdasV[i+1]
        end

        #  println("double_site_psi ", double_site_psi)
        #  println("gates[$i] ", gates[i])

        Res = gates[i] * double_site_psi

        # println("Res ", Res)

        indices = inds(Res)

        #println(indices)

        #find the indices for the SVD
        SVDLeftIndices = find_inds(indices, i)

        #NOTA: se uso mindim=mindim si rompre tutto!!
        A, S, B = length(SVDLeftIndices) == 1 ?
                  svd(Res, SVDLeftIndices[1]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim) :
                  svd(Res, SVDLeftIndices[1], SVDLeftIndices[2]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim)

        # return to Vidal form
        A = i != 1 ? LambdasInvEven[Int((i - 1) / 2)] * A : A
        B = i != N ? B * LambdasInvEven[Int((i + 1) / 2)] : B

        if (i != 1)
            push!(vOdd.lambdasV, v.lambdasV[i-1]) #in between the gates, leave the same lambdas
        end
        push!(vOdd.lambdasV, S)
        push!(vOdd.gammasV, A)
        push!(vOdd.gammasV, B)

    end


    if (N % 2 == 0) #N is even, so the number of sites is odd
        push!(vOdd.lambdasV, v.lambdasV[N])
        push!(vOdd.gammasV, v.gammasV[N+1])
    end

    #print_alternated_inds(GammasOdd, LambdasOdd)

    # unprime the indices
    noprime!.(vOdd.lambdasV)
    noprime!.(vOdd.gammasV)

    return v

end


# Application of EVEN GATES
function apply_even(gates, v; cutoff, mindim, maxBondDim)

    #println("Even apply")

    LambdasInvOdd = ITensor[]
    vEven = VidalGauge()
    N = length(v.lambdasV)
    SVDLeftIndices = Index[]

    # compute the inverse just for the right v.lambdasV (odd), 
    # the inversion is well defined as in the previous SVD there is a cutoff that
    # cuts the smallest singular values decreasing the bond dimension 
    for l in 1:2:N
        #LambdaInv = ITensor(inds(v.lambdasV[l]))
        #LambdaInv .= inv.(v.lambdasV[l])
        LambdaInv = diag_itensor(inds(v.lambdasV[l]))
        LambdaInv .= (v.lambdasV[l]) .^ (-1)
        push!(LambdasInvOdd, LambdaInv)
    end

    #println(LambdasInvOdd)

    #there is no first site to evolve, so the gamma is the same
    push!(vEven.gammasV, v.gammasV[1])

    for i in 2:2:N

        #println(i)

        empty!(SVDLeftIndices)

        #there is no first site to evolve
        if i == N
            double_site_psi = v.lambdasV[i-1] * v.gammasV[i] * v.lambdasV[i] * v.gammasV[i+1]
        else
            double_site_psi = v.lambdasV[i-1] * v.gammasV[i] * v.lambdasV[i] * v.gammasV[i+1] * v.lambdasV[i+1]
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
            push!(vEven.lambdasV, v.lambdasV[i-1]) #in between the gates, leave the same lambdas
        end
        push!(vEven.lambdasV, S)
        push!(vEven.gammasV, A)
        push!(vEven.gammasV, B)

    end


    if (N % 2 != 0) #N is odd, so the number of sites is even
        push!(vEven.lambdasV, v.lambdasV[N])
        push!(vEven.gammasV, v.gammasV[N+1])
    end

    #println("Even apply")
    #print_alternated_inds(GammasEven, LambdasEven)

    # unprime the indices
    noprime!.(vEven.lambdasV)
    noprime!.(vEven.gammasV)

    return vEven

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