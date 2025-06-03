
#Prepare gates for:
#-spinBoson: general spin-boson setting: interaction term is X on bath side
#-singleEx: number conserving interaction Hamiltonian

function spinBosonGates(ϵ::Float64,
Δ::Float64,
sysenvInt::String,
freqFile::String,
coupFile::String,
# Gammas::Array{ITensor},
# Lambdas::Array{ITensor},
s::Array{Index{Int64}},
tau::Float64;
isSingleExc::Bool=false
)::Union{Vector{ITensor},Nothing}

    ntot = length(s)
    ChainLength =ntot-1

    freqs = readdlm(freqFile,'\n' ,Float64); #chain fequencies
    coups = readdlm(coupFile,'\n',Float64);  #chain couplings; first element:sys-bath coupling
    println("Available freqs on file: ",length(freqs))
    println("Available coups on file: ", length(coups))

    if length(freqs)< ChainLength
    println("ERROR: Not enough frequency coefficients: required", ChainLength," available: ",length(freqs))
    return
    end

    if length(coups)< ChainLength-1
        println("ERROR: Not enough coupling coefficients: required", ChainLength-1," available: ",length(coups))
        return
    end




    
    gates = ITensor[]

    #Switch on isSingleExc:
    if isSingleExc
        h_start = 
                #Local energy system: ϵ(1+σz)/2 
                0.5 * ϵ * op("Id",s[1]) * op("Id",s[2])+
                0.5 * ϵ * op("Z",s[1]) * op("Id",s[2])+
                0.5 * freqs[1] * op("N", s[2]) * op("Id", s[1]) +
                #Here exchange interaction
                #(\sigmaMinus * Adag + sigmaPlus * A)  
                coups[1] * op("Splus", s[1]) * op("A", s[2]) +
                coups[1] * op("Adag", s[2]) * op("Sminus", s[1])
        push!(gates, exp(-im * tau / 2 * h_start)) #first one is always divided by 2
    else 
        h_start = 0.5 * ϵ * op("Z", s[1]) * op("Id", s[2]) +
        0.5 * Δ * op("X",s[1]) * op("Id",s[2]) +  
        0.5 * freqs[1] * op("N", s[2]) * op("Id", s[1]) +
        coups[1] * op(sysenvInt, s[1]) * op("A", s[2]) +
        coups[1] * op("Adag", s[2]) * op(sysenvInt, s[1])
        push!(gates, exp(-im * tau / 2 * h_start)) #first one is always divided by 2

    end

    #All the remaining gates are the same
        
    for j in 2:(ntot-2)
        t = isodd(j) ? tau / 2 : tau
        s1 = s[j]
        s2 = s[j+1]
        hj = coups[j] * op("Adag", s1) * op("A", s2) +
            coups[j] * op("A", s1) * op("Adag", s2) +
            0.5 * freqs[j-1] * op("N", s1) * op("Id", s2) +
            0.5 * freqs[j] * op("N", s2) * op("Id", s1)
        
        Gj = exp(-im * t * hj)
        push!(gates, Gj)
    end
    t = isodd(ntot) ? tau / 2 : tau
    h_end = freqs[ntot-1] * op("N", s[ntot]) * op("Id", s[ntot-1]) +
            0.5  * freqs[ntot-2] * op("N", s[ntot-1]) * op("Id", s[ntot]) +
            coups[ntot-1] * op("Adag", s[ntot]) * op("A", s[ntot-1]) +
            coups[ntot-1] * op("Adag", s[ntot-1]) * op("A", s[ntot])
    push!(gates, exp(-im * t * h_end))

    return gates

end