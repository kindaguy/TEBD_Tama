using DelimitedFiles
using ITensors

#import just for debugging purposes 
include("printing_functs.jl")
include("projectiveMeas.jl")
include("twoSiteMeas.jl")



function SpinBoson_evolution_TEBD(Gammas, Lambdas, s; 
    ϵ,
    Δ, 
    sysenvInt::String, 
    ChainLength, 
    tau, 
    ttotal,
    measStep,
    cutoff=1E-14, 
    minBondDim=5, 
    maxBondDim=100, 
    freqs, 
    coups,
    performProj=false,
    performProj2=false,
    performOcc = false,
    trackBondDim = false)

    
    #total_size
    ntot = ChainLength + 1

    println("ChainLength ", ntot)
    println("length lambda ", length(Lambdas))

   
    initmag = noprime!(op("Z",s[1])*Gammas[1])*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
    println("first value of exp Z: ", scalar(initmag))


    # # Define the gates as in the trotter formula e^[(F+G)t]=e^[F*t/2]*e^[G*t]*e^[F*t/2]
    # # F odd gate, G even gate

    gates = ITensor[]

   
    h_start = 0.5 * ϵ * op("Z", s[1]) * op("Id", s[2]) +
              0.5 * Δ * op("X",s[1]) * op("Id",s[2]) +  
              0.5 * freqs[1] * op("N", s[2]) * op("Id", s[1]) +
              coups[1] * op(sysenvInt, s[1]) * op("A", s[2]) +
              coups[1] * op("Adag", s[2]) * op(sysenvInt, s[1])
    push!(gates, exp(-im * tau / 2 * h_start)) #first one is always divided by 2

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

    # Open files to save data 
    ioMagMeas = open("Data/magMeas_TEBD.dat", "w")
    ioTv = open("Data/tv_TEBD.dat", "w")
    ioNormCheck = open("Data/normCheck_TEBD.dat", "w")
    ioPopMeas = open("Data/popMeas_TEBD.dat","w")
    
    #Projective measurement
    if performProj
        ioProjMeasUp=open("Data/projMeasUp_TEBD.dat","w")
        ioProjMeasDn=open("Data/projMeasDn_TEBD.dat","w")
    
    #Create basis for zero and single excitation subspace
    zeroLambdasUp,zeroGammasUp = envZeroBasisUp(s,ChainLength);
    zeroLambdasDn,zeroGammasDn = envZeroBasisDn(s,ChainLength);
    println("Creating single ex basis...")
    singleExUp =  createSingleExBasis(s,zeroGammasUp,ChainLength);
    singleExDn =  createSingleExBasis(s,zeroGammasDn,ChainLength);
    println("...done")

    end

    #ATTENTION: two-site (projections)measurements can be performed
    #only if single-site projections are enabled
    if performProj2 && performProj
        ioProjMeas2Up=open("Data/projMeas2Up_TEBD.dat","w")
        ioProjMeas2Dn=open("Data/projMeas2Dn_TEBD.dat","w")

        println("Creating doulble excitation basis...")
        twoExUp = createDoubleExBasis(s,zeroGammasUp,ChainLength)
        twoExDn = createDoubleExBasis(s,zeroGammasDn,ChainLength)
        println("...done!")
    end
     
    if performOcc
        ioOcc=open("Data/toOccData.dat","w")
    end

    if trackBondDim
        ioTrackBondDim=open("Data/bondDims.dat","w")
    end
    # If needed save them to be stored in arrays
    # magMeas = Vector{ComplexF64}()
    # tv = Vector{Float64}()
    # normCheck = Vector{ComplexF64}()
    # popMeas = Vector{Vector{ComplexF64}}

    #psiV = MPS[]
    #LambdasV = Vector{Vector{ITensor}}()
    #GammasV = Vector{Vector{ITensor}}()


    
    for (step,t) in enumerate(0.0:tau:ttotal)
        println("tempo: ", t)
        
        #psi = convert_to_MPS(Gammas, Lambdas, ntot)
        if (step-1) % measStep == 0
            println(ioTv, t)
            flush(ioTv)
            #push!(tv, t)
            
            
            #Magnetization measure
            spinMeasures=Vector{ComplexF64}([])
            appo =  noprime!(op("Z",s[1])*Gammas[1])*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
            push!(spinMeasures,scalar(appo))
            appo =  noprime!(op("X",s[1])*Gammas[1])*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
            push!(spinMeasures,scalar(appo))
            appo =  noprime!(op("Y",s[1])*Gammas[1])*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
            push!(spinMeasures,scalar(appo))

            writedlm(ioMagMeas, transpose(vcat(t, spinMeasures)), ',') 
            flush(ioMagMeas)
            #push!(magMeas, scalar(appo))

            #Chain occupation measure
            #this is delicate: we need to measure all the chain sites
            #This will be painful, but needed.
            #Here I implement full measure, but I suspect that it 
            #can be otimized somehow:
            #1) Far sites are perturbed at later Times
            #2) This can be parallelized over threads: instantiate whole array
            occMeasures = Vector{ComplexF64}([])
            for i in 2:ntot-1
                appoOcc = Lambdas[i-1]*noprime!(op("N",s[i])*Gammas[i])*Lambdas[i]*dag(Lambdas[i-1]*Gammas[i]*Lambdas[i])
                push!(occMeasures, scalar(appoOcc))
            end
            appoOcc = Lambdas[ntot-1]*noprime!(op("N",s[ntot])*Gammas[ntot])*dag(Lambdas[ntot-1]*Gammas[ntot])
            push!(occMeasures, scalar(appoOcc))
            writedlm(ioPopMeas, transpose(vcat(t, occMeasures)), ',')
            flush(ioPopMeas)

            #Norm measure
            appoNorm = Gammas[1]*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
            writedlm(ioNormCheck, [t scalar(appoNorm)],',')
            flush(ioNormCheck)
            #push!(normCheck, scalar(appoNorm))

            if performProj
                projUp = zeros(ComplexF64,ChainLength+1)
                projDn = zeros(ComplexF64,ChainLength+1)
                
                #Projection on zero exc subspace
                projUp[1] = project(Lambdas,Gammas,zeroLambdasUp,zeroGammasUp)
                projDn[1] = project(Lambdas,Gammas,zeroLambdasDn,zeroGammasDn)

                #Projection on single exc subspace
                Threads.@threads for i in 2:ChainLength+1
                    projUp[i] = project(Lambdas,Gammas,zeroLambdasUp,singleExUp[i-1])
                end

                Threads.@threads for i in 2:ChainLength+1
                    projDn[i] = project(Lambdas,Gammas,zeroLambdasDn,singleExDn[i-1])
                end

                writedlm(ioProjMeasUp, transpose(vcat(t, projUp)), ',')
                flush(ioProjMeasUp)
                writedlm(ioProjMeasDn, transpose(vcat(t, projDn)), ',')
                flush(ioProjMeasDn)
            end

            if performProj2 && performProj

                numOfProjs = size(twoExUp,1)
                projUp = zeros(ComplexF64,numOfProjs)
                projDn = zeros(ComplexF64,numOfProjs)
                
                #Projection on double exc subspace
                #we apply lightconing: we know that excitation have a 
                #propagation speed which is < 2 * k∞. This defines a light-cone
                #It is thus useless to measure the projection if
                #either of the sites are outside the light cone
                lcone = (step-1)* tau * 2 * coups[end]
                Threads.@threads for i=1:numOfProjs
                    p1 = div(i-1,ChainLength)   #i = p1+1
                    p2 = (i-1)%ChainLength      #j = p2+1
                    #...but we are conservative so...
                    if p1 <= lcone +2 && p2 <= lcone + 2
                        projUp[i] = project(Lambdas,Gammas,zeroLambdasUp,twoExUp[i])
                        projDn[i] = project(Lambdas,Gammas,zeroLambdasDn,twoExDn[i])
                    end
                end

                # Threads.@threads for i=1:numOfProjs
                #     projDn[i] = project(Lambdas,Gammas,zeroLambdasDn,twoExDn[i])
                # end

                writedlm(ioProjMeas2Up, transpose(vcat(t, projUp)), ',')
                flush(ioProjMeasUp)
                writedlm(ioProjMeas2Dn, transpose(vcat(t, projDn)), ',')
                flush(ioProjMeasDn)
            end
            

            if performOcc
                indTSM = indDoubleMeas(ChainLength)
                numTSM = length(indTSM)
                TSM = zeros(ComplexF64,numTSM)
                
                
                #Projection on double exc subspace
                #we apply lightconing: we know that excitation have a 
                #propagation speed which is < 2 * k∞. This defines a light-cone
                #It is thus useless to measure the projection if
                #either of the sites are outside the light cone
                lcone = (step-1)* tau * 2 * coups[end]
                Threads.@threads for i=1:numTSM
                    p1 = indTSM[i][1] 
                    p2 = indTSM[i][2]
                    #...but we are conservative so...
                    if  p2 <= max(20,lcone + 2)
                        #Here we create the projecting state
                        projGammas = applyTwoSiteOp(sysenv,Gammas, p2, "A",p1, "Adag")
                        projLambdas = Lambdas;
                        TSM[i] = projectBetween(Lambdas,Gammas,projLambdas,projGammas,p2,p1)
                    end
                end


                writedlm(ioOcc, transpose(vcat(t, TSM)), ',')
                flush(ioOcc)
            end



            if trackBondDim

                bondDim = [dim(inds(a)[1]) for a in Lambdas]
                writedlm(ioTrackBondDim, transpose(vcat(t,bondDim)), ',')
                flush(ioTrackBondDim)
            
            end

        
        
        end
        
        
        #push!(psiV, psi)
        #push!(GammasV, Gammas)
        #push!(LambdasV, Lambdas)

        Gammas, Lambdas = apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim=maxBondDim)

    end

    close(ioTv)
    close(ioMagMeas)
    close(ioNormCheck)
    close(ioPopMeas)

    if performProj
        close(ioProjMeasUp)
        close(ioProjMeasDn)
    end

    if performProj2 && performProj
        close(ioProjMeas2Up)
        close(ioProjMeas2Dn)
    end
    if performOcc
        close(ioOcc)
    end
    if trackBondDim
        close(ioTrackBondDim)
    end

    return Gammas, Lambdas

end

 function SingleEx_evolution_TEBD(Gammas, Lambdas, s; 
    ϵ,
    Δ, 
    sysenvInt::String, 
    ChainLength, 
    tau, 
    ttotal,
    measStep,
    cutoff=1E-14, 
    minBondDim=5, 
    maxBondDim=100, 
    freqs, 
    coups)

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
    initmag = noprime!(op("Z",s[1])*Gammas[1])*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
    println("first value of exp Z: ", scalar(initmag))


    # # Define the gates as in the trotter formula e^[(F+G)t]=e^[F*t/2]*e^[G*t]*e^[F*t/2]
    # # F odd gate, G even gate

    gates = ITensor[]

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

    # Open files to save data 
    ioMagMeas = open("Data/magMeas_TEBD.dat", "w")
    ioTv = open("Data/tv_TEBD.dat", "w")
    ioNormCheck = open("Data/normCheck_TEBD.dat", "w")
    ioPopMeas = open("Data/popMeas_TEBD.dat","w")

    # If needed save them to be stored in arrays
    #magMeas = Vector{ComplexF64}()
    #tv = Vector{Float64}()
    #normCheck = Vector{ComplexF64}()
    #popMeas = Vector{Vector{ComplexF64}}

    #psiV = MPS[]
    #LambdasV = Vector{Vector{ITensor}}()
    #GammasV = Vector{Vector{ITensor}}()


    
    for (step,t) in enumerate(0.0:tau:ttotal)
        println("tempo: ", t)
        
        #psi = convert_to_MPS(Gammas, Lambdas, ntot)
        if (step-1) % measStep == 0
            println(ioTv, t)
            #Immediately write to file without buffering
            flush(ioTv)
            #push!(tv, t)
            
            
            #Magnetization measure
            
            #Anna's approach
            #psiAppo = convert_to_MPS(Gammas, Lambdas, ntot)
            #writedlm(ioMagMeas, [t expect(psiAppo, "Z",sites=1:1)], ',')
            #appo =  noprime!(op("Z",s[1])*Gammas[1])*Lambdas[1]*Lambdas[1]*conj(Gammas[1])
            appo =  noprime!(op("Z",s[1])*Gammas[1])*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
            
            # #appo =  noprime!(op("Z",s[1])*Gammas[1])*conj(Gammas[1])
            writedlm(ioMagMeas, [t scalar(appo)], ',')
            #Immediately write to file without buffering
            flush(ioMagMeas)
            # push!(magMeas, scalar(appo))

            #Chain occupation measure
            #this is delicate: we need to measure all the chain sites
            #This will be painful, but needed.
            #Here I implement full measure, but I suspect that it 
            #can be otimized somehow:
            #1) Far sites are perturbed at later Times
            #2) This can be parallelized over threads: instantiate whole array
            occMeasures = Vector{ComplexF64}([])
            for i in 2:ntot-1
                appoOcc = Lambdas[i-1]*noprime!(op("N",s[i])*Gammas[i])*Lambdas[i]*dag(Lambdas[i-1]*Gammas[i]*Lambdas[i])
                #appoOcc = Lambdas[i-1]*Lambdas[i-1]*noprime!(op("N",s[i])*Gammas[i])*Lambdas[i]*Lambdas[i]*conj(Gammas[i])
                push!(occMeasures, scalar(appoOcc))
            end
            #Measure last chain site population
            appoOcc = Lambdas[ntot-1]*noprime!(op("N",s[ntot])*Gammas[ntot])*dag(Lambdas[ntot-1]*Gammas[ntot])
            push!(occMeasures, scalar(appoOcc))
            writedlm(ioPopMeas, transpose(vcat(t, occMeasures)), ',')
            #Immediately write to file without buffering
            flush(ioPopMeas)


            
            #writedlm(ioPopMeas, [t expect(psiAppo, "N",sites=2:2)], ',')

            #Norm measure
            #Maybe Lambdas[1]*Lambdas[1] is enough since > b < * dag(> b < ) = b * b 
            appoNorm = Gammas[1]*Lambdas[1]*dag(Gammas[1]*Lambdas[1])
            writedlm(ioNormCheck, [t sqrt(scalar(appoNorm))],',')
            #Immediately write to file without buffering
            flush(ioNormCheck)
            #push!(normCheck, scalar(appoNorm))
        end
        
        #push!(psiV, psi)
        #push!(GammasV, Gammas)
        #push!(LambdasV, Lambdas)

        Gammas, Lambdas = apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim=maxBondDim)

    end

    close(ioTv)
    close(ioMagMeas)
    close(ioNormCheck)
    close(ioPopMeas)

    return Gammas, Lambdas

end

function apply_TEBD(gates, Gammas, Lambdas; cutoff, mindim=5, maxBondDim)

    #println("first odd")
    Gammas, Lambdas = apply_odd(gates, Gammas, Lambdas; cutoff, mindim=1, maxBondDim=maxBondDim)
    #println("Gammas[1]",Gammas[1])
    #println("Gammas[2]",Gammas[2])
    #println("After first gate:")
    #println(Lambdas)
    #print_alternated_inds(Gammas, Lambdas)
    #pippo = readline()
    #println("the even step")
    Gammas, Lambdas = apply_even(gates, Gammas, Lambdas; cutoff, mindim=1, maxBondDim=maxBondDim)

    #println("second odd")
    #println("Gammas[1]: ", Gammas[1])
    Gammas, Lambdas = apply_odd(gates, Gammas, Lambdas; cutoff, mindim=1, maxBondDim=maxBondDim)    
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
    #sequentially and in multi-threading:
    
    #########ICE###############
    #BACKUP in TEDB_alt copy.jl
    #########ICE###############

    #Loop index arrays
    LambdasOdd = Array{ITensor}(undef,N) 
    #Site tensors are 1 more than link tensors
    GammasOdd = Array{ITensor}(undef,N+1)
    
    
    
    #print_inds(gates)

    # compute the inverse just for the right Lambdas (even), 
    # the inversion is well defined as in the previous SVD there is a cutoff that
    # cuts the smallest singular values decreasing the bond dimension
    # Q: perdo tanta informazione troncando?
    # Q: mindim non è compatibile, ma non capisco perchè
    LambdasInvEven = ITensor[]
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
    
    #println(Threads.threadid())
        
        #println("odd", i)
        #Temp array "local" to  for scope => one for each thread
        SVDLeftIndices = Index[]
        
        #empty!(SVDLeftIndices)
    
        #Apply the two site gate

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

        #NOTA: se uso mindim=mindim si rompre tutto!!
        A, S, B = length(SVDLeftIndices) == 1 ?
                  svd(Res, SVDLeftIndices[1]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim) :
                  svd(Res, SVDLeftIndices[1], SVDLeftIndices[2]; lefttags="u=$(i)", righttags="v=$(i)", cutoff=cutoff, use_absolute_cutoff=true, maxdim=maxBondDim)

        # return to Vidal form
        A = i != 1 ? LambdasInvEven[Int((i - 1) / 2)] * A : A
        B = i != N ? B * LambdasInvEven[Int((i + 1) / 2)] : B

        #All computation done.
        #Now we store the result

        if (i != 1)
            LambdasOdd[i-1] = Lambdas[i-1]
            #push!(LambdasOdd, Lambdas[i-1]) #in between the gates, leave the same lambdas
        end
        LambdasOdd[i] = S
        GammasOdd[i] = A
        GammasOdd[i+1] = B
        # push!(LambdasOdd, S)
        # push!(GammasOdd, A)
        # push!(GammasOdd, B)

    end


    if (N % 2 == 0) #N is even, so the number of sites is odd
        LambdasOdd[N] = Lambdas[N]
        GammasOdd[N+1] = Gammas[N+1]
        # push!(LambdasOdd, Lambdas[N])
        # push!(GammasOdd, Gammas[N+1])
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

    
    N = length(Lambdas)
    
    #println("Even apply")

    LambdasInvOdd = ITensor[]
    LambdasEven = Array{ITensor}(undef,N)
    GammasEven = Array{ITensor}(undef,N+1)
    
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
    GammasEven[1] = Gammas[1]
    #push!(GammasEven, Gammas[1])

Threads.@threads for i in 2:2:N

        #println("Even ",i)

        
        SVDLeftIndices = Index[]


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
            LambdasEven[i-1] = Lambdas[i-1]
            #push!(LambdasEven, Lambdas[i-1]) #in between the gates, leave the same lambdas
        end
        LambdasEven[i] = S
        GammasEven[i] = A
        GammasEven[i+1] = B
        # push!(LambdasEven, S)
        # push!(GammasEven, A)
        # push!(GammasEven, B)

    end


    if (N % 2 != 0) #N is odd, so the number of sites is even
        LambdasEven[N] = Lambdas[N]
        GammasEven[N+1] = Gammas[N+1]
        # push!(LambdasEven, Lambdas[N])
        # push!(GammasEven, Gammas[N+1])
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