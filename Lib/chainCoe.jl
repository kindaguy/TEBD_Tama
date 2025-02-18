using PolyChaos
using Integrals

function generateChainCoefficients(; sysEnergy, SpectralDensityInterval)

    e = sysEnergy

    # define the spectral density interval 
    domain = SpectralDensityInterval

    # define the spectral density function to use as weight
    weight(x) = e^2 * x * exp(-x / e)

    # defining the measure, i.e. the d_mu(omega)
    my_meas = Measure("my_meas", weight, domain, false, Dict())

    # define the recurrence coefficients of the monic orthogonal polynomials (cfr Note 1)
    my_OP = OrthoPoly("my_op", 100, my_meas; Nquad=200)

    # the frequency of the transformed oscillators
    freqs = my_OP.α
    # the coupling is sqrt of J(w)
    coups = sqrt.(my_OP.β)

    io = open("Data/freqs.dat", "w")
    for f in freqs
        println(io, f)
    end
    close(io)

    io = open("Data/coups.dat", "w")
    for f in coups
        println(io, f)
    end
    close(io)

end