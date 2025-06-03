function print_alternated_inds(tensorG, tensorL)


    println("Vidal:")
    for k in 1:length(tensorG)
        if k <= length(tensorG)
            println("  Gamma: ", k)
            for idx in inds(tensorG[k])
                println("  Indice: ", idx)
            end
        end
        if k <= length(tensorL)
            println("  lambda: ", k)
            for idx in inds(tensorL[k])
                println("  Indice: ", idx)
            end
        end
    end

end

function print_inds(tensori)

    for (n, T) in enumerate(tensori)
        println("Tensore $(n):")
        for idx in inds(T)
            println("  Indice: ", idx)
        end
    end

end