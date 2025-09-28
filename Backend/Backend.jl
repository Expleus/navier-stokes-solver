using Random, Plots, Distributions, RadialBasisFunctions, NonlinearSolve, SparseArrays, DelimitedFiles, LinearAlgebra, CairoMakie, ForwardDiff

include("DataStructures.jl")
include("Puntos.jl")
include("Matrices.jl")
include("Solucionador.jl")
include("Graficas.jl")

BLAS.set_num_threads(Sys.CPU_THREADS)

#RBF

ϕ = TPS(3) #RBF de velocidad
Φ = TPS(3) #RBF de presion
ψ = TPS(3) #RBF de temperatura

global PuntosRBF = []

function RBF_eval(x,y,w)
    punto = [x,y]
    return dot(w, [ϕ(punto-c) for c in PuntosRBF])
end

function problema(prob, x, select)
    
    global PuntosRBF

    if x == 1
        println(">> Problema con velocidades y presiones desconocidas")

        total_time = @elapsed begin

            if select == 1
                t = @elapsed centro, max, min = generar_fronteras_rng(prob)
                println("Generar puntos tomó $(round(t, digits=3)) s")
            elseif select == 2
                t = @elapsed centro, max, min = generar_fronteras_flujo(prob)
                println("Generar puntos tomó $(round(t, digits=3)) s")
            end

            t = @elapsed A, b, evaluacion = matriz_dimensional_vp(prob, centro, ϕ, Φ)
            println("Llenar la matriz tomó $(round(t, digits=3)) s")
            
            println("Puntos: $(length(evaluacion))")

            t = @elapsed w = solucionador_vp(A,b,evaluacion)
            println("Solucionar el sistema tomó $(round(t, digits=3)) s")

            PuntosRBF = evaluacion

            t = @elapsed fig = grafica_vp_dimensional(prob, w, min, max)
            println("Graficar tomó $(round(t, digits=3)) s")
        end

        println(">> Tiempo total: $(round(total_time, digits=3)) s")

        return w, evaluacion, fig

    elseif x == 2
        println(">> Problema con velocidades, presiones y temperaturas desconocidas")

        total_time = @elapsed begin

            if select == 1
                t = @elapsed centro, max, min = generar_fronteras_rng(prob)
                println("Generar puntos tomó $(round(t, digits=3)) s")
            elseif select == 2
                t = @elapsed centro, max, min = generar_fronteras_flujo(prob)
                println("Generar puntos tomó $(round(t, digits=3)) s")
            end

            t = @elapsed A, b, evaluacion = matriz_dimensional_vpt(prob, centro, ϕ, Φ, ψ)
            println("Llenar la matriz tomó $(round(t, digits=3)) s")
            
            println("Puntos: $(length(evaluacion))")

            t = @elapsed w = solucionador_vpt(A,b,evaluacion)
            println("Solucionar el sistema tomó $(round(t, digits=3)) s")

            PuntosRBF = evaluacion

            t = @elapsed fig = grafica_vpt_dimensional(prob, w, min, max)
            println("Graficar tomó $(round(t, digits=3)) s")
        end

        println(">> Tiempo total: $(round(total_time, digits=3)) s")

        return w, evaluacion, fig
    elseif x == 3
        println("Coming soon...")
    else
        println("Coming soon...")
    end

end