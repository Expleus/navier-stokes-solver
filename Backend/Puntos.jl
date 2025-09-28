function generar_fronteras_rng(prob)

    d = Uniform(0,prob.Lx)

    dominio_frontera_x = LinRange(0,prob.Lx,prob.p_fx)
    dominio_frontera_izq = LinRange(prob.inferior(0),prob.superior(0), prob.p_fy+2)
    dominio_frontera_der = LinRange(prob.inferior(prob.Lx),prob.superior(prob.Lx), prob.p_fy+2)

    frontera_superior = [[u, prob.superior(u)] for u in dominio_frontera_x]
    frontera_inferior = [[u, prob.inferior(u)] for u in dominio_frontera_x]

    max = maximum(p[2] for p in frontera_superior)
    min = minimum(p[2] for p in frontera_inferior)
    M = maximum(prob.superior(x) - prob.inferior(x) for x in LinRange(0,prob.Lx,1000))

    frontera_izquierda = [[0, u] for u in dominio_frontera_izq[2:end-1]]
    frontera_derecha   = [[prob.Lx, u] for u in dominio_frontera_der[2:end-1]]

    interior = vcat(frontera_izquierda, frontera_derecha)

    for k in 1:prob.p_interior
        while true
            x = rand(d)             # candidato en [a,b]
            u = rand() * M          # candidato en [0,M]
            ancho = S(x) - I(x)
            if u < ancho
                y = I(x) + ancho * rand()
                if minimum(norm.(vcat(interior, frontera_superior, frontera_inferior) .- [[x,y]])) > 0.1
                    interior = vcat([[x,y]], interior) 
                    break
                end
            end
        end
    end

    p = Plots.plot()

    for pts in (frontera_superior, frontera_inferior, frontera_izquierda, frontera_derecha, interior)
        xs = [p[1] for p in pts]
        ys = [p[2] for p in pts]
        Plots.scatter!(xs, ys, legend=false)
    end

    savefig("puntos.png")

    Plots.display(p)

    centro = centros(frontera_superior, frontera_inferior, frontera_izquierda, frontera_derecha, interior)

    return centro, max, min

end

function generar_fronteras_flujo(prob)

    d = Uniform(0,prob.Lx)

    dominio_frontera_x = LinRange(0,prob.Lx,prob.p_fx)
    dominio_frontera_izq = LinRange(prob.inferior(0),prob.superior(0), prob.p_fy+2)
    dominio_frontera_der = LinRange(prob.inferior(prob.Lx),prob.superior(prob.Lx), prob.p_fy+2)

    frontera_superior = [[u, prob.superior(u)] for u in dominio_frontera_x]
    frontera_inferior = [[u, I(u)] for u in dominio_frontera_x]

    max = maximum(p[2] for p in frontera_superior)
    min = minimum(p[2] for p in frontera_inferior)
    M = maximum(S(x) - I(x) for x in LinRange(0,prob.Lx,1000))

    frontera_izquierda = [[0, u] for u in dominio_frontera_izq[2:end-1]]
    frontera_derecha   = [[prob.Lx, u] for u in dominio_frontera_der[2:end-1]]

    t = Array(LinRange(0,1,prob.p_fy+2))[2:end-1]
    x_in = LinRange(frontera_superior[2][1],frontera_superior[end-1][1], prob.p_interior÷prob.p_fy)

    interior = vcat(frontera_izquierda, frontera_derecha)

    flujo(a,x) = (1-a)*prob.inferior(x)+a*prob.superior(x)

    for i in t
        v = [[u,flujo(i,u)] for u in x_in]
        interior = vcat(v, interior) 
    end

    p = Plots.plot()

    for pts in (frontera_superior, frontera_inferior, frontera_izquierda, frontera_derecha, interior)
        xs = [p[1] for p in pts]
        ys = [p[2] for p in pts]
        Plots.scatter!(xs, ys, legend=false)
    end

    savefig("puntos.png")

    Plots.display(p)

    centro = centros(frontera_superior, frontera_inferior, frontera_izquierda, frontera_derecha, interior)

    return centro, max, min

end
