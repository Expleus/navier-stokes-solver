

function grafica_vp_dimensional(prob, w, min, max)
    
    Lx = prob.Lx

    x_eval = LinRange(0, Lx, 100)
    y_eval = LinRange(min, max, 100)

    quiver_dense_x = 21
    quiver_dense_y = 21
    escala = 5
    δ = 0

    quiver_lines_x = Array(LinRange(δ, Lx-δ, quiver_dense_x+2))
    quiver_lines_y = Array(LinRange(min, max, quiver_dense_y))

    fig = Figure(size = (1200, 1400))

    # ============================================================
    # Panel 1: Quiver con magnitud
    # ============================================================
    ax1 = Axis(fig[1, 1],
        limits = ((0, Lx+1), (min-0.1, max+0.1)),
        xlabel = "x", ylabel = "y",
        title = "Perfil de velocidades"
    )

    lines!(ax1, x_eval, S.(x_eval), color = :black, linewidth = 10)
    lines!(ax1, x_eval, I.(x_eval), color = :black, linewidth = 10)

    all_x = Float64[]
    all_y = Float64[]
    all_u = Float64[]
    all_v = Float64[]

    for i in quiver_lines_x
        x = i .* ones(quiver_dense_y)
        y = Array(LinRange(I(i), S(i), quiver_dense_y))

        u = escala .* RBF_eval.(x, y, Ref(w.x))
        v = escala .* RBF_eval.(x, y, Ref(w.y))

        append!(all_x, x)
        append!(all_y, y)
        append!(all_u, u)
        append!(all_v, v)
    end

    mag_global = sqrt.(all_u.^2 .+ all_v.^2)

    arrows2d!(ax1, all_x, all_y, all_u, all_v;
        lengthscale = 1,
        tipwidth = 10,
        color = mag_global,
        colormap = :imola
    )

    Colorbar(fig[1, 2],
        colormap = :imola,
        label = "Magnitud (m/s)",
        limits = (minimum(mag_global), maximum(mag_global)),
    )

    # ============================================================
    # Panel 2: Heatmap con magnitud
    # ============================================================
    ax2 = Axis(fig[2, 1],
        limits = ((0, Lx+1), (min-0.1, max+0.1)),
        xlabel = "x", ylabel = "y",
        title = "Presión"
    )

    nx, ny = 200, 200
    xgrid = range(0, Lx, length=nx)
    ygrid = range(min, max, length=ny)

    M = [ (I(xi) ≤ yi ≤ S(xi)) ? RBF_eval(xi, yi, w.p) : NaN
        for xi in xgrid, yi in ygrid ]

    CairoMakie.heatmap!(ax2, xgrid, ygrid, M;
        colormap = :imola,
        transparency = true
    )

    lines!(ax2, x_eval, S.(x_eval), color = :black, linewidth = 10)
    lines!(ax2, x_eval, I.(x_eval), color = :black, linewidth = 10)

    Colorbar(fig[2, 2],
        colormap = :imola,
        label = "Presión (Pa)"
    )

    # ============================================================
    # Guardar figura combinada
    # ============================================================
    save("C:\\Users\\exple\\OneDrive\\Desktop\\Compartir\\figuras_combinadas.png", fig)
    save("Perfil_de_velocidades_combinado.png", fig)

    return fig

end

function grafica_vpt_dimensional(prob, w, min, max)
    
    Lx = prob.Lx

    x_eval = LinRange(0, Lx, 100)
    y_eval = LinRange(min, max, 100)

    quiver_dense_x = 21
    quiver_dense_y = 21
    escala = 5
    δ = 0

    quiver_lines_x = Array(LinRange(δ, Lx-δ, quiver_dense_x+2))
    quiver_lines_y = Array(LinRange(min, max, quiver_dense_y))

    fig = Figure(size = (1200, 2100))

    # ============================================================
    # Panel 1: Quiver con magnitud
    # ============================================================
    ax1 = Axis(fig[1, 1],
        limits = ((0, Lx+1), (min-0.1, max+0.1)),
        xlabel = "x", ylabel = "y",
        title = "Perfil de velocidades"
    )

    lines!(ax1, x_eval, S.(x_eval), color = :black, linewidth = 10)
    lines!(ax1, x_eval, I.(x_eval), color = :black, linewidth = 10)

    all_x = Float64[]
    all_y = Float64[]
    all_u = Float64[]
    all_v = Float64[]

    for i in quiver_lines_x
        x = i .* ones(quiver_dense_y)
        y = Array(LinRange(I(i), S(i), quiver_dense_y))

        u = escala .* RBF_eval.(x, y, Ref(w.x))
        v = escala .* RBF_eval.(x, y, Ref(w.y))

        append!(all_x, x)
        append!(all_y, y)
        append!(all_u, u)
        append!(all_v, v)
    end

    mag_global = sqrt.(all_u.^2 .+ all_v.^2) ./ escala

    arrows2d!(ax1, all_x, all_y, all_u, all_v;
        lengthscale = 1,
        tipwidth = 10,
        color = mag_global,
        colormap = :imola
    )

    Colorbar(fig[1, 2],
        colormap = :imola,
        label = "Magnitud (m/s)",
        limits = (minimum(mag_global), maximum(mag_global)),
    )

    # ============================================================
    # Panel 2: Heatmap con magnitud
    # ============================================================
    ax2 = Axis(fig[2, 1],
        limits = ((0, Lx+1), (min-0.1, max+0.1)),
        xlabel = "x", ylabel = "y",
        title = "Presión"
    )

    nx, ny = 200, 200
    xgrid = range(0, Lx, length=nx)
    ygrid = range(min, max, length=ny)

    M = [ (I(xi) ≤ yi ≤ S(xi)) ? RBF_eval(xi, yi, w.p) : NaN
        for xi in xgrid, yi in ygrid ]

    CairoMakie.heatmap!(ax2, xgrid, ygrid, M;
        colormap = :imola,
        transparency = true
    )

    lines!(ax2, x_eval, S.(x_eval), color = :black, linewidth = 10)
    lines!(ax2, x_eval, I.(x_eval), color = :black, linewidth = 10)

    Colorbar(fig[2, 2],
        colormap = :imola,
        label = "Presion (Pa)",
        limits = (minimum(M), maximum(M))
    )

    ax3 = Axis(fig[3, 1],
        limits = ((0, Lx+1), (min-0.1, max+0.1)),
        xlabel = "x", ylabel = "y",
        title = "Temperatura"
    )  

    Mt = [ (I(xi) ≤ yi ≤ S(xi)) ? RBF_eval(xi, yi, w.t) : NaN
        for xi in xgrid, yi in ygrid ]


    CairoMakie.heatmap!(ax3, xgrid, ygrid, Mt;
        colormap = :viridis,
        transparency = true
    )

    Colorbar(fig[3, 2],
        colormap = :viridis,
        label = "Temperatura (K)",
        limits = (minimum(Mt), maximum(Mt))
    )

    # ============================================================
    # Guardar figura combinada
    # ============================================================
    save("C:\\Users\\exple\\OneDrive\\Desktop\\Compartir\\figuras_combinadas.png", fig)
    save("Perfil_de_velocidades_combinado.png", fig)

    return fig

end