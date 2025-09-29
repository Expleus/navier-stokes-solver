

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
        label = "Presion (Pa)",
        limits = (minimum(filter(!isnan, M)), maximum(filter(!isnan, M))),
    )

    # ============================================================
    # Guardar figura combinada
    # ============================================================
    save("C:\\Users\\exple\\OneDrive\\Desktop\\Compartir\\figuras_combinadas.png", fig)
    save("Graficos/Perfil_de_velocidades_combinado.png", fig)

    return fig

end

function grafica_vpt_dimensional(prob, w, min, max)
   # Parámetros de evaluación y visualización
   Lx = prob.Lx
   x_eval = LinRange(0, Lx, 100)
   y_eval = LinRange(min, max, 100)
   quiver_dense_x = 21
   quiver_dense_y = 21
   escala = 3
   δ = 0
   quiver_lines_x = Array(LinRange(δ, Lx-δ, quiver_dense_x+2))
   quiver_lines_y = Array(LinRange(min, max, quiver_dense_y))

   # =====================
   # Inicializar figura y layout
   # =====================
   fig = Figure(size = (1600, 900))
   grid = fig[1, 1] = GridLayout()

   # =====================
   # Panel 1: Perfil de velocidades (izquierda)
   # =====================
   ax1 = Axis(grid[1:2, 1], limits = ((0, Lx+1), (min-0.1, max+0.1)),
      xlabel = "x", ylabel = "y", title = "Perfil de velocidades",
      titlesize = 20)
   # Fronteras superior e inferior
   lines!(ax1, x_eval, S.(x_eval), color = :black, linewidth = 10)
   lines!(ax1, x_eval, I.(x_eval), color = :black, linewidth = 10)

   # Calcular campo de velocidades
   all_x, all_y, all_u, all_v = Float64[], Float64[], Float64[], Float64[]
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
   arrows2d!(ax1, all_x, all_y, all_u, all_v; lengthscale = 1, tipwidth = 10, color = mag_global, colormap = :imola)
   cb1 = Colorbar(grid[1:2, 2], colormap = :imola, label = "Magnitud (m/s)", limits = (minimum(mag_global), maximum(mag_global)))

   # =====================
   # Panel 2: Presión (arriba derecha)
   # =====================
   ax2 = Axis(grid[1, 3], limits = ((0, Lx+1), (min-0.1, max+0.1)),
      xlabel = "x", ylabel = "y", title = "Presión",
      titlesize = 20)
   nx, ny = 400, 400
   xgrid = range(0, Lx, length=nx)
   ygrid = range(min, max, length=ny)
   M = [ (I(xi) ≤ yi ≤ S(xi)) ? RBF_eval(xi, yi, w.p) : NaN for xi in xgrid, yi in ygrid ]
   CairoMakie.heatmap!(ax2, xgrid, ygrid, M; colormap = :imola, transparency = true)
   lines!(ax2, x_eval, S.(x_eval), color = :black, linewidth = 10)
   lines!(ax2, x_eval, I.(x_eval), color = :black, linewidth = 10)
   cb2 = Colorbar(grid[1, 4], colormap = :imola, label = "Presion (Pa)", limits = (minimum(filter(!isnan, M)), maximum(filter(!isnan, M))))

   # =====================
   # Panel 3: Temperatura (abajo derecha)
   # =====================
   ax3 = Axis(grid[2, 3], limits = ((0, Lx+1), (min-0.1, max+0.1)),
      xlabel = "x", ylabel = "y", title = "Temperatura",
      titlesize = 20)
   Mt = [ (I(xi) ≤ yi ≤ S(xi)) ? RBF_eval(xi, yi, w.t) : NaN for xi in xgrid, yi in ygrid ]
   CairoMakie.heatmap!(ax3, xgrid, ygrid, Mt; colormap = :viridis, transparency = true)
   cb3 = Colorbar(grid[2, 4], colormap = :viridis, label = "Temperatura (K)", limits = (minimum(filter(!isnan, Mt)), maximum(filter(!isnan, Mt))))

   # =====================
   # Ajustar tamaño de columnas para layout
   # =====================
   colsize!(grid, 1, Relative(0.5))   # Panel de velocidad
   colsize!(grid, 2, Relative(0.05))  # Colorbar velocidad
   colsize!(grid, 3, Relative(0.40))  # Paneles presión/temperatura
   colsize!(grid, 4, Relative(0.05))  # Colorbars presión/temperatura

   # =====================
   # Guardar figura combinada
   # =====================
   save(raw"C:\Users\exple\OneDrive\Desktop\Compartir\figuras_combinadas.png", fig)
   save("Graficos/Perfil_de_velocidades_combinado.png", fig)

   return fig
end