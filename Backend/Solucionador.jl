function solucionador_vp(A,b, evaluacion)

    function mixer(u, v, w)
        n = length(u)
        m = length(w)
        T = promote_type(eltype(u), eltype(v), eltype(w))

        # Prealocar resultado con todo el espacio
        result = Vector{T}(undef, 2n + 3n^2 + m)

        # Copiar u
        copyto!(result, 1, u, 1, n)

        # uu = u ⊗ u
        idx0 = n + 1
        Threads.@threads for j in 1:n
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = u[i] * u[j]
            end
        end

        # vu = v ⊗ u
        idx0 = n + n^2 + 1
        Threads.@threads for j in 1:n
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = v[i] * u[j]
            end
        end

        # vv = v ⊗ v
        idx0 = n + 2n^2 + 1
        Threads.@threads for j in 1:n
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = v[i] * v[j]
            end
        end

        # Copiar v
        copyto!(result, n + 3n^2 + 1, v, 1, n)

        # Copiar w
        copyto!(result, 2n + 3n^2 + 1, w, 1, m)

        return result
    end

    numero_puntos = length(evaluacion)

    w0 = vcat(0.01*randn(3*numero_puntos), zeros(length(b)-3*numero_puntos))
    F = zeros(length(b))

    r! = function(F,w,p)

        residuo(F,w,A,b,numero_puntos)

    end
    
    function residuo(F,w,A,b,numero_puntos)
        
        wx = @view w[1:numero_puntos]
        wy = @view w[numero_puntos+1:2*numero_puntos]
        wp = @view w[2*numero_puntos+1:3*numero_puntos]

        x = mixer(wx,wy,wp)

        mul!(F,A,x)

        F .-= b

    end

    prob = NonlinearProblem(r!, w0)
    sol = solve(prob, LevenbergMarquardt();
                     maxiters = 1000,
                     abstol = 1e-8,
                     reltol = 1e-8)

    wx = sol.u[1:numero_puntos]
    wy = sol.u[numero_puntos+1:2*numero_puntos]
    wp = sol.u[2*numero_puntos+1:3*numero_puntos]

    peso = pesos(wx,wy,wp,wp)

    println("Residuo final: ", norm(A*mixer(wx, wy, wp)-b))

    return peso

end

function solucionador_vpt(A,b, evaluacion)

    function mixer(u, v, w, s)
        n = length(u)
        m = length(w)
        k = length(s)
        T = promote_type(eltype(u), eltype(v), eltype(w), eltype(s))

        # Tamaño total: u + uu + vu + vv + v + w + s + us + vs
        total_len = 2n + 3n^2 + m + k + 2n*k
        result = Vector{T}(undef, total_len)

        # --- Copiar u ---
        copyto!(result, 1, u, 1, n)

        # --- uu = u ⊗ u ---
        idx0 = n + 1
        Threads.@threads for j in 1:n
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = u[i] * u[j]
            end
        end

        # --- vu = v ⊗ u ---
        idx0 = n + n^2 + 1
        Threads.@threads for j in 1:n
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = v[i] * u[j]
            end
        end

        # --- vv = v ⊗ v ---
        idx0 = n + 2n^2 + 1
        Threads.@threads for j in 1:n
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = v[i] * v[j]
            end
        end

        # --- Copiar v ---
        copyto!(result, n + 3n^2 + 1, v, 1, n)

        # --- Copiar w ---
        copyto!(result, 2n + 3n^2 + 1, w, 1, m)

        # --- Copiar s ---
        copyto!(result, 2n + 3n^2 + m + 1, s, 1, k)

        # --- us = u ⊗ s ---
        idx0 = 2n + 3n^2 + m + k + 1
        Threads.@threads for j in 1:k
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = u[i] * s[j]
            end
        end

        # --- vs = v ⊗ s ---
        idx0 += n*k
        Threads.@threads for j in 1:k
            offset = idx0 + (j-1)*n
            for i in 1:n
                @inbounds result[offset + i - 1] = v[i] * s[j]
            end
        end

        return result
    end

    numero_puntos = length(evaluacion)

    w0 = vcat(0.01*randn(4*numero_puntos), zeros(length(b)-4*numero_puntos))
    F = zeros(length(b))

    r! = function(F,w,p)

        residuo(F,w,A,b,numero_puntos)

    end
    
    function residuo(F,w,A,b,numero_puntos)
        
        wx = @view w[1:numero_puntos]
        wy = @view w[numero_puntos+1:2*numero_puntos]
        wp = @view w[2*numero_puntos+1:3*numero_puntos]
        wt = @view w[3*numero_puntos+1:4*numero_puntos]

        x = mixer(wx,wy,wp,wt)

        mul!(F,A,x)

        F .-= b

    end

    prob = NonlinearProblem(r!, w0)
    sol = solve(prob, LevenbergMarquardt();
                     maxiters = 1000,
                     abstol = 1e-8,
                     reltol = 1e-8)

    wx = sol.u[1:numero_puntos]
    wy = sol.u[numero_puntos+1:2*numero_puntos]
    wp = sol.u[2*numero_puntos+1:3*numero_puntos]
    wt = sol.u[3*numero_puntos+1:4*numero_puntos]

    peso = pesos(wx,wy,wp,wt)

    println("Residuo final: ", norm(A*mixer(wx, wy, wp,wt)-b))

    return peso

end
