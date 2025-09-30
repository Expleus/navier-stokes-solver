function matriz_dimensional_vp(prob, centro, ϕ, Φ)

    S = prob.superior
    I = prob.inferior

    frontera_superior = centro.superior
    frontera_inferior = centro.inferior
    interior = centro.interior
    frontera_izquierda = centro.izquierda
    frontera_derecha = centro.derecha

    dS(x) = ForwardDiff.derivative(y -> S(y), x)
    dI(x) = ForwardDiff.derivative(y -> I(y), x)

    evaluacion = vcat(frontera_superior, frontera_inferior, interior)
    
    puntos_evaluacion = length(evaluacion)
    puntos_frontera = length(frontera_superior)
    puntos_interior = length(interior)
    puntos_presion = length(frontera_derecha)
    
    #Matrices con datos

    matriz_laplaciano = zeros(puntos_interior, puntos_evaluacion)
    matriz_frontera_superior = zeros(puntos_frontera, puntos_evaluacion)
    matriz_frontera_inferior = zeros(puntos_frontera, puntos_evaluacion)
    matriz_continuidad_x = zeros(puntos_interior, puntos_evaluacion)
    matriz_continuidad_y = zeros(puntos_interior, puntos_evaluacion)
    matriz_Wxx_Wyx = zeros(puntos_interior, puntos_evaluacion^2)
    matriz_Wyy_Wxy = zeros(puntos_interior, puntos_evaluacion^2)

    matriz_px = zeros(puntos_interior,puntos_evaluacion)
    matriz_py = zeros(puntos_interior,puntos_evaluacion)
    matriz_frontera_izquierda = zeros(puntos_presion,puntos_evaluacion)
    matriz_frontera_derecha = zeros(puntos_presion,puntos_evaluacion)
    matriz_laplaciano_f1_x = zeros(puntos_frontera, puntos_evaluacion)
    matriz_laplaciano_f1_y = zeros(puntos_frontera, puntos_evaluacion)
    matriz_laplaciano_f2_x = zeros(puntos_frontera, puntos_evaluacion)
    matriz_laplaciano_f2_y = zeros(puntos_frontera, puntos_evaluacion)
    matriz_f1_pn = zeros(puntos_frontera, puntos_evaluacion)
    matriz_f2_pn = zeros(puntos_frontera, puntos_evaluacion)

    #Matrices de relleno

    bloque_A = zeros(2*puntos_frontera,2*puntos_evaluacion^2)
    bloque_B = zeros(puntos_interior+2*puntos_frontera, puntos_evaluacion+puntos_evaluacion^2)
    bloque_C = zeros(puntos_interior,3*puntos_evaluacion^2)
    Bloque_D = zeros(puntos_evaluacion,puntos_evaluacion)
    Bloque_E = zeros(2*puntos_frontera,puntos_evaluacion)
    Bloque_F = zeros(2*puntos_presion,3*puntos_evaluacion^2+2*puntos_evaluacion)
    Bloque_G = zeros(2*puntos_frontera,3*puntos_evaluacion^2)

    #Matriz laplaciano 

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)#i:puntos_evaluacion
            matriz_laplaciano[i,j] = prob.μ*∇²(ϕ,interior[i]-evaluacion[j])
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            matriz_px[i,j] = -dX(Φ,interior[i]-evaluacion[j],1)
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            matriz_py[i,j] = -dX(Φ,interior[i]-evaluacion[j],2)
        end
    end

    #Condiciones de frontera
                                                                                
    #Velocidad inferior
    Threads.@threads for i in eachindex(frontera_superior)
        for j in eachindex(evaluacion)
            matriz_frontera_superior[i,j] = ϕ(frontera_superior[i]-evaluacion[j])
        end
    end
    #Velocidad superior
    Threads.@threads for i in eachindex(frontera_inferior)
        for j in eachindex(evaluacion)
            matriz_frontera_inferior[i,j] = ϕ(frontera_inferior[i]-evaluacion[j])
        end
    end

    #Condicion de continuidad

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            matriz_continuidad_x[i,j] = dX(ϕ,interior[i]-evaluacion[j],1)
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            matriz_continuidad_y[i,j] = dX(ϕ,interior[i]-evaluacion[j],2)
        end
    end

    #Condiciones a la entrada y salida

    Threads.@threads for i in eachindex(frontera_izquierda)
        for j in eachindex(evaluacion)
            matriz_frontera_izquierda[i,j] = Φ(frontera_izquierda[i]-evaluacion[j])
        end
    end    

    Threads.@threads for i in eachindex(frontera_derecha)
        for j in eachindex(evaluacion)
            matriz_frontera_derecha[i,j] = Φ(frontera_derecha[i]-evaluacion[j])
        end
    end

    #Balance normal

    Threads.@threads for i in eachindex(frontera_superior)
        for j in eachindex(evaluacion)
            n = [-dS(frontera_superior[i][1]),1]
            n = n ./ norm(n)
            matriz_laplaciano_f1_x[i,j] = n[1]*μ*∇²(ϕ,frontera_superior[i]-evaluacion[j])
            matriz_laplaciano_f1_y[i,j] = n[2]*μ*∇²(ϕ,frontera_superior[i]-evaluacion[j])
            matriz_f1_pn[i,j] = -(n[1]*dX(ϕ,frontera_superior[i]-evaluacion[j],1) + n[2]*dX(ϕ,frontera_superior[i]-evaluacion[j],2))
        end
    end

    Threads.@threads for i in eachindex(frontera_inferior)
        for j in eachindex(evaluacion)
            n = [dI(frontera_inferior[i][1]),-1]
            n = n ./ norm(n)
            matriz_laplaciano_f2_x[i,j] = n[1]*μ*∇²(ϕ,frontera_inferior[i]-evaluacion[j])
            matriz_laplaciano_f2_y[i,j] = n[2]*μ*∇²(ϕ,frontera_inferior[i]-evaluacion[j])
            matriz_f2_pn[i,j] = -(n[1]*dX(ϕ,frontera_inferior[i]-evaluacion[j],1) + n[2]*dX(ϕ,frontera_inferior[i]-evaluacion[j],2))
        end
    end

    #Matrices mixtas

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            for k in eachindex(evaluacion)
                #Los pesos en x tienen prioridad en la doble suma cruzada: wx1wy1+wx1+wy2+...
                matriz_Wxx_Wyx[i,length(evaluacion)*(j-1)+k] = -1*ρ*ϕ(interior[i]-evaluacion[j])*dX(ϕ,interior[i]-evaluacion[k],1)
            end
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            for k in eachindex(evaluacion)
                #Los pesos en x tienen prioridad en la doble suma cruzada: wx1wy1+wx1+wy2+...
                matriz_Wyy_Wxy[i,length(evaluacion)*(j-1)+k] = -1*ρ*ϕ(interior[i]-evaluacion[k])*dX(ϕ,interior[i]-evaluacion[j],2)
            end
        end
    end

    #Armado de la matriz

    Matriz_sup = [[matriz_laplaciano; matriz_frontera_superior; matriz_frontera_inferior] [[matriz_Wxx_Wyx matriz_Wyy_Wxy]; bloque_A] bloque_B]
    Matriz_continuidad = [matriz_continuidad_x bloque_C matriz_continuidad_y]
    Matriz_combinada = [[Matriz_sup; Matriz_continuidad] [matriz_px; Bloque_D]]
    Matriz_normal = [[matriz_laplaciano_f1_x ; matriz_laplaciano_f2_x] Bloque_G [[matriz_laplaciano_f1_y ; matriz_laplaciano_f2_y] [matriz_f1_pn ; matriz_f2_pn]]]
    Matriz_inf = [bloque_B [[matriz_Wxx_Wyx matriz_Wyy_Wxy];bloque_A] [matriz_laplaciano; matriz_frontera_superior; matriz_frontera_inferior] [matriz_py; Bloque_E]]
    Matriz_presion = [Bloque_F [matriz_frontera_izquierda; matriz_frontera_derecha]]

    A_temp = [Matriz_combinada; Matriz_normal; Matriz_inf]

    temp = size(A_temp,1)

    A_densa = [A_temp; Matriz_presion]

    A = sparse(A_densa)

    b = zeros(size(A,1))

    for i in (temp+1):(temp+puntos_presion)
        b[i] = prob.P
    end

    return A, b, evaluacion
end

function matriz_dimensional_vpt(prob, centro, ϕ, Φ, ψ)
    
    A, b, evaluacion = matriz_dimensional_vp(prob, centro, ϕ, Φ)

    frontera_superior = centro.superior
    frontera_inferior = centro.inferior
    interior = centro.interior
    frontera_izquierda = centro.izquierda
    frontera_derecha = centro.derecha

    puntos_evaluacion = length(evaluacion)
    puntos_frontera = length(frontera_superior)
    puntos_interior = length(interior)
    puntos_presion = length(frontera_derecha)

    #Matrices de relleno

    bloque_H = zeros(size(A,1), puntos_evaluacion+2*puntos_evaluacion^2)
    bloque_I1 = zeros(puntos_interior, puntos_evaluacion)
    bloque_I2 = zeros(puntos_interior, 2*puntos_evaluacion)
    bloque_I = zeros(2*puntos_frontera+2*puntos_presion, size(A,2))
    bloque_J = zeros(2*puntos_frontera+2*puntos_presion, 2*puntos_evaluacion^2)

    #Matrices con datos 

    matriz_laplaciano = zeros(puntos_interior, puntos_evaluacion)
    matriz_frontera_superior = zeros(puntos_frontera, puntos_evaluacion)
    matriz_frontera_inferior = zeros(puntos_frontera, puntos_evaluacion)
    matriz_frontera_izquierda = zeros(puntos_presion,puntos_evaluacion)
    matriz_frontera_derecha = zeros(puntos_presion,puntos_evaluacion)
    matriz_Wxx = zeros(puntos_interior, puntos_evaluacion^2)
    matriz_Wyy = zeros(puntos_interior, puntos_evaluacion^2)
    matriz_Wxy = zeros(puntos_interior, puntos_evaluacion^2)
    matriz_Wxt = zeros(puntos_interior, puntos_evaluacion^2)
    matriz_Wyt = zeros(puntos_interior, puntos_evaluacion^2)

    #Matriz laplaciano 

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)#i:puntos_evaluacion
            matriz_laplaciano[i,j] = prob.k*∇²(ψ,interior[i]-evaluacion[j])
        end
    end

    #Matrices mixtas

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            for k in eachindex(evaluacion)
                #Los pesos en x tienen prioridad en la doble suma cruzada: wx1wy1+wx1+wy2+...
                matriz_Wxx[i,length(evaluacion)*(j-1)+k] = prob.μ*(4*dX(ϕ,interior[i]-evaluacion[j],1)*dX(ϕ,interior[i]-evaluacion[k],1)
                                                                + dX(ϕ,interior[i]-evaluacion[j],2)*dX(ϕ,interior[i]-evaluacion[k],2))
            end
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            for k in eachindex(evaluacion)
                #Los pesos en x tienen prioridad en la doble suma cruzada: wx1wy1+wx1+wy2+...
                matriz_Wyy[i,length(evaluacion)*(j-1)+k] = prob.μ*dX(ϕ,interior[i]-evaluacion[j],1)*dX(ϕ,interior[i]-evaluacion[k],1)
            end
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            for k in eachindex(evaluacion)
                #Los pesos en x tienen prioridad en la doble suma cruzada: wx1wy1+wx1+wy2+...
                matriz_Wxy[i,length(evaluacion)*(j-1)+k] = 2*prob.μ*dX(ϕ,interior[i]-evaluacion[j],2)*dX(ϕ,interior[i]-evaluacion[k],1)
            end
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            for k in eachindex(evaluacion)
                #Los pesos en x tienen prioridad en la doble suma cruzada: wx1wy1+wx1+wy2+...
                matriz_Wxt[i,length(evaluacion)*(j-1)+k] = -prob.ρ*prob.Cp*ϕ(interior[i]-evaluacion[j])*dX(ψ,interior[i]-evaluacion[k],1)
            end
        end
    end

    Threads.@threads for i in eachindex(interior)
        for j in eachindex(evaluacion)
            for k in eachindex(evaluacion)
                #Los pesos en x tienen prioridad en la doble suma cruzada: wx1wy1+wx1+wy2+...
                matriz_Wyt[i,length(evaluacion)*(j-1)+k] = -prob.ρ*prob.Cp*ϕ(interior[i]-evaluacion[j])*dX(ψ,interior[i]-evaluacion[k],2)
            end
        end
    end

    #Condiciones de frontera

    Threads.@threads for i in eachindex(frontera_izquierda)
        for j in eachindex(evaluacion)
            matriz_frontera_izquierda[i,j] = dX(ψ,frontera_izquierda[i]-evaluacion[j],1)
        end
    end    

    Threads.@threads for i in eachindex(frontera_derecha)
        for j in eachindex(evaluacion)
            matriz_frontera_derecha[i,j] = dX(ψ,frontera_derecha[i]-evaluacion[j],1)
        end
    end

    Threads.@threads for i in eachindex(frontera_superior)
        for j in eachindex(evaluacion)
            matriz_frontera_superior[i,j] = ψ(frontera_superior[i]-evaluacion[j])
        end
    end    

    Threads.@threads for i in eachindex(frontera_inferior)
        for j in eachindex(evaluacion)
            matriz_frontera_inferior[i,j] = ψ(frontera_inferior[i]-evaluacion[j])
        end
    end

    matriz_interior = [bloque_I1 matriz_Wxx matriz_Wxy matriz_Wyy bloque_I2 matriz_laplaciano matriz_Wxt matriz_Wyt]
    matriz_columna_fronteras = [matriz_frontera_superior; matriz_frontera_inferior; matriz_frontera_izquierda; matriz_frontera_derecha]

    A_dense = [[A bloque_H]; matriz_interior; [bloque_I matriz_columna_fronteras bloque_J]]

    A = sparse(A_dense)

    bt = zeros(puntos_interior+size(matriz_columna_fronteras,1))

    bt[puntos_interior+1:puntos_interior+puntos_frontera] .= prob.Tsup
    bt[puntos_interior+puntos_frontera+1:puntos_interior+2*puntos_frontera] .= prob.Tinf

    b = vcat(b,bt)

    return A, b, evaluacion

end