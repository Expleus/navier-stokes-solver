include("Backend/Backend.jl")

#Parametros fisicos

μ = 1 #Viscosidad (Pa)
ρ = 1 #Densidad (Kg/m^3)
k = 1
Cp = 1
Pin = 1 #Presion a la entrada (Pa)
Tsup = 0
Tinf = 0

#Parametros geometricos

Lx = 5 #Longitud del tubo (m)
Ly = 1 #Longitud de referencia vertical (m)
S(x) = 1 #Parte superior del tubo
I(x) = 0 #Parte inferior del tubo

#Parametros de la simulacion

puntos_frontera_x = 20
puntos_frontera_y = 7
puntos_interior = 130
puntos = 1 # Opciones: rng (1), flux (2)

#Definicion del problema

prob= problema(μ, ρ, k, Cp, Pin, 
               Lx, Ly, Tsup, Tinf, S, I, 
               puntos_frontera_x, puntos_frontera_y, puntos_interior)

#Generar los puntos de entrenamiento

#Segunda entrada problema: velocidad y presion desconocida (1), 

w,evaluacion,fig = problema(prob, 2, puntos) 
            
fig

adjkzhdj(y) = ((Pin/Lx)^2/192)*(Ly^4-(Ly-2*y)^4)

Plots.plot(RBF_eval.(Ref(2),LinRange(0,1,100),Ref(w.t)), LinRange(0,1,100))
Plots.plot!(adjkzhdj.(LinRange(0,1,100)), LinRange(0,1,100))

y = LinRange(0,1,100)

num  = RBF_eval.(Ref(5), y, Ref(w.t))
anal = adjkzhdj.(y)

err_pct = (norm(num - anal) / norm(anal)) * 100

fig = grafica_vpt_dimensional(prob, w, 0, 1)

fig