# =====================
# Inclusión de módulos principales
# =====================

include("Backend/Backend.jl") # Importa toda la funcionalidad del backend


# =====================
# Parámetros físicos del problema
# =====================

μ = 1   # Viscosidad dinámica del fluido (Pa·s)
ρ = 1   # Densidad del fluido (kg/m^3)
k = 1   # Conductividad térmica (W/m·K)
Cp = 1  # Calor específico a presión constante (J/kg·K)
Pin = 1 # Presión en la entrada del tubo (Pa)
Tsup = 0 # Temperatura en la frontera superior (°C)
Tinf = 0 # Temperatura en la frontera inferior (°C)


# =====================
# Parámetros geométricos del dominio
# =====================

Lx = 5 # Longitud del tubo (m)
Ly = 1 # Altura de referencia vertical (m)
S(x) = log(0.05*(x+2))+4  # Función que define la frontera superior del tubo en x
I(x) = (sin(x))^2 # Función que define la frontera inferior del tubo en x


# =====================
# Parámetros de la simulación numérica
# ====================

puntos_frontera_x = 20 # Número de puntos en la frontera horizontal (entrada/salida)
puntos_frontera_y = 7  # Número de puntos en la frontera vertical (superior/inferior)
puntos_interior = 140  # Número de puntos interiores para el dominio
puntos = 1 # Tipo de generación de puntos: 1=rng (aleatorio), 2=flux (basado en la geometria)


# =====================
# Definición de la estructura del problema físico
# =====================

prob = problema(
    μ, ρ, k, Cp, Pin, 
    Lx, Ly, Tsup, Tinf, S, I, 
    puntos_frontera_x, puntos_frontera_y, puntos_interior
)


# =====================
# Ejecución de la simulación principal
# =====================


# Segunda y tercera entrada del problema:
#   - x=1: velocidad y presión desconocidas
#   - x=2: velocidad, presión y temperatura desconocidas
#   - puntos: tipo de generación de puntos

w, evaluacion, fig = problema(prob, 2, puntos)

# Muestra la figura generada por la simulación
fig