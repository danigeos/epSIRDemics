# MODELO SIRD (modificado del de JMO)
# S (susceptibles), 
# I (infecciosos), 
# R (recuperados) 
# D (difuntos) 
# C (casos detectados acumulados)
#
# tiempo en días y las demás variables en millares de personas.

library(deSolve) # para usar la funci?n "ode" 

# Paso 1: Sistemas de ecuaciones diferenciales en R
# Ojo porque este método ignora el retraso o latencia entre  la infección y la infeciosidad.
# Podemos resolver numericamente estos sistemas en R gracias a la función ode() del package deSolve
# Instalar con:
##  install.packages("deSolve")
# Se usa con la función with() function abajo:
sird_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- -alpha(time) * I * S + gamma * R
    dI <-  alpha(time) * I * S - beta * I - delta*I
    dR <-  beta * I
    dD <-  delta * I
    dC <-  alpha(time) * I * S * q
    return(list(c(dS, dI, dR, dD, dC)))
  })
}

#PARÁMETROS
poblacion <- 47000 #en miles de individuos
p <- 0.3 #proporción inicial de susceptibles
q <- .016  #proporción de infecciosos que son detectados (a partir del estudio islandés)
S0 <- p*poblacion #susceptibles iniciales [miles]
I0 <- .45  #infecciosos iniciales [miles]
R0 <- 0   #recuperados iniciales [miles]
D0 <- 0   #difuntos inciales [miles]
C0 <- I0*q  #casos iniciales

time_first_data <- 15 #día equivalente al día del primer dato (ver abajo)
time_confinement <- time_first_data + 124 #onset of alpha2, days after initial time=0
alpha1 <- 22.8e-6 #tasa capacidad de infeccion-transmision (/persona^2/dia)
alpha2 <- 1.3e-6 #alpha1/1.5
alpha <- function(time){if (time<time_confinement){alpha1} else{alpha2}} 
beta <- .07      #tasa diaria de recuperación de infectados [/día]. Basado en los 14 días de recuperación media.
gamma <- 0e-3    #tasa diaria de recuperados que devienen susceptibles de nuevo (/día) 
delta <- 1.5e-4  #tasa diaria de muerte de un infeccioso [/día]. Basado en la letalidad del 0.2% del estudio islandés.

parameters_values <- c(
  beta ,
  delta    
)

#valores iniciales de las variables en un vector:
initial_values <- c(S = S0, I = I0, R = R0, D = D0, C = C0)

#Tiempos en que calcular las variables
time_values <- seq(0, 365, by = .1) # los 30 primeros días,  paso intermedios de 0.1


#Resolver numericamente el sistema:
sird_values <- ode(
  y = initial_values,
  times = time_values,
  func = sird_equations,
  parms = parameters_values 
)


# Podemos ver los valores calculados:
sird_values
# Podremos usarlos posteriormente. 
# De momento coloqu?moslos en un data frame
# despu?s los podremos usar con la funci?n with()
sird_values <- as.data.frame(sird_values)
  
# Datos reales en España de coronavirus 
# primer dato: 29 de febrero 2020)

dia <- seq(0, 30, by = 1)
dia <- dia + time_first_data 

#Casos detectados en miles, de https://ourworldindata.org/coronavirus#confirmed-covid-19-cases-by-country
casos <- c(.083,.114,.151,.200,.261,.374,.430,.589,1.204,1.639,2.140,3.004,4.231,5.753,7.753,9.191,11.178,13.716,17.147,19.980,24.926,28.572,33.089,39.673,47.610,56.188,64.059,72.248,78.797,
           85.195,94.417) 

#Muertes en miles, https://ourworldindata.org/coronavirus#confirmed-covid-19-deaths-by-country
#Primera muerte el 3 de marzo.
dead <- c(0.000,0.000,0.000,0.001,0.003,0.005,0.009,
          0.019,0.028,0.036,0.050,0.084,0.121,0.136,0.288,
          0.309,0.491,0.598,0.767,1.002,1.326,1.720,2.182,
          2.696,3.434,4.089,4.858,5.690,6.528, 
          7.340, 8.189)
datos <- data.frame (dia,casos,dead)    
 
#modifiquemos los márgenes para permitir los dibujos (según pantalla)  
par(mar=c(4,4,2,1))
#gráfica general
with(sird_values, {
  plot(
    time, C, lty = 3, type = "l", col = "dark red",
    xlim = c(0,65), cex.axis=0.9,
    ylim = c(-10,122),
    xlab = "tiempo (días)", 
    ylab = "personas (x1000)",
    main = "Evolución epidemia - modelo SIRD"
  )
  lines(time, S, col = "blue")
  lines(time, I, col = "red")
  lines(time, R, col = "green")
  lines(time, D, col = "black")
  #abline(a=50,b=0,col="blue")
  #abline(v=33.5,col="blue")
  text(time_confinement,-11,"conf",cex=0.6,col="black")
  text(0,-10,"14 feb",cex=0.6,col="black")
  #text(80,32,"aprox. 45,5 mil",cex=0.6,col="black")
  with(datos,points(dia, casos, pch = 20, cex=0.5,col = "dark red"))
  with(datos,points(dia, dead, pch = 20, cex=0.5,col = "black"))
})

#Leyenda:
legend("right", c("Susceptibles", "Infecciosos activos", "Positivos acumulados", "Recuperados","Difuntos"),
  col = c( "blue", "red", "dark red", "green", "black"), cex=0.6, lty = 1, bty = "n")
