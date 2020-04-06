# MODELO SIRD (modificado del de JMO)
# S (susceptibles), 
# I (infecciosos), 
# R (recuperados) 
# D (difuntos) 
# C (casos detectados acumulados)
#
# tiempo en días y las demás variables en millares de personas.

library(deSolve) # para usar la función "ode" http://biom300.weebly.com/simulating-disease-models-examples.html

# Paso 1: Sistemas de ecuaciones diferenciales en R
# Ojo porque este método ignora el retraso o latencia entre  la infección y la infeciosidad.
# Podemos resolver numericamente estos sistemas en R gracias a la función ode() del package deSolve
# Instalar con:
##  install.packages("deSolve")
# Se usa con la función with() function abajo:
sird_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- -alpha(time) * I * S + gamma * R
    dI <- +alpha(time) * I * S - beta * I - delta*I
    dR <- +beta * I - gamma * R
    dD <- +delta * I
    dC <- +alpha(time) * I * S
    return(list(c(dS, dI, dR, dD, dC)))
  })
}

día_final_modelo=365
dia_cero_lbl = "14-Feb" #día correspondiente al inicio del modelo.
dia_primer_dato <- 15 #29 de febrero 2020, día del modelo correspondiente al día del primer dato (ver abajo)
dia_primera_muerte <- dia_primer_dato + 3
dia_primera_muerte_lbl="3-Mar"
dia_confinamiento <- dia_primer_dato + 20 #13-April onset of confinement and alpha2
dia_confinamiento_lbl="Conf."
# Datos reales de Covid-19 en España 
dias <- dia_primer_dato + seq(0, 34, by = 1) 
#Casos detectados en miles, de https://ourworldindata.org/coronavirus#confirmed-covid-19-cases-by-country
casos <- c(.083, #29 febrero
           .114,.151,.200,.261,.374,.430,.589,1.204,1.639,2.140,3.004,4.231,5.753,7.753, #hasta el 14 de marzo
           9.191,11.178,13.716,17.147,19.980,24.926,28.572,33.089,39.673,47.610,56.188,64.059,72.248,78.797, 85.195,
           94.417, 102.136, 110.238, 117.710, 124.136) #from April 1st
#Muertes en miles, https://ourworldindata.org/coronavirus#confirmed-covid-19-deaths-by-country
dead <- c(0.000, #29 febrero
          0.000,0.000,0.001,0.003,0.005,0.009, 0.019,0.028,0.036,0.050,0.084,0.121,0.136,0.288,
          0.309,0.491,0.598,0.767,1.002,1.326,1.720,2.182, 2.696,3.434,4.089,4.858,5.690,6.528, 7.340, 
          8.189, 9.053, 10.003, 10.935, 11.744) #from April 1st
datos <- data.frame (dias,casos,dead)

#Tiempos en que calcular las variables
tiempos <- seq(0, día_final_modelo, by = .1) #paso de 0.1 días


#PARÁMETROS
población <- 47000 #en miles de individuos
p <- 0.5 #proporción inicial de susceptibles
q <- .019 #.17 proporción de infecciosos que son detectados. Vò study: 89 positives nearly none knew.
S0 <- p*población #susceptibles iniciales [miles]
I0 <- 2.0  #infecciosos iniciales [miles]
R0 <- 0   #recuperados iniciales [miles]
D0 <- 0   #difuntos inciales [miles]
C0 <- I0*q  #casos iniciales


alpha1 <- 9.9e-6 #tasa capacidad de infeccion-transmision`` [/persona/dia]
alpha2 <- 8e-6 #tras el confinamiento
alpha <- function(time){if (time<dia_confinamiento){alpha1} else{alpha2}} 
beta <- .05      #tasa diaria de recuperación de infectados [/día]. Basado en los 14 días de recuperación media. Vò (Veneto): onl6 6 out of 89 remain infected after 9 days. Most were assymptomatic: https://www.theguardian.com/commentisfree/2020/mar/20/eradicated-coronavirus-mass-testing-covid-19-italy-vo?CMP=share_btn_tw
gamma <- 0e-3    #tasa diaria de recuperados que devienen susceptibles de nuevo [/día]
delta <- 3e-4  #tasa diaria de muerte de un infeccioso [/día]. 1.5e-4 basado en la letalidad del 0.2% https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/

parameters_values <- c(
  beta,
  delta    
)
#valores iniciales de las variables en un vector:
initial_values <- c(S = S0, I = I0, R = R0, D = D0, C = C0)



#Resolver numericamente el sistema:
sird_values <- ode(
  y = initial_values,
  times = tiempos,
  func = sird_equations,
  parms = parameters_values 
)


# Para ver los valores calculados:
sird_values
# Los ponemos en un data frame para usarlos con la función with()
sird_values <- as.data.frame(sird_values)
  

 
#Márgenes del gráfico  
par(mar=c(4,4,2,1))
#Gráfica general
with(sird_values, {
  plot(
    time, C*q, lty = 3, type = "l", col = "dark red",
    xlim = c(0,60), cex.axis=0.9,
    ylim = c(-10,140),
    xlab = "tiempo (días)", 
    ylab = "personas (x1000)",
    main = "Evolución epidemia - modelo SIRD"
  )
  lines(time, S, col = "blue")
  lines(time, I, col = "red")
  lines(time, R, col = "green")
  lines(time, D, col = "black")
  abline(a=0,b=0,lty = 3, col="black")
  #abline(v=33.5,col="blue")
  text(0,-10,dia_cero_lbl,cex=0.6,col="black")
  text(dia_primera_muerte,-10,dia_primera_muerte_lbl,cex=0.6,col="black")
  text(dia_confinamiento,-10,dia_confinamiento_lbl,cex=0.6,col="black")
  with(datos,points(dias, casos, pch = 20, cex=0.5,col = "dark red"))
  with(datos,points(dias, dead, pch = 20, cex=0.5,col = "black"))
})

#Leyenda:
legend("right", c("Susceptibles", "Infecciosos activos", "Positivos acumulados", "Recuperados","Difuntos"),
  col = c( "blue", "red", "dark red", "green", "black"), cex=0.6, lty = 1, bty = "n", pch=20)

cat("Final deaths:", max(sird_values$D), "k pp")


#REFS:
# https://reason.com/2020/04/03/what-we-should-have-learned-from-icelands-response-to-covid-19/
# https://www.statista.com/statistics/1106855/tested-and-confirmed-coronavirus-cases-in-iceland/
# https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/
# https://cleantechnica.com/2020/03/21/iceland-is-doing-science-50-of-people-with-covid-19-not-showing-symptoms-50-have-very-moderate-cold-symptoms/
# https://www.government.is/news/article/2020/03/15/Large-scale-testing-of-general-population-in-Iceland-underway/
# https://www.theguardian.com/commentisfree/2020/mar/20/eradicated-coronavirus-mass-testing-covid-19-italy-vo?CMP=share_btn_tw


