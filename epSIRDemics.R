# MODELO SIERD (modificado del de JMO)
# Modelo de cajas para propagación de epidemias. 
# Cajas:
# S susceptibles
# I infecciosos 
# E Enfermos
# R recuperados 
# D difuntos
# C casos detectados acumulados
#
# Procesos que este modelo no tiene en cuenta:
# Compartimentalización espacial de la población.
# Distribución por edades de la población.
# Retraso o latencia entre la infección y la infeciosidad.
# 
# Unidades: tiempo en días y las demás variables en millares de personas.

library(deSolve) #Para usar la función "ode" http://biom300.weebly.com/simulating-disease-models-examples.html

# Sistemas de ecuaciones diferenciales en R
# Podemos resolver numericamente estos sistemas en R gracias a la función ode() del package deSolve
# Instalar con:
##  install.packages("deSolve")
# Se usa con la función with() function abajo:
sird_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- - alpha(time)*I*S + gamma*R
    dI <- + alpha(time)*I*S - epsilon*(I-E)
    dE <- + epsilon*(I-E) - beta*E - delta*E
    dR <- + beta*E - gamma*R
    dD <- + delta*E
    dC <- + alpha(time)*E*S
    return(list(c(dS, dI, dR, dD, dC, dE)))
  })
}

#Salida gráfica:
ymin=-2
ymax=020

#Tiempos del modelo:
día_final_modelo=55
dia_cero_lbl = "14-Feb" #día correspondiente al inicio del modelo.
dia_primer_dato <- 15 #29 de febrero 2020, día del modelo correspondiente al día del primer dato (ver abajo)
dia_primera_muerte <- dia_primer_dato + 3
dia_primera_muerte_lbl="3-Mar"
dia_confinamiento <- dia_primer_dato + 13 #13-April onset of confinement and alpha2
dia_confinamiento_lbl="Conf."
dia_confinamiento2 <- dia_primer_dato + 28 #Lunes 30-April onset of confinement and alpha2
dia_confinamiento2_lbl="Conf2"
# Datos reales de Covid-19 en España 
numero_datos = 39
dias <- dia_primer_dato + seq(0, numero_datos, by = 1) 
#Casos detectados en miles, de https://ourworldindata.org/coronavirus#confirmed-covid-19-cases-by-country
casos <- c(.083, #29 febrero
           .114,.151,.200,.261,.374,.430,.589,1.204,1.639,2.140,3.004,4.231,5.753,7.753, #hasta el 14 de marzo
           9.191,11.178,13.716,17.147,19.980,24.926,28.572,33.089,39.673,47.610,56.188,64.059,72.248,78.797, 85.195,
           94.417, 102.136, 110.238, 117.710, 124.136,  #from April 1st
          131.646, 136.675, 141.942, 148.220, 153.222 ) #from April 5th, data from https://www.worldometers.info/coronavirus/country/spain/
#Muertes en miles, https://ourworldindata.org/coronavirus#confirmed-covid-19-deaths-by-country
dead <- c(0.000, #29 febrero
          0.000,0.000,0.001,0.003,0.005,0.009, 0.019,0.028,0.036,0.050,0.084,0.121,0.136,0.288,0.309,0.491,0.598,0.767,1.002,1.326,1.720,2.182, 2.696,3.434,4.089,4.858,5.690,6.528, 7.340, 
          8.189, 9.053, 10.003, 10.935, 11.744, #from April 1st
          12.641, 13.341, 14.045, 14.792, 15.447 ) #from April 5th, data from https://www.worldometers.info/coronavirus/country/spain/
datos <- data.frame (dias,casos,dead)

#Tiempos en que calcular las variables
tiempos <- seq(0, día_final_modelo, by = .1) #paso de 0.1 días


#PARÁMETROS
población <- 47000 #en miles de individuos
p <- 1 #1 #proporción inicial de susceptibles
q <- .007 #.17 proporción de infecciosos que son detectados. Vò study: 89 positives nearly none knew.
#valores iniciales [miles]:
S0 <- p*población 
I0 <- 1.5 #el 14-F había en España 2 casos confirmados
R0 <- 0
D0 <- 0
C0 <- I0*q
E0 <- 0


alpha1 <- 8.3e-6 #tasa capacidad de infeccion de Susceptibles`` [%/100/persona/dia]
alpha2 <- alpha1/1  #0e-6 #tras el confinamiento
alpha3 <- alpha2/7 #0e-6 #tras el cierre salvo esenciales
alpha <- function(time){if (time<dia_confinamiento2){if (time<dia_confinamiento){alpha1} else{alpha2}} else{alpha3}} 
beta <- .2      #tasa diaria de recuperación de un Enfermo [%/100/día]. Basado en los 14 días de recuperación media. Vò (Veneto): onl6 6 out of 89 remain infected after 9 days. Most were assymptomatic: https://www.theguardian.com/commentisfree/2020/mar/20/eradicated-coronavirus-mass-testing-covid-19-italy-vo?CMP=share_btn_tw
gamma <- 0e-3    #tasa diaria de Recuperados que devienen Susceptibles de nuevo [%/100/día]
delta <- 1.5e-4  #tasa diaria de muerte de un Enfermo [%/100/día]. 1.5e-4 basado en la letalidad del 0.2% en dos semanas https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/
epsilon <- .2   #tasa diaria de Infectados que devienen Enfermos [%/100/día]. 1.5e-4 basado en el periodo de incubación en 5 días


parameters_values <- c(
  beta,
  delta    
)
#valores iniciales de las variables en un vector:
initial_values <- c(S = S0, I = I0, R = R0, D = D0, C = C0, E=E0)



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
    time, C*q, lty = 1, type = "l", col = "dark red",
    xlim = c(0,día_final_modelo), cex.axis=0.9,
    ylim = c(ymin,ymax),
    xlab = "tiempo (días)", 
    ylab = "personas (x1000)",
    main = "Evolución epidemia - modelo SIRD"
  )
  lines(time, S, col = "blue")
  lines(time, I, col = "red")
  lines(time, E, col = "orange")
  lines(time, R, col = "green")
  lines(time, D, col = "black")
  lines(time, S+I+R+E+D, col = "pink")
  abline(a=0,b=0,lty = 3, col="black")
  abline(v=numero_datos+dia_primer_dato,lty=3,col="black")
  text(0,ymin,dia_cero_lbl,cex=0.6,col="black")
  text(dia_primera_muerte,ymin,dia_primera_muerte_lbl,cex=0.6,col="black")
  text(dia_confinamiento ,ymin,dia_confinamiento_lbl, cex=0.6,col="black")
  text(dia_confinamiento2,ymin,dia_confinamiento2_lbl,cex=0.6,col="black")
  text(numero_datos+dia_primer_dato,ymin,"ayer",cex=0.6,col="black")
  with(datos,points(dias, casos, pch = 20, cex=0.5,col = "dark red"))
  with(datos,points(dias, dead, pch = 20, cex=0.5,col = "black"))
})

#Leyenda:
legend("topleft", c("Susceptibles","Infecciosos","Positivos acumulados","Enfermos","Recuperados","Difuntos", "Pobllación total"),
  col = c( "blue", "red", "dark red", "orange", "green", "black", "pink"), cex=0.6, lty = 1, bty = "n", pch=20)

cat(  "S Min. Susceptibles:", min(sird_values$S), "k pp")
cat("\nI Max. Infecciosos: ", max(sird_values$I), "k pp")
cat("\nE Max. Enfermos:    ", max(sird_values$E), "k pp")
cat("\nR Final Recuperados:", max(sird_values$R), "k pp")
cat("\nD Final Deaths:     ", max(sird_values$D), "k pp")


#REFS:
# https://reason.com/2020/04/03/what-we-should-have-learned-from-icelands-response-to-covid-19/
# https://www.statista.com/statistics/1106855/tested-and-confirmed-coronavirus-cases-in-iceland/
# https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/
# https://cleantechnica.com/2020/03/21/iceland-is-doing-science-50-of-people-with-covid-19-not-showing-symptoms-50-have-very-moderate-cold-symptoms/
# https://www.government.is/news/article/2020/03/15/Large-scale-testing-of-general-population-in-Iceland-underway/
# https://www.theguardian.com/commentisfree/2020/mar/20/eradicated-coronavirus-mass-testing-covid-19-italy-vo?CMP=share_btn_tw
# I=0.15 https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-13-europe-npi-impact/
# I=0.15 https://catalunyaplural.cat/ca/la-prudencia-en-epidemiologia/
