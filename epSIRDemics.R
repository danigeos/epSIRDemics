# MODELO SIERD 
# Modelo de cajas para propagación de epidemias. 
# Cajas:
# S susceptibles
# I infectados
# E infEctious
# R recuperados 
# D difuntos
# C infectados detectados acumulados
#
# Procesos que este modelo no tiene en cuenta:
# Compartimentalización espacial de la población.
# Comportamiento por edades de la población.  
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
    dI <- + alpha(time)*I*S - epsilon*I
    dE <- + epsilon*I - beta*E - delta*E
    dR <- + beta*E - gamma*R
    dD <- + delta*E
    dC <- + epsilon*I
    return(list(c(dS, dI, dR, dD, dC, dE)))
  })
}


#Tiempos del modelo:
día_final_modelo=365
dia_cero_lbl = "14-Feb" #día correspondiente al inicio del modelo.
dia_primer_dato <- 15 #29 de febrero 2020, día del modelo correspondiente al día del primer dato (ver abajo)
dia_primera_muerte <- dia_primer_dato + 3; dia_primera_muerte_lbl="3-Mar"
dia_conf1 <- dia_primer_dato + 13; dia_conf1_lbl="C1" #13-March (oficialmente lunes 16-M) onset of confinement and alpha1. Caída movilidad según Google. 
dia_conf2 <- dia_primer_dato + 23; dia_conf2_lbl="C2" #Necesito adelantarlo! 30-March, lunes, cierre producción and onset alpha2
dia_conf3 <- dia_primer_dato + 44; dia_conf3_lbl="C3" #13-April, lunes, reapertura producción and alpha3
dia_conf4 <- dia_primer_dato + 57; dia_conf4_lbl="C4" #26-April, domingo, paseos permitidos
dia_conf5 <- dia_primer_dato + 72; dia_conf5_lbl="C5" #11-Mayo, lunes, Fase 1 en media España

# Datos reales de Covid-19 en España 
numero_datos = 74
dias <- dia_primer_dato + seq(0, numero_datos, by = 1) 
#Casos detectados en miles, de https://ourworldindata.org/coronavirus#confirmed-covid-19-cases-by-country
casos <- c(.083, #29 febrero
           .114,.151,.200,.261,.374,.430,.589,1.204,1.639,2.140,3.004,4.231,5.753,7.753, #1-14 de marzo
           9.191,11.178,13.716,17.147,19.980,24.926,28.572,33.089,39.673,47.610,56.188,64.059,72.248,78.797, 85.195, 87.956, 95.923, #15-31 Mar
           104.118, 112.065, 119.199, 126.168,  131.646, 136.675, 141.942, 148.220, 153.222, 158.273, 163.027, 166.831, 170.099, 174.060, #1-14 April, from now on data from https://www.worldometers.info/coronavirus/country/spain/
           180.659, 184.948, 190.839, 194.416, 198.674, 200.210, 204.178, 208.389, 213.024, 219.764, 223.759, 226.629, 229.422, 232.128, 236.899, 239.340, #April 15-
           242.974, 245.589, 247.122, 248.301, 250.561, 253.682, 256.855, 260.117, 262.783, 264.663, 268.143, 269.520, 271.095 ) #May 1 - 
#Muertes en miles, https://ourworldindata.org/coronavirus#confirmed-covid-19-deaths-by-country
dead <- c(0.000, #29 febrero
          0.000,0.000,0.001,0.003,0.005,0.009, 0.019,0.028,0.036,0.050,0.084,0.121,0.136,0.288, #1-14 marzo
          0.309,0.491,0.598,0.767,1.002,1.326,1.720,2.182, 2.696,3.434,4.089,4.858,5.690,6.528, 6.803, 7.716, 8.464, #15-31 Mar
          9.387, 10.348, 11.198, 11.947, 12.641, 13.341, 14.045, 14.792, 15.447, 16.081, 16.606, 17.209, 17.756, 18.255, #1-14 April, from now on data from https://www.worldometers.info/coronavirus/country/spain/
          18.812, 19.130, 19.478, 20.043, 20.453, 20.852, 21.282, 21.717, 22.157, 22.524, 22.902, 23.190, 23.521, 23.822, 24.275, 24.543, #April 15-30
          24.824, 25.100, 25.264, 25.428, 25.613, 25.857, 26.070, 26.299, 26.478, 26.621, 26.744, 26.920, 27.104 ) #May-1
          
datos <- data.frame (dias,casos,dead)

datos

#Tiempos en que calcular las variables
time_step = 0.1 #[days]
tiempos <- seq(0, día_final_modelo, by = time_step) 

#PARÁMETROS
población <- 47000 #en miles de individuos
p <- 1 #1 #proporción inicial de susceptibles
q <- .08 #.17 proporción de Enfermos que son detectados. Vò study: 89 positives nearly none knew.
#valores iniciales [miles]:
S0 <- p*población 
I0 <- .005 #el 14-F había en España 2 casos confirmados
R0 <- 0
D0 <- 0
C0 <- I0*q
E0 <- 0


alpha0 <- 10.2e-6 #tasa de infección de Susceptibles`` [%/100/persona/dia]
alpha1 <- alpha0*.6 #tras el confinamiento
alpha2 <- alpha0*.15 #tras el cierre producción salvo esenciales
alpha3 <- alpha0*.12 #tras levantamiento del cierre producción
alpha4 <- alpha0*.13 #paseos, deporte permitido
alpha5 <- alpha0*.35  #Fase 1 en media España
alpha <- function(time){if (time<dia_conf5){if (time<dia_conf4){if (time<dia_conf3){if (time<dia_conf2){if (time<dia_conf1){alpha0} else{alpha1}} else{alpha2}}else{alpha3}}else{alpha4}}else{alpha5}} 
beta <- .26     #tasa diaria de recuperación de Enfermos [%/100/día]. Basado en los 14 días de recuperación media. Vò (Veneto): onl6 6 out of 89 remain infected after 9 days (89*.74^9=5.9). Most were assymptomatic: https://www.theguardian.com/commentisfree/2020/mar/20/eradicated-coronavirus-mass-testing-covid-19-italy-vo?CMP=share_btn_tw
gamma <- 0e-3   #tasa diaria de Recuperados que devienen Susceptibles de nuevo [%/100/día]
delta <- 2e-3   #tasa diaria de muerte de los Enfermos [%/100/día]. Letalidad actualmente en 1% (estudio de prevalencia España mayo). OLD: 1.5e-4 basado en la letalidad del 0.2% en dos semanas https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/
epsilon <- .1   #tasa diaria de Infectados que devienen Enfermos [%/100/día]. 1.5e-4 basado en el periodo de incubación en 5 días



#Salida gráfica:
día_max_plot=190
logaxis = ""  #"y" for a logarithmic vertical scale
ymax=población*p*.1001
ymin=ymax/5000
y_label=-ymax/50
escale_lbl = .7

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

daily_death_max =0
daily_death <- seq(0, día_final_modelo, by = time_step) 
for(i in 1:((día_final_modelo/time_step)-2)) {
  daily_death[i] = ((sird_values$D[i+1]-sird_values$D[i])/time_step)
  if (daily_death[i] > daily_death_max) {daily_death_max = daily_death[i]; daily_death_max_day=i*time_step}
}

#Márgenes del gráfico  
par(mar=c(4,4,2,1))
#Gráfica general
with(sird_values, {
  plot(
    time, S+I+R+E+D, lty = 1, type = "l", col = "pink",
    xlim = c(0,día_max_plot), cex.axis=0.9,
    ylim = c(ymin,ymax),
    xlab = "days since 14-F", 
    ylab = "people (x1000)",
    main = "SIRD model",
    log = logaxis
  )
  lines(time, S, col = "blue")
  lines(time, I, col = "red")
  lines(time, E, col = "orange")
  lines(time, R, col = "green")
  lines(time, D, col = "black")
  lines(time, C, col = "dark red")
  lines(time, C*q, col = "dark red", lty = 2)
  abline(a=0,b=0,col="black", lty = 3)
  abline(v=numero_datos+dia_primer_dato+1,lty=3,col="black")
  text(0,y_label,dia_cero_lbl,cex=escale_lbl,col="black")
  text(dia_primera_muerte,y_label,dia_primera_muerte_lbl,cex=escale_lbl,col="black")
  text(dia_conf1,y_label,dia_conf1_lbl,cex=escale_lbl,col="black")
  text(dia_conf2,y_label,dia_conf2_lbl,cex=escale_lbl,col="black")
  text(dia_conf3,y_label,dia_conf3_lbl,cex=escale_lbl,col="black")
  text(dia_conf4,y_label,dia_conf4_lbl,cex=escale_lbl,col="black")
  text(dia_conf5,y_label,dia_conf5_lbl,cex=escale_lbl,col="black")
  text(dia_primer_dato+33,y_label,"pico",cex=escale_lbl,col="black") #2-Apr
  text(dia_primer_dato+65,y_label,"4-M",cex=escale_lbl,col="black") #4-Mayo, fecha media estudio prevalencia
  text(daily_death_max_day,y_label,"Peak",cex=escale_lbl,col="black") #pico del modelo
  text(numero_datos+dia_primer_dato+1,y_label,"now",cex=escale_lbl,col="black")
  with(datos,points(dias, casos, pch = 20, cex=0.5,col = "dark red"))
  with(datos,points(dias, dead, pch = 20, cex=0.5,col = "black"))
})

#Leyenda:
legend("topleft", c("S-Susceptibles","I-Infected","Infected (cumul.)","E-infEctious","R-Recovered","D-Dead", "Total"),
       col = c( "blue", "red", "dark red", "orange", "green", "black", "pink"), cex=0.6, lty = 1, bty = "n", pch=20)

cat(  "S Min. Susceptibles: ", min(sird_values$S), "k pp")
cat("\nI Max. Infected:   ", max(sird_values$I), "k pp")
cat("\nE Max. Enfermos:     ", max(sird_values$E), "k pp")
cat("\nR Final Recovered:   ", max(sird_values$R), "k pp")
cat("\nR Cumul cases:       ", max(sird_values$C), "k pp")
cat("\nD Final Deaths:      ", max(sird_values$D), "k pp")

cat("\n\nRatio infectados 5 mayo: ", sird_values$C[65/time_step]/población*100, "%")
cat("\nDaily death peak: ", daily_death_max, "k pp/day on day", daily_death_max_day)

#REFS:
# Kermack and McKendrick [Proc. R. Soc. A, 115, 772 (1927)]
# https://reason.com/2020/04/03/what-we-should-have-learned-from-icelands-response-to-covid-19/
# https://www.statista.com/statistics/1106855/tested-and-confirmed-coronavirus-cases-in-iceland/
# https://www.cebm.net/covid-19/global-covid-19-case-fatality-rates/
# https://cleantechnica.com/2020/03/21/iceland-is-doing-science-50-of-people-with-covid-19-not-showing-symptoms-50-have-very-moderate-cold-symptoms/
# https://www.government.is/news/article/2020/03/15/Large-scale-testing-of-general-population-in-Iceland-underway/
# https://www.theguardian.com/commentisfree/2020/mar/20/eradicated-coronavirus-mass-testing-covid-19-italy-vo?CMP=share_btn_tw
# I=0.15 https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-13-europe-npi-impact/
# I=0.15 https://catalunyaplural.cat/ca/la-prudencia-en-epidemiologia/
# Death rates from Covid reliable: https://elpais.com/sociedad/2020/04/09/actualidad/1586435286_092135.html
# Estudio seroprevalencia España Mayo: https://www.rtve.es/noticias/20200513/5-poblacion-espanola-se-contagiado-del-virus-segun-primera-oleada-del-estudio-seroprevalencia/2013992.shtml