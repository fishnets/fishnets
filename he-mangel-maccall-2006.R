#' A prior for steepness based on natural mortality and recruitment variability 
#' [He et al 2006](http://fishbull.noaa.gov/1043/he.pdf)

HeMangelMaccall2006 <- function(){

  self <- object('SteepnessHe2006')

  #' Parameters of the prior from Table 1 of He et al (2006).
  self$table <- data.frame(
    # Create an iteration over M, parameter and sigmar representing order of values
    expand.grid(
      m = c(0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7),
      # Note that in the table caption the authors say that the parameters are in the order omega1,omega2,omega3
      # However, if one assumes this order you don't get sensible results. If one assumes the order omega1,omega3,omega2 do get sensible priors
      par = c('o1','o3','o2'),
      sigmar = seq(0.2,1.6,0.2)
    ),
    # Note that the value with "#Iterpolated" next to it was missing in the original table and is interpolated from the omega1 values 'above' and 'below' it in the table
    # Note that the triptych marked with "#As for <-" was marked as "not fitted" in the table and I used the values for Sigma=1.6, M=0.6
    value = c(
    -1.0000000,-1.0000000,-1.0000000,-1.0000000,-1.0000000,1.00000e+0,1.00000e+0,1.00000e+0,1.00000e+0,1.00000e+0,
    -1.0000000,-1.0000000,-1.0000000,-1.0000000,-1.0000000,6.69600e+2,4.43000e+2,3.06130e+2,2.82660e+2,2.78160e+2,
    -1.0000000,-1.0000000,-1.0000000,-1.0000000,-1.0000000,9.28900e-1,8.04600e-1,6.66200e-1,5.01500e-1,3.72400e-1,
    
    -1.0000000,9.99400e-1,1.00000e+0,1.00000e+0,1.00000e+0,1.00000e+0,9.99900e-1,9.99900e-1,9.99900e-1,9.99700e-1,
    -1.0000000,1.39900e+2,2.48400e+2,2.23600e+2,1.96430e+2,1.99300e+2,1.82220e+2,1.57820e+2,1.48630e+2,1.45480e+2,
    -1.0000000,9.80000e-1,7.97300e-1,5.58500e-1,3.73900e-1,2.14600e-1,8.61900e-2,4.45400e-2,2.26400e-2,1.15400e-2,
    
    -1.0000000,1.00000e+0,1.00000e+0,1.00000e+0,9.99800e-1,9.99500e-1,9.99200e-1,9.98800e-1,9.98500e-1,9.97000e-1,#Iterpolated
    -1.0000000,1.49700e+2,1.33600e+2,1.16740e+2,1.09860e+2,1.05320e+2,8.95500e+1,7.87900e+1,6.92900e+1,6.14800e+1,
    -1.0000000,7.55200e-1,3.89500e-1,1.79100e-1,8.06900e-2,3.80500e-2,1.54300e-2,7.60200e-3,5.07400e-3,3.86100e-3,
    
    1.00000e+0,1.00000e+0,1.00000e+0,9.99700e-1,9.99200e-1,9.99000e-1,9.98200e-1,9.97200e-1,9.95900e-1,9.94300e-1,
    8.30000e+1,9.31300e+1,8.42500e+1,7.62900e+1,6.87300e+1,6.55300e+1,4.88400e+1,4.05100e+1,3.40300e+1,2.91200e+1,
    9.25900e-1,4.31700e-1,1.37600e-1,4.54000e-2,2.02200e-2,1.19700e-2,5.79600e-3,3.74500e-3,3.04100e-3,2.72200e-3,
    
    1.00000e+0,1.00000e+0,9.99500e-1,9.98700e-1,9.98000e-1,9.97000e-1,9.94400e-1,9.90200e-1,9.83000e-1,9.70900e-1,
    7.53000e+1,6.44200e+1,5.75200e+1,4.87700e+1,4.14400e+1,3.55700e+1,2.71500e+1,2.15800e+1,1.77100e+1,1.49700e+1,
    9.52600e-1,1.92900e-1,4.46400e-2,1.70000e-1,9.24000e-3,6.19900e-3,3.90600e-3,3.07600e-3,2.68000e-3,2.40200e-3,
    
    1.00000e+0,9.99500e-1,9.98300e-1,9.96700e-1,9.94600e-1,9.91500e-1,9.80100e-1,9.54300e-1,9.01100e-1,8.08500e-1,
    6.93600e+1,4.81600e+1,3.87100e+1,3.09000e+1,2.53800e+1,2.12600e+1,1.57500e+1,1.26400e+1,1.10100e+1,1.03400e+1,
    5.92000e-1,7.18200e-2,1.84400e-2,9.03600e-3,5.81900e-3,4.42000e-3,3.13400e-3,2.24300e-3,1.32000e-3,5.88200e-4,
    
    1.00000e+0,9.98300e-1,9.95600e-1,9.91300e-1,9.83500e-1,9.69700e-1,9.02500e-1,7.42800e-1,4.97300e-1,2.42800e-1,
    4.77900e+1,3.50600e+1,2.58100e+1,1.98700e+1,1.59400e+1,1.33900e+1,1.07400e+1,1.02500e+1,1.08800e+1,1.26200e+1,
    2.02400e-1,2.50400e-2,9.16300e-3,5.64700e-3,4.01000e-3,3.08700e-3,1.37300e-3,3.18300e-4,3.68300e-5,1.60000e-6,
    
    9.99200e-1,9.95500e-1,9.88300e-1,9.71600e-1,9.33800e-1,8.58500e-1,5.46300e-1,2.02700e-1,2.79000e-2,2.79000e-2,#As for <-
    3.90900e+1,2.50300e+1,1.72200e+1,1.32200e+1,1.10800e+1,1.02500e+1,1.07500e+1,1.29000e+1,1.68500e+1,1.68500e+1,#As for <-
    3.00100e-2,8.90000e-3,4.92300e-3,3.34000e-3,1.97100e-3,8.68600e-4,5.31500e-5,7.90000e-7,1.00000e-8,1.00000e-8 #As for <-
  ))
  self$table <- cast(self$table,sigmar+m~par)
  
  #' Get parameter values
  self$parameters <- function(m,sigmar){
    m_ <- self$table$m[which.min(abs(self$table$m-m))]
    sigmar_ <- self$table$sigmar[which.min(abs(self$table$sigmar-sigmar))]
    subset(self$table,m==m_ & sigmar==sigmar_)
  }  
  
  #' Get the *relative* probability of a given steepness, given the parameters of the prior
  self$relative <- function(steepness, parameters){
    with(parameters,{
      o1/(1+(o1-o2)/o2*exp(-o3*(steepness-0.2)/0.8))
    })
  }
  
  #' Get probability densities for a given m and sigmar
  self$densities <- function(m,sigmar){
    parameters <- self$parameters(m,sigmar)
    relatives <- self$relative(seq(0.2,1,0.01),parameters)
    relatives/sum(relatives)
  }
  
  self$sample <- random(m,sigmar) {
  }

  self
}





