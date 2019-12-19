# Decode combination code to get simulation settings
# (N)(tau)(cop.int)
# N: 1=500, 2=1000
# tau: tau/10 .1, .3, .6
# cop.int: 1=normal, 2=clayton, 3=gumbel, 4=frank

decodeCombo <- function(combostr) {
  N <- switch(as.integer(substr(combostr, 1, 1)), 500, 1000) 
  tau <- as.integer(substr(combostr, 2, 2)) / 10
  cop.int <- as.integer(substr(combostr, 3, 3))
  return(list(N = N, tau = tau, cop.int = cop.int))
}

# Decode combination code to get simulation settings: Model selection simulation
# (cop.gen.int)(tau)
# tau: tau/10 .1, .3, .6
# cop.int: 1=normal, 2=clayton, 3=gumbel, 4=frank  Generation, not estimation
decodeCombo2 <- function(combostr) {
  cop.int <- as.integer(substr(combostr, 1, 1))
  tau <- as.integer(substr(combostr, 2, 2)) / 10
  return(list(tau = tau, cop.int = cop.int))
}