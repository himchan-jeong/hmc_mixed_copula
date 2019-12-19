# cop.int code:
# 1=normal, 2=clayton, 3=gumbel, 4=frank

link_normal <- function(x) atanh(x)

link_clayton <- function(x) log(x + 1)

link_gumbel <- function(x) log(x - 1)

link_frank <- function(x) x

invlink_normal <- function(x) tanh(x)

invlink_clayton <- function(x) exp(x) - 1

invlink_gumbel <- function(x) 1 + exp(x)

invlink_frank <- function(x) x

# transform copula from original support to real line support
copula_link <- function(x, cop.int) {
  switch(cop.int, link_normal(x), link_clayton(x), link_gumbel(x), 
         link_frank(x))
}

# untransform back to original support
copula_invlink <- function(x, cop.int) {
  switch(cop.int, invlink_normal(x), invlink_clayton(x), invlink_gumbel(x), 
         invlink_frank(x))
}