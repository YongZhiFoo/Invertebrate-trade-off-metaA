# effect size functions for lnRR and SMD

# lnRR
# this includes the second order correction - see 
lnRR_arcsin <- function(m1, m2, n1, n2) { # m1 and m2 should be between 0 and 1
  # arcsine transformation
  asin_trans <- function(p) { asin(sqrt(p)) }
  # varaince for arcsine 
  var1 <- 1/(8)
  var2 <- 1/(8)
  # lnRR - with 2nd order correction (see Lajeunesse 2015)
  lnrr <- log(asin_trans(m1)/asin_trans(m2)) +
    0.5 * ((var1 / (n1 * asin_trans(m1)^2)) - (var2 / (n2 * asin_trans(m2)^2)))	
  var <- var1 / (n1 * asin_trans(m1)^2) + var1^2 / (2 * n1^2 * asin_trans(m1)^4)  +
    var2 / (n2 * asin_trans(m2)^2) + var2^2 / (2 * n2^2 * asin_trans(m2)^4)
  effect <- data.frame(yi = lnrr ,vi = var)
  effect
}

# SMD

SMD_arcsin <- function(m1, m2, n1, n2) { # m1 and m2 should be between 0 and 1
  # arcsine transformation
  asin_trans <- function(p) { asin(sqrt(p)) }
  
  # variance for arcsine 
  var1 <- 1/(8)
  var2 <- 1/(8)
  # SMD
  s_pooled <- sqrt(((n1 - 1)*var1 + (n2- 1)*var2) / (n1 + n2 - 2))	
  d <- (asin_trans(m1) - asin_trans(m2)) / s_pooled
  # with small sample correction
  j <-  (1 - (3 / (4 * (n1 + n2 - 2) - 1)))
  SMD <- d * j
  var <- ((n1 + n2) / (n1 * n2)) + ((SMD^2) / (2 * (n1 + n2 - 2)))*j^2
  effect <- data.frame(yi = SMD , vi = var)
  effect
}

