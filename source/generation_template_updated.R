## This file prepares a template of all possible sequences of
## viral introductions and subsequent transmission events within household
## for various final sizes j (aka total # infected in a HH) and household sizes k

load.generation.comb <- function(){
  library(gtools)
  library(dplyr)
  
  ## a function that generates the one combination from 0/k households
  get_m0k <- function(k){
    tmp <- t(matrix(rep(0, 3*k), ncol = k))
    m0k <- cbind(9999, tmp, 0, k, 1:k, 0, 1)
    colnames(m0k) <- c("gen", "inf_out",	"inf_in_n", "avoid_in_n", "hh_inf_n", "hh_size", "ind_id", "pos", "set")
    return(m0k)
  }
  
  ## a function that generates all combinations in j/k households based on 
  ## combinations from j/j households
  
  get_mjk <- function(j, k, mjj){ # j infected in k person HH
    # "inf_out,	inf_in_n, avoid_in_n" of those NOT infected
    tmp <- t(matrix(rep(c(0,0,j),(k-j)*max(mjj[,"set"])),
                    ncol = (k-j)*max(mjj[,"set"])))
    
    # all columns "gen,  inf_out,	inf_in_n, avoid_in_n, hh_inf_n, hh_size, ind_id, pos, set" 
    # of those NOT infected
    tmp_not_infected <- cbind(9999, tmp, j, k, (j+1):k, 0, rep(1:max(mjj[,"set"]), each=k-j)) 
    
    # combine with that from those infected
    mjk <- rbind(tmp_not_infected, mjj)
    colnames(mjk) <- c("gen", "inf_out",	"inf_in_n", "avoid_in_n", "hh_inf_n", "hh_size", "ind_id", "pos", "set")
    mjk[,"hh_size"] <- k
    return(mjk)
  }
  
  
  # to_assign = 1, output shows which generation each person is infected at
  
  assign_gens <- function (to_assign, assignment, n_inf) {
    if (to_assign>length(assignment)) {
      #make sure it is legal
      assign_legal <- TRUE
      for(i in 0:max(assignment)) {
        assign_legal <- assign_legal * (i%in%assignment)
      }
      if (assign_legal) {
        return(assignment)
      } else{
        return(NULL)
      }
    }
    rc<-NULL
    ##iterate over assignments
    for(i in 0:(n_inf-1)) {
      assignment[to_assign] <- i
      rc <- rbind(rc,assign_gens(to_assign+1,assignment, n_inf))
    }
    return(rc)
  }
  
  # create input to stan based on output of assign_gens
  create_stan_input_kk_hh <- function(dat){
    gen <- as.vector(t(dat))
    inf_out <- as.vector(t(dat) == 0)
    inf_in_n <- c() # total infected the generation before hh member "i" was infected
    for(i in 1:nrow(dat)){
      for (j in 1:ncol(dat)){
        inf_in_n <- c(inf_in_n, sum(dat[i,j] == dat[i,]+1))
      }
    }
    avoid_in_n <- c() # total number of infected household members that individual "i" avoided infection up to the generation before "i" was infected in
    for(i in 1:nrow(dat)){
      for (j in 1:ncol(dat)){
        avoid_in_n <- c(avoid_in_n, sum(dat[i,j] > dat[i,]+1))
      }
    }
    first4col <- cbind(gen, inf_out, inf_in_n, avoid_in_n)
    # add "hh_inf_n", "hh_size", "ind_id", "pos", "set"
    k = ncol(dat)
    output <- cbind(first4col, k, k, 1:k, 1, rep(1:nrow(dat), each = k))
    # add "set"
    colnames(output) <- c("gen", "inf_out",	"inf_in_n", "avoid_in_n", "hh_inf_n", "hh_size", "ind_id", "pos", "set")
    return(output)
  }
  
  ## generate all possible sequence of viral intro and subsequent transmission events within hh from j/k households
  ## j = number infected, k = hhsize
  
  m01 <- get_m0k(1)
  m02 <- get_m0k(2)
  m03 <- get_m0k(3)
  m04 <- get_m0k(4)
  m05 <- get_m0k(5)
  m06 <- get_m0k(6)
  m07 <- get_m0k(7)
  
  m11 <- create_stan_input_kk_hh(assign_gens(1, c(NA), n_inf=1))
  m12 <- get_mjk(1,2,m11)
  m13 <- get_mjk(1,3,m11)
  m14 <- get_mjk(1,4,m11)
  m15 <- get_mjk(1,5,m11)
  
  m22 <- create_stan_input_kk_hh(assign_gens(1, c(NA,NA), n_inf=2))
  m23 <- get_mjk(2,3,m22)
  m24 <- get_mjk(2,4,m22)
  m25 <- get_mjk(2,5,m22)
  
  m33 <- create_stan_input_kk_hh(assign_gens(1,c(NA,NA,NA), n_inf = 3))
  m34 <- get_mjk(3,4,m33)
  m35 <- get_mjk(3,5,m33)
  
  m44 <- create_stan_input_kk_hh(assign_gens(1,c(NA,NA,NA,NA), n_inf = 4))
  m45 <- get_mjk(4,5,m44)
  
  m55 <- create_stan_input_kk_hh(assign_gens(1,c(NA,NA,NA,NA,NA), n_inf = 5))
  
  ## combine everything to one dataset
  dat_comb <- rbind(m01,m11,
                    m02,m12,m22,
                    m03,m13,m23,m33,
                    m04,m14,m24,m34,m44,
                    m05,m15,m25,m35,m45,m55, #no m45 in gva
                    m06, #only m06 exists in gva
                    m07) #only m07 exists gva
  return(dat_comb)
}

