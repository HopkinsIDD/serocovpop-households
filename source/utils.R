## reloads packages and source file
reload_source <- function(){
  if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
  if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
  if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
  if (!require('purrr')) install.packages('purrr'); library('purrr')
  if (!require('readr')) install.packages('readr'); library('readr')
  if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
  if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
  if (!require('knitr')) install.packages('knitr'); library('knitr')

  if (!require('rstan')) install.packages('rstan'); library('rstan')
  if (!require('kableExtra')) install.packages('kableExtra'); library('kableExtra')
  if (!require('boot')) install.packages('boot'); library('boot')
  if (!require('loo')) install.packages('loo'); library('loo')
  if (!require('scales')) install.packages('scales'); library('scales')
  if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
  if (!require('Hmisc')) install.packages('Hmisc'); library('Hmisc')
  if (!require('forcats')) install.packages('forcats'); library('forcats')
  if (!require('cowplot')) install.packages('cowplot'); library('cowplot')

  #library(tidyverse)
  #library(janitor)
  source(here("source/utils.R"))
}


run_analysis_bq <- function(bq_model,
                            dat_expand,
                            coef_eqn,
                            save_stan=F,
                            ...){
  
  ## aggregate data to household

  ana_dat_expand <- dat_expand %>% 
    ## arrange data by household, set, and ind_id
    arrange(Household_codbar,set,ind_id) %>%
    # create individual id (1 to 4534), household id (1 to 2267), set id (1 to 3287)
    mutate(household_id=as.numeric(as.factor(Household_codbar)),
           individual_id=as.numeric(as.factor(ind_id/10+Household_codbar)),
           set_id=as.numeric(as.factor(set/1000+Household_codbar)),
           row_id = row_number()) 
  
  # create model matrix based on individual level characteristics that Q and B are based on
  X <- model.matrix(as.formula(paste("pos ~", coef_eqn)), data=ana_dat_expand) 
  
  ## run stan model
  if(save_stan){
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   X = X,# model matrix,
                                   inf_out = ana_dat_expand$inf_out, # infected? 0,1 
                                   inf_in_n = ana_dat_expand$inf_in_n, # total symptomatically infected the generation before "i" was infected
                                   avoid_in_n = ana_dat_expand$avoid_in_n,
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp)# row number of the start of each set
                         ),
                         ...)
  }else{
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   X = X,# model matrix,
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   inf_out = ana_dat_expand$inf_out, # infected?0,1
                                   inf_in_n = ana_dat_expand$inf_in_n,# total symptomatically infected the generation before "i" was infected
                                   avoid_in_n = ana_dat_expand$avoid_in_n,
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp)# row number of the start of each set
                         ),
                         ...)
  }
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  if(save_stan){
    return(list(model_mtx = X,
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                stan_ll=stan_ll,
                stan_est=stan_est))
  } else{
    return(list(beta = extract(stan_est, pars="beta_B")[[1]],
                theta = extract(stan_est, pars="beta_Q")[[1]],
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                model_mtx = X,
                obs = nrow(ana_dat),
                stan_ll=stan_ll))
  }
}


run_analysis_bq_contact <- function(bq_model,
                            dat_expand,
                            coef_eqn1,
                            coef_eqn,
                            save_stan=F,
                            ...){
  
  ## aggregate data to household
  
  ana_dat_expand <- dat_expand %>% 
    ## arrange data by household, set, and ind_id
    arrange(Household_codbar,set,ind_id) %>%
    # create individual id (1 to 4534), household id (1 to 2267), set id (1 to 3287)
    mutate(household_id=as.numeric(as.factor(Household_codbar)),
           individual_id=as.numeric(as.factor(ind_id/10+Household_codbar)),
           set_id=as.numeric(as.factor(set/1000+Household_codbar)),
           row_id = row_number()) 
  
  # create model matrix based on individual level characteristics that Q and B are based on
  X <- model.matrix(as.formula(paste("pos ~", coef_eqn)), data=ana_dat_expand) 
  X1 <- model.matrix(as.formula(paste("pos ~", coef_eqn1)), data=ana_dat_expand) 
  
  ## run stan model
  if(save_stan){
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   p_vars1 = ncol(X1), 
                                   X = X,# model matrix
                                   X1 = X1,
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   inf_out = ana_dat_expand$inf_out, # infected outside? 0,1
                                   inf_in_n = ana_dat_expand$inf_in_n,# total infected the generation before "i" was infected
                                   avoid_in_n = ana_dat_expand$avoid_in_n, # total infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp)# row number of the start of each set
                         ),
                         ...)
  }else{
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   X = X,# model matrix
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   inf_out = ana_dat_expand$inf_out,# infected outside? 0,1
                                   inf_in_n = ana_dat_expand$inf_in_n, # total infected the generation before "i" was infected
                                   avoid_in_n = ana_dat_expand$avoid_in_n,# total infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp)# row number of the start of each set
                         ),
                         ...)
  }
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  if(save_stan){
    return(list(model_mtx = X,
                model_mtx = X1,
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                stan_ll=stan_ll,
                stan_est=stan_est))
  } else{
    return(list(beta = extract(stan_est, pars="beta_B")[[1]],
                theta = extract(stan_est, pars="beta_Q")[[1]],
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                model_mtx = X,
                model_mtx = X1,
                obs = nrow(ana_dat),
                stan_ll=stan_ll))
  }
}



run_analysis_bq_infector <- function(bq_model,
                                         dat_expand,
                                         coef_eqn,
                                         save_stan=F,
                                         ...){
  
  ## aggregate data to household
  
  ana_dat_expand <- dat_expand %>% 
    ## arrange data by household, set, and ind_id
    arrange(Household_codbar,set,ind_id) %>%
    # create individual id (1 to 4534), household id (1 to 2267), set id (1 to 3287)
    mutate(household_id=as.numeric(as.factor(Household_codbar)),
           individual_id=as.numeric(as.factor(ind_id/10+Household_codbar)),
           set_id=as.numeric(as.factor(set/1000+Household_codbar)),
           row_id = row_number()) 
  
  # create model matrix based on individual level characteristics that Q and B are based on
  X <- model.matrix(as.formula(paste("pos ~", coef_eqn)), data=ana_dat_expand) 
  
  ## run stan model
  if(save_stan){
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   X = X,# model matrix
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   inf_out = ana_dat_expand$inf_out,# infected outside? 0,1
                                   inf_in_n = ana_dat_expand$inf_in_n, # total infected the generation before "i" was infected
                                   inf_in_n_0509 = ana_dat_expand$inf_in_n_sym_0509 + ana_dat_expand$inf_in_n_asym_0509, # total in this age group infected the generation before "i" was infected
                                   avoid_in_n_0509 = ana_dat_expand$avoid_in_n_sym_0509 + ana_dat_expand$avoid_in_n_asym_0509, # total infected in this age group that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_1019 = ana_dat_expand$inf_in_n_sym_1019 + ana_dat_expand$inf_in_n_asym_1019,
                                   avoid_in_n_1019 = ana_dat_expand$avoid_in_n_sym_1019 + ana_dat_expand$avoid_in_n_asym_1019,
                                   inf_in_n_2049 = ana_dat_expand$inf_in_n_sym_2049 + ana_dat_expand$inf_in_n_asym_2049,
                                   avoid_in_n_2049 = ana_dat_expand$avoid_in_n_sym_2049 + ana_dat_expand$avoid_in_n_asym_2049,
                                   inf_in_n_5064 = ana_dat_expand$inf_in_n_sym_5064 + ana_dat_expand$inf_in_n_asym_5064,
                                   avoid_in_n_5064 = ana_dat_expand$avoid_in_n_sym_5064 + ana_dat_expand$avoid_in_n_asym_5064,
                                   inf_in_n_65up = ana_dat_expand$inf_in_n_sym_65up + ana_dat_expand$inf_in_n_asym_65up,
                                   avoid_in_n_65up = ana_dat_expand$avoid_in_n_sym_65up + ana_dat_expand$avoid_in_n_asym_65up,
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id), #hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp)# row number of the start of each set
                         ),
                         ...)
  }else{
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   X = X,# model matrix
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   inf_out = ana_dat_expand$inf_out,# infected outside? 0,1
                                   inf_in_n = ana_dat_expand$inf_in_n,# total  infected the generation before "i" was infected
                                   inf_in_n_0509 = ana_dat_expand$inf_in_n_sym_0509 + ana_dat_expand$inf_in_n_asym_0509, # total in this age gruop infected the generation before "i" was infected
                                   avoid_in_n_0509 = ana_dat_expand$avoid_in_n_sym_0509 + ana_dat_expand$avoid_in_n_asym_0509,# total infected in this age group that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_1019 = ana_dat_expand$inf_in_n_sym_1019 + ana_dat_expand$inf_in_n_asym_1019,
                                   avoid_in_n_1019 = ana_dat_expand$avoid_in_n_sym_1019 + ana_dat_expand$avoid_in_n_asym_1019,
                                   inf_in_n_2049 = ana_dat_expand$inf_in_n_sym_2049 + ana_dat_expand$inf_in_n_asym_2049,
                                   avoid_in_n_2049 = ana_dat_expand$avoid_in_n_sym_2049 + ana_dat_expand$avoid_in_n_asym_2049,
                                   inf_in_n_5064 = ana_dat_expand$inf_in_n_sym_5064 + ana_dat_expand$inf_in_n_asym_5064,
                                   avoid_in_n_5064 = ana_dat_expand$avoid_in_n_sym_5064 + ana_dat_expand$avoid_in_n_asym_5064,
                                   inf_in_n_65up = ana_dat_expand$inf_in_n_sym_65up + ana_dat_expand$inf_in_n_asym_65up,
                                   avoid_in_n_65up = ana_dat_expand$avoid_in_n_sym_65up + ana_dat_expand$avoid_in_n_asym_65up,
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id), #hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp)# row number of the start of each set
                         ),
                         ...)
  }
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  if(save_stan){
    return(list(model_mtx = X,
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                stan_ll=stan_ll,
                stan_est=stan_est))
  } else{
    return(list(beta = extract(stan_est, pars="beta_B")[[1]],
                theta = extract(stan_est, pars="beta_Q")[[1]],
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                model_mtx = X,
                obs = nrow(ana_dat),
                stan_ll=stan_ll))
  }
}



run_analysis_bq_sym_infector <- function(bq_model,
                                dat_expand,
                                coef_eqn,
                                save_stan=F,
                                ...){
  
  ## aggregate data to household
  
  ana_dat_expand <- dat_expand %>% 
    ## arrange data by household, set, and ind_id
    arrange(Household_codbar,set,ind_id) %>%
    # create individual id (1 to 4534), household id (1 to 2267), set id (1 to 3287)
    mutate(household_id=as.numeric(as.factor(Household_codbar)),
           individual_id=as.numeric(as.factor(ind_id/10+Household_codbar)),
           set_id=as.numeric(as.factor(set/1000+Household_codbar)),
           row_id = row_number()) 
  
  # create model matrix based on individual level characteristics that Q and B are based on
  X <- model.matrix(as.formula(paste("pos ~", coef_eqn)), data=ana_dat_expand) 
  
  # create "infector_rowid_matrix", potential infector of each rowid, 0 if no infector
  infector_rowid_matrix <- matrix(0,nrow = nrow(ana_dat_expand),ncol=max(ana_dat_expand$inf_in_n))
  for (i in 1:nrow(ana_dat_expand)){
    tmp <- NULL
    tmp <- which((ana_dat_expand$gen[i] == ana_dat_expand$gen + 1) & (ana_dat_expand$household_id[i] == ana_dat_expand$household_id) & (ana_dat_expand$set_id[i] == ana_dat_expand$set_id))
    if (length(tmp)!=0) {infector_rowid_matrix[i,1:length(tmp)] <- tmp}
  }
  
  ## run stan model
  if(save_stan){
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   X = X,# model matrix
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   inf_out = ana_dat_expand$inf_out, # infected outside? 0,1
                                   inf_in_n = ana_dat_expand$inf_in_n, # total infected the generation before "i" was infected
                                   inf_in_n_sym_0509 = ana_dat_expand$inf_in_n_sym_0509,# total in this age group symptomatically infected the generation before "i" was infected
                                   avoid_in_n_sym_0509 = ana_dat_expand$avoid_in_n_sym_0509, # total in this age group symptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_asym_0509 = ana_dat_expand$inf_in_n_asym_0509,# total in this age group asymptomatically infected the generation before "i" was infected
                                   avoid_in_n_asym_0509 = ana_dat_expand$avoid_in_n_asym_0509, # total in this age group asymptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_sym_1019 = ana_dat_expand$inf_in_n_sym_1019,
                                   avoid_in_n_sym_1019 = ana_dat_expand$avoid_in_n_sym_1019,
                                   inf_in_n_asym_1019 = ana_dat_expand$inf_in_n_asym_1019,
                                   avoid_in_n_asym_1019 = ana_dat_expand$avoid_in_n_asym_1019,
                                   inf_in_n_sym_2049 = ana_dat_expand$inf_in_n_sym_2049,
                                   avoid_in_n_sym_2049 = ana_dat_expand$avoid_in_n_sym_2049,
                                   inf_in_n_asym_2049 = ana_dat_expand$inf_in_n_asym_2049,
                                   avoid_in_n_asym_2049 = ana_dat_expand$avoid_in_n_asym_2049,
                                   inf_in_n_sym_5064 = ana_dat_expand$inf_in_n_sym_5064,
                                   avoid_in_n_sym_5064 = ana_dat_expand$avoid_in_n_sym_5064,
                                   inf_in_n_asym_5064 = ana_dat_expand$inf_in_n_asym_5064,
                                   avoid_in_n_asym_5064 = ana_dat_expand$avoid_in_n_asym_5064,
                                   inf_in_n_sym_65up = ana_dat_expand$inf_in_n_sym_65up,
                                   avoid_in_n_sym_65up = ana_dat_expand$avoid_in_n_sym_65up,
                                   inf_in_n_asym_65up = ana_dat_expand$inf_in_n_asym_65up,
                                   avoid_in_n_asym_65up = ana_dat_expand$avoid_in_n_asym_65up,
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id), #hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp),# row number of the start of each set
                                   max_infectors = max(ana_dat_expand$inf_in_n),#max number of potential infectors among all
                                   infector_rowid_matrix = infector_rowid_matrix,#potential infectors' rowid per person
                                   sym = ana_dat_expand$symptom_any_mod, # symptom status of each person
                                   sym_age = 10*ana_dat_expand$symptom_any_mod + as.numeric(factor(ana_dat_expand$age_cat,levels=c("5-9","10-19","20-49","50-64","65+"))) # identifier of symptom and age category of each person
                         ),
                         ...)
  }else{
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand),#number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   p_vars = ncol(X),  # number of variables to adjust for
                                   X = X,# model matrix
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   inf_out = ana_dat_expand$inf_out, # infected outside? 0,1
                                   inf_in_n = ana_dat_expand$inf_in_n, # total infected the generation before "i" was infected
                                   inf_in_n_sym_0509 = ana_dat_expand$inf_in_n_sym_0509, # total in this age group symptomatically infected the generation before "i" was infected
                                   avoid_in_n_sym_0509 = ana_dat_expand$avoid_in_n_sym_0509, # total in this age group symptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_asym_0509 = ana_dat_expand$inf_in_n_asym_0509,# total in this age group asymptomatically infected the generation before "i" was infected
                                   avoid_in_n_asym_0509 = ana_dat_expand$avoid_in_n_asym_0509, # total in this age group asymptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_sym_1019 = ana_dat_expand$inf_in_n_sym_1019,
                                   avoid_in_n_sym_1019 = ana_dat_expand$avoid_in_n_sym_1019,
                                   inf_in_n_asym_1019 = ana_dat_expand$inf_in_n_asym_1019,
                                   avoid_in_n_asym_1019 = ana_dat_expand$avoid_in_n_asym_1019,
                                   inf_in_n_sym_2049 = ana_dat_expand$inf_in_n_sym_2049,
                                   avoid_in_n_sym_2049 = ana_dat_expand$avoid_in_n_sym_2049,
                                   inf_in_n_asym_2049 = ana_dat_expand$inf_in_n_asym_2049,
                                   avoid_in_n_asym_2049 = ana_dat_expand$avoid_in_n_asym_2049,
                                   inf_in_n_sym_5064 = ana_dat_expand$inf_in_n_sym_5064,
                                   avoid_in_n_sym_5064 = ana_dat_expand$avoid_in_n_sym_5064,
                                   inf_in_n_asym_5064 = ana_dat_expand$inf_in_n_asym_5064,
                                   avoid_in_n_asym_5064 = ana_dat_expand$avoid_in_n_asym_5064,
                                   inf_in_n_sym_65up = ana_dat_expand$inf_in_n_sym_65up,
                                   avoid_in_n_sym_65up = ana_dat_expand$avoid_in_n_sym_65up,
                                   inf_in_n_asym_65up = ana_dat_expand$inf_in_n_asym_65up,
                                   avoid_in_n_asym_65up = ana_dat_expand$avoid_in_n_asym_65up,
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id), #hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp),# row number of the start of each set
                                   max_infectors = max(ana_dat_expand$inf_in_n),#max number of potential infectors among all
                                   infector_rowid_matrix = infector_rowid_matrix,#potential infectors' rowid per person
                                   sym_age = 10*ana_dat_expand$symptom_any_mod + as.numeric(factor(ana_dat_expand$age_cat,levels=c("5-9","10-19","20-49","50-64","65+")))# identifier of symptom and age category of each person
                         ),
                         ...)
  }
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  if(save_stan){
    return(list(model_mtx = X,
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                stan_ll=stan_ll,
                stan_est=stan_est))
  } else{
    return(list(beta = extract(stan_est, pars="beta_B")[[1]],
                theta = extract(stan_est, pars="beta_Q")[[1]],
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                model_mtx = X,
                obs = nrow(ana_dat),
                stan_ll=stan_ll))
  }
}


run_analysis_bq_sym <- function(bq_model,
                                dat_expand,
                                coef_eqn,
                                save_stan=F,
                                ...){
  
  ana_dat_expand <- dat_expand %>% 
    ## arrange data by household, set, and ind_id
    arrange(Household_codbar,set,ind_id) %>%
    # create individual id (1 to 4534), household id (1 to 2267), set id (1 to 3287)
    mutate(household_id=as.numeric(as.factor(Household_codbar)),
           individual_id=as.numeric(as.factor(ind_id/10+Household_codbar)),
           set_id=as.numeric(as.factor(set/1000+Household_codbar)),
           row_id = row_number()) 
  
  # create model matrix based on individual level characteristics that Q and B are based on
  X <- model.matrix(as.formula(paste("pos ~", coef_eqn)), data=ana_dat_expand) 
  
  # create "infector_rowid_matrix", potential infector of each rowid, 0 if no infector
  infector_rowid_matrix <- matrix(0,nrow = nrow(ana_dat_expand),ncol=max(ana_dat_expand$inf_in_n))
  for (i in 1:nrow(ana_dat_expand)){
    tmp <- NULL
    tmp <- which((ana_dat_expand$gen[i] == ana_dat_expand$gen + 1) & (ana_dat_expand$household_id[i] == ana_dat_expand$household_id) & (ana_dat_expand$set_id[i] == ana_dat_expand$set_id))
    if (length(tmp)!=0) {infector_rowid_matrix[i,1:length(tmp)] <- tmp}
  }
  
  
  ## run stan model
  if(save_stan){
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand), #number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   p_vars = ncol(X), # number of variables to adjust for
                                   X = X,# model matrix
                                   inf_out = ana_dat_expand$inf_out, # infected outside? 0,1
                                   inf_in_n = ana_dat_expand$inf_in_n,  # total infected the generation before "i" was infected
                                   inf_in_n_sym = ana_dat_expand$inf_in_n_sym,# total symptomatically infected the generation before "i" was infected
                                   avoid_in_n_sym = ana_dat_expand$avoid_in_n_sym,  # total in this age group symptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_asym = ana_dat_expand$inf_in_n_asym,# total asymptomatically infected the generation before "i" was infected
                                   avoid_in_n_asym = ana_dat_expand$avoid_in_n_asym, # total in this age group asymptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   gen = ana_dat_expand$gen, # generation each is infected in
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp),# row number of the start of each set
                                   max_infectors = max(ana_dat_expand$inf_in_n),#max number of potential infectors among all
                                   infector_rowid_matrix = infector_rowid_matrix,#potential infectors' rowid per person
                                   sym = ana_dat_expand$symptom_any_mod # symptom status of each person
                         ),
                         ...)
  }else{
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand), #number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   p_vars = ncol(X), # number of variables to adjust for
                                   X = X,# model matrix
                                   inf_out = ana_dat_expand$inf_out,
                                   inf_in_n = ana_dat_expand$inf_in_n, # total infected the generation before "i" was infected
                                   inf_in_n_sym = ana_dat_expand$inf_in_n_sym,# total symptomatically infected the generation before "i" was infected
                                   avoid_in_n_sym = ana_dat_expand$avoid_in_n_sym,  # total in this age group symptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_asym = ana_dat_expand$inf_in_n_asym,# total asymptomatically infected the generation before "i" was infected
                                   avoid_in_n_asym = ana_dat_expand$avoid_in_n_asym, # total in this age group asymptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   gen = ana_dat_expand$gen,# generation each is infected in
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp),# row number of the start of each set
                                   max_infectors = max(ana_dat_expand$inf_in_n),#max number of potential infectors among all
                                   infector_rowid_matrix = infector_rowid_matrix,#potential infectors' rowid per person
                                   sym = ana_dat_expand$symptom_any_mod # symptom status of each person
                         ),
                         ...)
  }
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  if(save_stan){
    return(list(model_mtx = X, # save model matrix
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                stan_ll=stan_ll,
                stan_est=stan_est))
  } else{
    return(list(beta = extract(stan_est, pars="beta_B")[[1]],
                theta = extract(stan_est, pars="beta_Q")[[1]],
                model_mtx = X,
                obs = nrow(ana_dat),
                stan_ll=stan_ll))
  }
}



run_analysis_bq_sym_contact <- function(bq_model,
                                    dat_expand,
                                    coef_eqn,
                                    coef_eqn1,
                                    save_stan=F,
                                    ...){
  
  ana_dat_expand <- dat_expand %>% 
    ## arrange data by household, set, and ind_id
    arrange(Household_codbar,set,ind_id) %>%
    # create individual id (1 to 4534), household id (1 to 2267), set id (1 to 3287)
    mutate(household_id=as.numeric(as.factor(Household_codbar)),
           individual_id=as.numeric(as.factor(ind_id/10+Household_codbar)),
           set_id=as.numeric(as.factor(set/1000+Household_codbar)),
           row_id = row_number()) 
  
  # create model matrix based on individual level characteristics that Q and B are based on
  X <- model.matrix(as.formula(paste("pos ~", coef_eqn)), data=ana_dat_expand) 
  X1 <- model.matrix(as.formula(paste("pos ~", coef_eqn1)), data=ana_dat_expand) 
  
  # create "infector_rowid_matrix", potential infector of each rowid, 0 if no infector
  infector_rowid_matrix <- matrix(0,nrow = nrow(ana_dat_expand),ncol=max(ana_dat_expand$inf_in_n))
  for (i in 1:nrow(ana_dat_expand)){
    tmp <- NULL
    tmp <- which((ana_dat_expand$gen[i] == ana_dat_expand$gen + 1) & (ana_dat_expand$household_id[i] == ana_dat_expand$household_id) & (ana_dat_expand$set_id[i] == ana_dat_expand$set_id))
    if (length(tmp)!=0) {infector_rowid_matrix[i,1:length(tmp)] <- tmp}
  }
  
  
  ## run stan model
  if(save_stan){
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand), #number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   p_vars = ncol(X), # number of variables to adjust for for Q
                                   p_vars1 = ncol(X1), # number of variables to adjust for for B
                                   X = X, # model matrix for Q
                                   X1 = X1, # model matrix for B
                                   inf_out = ana_dat_expand$inf_out, # total infected outside of a hh
                                   inf_in_n = ana_dat_expand$inf_in_n, # total infected the generation before "i" was infected
                                   inf_in_n_sym = ana_dat_expand$inf_in_n_sym,# total symptomatically infected the generation before "i" was infected
                                   avoid_in_n_sym = ana_dat_expand$avoid_in_n_sym, # total symptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   inf_in_n_asym = ana_dat_expand$inf_in_n_asym, # total asymptomatically infected the generation before "i" was infected
                                   avoid_in_n_asym = ana_dat_expand$avoid_in_n_asym, # total asymptomatically infected that individual "i" avoided infection from up to the generation before "i" was infected in
                                   gen = ana_dat_expand$gen, # the generation each id is infected in
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp),# row number of the start of each set
                                   max_infectors = max(ana_dat_expand$inf_in_n),#max number of potential infectors among all
                                   infector_rowid_matrix = infector_rowid_matrix,#potential infectors' rowid per person
                                   sym = ana_dat_expand$symptom_any_mod  # symptom status of each person
                         ),
                         ...)
  }else{
    stan_est <- sampling(bq_model,
                         data=list(N_survey_expand = nrow(ana_dat_expand), #number of rows
                                   N_ind = length(unique(ana_dat_expand$codbar_entry)),# number of individuals
                                   N_set = max(ana_dat_expand$set_id),#number of sets 
                                   N_hh = max(ana_dat_expand$household_id),#number of HHs in the survey
                                   p_vars = ncol(X), # number of variables to adjust for
                                   p_vars1 = ncol(X1), # number of variables to adjust for
                                   X = X, # model matrix for Q
                                   X1 = X1, # model matrix for B
                                   inf_out = ana_dat_expand$inf_out,
                                   inf_in_n = ana_dat_expand$inf_in_n,
                                   inf_in_n_sym = ana_dat_expand$inf_in_n_sym,
                                   avoid_in_n_sym = ana_dat_expand$avoid_in_n_sym,
                                   inf_in_n_asym = ana_dat_expand$inf_in_n_asym,
                                   avoid_in_n_asym = ana_dat_expand$avoid_in_n_asym,
                                   gen = ana_dat_expand$gen,
                                   set_per_ind = ana_dat_expand$set_id,# set id of each row
                                   indid_per_row = ana_dat_expand$individual_id, # individual id of each row
                                   hhid_per_set = ana_dat_expand %>% distinct(household_id, set_id) %>% pull(household_id),#hh id of each set
                                   first_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=min(set_id))%>% pull(tmp),#first set id of a hh
                                   last_setid_per_hh = ana_dat_expand %>% group_by(household_id) %>% summarise(tmp=max(set_id))%>% pull(tmp),#last set id of a hh
                                   hhsize = ana_dat_expand %>% distinct(household_id, hh_size) %>% pull(hh_size),#hhsize per household
                                   row_num_first = ana_dat_expand %>% group_by(set_id) %>% summarise(tmp=min(row_id))%>% pull(tmp),# row number of the start of each set
                                   max_infectors = max(ana_dat_expand$inf_in_n),#max number of potential infectors among all
                                   infector_rowid_matrix = infector_rowid_matrix,#potential infectors' rowid per person
                                   sym = ana_dat_expand$symptom_any_mod   # symptom status of each person
                         ),
                         ...)
  }
  stan_ll <- loo::extract_log_lik(stan_est) %>% loo::loo()
  if(save_stan){
    return(list(model_mtx = X, # save model matrix
                model_mtx1 = X1, # save model matrix
                input_data = ana_dat_expand, # save input data with hh id, set id, indivdiual id and individual level characteristics
                stan_ll=stan_ll,
                stan_est=stan_est))
  } else{
    return(list(beta = extract(stan_est, pars="beta_B")[[1]],
                theta = extract(stan_est, pars="beta_Q")[[1]],
                model_mtx = X,
                model_mtx1 = X1,
                obs = nrow(ana_dat),
                stan_ll=stan_ll))
  }
}
