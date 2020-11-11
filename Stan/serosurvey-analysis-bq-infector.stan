//
// This Stan program .....
//
data {
  int<lower=1> N_survey_expand; //number of rows in the expanded dataset
  int<lower=1> N_ind; // number of individuals
  int<lower=1> p_vars; //number of variables to adjust for
  matrix[N_survey_expand, p_vars] X; //variable model matrix
  int<lower=1> N_set; //number of unique set 
  int<lower=1> N_hh; //number of HH in the survey
  vector<lower=0>[N_survey_expand] inf_out; // 1 if infected outside 0 otherwise
  vector<lower=0>[N_survey_expand] inf_in_n;// total infected the generation before hh member "i" was infected
  vector<lower=0>[N_survey_expand] inf_in_n_0509; // total in xx age group infected the generation before hh member "i" was infected
  vector<lower=0>[N_survey_expand] avoid_in_n_0509; // total infected in this age group that individual "i" avoided infection from up to the generation before "i" was infected in
  vector<lower=0>[N_survey_expand] inf_in_n_1019;
  vector<lower=0>[N_survey_expand] avoid_in_n_1019;
  vector<lower=0>[N_survey_expand] inf_in_n_2049;
  vector<lower=0>[N_survey_expand] avoid_in_n_2049;
  vector<lower=0>[N_survey_expand] inf_in_n_5064;
  vector<lower=0>[N_survey_expand] avoid_in_n_5064;
  vector<lower=0>[N_survey_expand] inf_in_n_65up;
  vector<lower=0>[N_survey_expand] avoid_in_n_65up;
  vector<lower=0>[N_survey_expand] gen; // generation of infection, 0 out, 9999 not infected
  int set_per_ind [N_survey_expand]; //set id of each row
  int indid_per_row[N_survey_expand]; // ind id of each row
  int hhid_per_set [N_set]; //hh id of each set
  int first_setid_per_hh [N_hh]; // first set id of a hh
  int last_setid_per_hh [N_hh]; // last set id of a hh
  int hhsize [N_hh]; // household size per hh
  int row_num_first [N_set]; // row number of the first element in each set
  
}

parameters {
  vector[p_vars] beta_B; // coefficients for community vars
  vector[p_vars] beta_Q; // coefficients for hh vars
  real alpha0509;// coefficient for infector's age group
  real alpha1019;// coefficient for infector's age group
  real alpha5064;// coefficient for infector's age group
  real alpha65up;// coefficient for infector's age group
}

transformed parameters {
  vector<lower=0, upper=1>[N_survey_expand] B; // probability of escaping infection outside of the HH
  vector<lower=0, upper=1>[N_survey_expand] Q0509; // probability of escaping infection from one symp infected person in this age range inside of the HH
  vector<lower=0, upper=1>[N_survey_expand] Q1019; 
  vector<lower=0, upper=1>[N_survey_expand] Q2049; 
  vector<lower=0, upper=1>[N_survey_expand] Q5064; 
  vector<lower=0, upper=1>[N_survey_expand] Q65up; 
  vector<lower=0, upper=1>[N_survey_expand] m_row; // inividual contribution to m_ij per HH
  vector<lower=0, upper=1>[N_set] set_probs; // product of inividual contribution to likelihood for every set in every HH
  vector<lower=0, upper=1>[N_hh] hh_probs; // sum of set contribution to likelihood in every HH
  
  // get Bs and Qs
  B = inv_logit(X * beta_B);
  Q0509 = inv_logit(X * beta_Q + alpha0509);
  Q1019 = inv_logit(X * beta_Q + alpha1019);
  Q2049 = inv_logit(X * beta_Q);
  Q5064 = inv_logit(X * beta_Q + alpha5064);
  Q65up = inv_logit(X * beta_Q + alpha65up);
    
  // individual contribution to likelihood given the
  // generation this person is infected in.
  for (i in 1:N_survey_expand){
      m_row[i] = pow(1 - B[i], inf_out[i]) * pow(B[i], 1-inf_out[i]) * 
      (1 - (1 - pow(0, inf_in_n[i]))  *  
      pow(Q0509[i], inf_in_n_0509[i]) * 
      pow(Q1019[i], inf_in_n_1019[i]) * 
      pow(Q2049[i], inf_in_n_2049[i]) *
      pow(Q5064[i], inf_in_n_5064[i]) * 
      pow(Q65up[i], inf_in_n_65up[i])) * 
      pow(Q0509[i], avoid_in_n_0509[i])* 
      pow(Q1019[i], avoid_in_n_1019[i])*
      pow(Q2049[i], avoid_in_n_2049[i])* 
      pow(Q5064[i], avoid_in_n_5064[i])*
      pow(Q65up[i], avoid_in_n_65up[i]);
  }
  
  // contribution to likelihood based on each sequence of viral intro and subsequent transmission events (a set) of a HH = Pr(HH_h,k)
  set_probs = rep_vector(1, N_set);

  for (i in 1:N_survey_expand){
      // then multiply prob of the same set
      set_probs[set_per_ind[i]] = set_probs[set_per_ind[i]] * m_row[i];
  }

  // contribution to likelihood of each HH,  Pr(HH_k)
  // this is essentially the probability that j of k initial susceptibles within a household 
  // are infected during the course of the epidemic
  hh_probs = rep_vector(0, N_hh);

  for (i in 1:N_set){
      // then add up prob in the same household
      hh_probs[hhid_per_set[i]] = hh_probs[hhid_per_set[i]]  + set_probs[i];
  }
}
  

model {
  target+= normal_lpdf(beta_B | 0, 1.5); // beta_B ~ normal(0,1); prior for beta_B
  target+= normal_lpdf(beta_Q | 0, 1.5); // beta_Q ~ normal(0,1); prior for beta_Q
  target+= normal_lpdf(alpha0509 | 0, 1.5); 
  target+= normal_lpdf(alpha1019 | 0, 1.5); 
  target+= normal_lpdf(alpha5064 | 0, 1.5); 
  target+= normal_lpdf(alpha65up | 0, 1.5); 
  
  // Likelihood = multiply Pr(HH_k) from all households 
  target += sum(log(hh_probs));
}

generated quantities {
 int set_id [N_hh];//set id of the one set drawn from each household
  int cat [N_hh];
  int cat_tmp;
  real bern_tmp;
  real probs_ind_tmp; // probabilities of being infected by a symptomatic person (input into Bernoulli)
  real probs_ind[N_ind]; // store the above in a vector, those not infected are assigned 9999
  vector[N_hh] log_lik;
  matrix[N_ind,5] infoutsymasymsize; 

  
  for (i in 1:N_hh){
    // how many sets per household
    int num_sets_per_hh = last_setid_per_hh[i] - first_setid_per_hh[i] + 1; 
    // placeholder of likelihood from all sets of a household
    vector [num_sets_per_hh] probs; 
    for (j in 1:num_sets_per_hh) {
      probs[j] = set_probs[first_setid_per_hh[i] + j - 1];
    }
    // set id of the one set drawn from each household
    cat_tmp = categorical_rng(softmax(probs));
    set_id[i] = first_setid_per_hh[i] + cat_tmp - 1;
    cat[i] = cat_tmp;
    // prob of being infected from a symptomatic person
    for (k in 1:hhsize[i]){ // for each person in that set
      int rowid = row_num_first[set_id[i]] + k - 1; // row number of that person
      int indid = indid_per_row[rowid]; // indiviudal id of that person
      
      if(gen[rowid]==9999){ // if not infected
        infoutsymasymsize[indid] = [0,0,0,0,hhsize[i]];
      } else if (gen[rowid]==0){ // if infected outside
        infoutsymasymsize[indid] = [1,1,0,0,hhsize[i]];
      } else { // if infected within a household
        infoutsymasymsize[indid] = [1,0,0,0,hhsize[i]];
      }
    }
  }
  log_lik = log(hh_probs);
}