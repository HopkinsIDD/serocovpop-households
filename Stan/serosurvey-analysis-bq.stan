//
// This Stan program .....
//
data {
  int<lower=1> N_survey_expand; //number of rows in the stan input data
  int<lower=1> N_ind; // number of individuals
  int<lower=1> N_set; //number of sets 
  int<lower=1> N_hh; //number of HH 
  int<lower=1> p_vars; //number of variables to adjust for
  matrix[N_survey_expand, p_vars] X; //variable model matrix
  vector<lower=0>[N_survey_expand] inf_out;// 1 if infected outside 0 otherwise
  vector<lower=0>[N_survey_expand] inf_in_n; // total infected the generation before hh member "i" was infected
  vector<lower=0>[N_survey_expand] avoid_in_n; // total number of infected household members that individual "i" avoided infection from up to the generation before "i" was infected in
  vector<lower=0>[N_survey_expand] gen; // generation of infection, 0 out, 9999 not infected
  int set_per_ind [N_survey_expand]; //set id of each row
  int indid_per_row[N_survey_expand]; // ind id of each row
  int hhid_per_set [N_set]; //hh id of each set
  int first_setid_per_hh [N_hh]; // first set id of a hh
  int last_setid_per_hh [N_hh]; // last set id of a hh
  int hhsize [N_hh]; // household size per hh
  int row_num_first [N_set]; // row number of the first element in each set
  int max_gen_per_set [N_set]; // maximum generation of each set
}

parameters {
  vector[p_vars] beta_B; // coefficients for community vars
  vector[p_vars] beta_Q; // coefficients for hh vars
}

transformed parameters {
  vector<lower=0, upper=1>[N_survey_expand] B; // probability of escaping infection outside of a HH
  vector<lower=0, upper=1>[N_survey_expand] Q; // probability of escaping infection from one infected person inside of a HH
  vector<lower=0, upper=1>[N_survey_expand] m_row; // inividual contribution to m_ij per HH
  vector<lower=0, upper=1>[N_set] set_probs; // product of inividual contribution to likelihood for every set in every HH
  vector<lower=0, upper=1>[N_hh] hh_probs; // sum of set contribution to likelihood in every HH
  
  // get Bs and Qs
  B = inv_logit(X * beta_B);
  Q = inv_logit(X * beta_Q);
    
  // individual contribution to likelihood given the
  // generation this person is infected in.
  for (i in 1:N_survey_expand){
      m_row[i] = pow(1 - B[i], inf_out[i]) * pow(B[i], 1-inf_out[i]) * (1 - (1 - pow(0, inf_in_n[i])) * pow(Q[i], inf_in_n[i])) * pow(Q[i], avoid_in_n[i]);
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
  
  // Likelihood = multiply Pr(HH_k) from all households 
  target += sum(log(hh_probs));
}

generated quantities {
  int set_id [N_hh];//set id of the one set drawn from each household
  int max_gen [N_hh]; //maximum generation of each simulated set
  int cat [N_hh];
  int cat_tmp;
  vector[N_hh] log_lik;
  matrix[N_ind,5] infoutsymasymsize; 

  
  for (i in 1:N_hh){
    // how many sets per household
    int num_sets_per_hh = last_setid_per_hh[i] - first_setid_per_hh[i] + 1; 
    // placeholder of likelihood from all sets of a household
    vector [num_sets_per_hh] probs; 
    vector [num_sets_per_hh] log_probs; 
    for (j in 1:num_sets_per_hh) {
      probs[j] = set_probs[first_setid_per_hh[i] + j - 1];
      log_probs[j] = log(probs[j]);
    }
    // set id of the one set drawn from each household
    cat_tmp = categorical_rng(softmax(log_probs));
    set_id[i] = first_setid_per_hh[i] + cat_tmp - 1;
    cat[i] = cat_tmp;
    // get maximum generation of the selected set
    max_gen[i] = max_gen_per_set[set_id[i]];
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