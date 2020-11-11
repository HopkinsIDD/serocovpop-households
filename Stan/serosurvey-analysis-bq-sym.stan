//
// This Stan program .....
//
  data {
    int<lower=1> N_survey_expand; // number of rows
    int<lower=1> N_ind; // number of individuals
    int<lower=1> N_set; // number of sets 
    int<lower=1> N_hh; // number of HHs in the survey
    int<lower=1> p_vars; // number of variables to adjust for
    matrix[N_survey_expand, p_vars] X; // variable model matrix
    int inf_out [N_survey_expand]; // 1 if infected outside 0 otherwise
    int inf_in_n [N_survey_expand]; // total infected the generation before hh member "i" was infected
    int inf_in_n_sym [N_survey_expand]; // total symptomatically infected the generation before hh member "i" was infected
    int avoid_in_n_sym [N_survey_expand];  // total number of asymptomatically infected household members that individual "i" avoided infection from up to the generation before "i" was infected in
    int inf_in_n_asym [N_survey_expand];   // total asymptomatically infected the generation before hh member "i" was infected
    int avoid_in_n_asym [N_survey_expand]; // total number of symptomatically infected household members that individual "i" avoided infection from up to the generation before "i" was infected in
    int gen [N_survey_expand]; // generation of infection, 0 out, 9999 not infected
    int set_per_ind [N_survey_expand]; // set id of each row
    int indid_per_row[N_survey_expand]; // ind id of each row
    int hhid_per_set [N_set]; // hh id of each set
    int first_setid_per_hh [N_hh]; // first set id of a hh
    int last_setid_per_hh [N_hh]; // last set id of a hh
    int hhsize [N_hh]; // household size per hh
    int row_num_first [N_set]; // row number of the first element in each set
    int<lower=1> max_infectors; //max number of potential infectors of our study participants
    int infector_rowid_matrix [N_survey_expand,max_infectors];// rowid of all potential infectors of a hh member, depending on the sequences of viral intro and trans events within
    int sym [N_survey_expand]; // symptom status of each ind ordered by row_id
  }

parameters {
  vector[p_vars] beta_B; // coefficients for community vars
  vector[p_vars] beta_Q; // coefficients for hh vars
  real alpha_sym; // coefficients for infector's symptom status
}

transformed parameters {
  vector<lower=0, upper=1>[N_survey_expand] B; // probability of escaping infection outside of the HH
  vector<lower=0, upper=1>[N_survey_expand] Q_sym; // probability of escaping infection from one symp infected person inside of the HH
  vector<lower=0, upper=1>[N_survey_expand] Q_asym; // probability of escaping infection from one asymp infected person inside of the HH
  vector<lower=0, upper=1>[N_survey_expand] m_row; // inividual contribution to the likelihood
  vector<lower=0, upper=1>[N_set] set_probs; // product of individual contribution to likelihood, one per set
  vector<lower=0, upper=1>[N_hh] hh_probs; // sum of set contribution to likelihood in every HH
  
  // get Bs and Qs
  B = inv_logit(X * beta_B);
  Q_sym = inv_logit(X * beta_Q + alpha_sym);
  Q_asym = inv_logit(X * beta_Q);
  
  // individual contribution to likelihood given the
  // generation this person is infected in.
  for (i in 1:N_survey_expand){
    m_row[i] = pow(1 - B[i], inf_out[i]) * pow(B[i], 1-inf_out[i]) * 
    (1 - (1 - pow(0, inf_in_n[i])) * pow(Q_sym[i], inf_in_n_sym[i]) * 
    pow(Q_asym[i], inf_in_n_asym[i])) * 
    pow(Q_sym[i], avoid_in_n_sym[i]) * 
    pow(Q_asym[i], avoid_in_n_asym[i]);
  }
  
  // contribution to likelihood based on each sequence of viral intro and subsequent transmission events (a set) of a HH = Pr(HH_h,k)
  set_probs = rep_vector(1, N_set);
  
  for (i in 1:N_survey_expand){
    // then multiply prob from the same set
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
  target+= normal_lpdf(alpha_sym | 0, 1.5); // beta_Q ~ normal(0,1); prior for alpha_sym
  
  // Likelihood = multiply Pr(HH_k) from all households 
  target += sum(log(hh_probs));
}

generated quantities {
  int set_id [N_hh];//set id of the one set drawn from each household
  int cat [N_hh];
  int cat_sym [N_ind];
  int cat_tmp;
  int cat_sym_tmp;
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
      } else {// if infected within a household
        // get rowid of this person's infectors
        // number of potential infectors per ind from the chosen set
        int num_infectors_per_ind = inf_in_n[rowid];
      
        // create placeholder
        // placeholder of 1-Q from all potential infectors
        vector [num_infectors_per_ind] probs_sym; 
        int all_rowid_infector [num_infectors_per_ind]; //this stores rowid of all infectors of this person
        
        for (n in 1:num_infectors_per_ind){ // for each potential infector
          // get potential infector's row id based on set_id and individual id
          int rowid_infector = infector_rowid_matrix[rowid,n];
          all_rowid_infector[n] = rowid_infector;
          //
          if(sym[rowid_infector]==1){
            probs_sym[n] = 1-Q_sym[rowid_infector];
          } else{
            probs_sym[n] = 1-Q_asym[rowid_infector];
          }
        }
        
        // draw from infector
        cat_sym_tmp = categorical_rng(softmax(probs_sym));
        cat_sym[indid] = cat_sym_tmp;
        infoutsymasymsize[indid] = [1,0,sym[all_rowid_infector[cat_sym_tmp]],1-sym[all_rowid_infector[cat_sym_tmp]],hhsize[i]];
      }
    }
  }
  log_lik = log(hh_probs);
}


