Code to reproduce results from _Household Transmission of SARS-COV-2: Insights from a Population-based Serological Survey_.

## Code Notes

In the `Stan/` directory are stan code for running a series of different models we compared in the analyses. Code that runs the stan models are in HHtrans_public.Rmd. _Note that generated data are saved using [git-lfs](https://www.atlassian.com/git/tutorials/git-lfs) so if you want them locally you need to use `git lfs pull`_.

The models are as follows:

| file name |model number | description |
|-----------| ---| ----|
| serosurvey-analysis-bq.stan | 0 | null model |
| serosurvey-analysis-bq.stan | 1 | age (infectee) |  
| serosurvey-analysis-bq.stan | 2 | age (infectee), sex (infectee) |
| serosurvey-analysis-bq.stan | 3 | age (infectee), sex (infectee), age and sex interaction term (infectee) |
| serosurvey-analysis-bq-sym.stan | 4 | age (infectee), sex (infectee), symptom (infector) |
| serosurvey-analysis-bq-sym-contact.stan | 5 | age (infectee), sex (infectee), reduced extra-hh contact (infectee), symptom (infector) |
| serosurvey-analysis-bq-sym-contact.stan | 6 | age (infectee), sex (infectee), reduced extra-hh contact (infectee), extra-hh contact freq (infectee), symptom (infector) |
| serosurvey-analysis-bq-sym-infector-2049ref.stan | 7 | age (infectee), sex (infectee), symptom (infector), age (infector) | 
| serosurvey-analysis-bq-sym-infector-2049ref.stan | 8 | sex (infectee), symptom (infector), age (infector) | 
| serosurvey-analysis-bq-sym-infector.stan | 9 | age (infectee), sex (infectee), age (infector) | 


## Data Notes
We have included dummy data for in order to help users test the code. Data will be made available upon request by contacting azman@jhu.edu. 

We have tried our best to make this code understandable and usable but please do get in touch if you spot an issue or have a question. 
