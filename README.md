# COVID-19-Policy-and-Changepoints
This repository contains the dataset and R code used for our study on the lag time between state-level policy interventions and changepoints in COVID-l9 outcomes in the United States.
## File description:
1. "**analysis_GA-RWD_cases_v1-gh.R**" is the R code for our changepoint method to COVID-19 daily confirmed cases using a stepwise drift random walk model.
2. "**analysis_GA-RWD_deaths_v1-gh.R**" is the R code for our changepoint method to COVID-19 daily deaths using a stepwise drift random walk model.
3. "**lib_ga-RWD_v1-gh.R**" is the R code that contains all required functions to run GA changepoint method and stepwise drift random walk model.
4. "**modeling_random_variables.R**" is the R code that performs tests for lag times using the log-normal distribution and compares the log-normal model with other models.
5. "**random_variable_confirmed_cases.xlsx**" contains the values of the lag time random variable (Y) for confirmed cases.
6. "**random_variable_deaths.xlsx**" contains the values of the lag time random variable (Y) for deaths.
7. "**updated_policy_and_rollback_changes_till_22Feb2021.xlsx**" is a database for state-level policies during the study period: March 8, 2020 - February 28, 2021.
8. "**us-states_20210303.csv**" is the dataset of the COVID-19 United State daily confrimed cases and deaths recorded up to March 3, 2021. The data was downloaded at approximately 4:00 MST on March 6, 2021 from the New York Times COVID-19 data repository at https://github.com/nytimes/covid-19-data.
