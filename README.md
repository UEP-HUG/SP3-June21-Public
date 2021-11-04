# Seroprevalence of anti-SARS-CoV-2 antibodies 6 months into the vaccination campaign in Geneva, Switzerland, 1 June to 7 July 2021

This repo includes R and Stan code used in the above paper which can be found at https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2021.26.43.2100830

## Folder structure
- `scripts_specchio/{R, Stan}` contains the R and Stan scripts for running the models
- `data/processed` contains Geneva population data (and is where input data ready for analysis lives)
- `output/results` is where output from Stan and/or R is stored
- `latex` Supplementary Methods info concerning the statistical modelling

## Input data format
To run the code, the input data on line 22 of the scripts_specchio/R/02_run_Stan_* files has the following format:

```r
R> head(input_dat)
# A tibble: 6 x 6
  sex     age N_interp S_interp hh_id  vaccinated
  <chr> <int> <chr>    <chr>    <chr>  <lgl>     
1 m        42 neg      pos      000001 TRUE      
2 f        37 neg      neg      000002 FALSE     
3 m        38 neg      neg      000002 FALSE     
4 f         7 neg      neg      000002 FALSE     
5 f        65 pos      pos      000003 TRUE      
6 m        79 neg      pos      000004 TRUE    
```

and for the models with education (only adults), the same as above but with an extra column indicating their education level.
```r
R> head(input_dat)
# A tibble: 6 x 7
  sex     age N_interp S_interp hh_id  vaccinated edu      
  <chr> <int> <chr>    <chr>    <chr>  <lgl>      <fct>    
1 m        42 neg      pos      000001 TRUE       Tertiary 
2 f        37 neg      neg      000002 FALSE      Secondary
3 m        38 neg      neg      000002 FALSE      Tertiary 
4 f        65 pos      pos      000003 TRUE       Secondary
5 m        79 neg      pos      000004 TRUE       Mandatory 
6 m        75 neg      pos      000004 TRUE       Secondary
```

## To run the models
Running through scripts_specchio/R/02_run_Stan_multinomial_SP3_pred.R will load the relevant helper functions file, as well as the corresponding Stan script, load the data and do some further data processing, assign the different tests' validation data, read in the population data, call Stan for inference, then process the results into a final table.

The same for scripts_specchio/R/02_run_Stan_multinomial_SP3_pred_with_edu.R
