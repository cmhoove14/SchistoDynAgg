# -----------------------------
# SchistoDynAgg Process ZEST data for analysis
# Chris Hoover
# -----------------------------

library(tidyverse)
devtools::load_all()

# Convert individual level datasets to community level summaries --------------------

# Load raw data and convert to rds
chld <- haven::read_dta(here::here("Data", "Raw", "Merged_datasets_Schools_Unguja_Pemba 2012-2017_v30_FINAL_plus_MDA_GPS_short_(7.2.2018)-CHooverv2.dta"))

  saveRDS(chld, file = here::here("Data", "Raw", "children_Unguja_Pemba_2012_2017.rds"))

adlt <- read_csv(here::here("Data", "Raw", "Merged_datasets_Adults_Unguja_Pemba_2012-2017_v15_FINAL_plus_MDA_short_(7.2.2018)-CHoover.csv"),
                 guess_max = 50000)
  saveRDS(adlt, file = here::here("Data", "Raw", "adults_Unguja_Pemba_2012_2017.rds"))

# Get rid of unneeded variables
chld_reduce <- chld %>% 
  dplyr::select(c("isl", "intervention", "year", "shehia_of_school",
           "id_num", "sex", "age", "vis_haem", "haemastix_haem",
           "uf_sh", "present_at_school_mda", "percent_total_pop_rxed_r1", "percent_total_pop_rxed_r2")) 

adlt_reduce <- adlt %>% 
  dplyr::select(c("Isl", "Intervention", "Year", "Shehia",
           "ID_num", "Sex", "Age", "Vis_haem", "Haemastix_haem",
           "UF_Sh", "Recpzq", "Percent_Total_Pop_Rxed_R1", "Percent_Total_Pop_Rxed_R2"))

colnames(chld_reduce) <- colnames(adlt_reduce)

zest_all <- chld_reduce %>% 
  bind_rows(adlt_reduce) %>% 
  mutate(age_group = cut(Age, breaks = c(0,5,8,12,15,20,100)))

comm_sums <- zest_all %>% 
  group_by(Isl, Shehia, Intervention, Year) %>% 
    summarise(n_ppl = n(),
              n_chld = sum(Age <= 14, na.rm = T),
              n_adlt = sum(Age > 14, na.rm = T),
              haemastix_test = sum(!is.na(Haemastix_haem)),
              haemastix_pos = sum(Haemastix_haem > 0, na.rm = T),
              haemastix_prev = haemastix_pos/haemastix_test,
              UF_test = sum(!is.na(UF_Sh)),
              UF_pos = sum(UF_Sh > 0, na.rm = T),
              UF_neg = sum(UF_Sh == 0, na.rm = T),
              UF_max = max(UF_Sh, na.rm = T),
              UF_prev = UF_pos/UF_test,
              UF_mean = mean(UF_Sh, na.rm = T),
              UF_var = var(UF_Sh, na.rm = T),
              UF_alpha_mle = mle_alpha(UF_Sh),
              UF_alpha_mle_se = mle_alpha_se(UF_Sh),
              UF_disp_mle = mle_disp_par(UF_Sh),
              UF_disp_mle_sd = mle_disp_par_sd(UF_Sh),
              UF_mean_mle = mle_mean_par(UF_Sh),
              UF_mean_mle_sd = mle_mean_par_sd(UF_Sh),
              UF_disp2 = if_else(UF_mean < UF_var | UF_mean == 0, mean_var_agg_par(UF_Sh), NA_real_),
              UF_disp = if_else(UF_mean < UF_var | UF_mean == 0, mean_var_agg_crctd(UF_Sh), NA_real_))

saveRDS(comm_sums, file = here::here("Data", "Derived", "adults&chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds"))

chld_sums <- chld_reduce %>% 
  group_by(Isl, Shehia, Intervention, Year) %>% 
    summarise(n_ppl = n(),
              haemastix_test = sum(!is.na(Haemastix_haem)),
              haemastix_pos = sum(Haemastix_haem > 0, na.rm = T),
              haemastix_prev = haemastix_pos/haemastix_test,
              UF_test = sum(!is.na(UF_Sh)),
              UF_pos = sum(UF_Sh > 0, na.rm = T),
              UF_neg = sum(UF_Sh == 0, na.rm = T),
              UF_max = max(UF_Sh, na.rm = T),
              UF_prev = UF_pos/UF_test,
              UF_mean = mean(UF_Sh, na.rm = T),
              UF_var = var(UF_Sh, na.rm = T),
              UF_alpha_mle = mle_alpha(UF_Sh),
              UF_alpha_mle_se = mle_alpha_se(UF_Sh),
              UF_disp_mle = mle_disp_par(UF_Sh),
              UF_disp_mle_sd = mle_disp_par_sd(UF_Sh),
              UF_mean_mle = mle_mean_par(UF_Sh),
              UF_mean_mle_sd = mle_mean_par_sd(UF_Sh),
              UF_disp2 = if_else(UF_mean < UF_var | UF_mean == 0, mean_var_agg_par(UF_Sh), NA_real_),
              UF_disp = if_else(UF_mean < UF_var | UF_mean == 0, mean_var_agg_crctd(UF_Sh), NA_real_),
              n_asked_pzq = sum(Recpzq != ""),
              pct_asked = n_asked_pzq/n_ppl,
              n_recpzq = sum(Recpzq == 1),
              cvrg_pzq = n_recpzq / n_asked_pzq,
              cvrg_R1 = first(Percent_Total_Pop_Rxed_R1),
              cvrg_R2 = first(Percent_Total_Pop_Rxed_R2),
              cvrg_R12 = mean(c(cvrg_R1, cvrg_R2), na.rm = T))

saveRDS(chld_sums, file = here::here("Data", "Derived", "chldrn_shehia_sums_Unguja_Pemba_2012_2017.rds"))

adlt_sums <- adlt_reduce %>% 
  group_by(Isl, Shehia, Intervention, Year) %>% 
    summarise(n_ppl = n(),
              haemastix_test = sum(!is.na(Haemastix_haem)),
              haemastix_pos = sum(Haemastix_haem > 0, na.rm = T),
              haemastix_prev = haemastix_pos/haemastix_test,
              UF_test = sum(!is.na(UF_Sh)),
              UF_pos = sum(UF_Sh > 0, na.rm = T),
              UF_neg = sum(UF_Sh == 0, na.rm = T),
              UF_max = max(UF_Sh, na.rm = T),
              UF_prev = UF_pos/UF_test,
              UF_mean = mean(UF_Sh, na.rm = T),
              UF_var = var(UF_Sh, na.rm = T),
              UF_alpha_mle = mle_alpha(UF_Sh),
              UF_alpha_mle_se = mle_alpha_se(UF_Sh),
              UF_disp_mle = mle_disp_par(UF_Sh),
              UF_disp_mle_sd = mle_disp_par_sd(UF_Sh),
              UF_mean_mle = mle_mean_par(UF_Sh),
              UF_mean_mle_sd = mle_mean_par_sd(UF_Sh),
              UF_disp2 = if_else(UF_mean < UF_var | UF_mean == 0, mean_var_agg_par(UF_Sh), NA_real_),
              UF_disp = if_else(UF_mean < UF_var | UF_mean == 0, mean_var_agg_crctd(UF_Sh), NA_real_),
              n_askpzq = sum(!is.na(Recpzq)),
              n_recpzq = sum(Recpzq == "Yes"),
              pzq_cvrg = n_recpzq/n_askpzq,
              cvrg_R1 = first(Percent_Total_Pop_Rxed_R1),
              cvrg_R2 = first(Percent_Total_Pop_Rxed_R2),
              cvrg_R12 = mean(c(cvrg_R1, cvrg_R2), na.rm = T))

saveRDS(adlt_sums, file = here::here("Data", "Derived", "adults_shehia_sums_Unguja_Pemba_2012_2017.rds"))







# Shehia level summaries from individual data------------------
chld_pars <- chld %>% 
  group_by(isl, shehia_of_school, intervention, year) %>% 
    summarise(C_date = first(date_new),
              C_obs = n(),
              C_samps = sum(!is.na(uf_sh)),
              C_haems = sum(!is.na(haemastix_haem)),
              C_haem_pos = sum(haemastix_haem > 0, na.rm = T),
              C_haem_prev = C_haem_pos/C_haems,
              C_UFs = sum(!is.na(uf_sh)),
              C_UF_pos = sum(uf_sh > 0, na.rm = T),
              C_UF_neg = sum(uf_sh == 0, na.rm = T),
              C_UF_prev = C_UF_pos/C_UFs,
              C_UF_mean = mean(uf_sh, na.rm = T),
              C_UF_var = var(uf_sh, na.rm = T),
              C_UF_alpha = mle_alpha(uf_sh),
              C_UF_alpha_se = mle_alpha_se(uf_sh),
              Census_Pop2012 = first(total_pop_census_2012),
            #Round 1 treatment summaries  
              MDA_Date_Rx1 = first(date_cwt_r1),
              CWT_Pop_Rx1 = first(as.numeric(total_round_1)),
              CWT_Elig_Rx1 = first(as.numeric(eligible_round_1)),
              CWT_Treated_Rx1 = first(as.numeric(treated_round_1)),
              CWT_Cvrg_Rx1 = CWT_Treated_Rx1/CWT_Pop_Rx1,
              CWT_RprtCvrg_Rx1 = first(percent_total_pop_rxed_r1),
              SBT_Pop_Rx1 = first(as.numeric(f_registered_r1)+as.numeric(m_registered_r1)),
              SBT_Treated_Rx1 = first(as.numeric(f_treated_r1)+as.numeric(m_treated_r1)),
              SBT_Cvrg_Rx1 = SBT_Treated_Rx1/SBT_Pop_Rx1,
              SBT_RprtCvrg_Rx1 = first(percent_sac_rxed_r1),
            #Round2 treatment summaries  
              MDA_Date_Rx2 = first(date_cwt_r2),
              CWT_Pop_Rx2 = first(as.numeric(total_round_2)),
              CWT_Elig_Rx2 = first(as.numeric(eligible_round_2)),
              CWT_Treated_Rx2 = first(as.numeric(treated_round_2)),
              CWT_Cvrg_Rx2 = CWT_Treated_Rx2/CWT_Pop_Rx2,
              CWT_RprtCvrg_Rx2 = first(percent_total_pop_rxed_r2),
              SBT_Pop_Rx2 = first(as.numeric(f_registered_r2)+as.numeric(m_registered_r2)),
              SBT_Treated_Rx2 = first(as.numeric(f_treated_r2)+as.numeric(m_treated_r2)),
              SBT_Cvrg_Rx2 = SBT_Treated_Rx2/SBT_Pop_Rx2,
              SBT_RprtCvrg_Rx2 = first(percent_sac_rxed_r2))

adlt_pars <- adlt %>% 
  mutate(Year = if_else(Year == 2011, 2012, Year)) %>%  # Change to 2012 to match baseline year of children
  group_by(Isl, Shehia, Intervention, Year) %>% 
    summarise(A_date = first(Date_new),
              A_obs = n(),
              A_samps = sum(!is.na(UF_Sh)),
              A_haems = sum(!is.na(Haemastix_haem)),
              A_haem_pos = sum(Haemastix_haem > 0, na.rm = T),
              A_haem_prev = A_haem_pos/A_haems,
              A_UFs = sum(!is.na(UF_Sh)),
              A_UF_pos = sum(UF_Sh > 0, na.rm = T),
              A_UF_neg = sum(UF_Sh == 0, na.rm = T),
              A_UF_prev = A_UF_pos/A_UFs,
              A_UF_mean = mean(UF_Sh, na.rm = T),
              A_UF_var = var(UF_Sh, na.rm = T),
              A_UF_alpha = mle_alpha(UF_Sh),
              A_UF_alpha_se = mle_alpha_se(UF_Sh))

all_pars <- adlt_pars %>% 
  full_join(chld_pars,
            by = c("Isl" = "isl", "Shehia" = "shehia_of_school", "Intervention" = "intervention", "Year" = "year")) %>% 
  mutate(A_date = as.Date(A_date, format = "%m/%d/%Y"))

saveRDS(all_pars, here::here("Data", "Derived", "Obs_to_model_pars.rds"))

# Island-intervention level summaries ------------------
chld_int <- chld %>% 
  group_by(isl, intervention, year) %>% 
    summarise(C_date = mean(date_new),
              C_obs = n(),
              C_samps = sum(!is.na(uf_sh)),
              C_haems = sum(!is.na(haemastix_haem)),
              C_haem_pos = sum(haemastix_haem > 0, na.rm = T),
              C_haem_prev = C_haem_pos/C_haems,
              C_UFs = sum(!is.na(uf_sh)),
              C_UF_pos = sum(uf_sh > 0, na.rm = T),
              C_UF_neg = sum(uf_sh == 0, na.rm = T),
              C_UF_prev = C_UF_pos/C_UFs,
              C_UF_mean = mean(uf_sh, na.rm = T),
              C_UF_var = var(uf_sh, na.rm = T),
              C_UF_alpha = mle_alpha(uf_sh),
              C_UF_alpha_se = mle_alpha_se(uf_sh),
              Census_Pop2012 = sum(as.numeric(unique(total_pop_census_2012))),
            #Round 1 treatment summaries  
              MDA_Date_Rx1 = mean(date_cwt_r1, na.rm = T),
              CWT_Pop_Rx1 = sum(as.numeric(unique(total_round_1)), na.rm = T),
              CWT_Elig_Rx1 = sum(as.numeric(unique(eligible_round_1)), na.rm = T),
              CWT_Treated_Rx1 = sum(as.numeric(unique(treated_round_1)), na.rm = T),
              CWT_Cvrg_Rx1 = CWT_Treated_Rx1/CWT_Pop_Rx1,
              CWT_EligCvrg_Rx1 = CWT_Treated_Rx1/CWT_Elig_Rx1,
              SBT_Pop_Rx1 = sum(as.numeric(unique(f_registered_r1))+as.numeric(unique(m_registered_r1)), na.rm = T),
              SBT_Treated_Rx1 = sum(as.numeric(unique(f_treated_r1))+as.numeric(unique(m_treated_r1)), na.rm = T),
              SBT_Cvrg_Rx1 = SBT_Treated_Rx1/SBT_Pop_Rx1,
            #Round2 treatment summaries  
              MDA_Date_Rx2 = mean(date_cwt_r2, na.rm = T),
              CWT_Pop_Rx2 = sum(as.numeric(unique(total_round_2)), na.rm = T),
              CWT_Elig_Rx2 = sum(as.numeric(unique(eligible_round_2)), na.rm = T),
              CWT_Treated_Rx2 = sum(as.numeric(unique(treated_round_2)), na.rm = T),
              CWT_Cvrg_Rx2 = CWT_Treated_Rx2/CWT_Pop_Rx2,
              CWT_EligCvrg_Rx2 = CWT_Treated_Rx2/CWT_Elig_Rx2,
              SBT_Pop_Rx2 = sum(as.numeric(unique(f_registered_r2))+as.numeric(unique(m_registered_r2)), na.rm = T),
              SBT_Treated_Rx2 = sum(as.numeric(unique(f_treated_r2))+as.numeric(unique(m_treated_r2)), na.rm = T),
              SBT_Cvrg_Rx2 = SBT_Treated_Rx2/SBT_Pop_Rx2)

# 2014 community populations for round 2 were missing for some reason, so assume same as round 1
chld_int$CWT_Pop_Rx2[which(chld_int$CWT_Pop_Rx2 == 0)] <- chld_int$CWT_Pop_Rx1[which(chld_int$CWT_Pop_Rx2 == 0)]
chld_int$CWT_Cvrg_Rx2 <- chld_int$CWT_Treated_Rx2/chld_int$CWT_Pop_Rx2

# Same for adults
adlt_int <- adlt %>% 
  mutate(Year = if_else(Year == 2011, 2012, Year)) %>%  # Change to 2012 to match baseline year of children
  group_by(Isl, Intervention, Year) %>% 
    summarise(A_date = first(Date_new),
              A_obs = n(),
              A_samps = sum(!is.na(UF_Sh)),
              A_haems = sum(!is.na(Haemastix_haem)),
              A_haem_pos = sum(Haemastix_haem > 0, na.rm = T),
              A_haem_prev = A_haem_pos/A_haems,
              A_UFs = sum(!is.na(UF_Sh)),
              A_UF_pos = sum(UF_Sh > 0, na.rm = T),
              A_UF_neg = sum(UF_Sh == 0, na.rm = T),
              A_UF_prev = A_UF_pos/A_UFs,
              A_UF_mean = mean(UF_Sh, na.rm = T),
              A_UF_var = var(UF_Sh, na.rm = T),
              A_UF_alpha = mle_alpha(UF_Sh),
              A_UF_alpha_se = mle_alpha_se(UF_Sh))

all_int <- adlt_int %>% 
  full_join(chld_int,
            by = c("Isl" = "isl","Intervention" = "intervention", "Year" = "year")) %>% 
  mutate(A_date = as.Date(A_date, format = "%m/%d/%Y"))

saveRDS(all_int, here::here("Data", "Derived", "Int_Isl_Year_Sums.rds"))
