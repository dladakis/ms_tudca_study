library(tidyverse)
library(readxl)
library(anytime)
library(nlme)
library(janitor)
library(rlang)
library(ggbeeswarm)
library(DescTools)
library(patchwork)
library(gt)
library(gtsummary)
library(effects)
library(vegan)
library(microbiome)
library(NBZIMM)
library(pheatmap)
library(RColorBrewer)
library(pairwiseAdonis)
library(viridis)
set.seed(10271996)
theme_gtsummary_journal("jama")

write_csv(final %>% select(-c(full_name, dob)), "./ba_paper/submission documents/code & data/final_ba_df.csv")
df <- read.csv("final_ba_df.csv")

#Info for flow chart----
length(unique((df %>% filter(visit == 1))$baid))
table(df$visit)

baseline_all_59 <- df %>% arrange(baid, visit) %>% filter(visit == first(visit))
sum(is.na(baseline_all_59$tx_group)); (baseline_all_59 %>% filter(is.na(tx_group)))$baid #people that were not randomized,BA11 Excluded (antibiotics), rest Unknown
randomized_54 <- df %>% filter(!is.na(tx_group))
table((randomized_54 %>% filter(visit == 1))$tx_group)
randomized_no_fup <- randomized_54 %>% filter(n() == 1)
table(randomized_no_fup$tx_group); table(randomized_no_fup$event_name)
randomized_fup <- randomized_54 %>% filter(!(baid %in% randomized_no_fup$baid))
randomized_3rd <- randomized_fup %>% filter(event_name=="End of Study Visit") %>% distinct(baid, .keep_all = T)
table(randomized_3rd$tx_group)
randomized_2_total_visits <- randomized_fup %>% filter(n() == 2 & any(event_name == "Mid-Study Visit")) 
table(randomized_2_total_visits$event_name); table((randomized_2_total_visits %>% distinct(baid, .keep_all = T))$tx_group)
drop_outs <- randomized_no_fup$baid
drop_outs_2 <- (randomized_2_total_visits %>% distinct(baid, .keep_all = T))$baid




table((ba_levels %>% filter(analyte == "Tauroursodeoxycholic Acid"))$visit)
table(df$visit)
ba_levels <- ba_levels %>% mutate(unique_id = paste0(baid, "_", visit)) 
df <- df %>% mutate(unique_id = paste0(baid, "_", visit))  


randomized_54 <- randomized_54 %>% mutate(first_sx_year = if_else(baid == "BA046" & visit == 1, 2005, first_sx_year),
                                          dx_year = if_else(baid == "BA046" & visit == 1, 2015, dx_year))

randomized_54_base <- randomized_54 %>% filter(visit == 1) %>% mutate(symptoms_dur = as.numeric(format(visit_date, "%Y")) - first_sx_year,
                                                                      disease_dur = as.numeric(format(visit_date, "%Y")) - dx_year) %>% 
        left_join(ba_levels %>% filter(analyte == "Tauroursodeoxycholic Acid") %>% select(baid, visit, analyte, result) %>% 
                          rename(tudca_levels = result))
randomized_54_base <- randomized_54_base %>% mutate(ms_type = ifelse(full_name == "Terrence McCarthy", "PPMS", ms_type))

#Table 1 ----
table_1 <- randomized_54_base %>% ungroup() %>%
        dplyr::select(gender, race, age, tx_group, ms_type, current_ms_dmt,edss, symptoms_dur, disease_dur, tudca_levels) %>%
        gtsummary::tbl_summary(by = tx_group,
                               missing = "ifany",
                               statistic = list(age ~ "{mean} ({SD})",
                                                c(edss, symptoms_dur, disease_dur, tudca_levels) ~ "{median} ({p25} to {p75})"),
                               digits = c(age,edss,tudca_levels) ~ 1) %>%
        add_p(test =list(age ~ "t.test",
                         c(gender, race, ms_type) ~ "chisq.test",
                         current_ms_dmt ~ "fisher.test",
                         c(edss, symptoms_dur, disease_dur, tudca_levels) ~ "wilcox.test"),
              pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
        bold_p() %>%
        bold_labels()%>%
        as_flex_table() %>%
        flextable::save_as_docx(path = "C:/Users/dladaki1/OneDrive - Johns Hopkins/Documents/Bile Acid Local/ba_analysis/table_1_final.docx")


df <- df %>% left_join((ba_levels %>% filter(analyte == "Tauroursodeoxycholic Acid"))[c("baid", "visit", "analyte", "result")]) %>%
        rename(tudca_levels = result)

#NutritionQuest
food_q_year <- read.csv("./Food Questionnaires/Block 2014 - past year/Bhargava1_B14-Year-raw-data-and-results.csv", na = "M")
replace_missing <- function(x) {str_replace_all(x, "M", "")}
food_q_year <- apply(food_q_year, 2, replace_missing)
food_q_year[food_q_year == ""] <- NA
food_q_year <- apply(food_q_year, 2, FUN = as.numeric)
food_q_year <- food_q_year %>% clean_names() %>% as.data.frame() %>% select(respondentid, f_total, f_whole, v_total, v_drkgr, v_legumes, g_whole, g_refined, pf_total, pf_legumes, pf_seafd_hi, pf_seafd_low, pf_soy, pf_nutsds, d_total, add_sugars, adsug_na, dt_sfat, dt_mfat, dt_pfat, dt_sodi, dt_kcal)
food_q_year <- food_q_year %>% mutate(total_fruits = f_total*1000/dt_kcal,
                                      whole_fruits = f_whole*1000/dt_kcal,
                                      total_vegetables = (v_total + v_legumes)*1000/dt_kcal,
                                      greens_and_beans = (v_drkgr + v_legumes)*1000/dt_kcal,
                                      whole_grains = g_whole*1000/dt_kcal,
                                      dairy = d_total*1000/dt_kcal,
                                      total_protein = (pf_total + pf_legumes)*1000/dt_kcal,
                                      seafood_plant_protein = (pf_seafd_hi + pf_seafd_low + pf_nutsds + pf_soy + pf_legumes)*1000/dt_kcal,
                                      fas = (dt_mfat + dt_pfat)/dt_sfat,
                                      refined_grains = g_refined*1000/dt_kcal,
                                      sodium = dt_sodi/dt_kcal, #sodium was in mg so we also divide by 1000 
                                      added_sugars = add_sugars*16*100/dt_kcal,
                                      sat_fats = dt_sfat*9*100/dt_kcal,
                                      hei_tf = ifelse(total_fruits >= 0.8, 5, total_fruits*5/0.8),
                                      hei_wf = ifelse(whole_fruits >= 0.4, 5, whole_fruits*5/0.4),
                                      hei_tv = ifelse(total_vegetables >= 1.1, 5, total_vegetables*5/1.1),
                                      hei_gb = ifelse(greens_and_beans >= 0.2, 5, greens_and_beans*5/0.2),
                                      hei_wg = ifelse(whole_grains >= 1.5, 10, whole_grains*10/1.5),
                                      hei_d = ifelse(dairy >= 1.3, 10, dairy*10/1.3),
                                      hei_tp = ifelse(total_protein >= 2.5, 5, total_protein*2),
                                      hei_spp = ifelse(seafood_plant_protein >= 0.8, 5, seafood_plant_protein*5/0.8),
                                      hei_fa = ifelse(fas >= 2.5, 10, 
                                                      ifelse(fas <= 1.2, 0, fas*10/(2.5-1.2))),
                                      hei_rg = ifelse(refined_grains <= 1.8, 10, 
                                                      ifelse(refined_grains >= 4.3, 0, (4.3-refined_grains)*10/(4.3-1.8))),
                                      hei_Na = ifelse(sodium <= 1.1, 10, 
                                                      ifelse(sodium >= 2, 0, (2-sodium)*10/(2-1.1))),
                                      hei_sug = ifelse(added_sugars <= 6.5, 10, 
                                                       ifelse(added_sugars >= 26, 0, (26-added_sugars)*10/(26-6.5))),
                                      hei_satf = ifelse(sat_fats <= 8, 10, 
                                                        ifelse(sat_fats >= 16, 0, (16-sat_fats)*10/(16-8))))
food_q_year <- food_q_year %>% rowwise() %>% mutate(hei = sum(c_across(hei_tf:hei_satf)))
food_q_year <- food_q_year %>% mutate(baid = paste0("BA", str_pad(str_replace(respondentid, "[:digit:]$", ""), width = 3, side = "left", pad = "0")),
                                      visit = str_extract(respondentid, "[:digit:]$")) %>%
        relocate(baid, visit, .after = 1) %>% 
        mutate(duration = "year") %>%
        relocate(duration, .after = visit)

food_q_month <- read.csv("./Food Questionnaires/Block 2014 - past month/Bhargava1_B14-month-raw-data-and-results.csv", na = "M")
food_q_month <- apply(food_q_month, 2, replace_missing)
food_q_month[food_q_month == ""] <- NA
food_q_month <- apply(food_q_month, 2, FUN = as.numeric)
food_q_month <- food_q_month %>% clean_names() %>% as.data.frame() %>% select(respondentid, f_total, f_whole, v_total, v_drkgr, v_legumes, g_whole, g_refined, pf_total, pf_legumes, pf_seafd_hi, pf_seafd_low, pf_soy, pf_nutsds, d_total, add_sugars, adsug_na, dt_sfat, dt_mfat, dt_pfat, dt_sodi, dt_kcal)
food_q_month <- food_q_month %>% mutate(total_fruits = f_total*1000/dt_kcal,
                                        whole_fruits = f_whole*1000/dt_kcal,
                                        total_vegetables = (v_total + v_legumes)*1000/dt_kcal,
                                        greens_and_beans = (v_drkgr + v_legumes)*1000/dt_kcal,
                                        whole_grains = g_whole*1000/dt_kcal,
                                        dairy = d_total*1000/dt_kcal,
                                        total_protein = (pf_total + pf_legumes)*1000/dt_kcal,
                                        seafood_plant_protein = (pf_seafd_hi + pf_seafd_low + pf_nutsds + pf_soy + pf_legumes)*1000/dt_kcal,
                                        fas = (dt_mfat + dt_pfat)/dt_sfat,
                                        refined_grains = g_refined*1000/dt_kcal,
                                        sodium = dt_sodi/dt_kcal, #sodium was in mg so we also divide by 1000 
                                        added_sugars = add_sugars*16*100/dt_kcal,
                                        sat_fats = dt_sfat*9*100/dt_kcal,
                                        hei_tf = ifelse(total_fruits >= 0.8, 5, total_fruits*5/0.8),
                                        hei_wf = ifelse(whole_fruits >= 0.4, 5, whole_fruits*5/0.4),
                                        hei_tv = ifelse(total_vegetables >= 1.1, 5, total_vegetables*5/1.1),
                                        hei_gb = ifelse(greens_and_beans >= 0.2, 5, greens_and_beans*5/0.2),
                                        hei_wg = ifelse(whole_grains >= 1.5, 10, whole_grains*10/1.5),
                                        hei_d = ifelse(dairy >= 1.3, 10, dairy*10/1.3),
                                        hei_tp = ifelse(total_protein >= 2.5, 5, total_protein*2),
                                        hei_spp = ifelse(seafood_plant_protein >= 0.8, 5, seafood_plant_protein*5/0.8),
                                        hei_fa = ifelse(fas >= 2.5, 10, 
                                                        ifelse(fas <= 1.2, 0, fas*10/(2.5-1.2))),
                                        hei_rg = ifelse(refined_grains <= 1.8, 10, 
                                                        ifelse(refined_grains >= 4.3, 0, (4.3-refined_grains)*10/(4.3-1.8))),
                                        hei_Na = ifelse(sodium <= 1.1, 10, 
                                                        ifelse(sodium >= 2, 0, (2-sodium)*10/(2-1.1))),
                                        hei_sug = ifelse(added_sugars <= 6.5, 10, 
                                                         ifelse(added_sugars >= 26, 0, (26-added_sugars)*10/(26-6.5))),
                                        hei_satf = ifelse(sat_fats <= 8, 10, 
                                                          ifelse(sat_fats >= 16, 0, (16-sat_fats)*10/(16-8))))
food_q_month <- food_q_month %>% rowwise() %>% mutate(hei = sum(c_across(hei_tf:hei_satf)))
food_q_month <- food_q_month %>% mutate(baid = paste0("BA", str_pad(str_replace(respondentid, "[:digit:]$", ""), width = 3, side = "left", pad = "0")),
                                        visit = str_extract(respondentid, "[:digit:]$")) %>%
        relocate(baid, visit, .after = 1) %>% 
        mutate(duration = "month") %>%
        relocate(duration, .after = visit)

food_q <- bind_rows(food_q_year, food_q_month) %>% arrange(baid, visit) %>% left_join(df %>% distinct(baid, .keep_all=T) %>% select(baid, tx_group)) %>% relocate(tx_group, .after = duration)
food_q$visit <- as.numeric(food_q$visit)
summary(lme(hei ~ visit*tx_group, 
            random = ~ visit | baid,
            data = food_q %>% group_by(baid) %>% filter(n() >= 2),
            method = "REML",
            na.action = na.omit))





# 
# test_base <- read.csv("C:/Users/dladaki1/OneDrive - Johns Hopkins/Documents/Bile Acid Local/ba_analysis/ba_base_df_demo.csv")
# # theme_gtsummary_journal("jama")
# #
# test_base %>% dplyr::select(gender, race, age, tx_group, ms_type, current_ms_dmt, edss, symptoms_dur, disease_dur, tudca_levels) %>%
#         gtsummary::tbl_summary(by = tx_group,
#                                missing = "ifany",
#                                statistic = list(age ~ "{mean} ({SD})",
#                                                 c(edss, symptoms_dur, disease_dur, tudca_levels) ~ "{median} ({p25} to {p75})"),
#                                digits = c(age,edss,tudca_levels) ~ 1) %>%
#         add_p(test =list(age ~ "t.test",
#                          c(gender, race, ms_type) ~ "chisq.test",
#                          current_ms_dmt ~ "fisher.test",
#                          c(edss, symptoms_dur, disease_dur, tudca_levels) ~ "wilcox.test"),
#               pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
#         bold_p() %>%
#         bold_labels() %>%
#         as_flex_table() %>%
#         flextable::save_as_docx(path = "C:/Users/dladaki1/OneDrive - Johns Hopkins/Documents/Bile Acid Local/ba_analysis/out/table_1_final.docx")
# 
# 

#Clinical variables over time ----
df <- df %>%group_by(baid) %>% rowwise() %>% mutate(t25f = mean(c(t25f_1, t25f_2), na.rm = T),
                                                          d9hp = mean(c(d9hp_1, d9hp_2), na.rm = T),
                                                          nd9hp = mean(c(nd9hp_1, nd9hp_2), na.rm = T))
df <- df %>% mutate(t25f = if_else(str_detect(non_completion_comments, "wheelchair"), max(df$t25f, na.rm = T), t25f, missing = t25f),
                          d9hp = ifelse(is.na(d9hp) & !is.na(dominant_hand), 777, d9hp),
                          nd9hp = ifelse(is.na(nd9hp) & !is.na(dominant_hand), 777, nd9hp),
                          hp_inv = (1/d9hp + 1/nd9hp)/2,
                          hp = mean(c(d9hp, nd9hp), na.rm = T))

baseline <- df %>% filter(visit == 1)


df <- df %>% group_by(baid) %>% arrange(baid, visit_date) %>% mutate(t25f_z = -(t25f - mean(baseline$t25f, na.rm = T))/sd(baseline$t25f, na.rm = T),
                                                                           pasat_z = (pasat - mean(baseline$pasat, na.rm = T))/sd(baseline$pasat, na.rm = T),
                                                                           hp_z = (hp_inv - mean(baseline$hp_inv, na.rm = T))/sd(baseline$hp_inv, na.rm = T),
                                                                           msfc = (t25f_z + hp_z + pasat_z)/3) %>%
        mutate(across(t25f_z:msfc, ~ first(.x), .names = "{.col}_base"),
               edss_base = first(edss),
               across(msqol_physical:msqol_mental_composite, ~ first(.x), .names = "{.col}_base"))

df <- df %>% left_join(ba_main[c("baid", "jhhmrn")]) %>% relocate(jhhmrn, .after = os_id) %>% arrange(baid, visit_date) %>%
        mutate(time = (time_length(difftime(visit_date, first(visit_date)), unit = "week"))/16) #time of 1 = 16 weeks

baid_3_visits <- (unique((df %>% group_by(baid) %>% filter(n()==3))$baid))
baid_2_visits <- (unique((df %>% group_by(baid) %>% filter(n()==2))$baid))


baseline <- df %>% filter(visit == 1)
baseline <- baseline %>% mutate(symptoms_dur = as.numeric(format(visit_date, "%Y")) - first_sx_year,
                                disease_dur = as.numeric(format(visit_date, "%Y")) - dx_year)


edss_ids <- unique((df %>% filter(!is.na(edss)) %>% group_by(baid) %>% filter(n() >= 2))$baid)
t25f_ids <- unique((df %>% filter(!is.na(t25f)) %>% group_by(baid) %>% filter(n() >= 2))$baid)
msqol_physical_ids <- unique((df %>% filter(!is.na(msqol_physical)) %>% group_by(baid) %>% filter(n() >= 2))$baid)
hp_ids <- unique((df %>% filter(!is.na(hp_z)) %>% group_by(baid) %>% filter(n() >= 2))$baid)
msqol_mental_ids <- unique((df %>% filter(!is.na(msqol_mental_composite)) %>% group_by(baid) %>% filter(n() >= 2))$baid)


baseline <- baseline %>% 
        mutate(dmt_class = ifelse(current_ms_dmt %in% c("Ocrevus", "Tysabri", "Rituximab"), "Infusion", current_ms_dmt))

df <- df %>% left_join(baseline[c("baid", "dmt_class")])

df <- df %>% ungroup() %>% group_by(baid) %>% arrange(baid, visit_date) %>% mutate(age_base = first(age))



adherence_distinct <- adherence_distinct %>% mutate(adh = 1-ratio_total)
df <- df %>% left_join(adherence_distinct[c("baid", "adh")])
df <- df %>% group_by(baid) %>% filter(!is.na(tx_group))
df_p <- df %>% filter(tx_group == 0)
df_tx <- df %>% filter(tx_group == 1)

##Mixed effects models ----

clinical_results <- c()
for(i in c("edss", "t25f_z", "hp_z", "pasat_z", "msfc", "msqol_physical_composite", "msqol_mental_composite")){
        df_temp <- df[c("baid", "visit", "tx_group", i, paste0(i,"_base"), "time")] %>% filter(!is.na(!!as.symbol(i))) %>%
                group_by(baid) %>% filter(n()>=2)
        formula <- as.formula(paste0(i, " ~ time*factor(tx_group) + time*",i,"_base"))
        model <- lme(formula,
                     random = ~ time | baid,
                     data = df_temp,
                     control = lmeControl(maxIter = 5000, opt = "optim"),
                     na.action = na.omit)
        coef <- cbind(outcome = i, n = model$dims$ngrps[[1]], as.data.frame(summary(model)$tTable)[5,], as.data.frame(unclass(intervals(model, which = "fixed"))$fixed)[5,])
        clinical_results <- rbind(clinical_results, coef)
}

clinical_results <- clinical_results %>% mutate(dif_beta_ci = paste0(round(Value, 2), " (", round(lower, 2), " to ", round(upper, 2), ")")) %>% rename(dif_pval = "p-value") %>% relocate(dif_pval, .after = dif_beta_ci)

clinical_results_tx <- c()
for(i in c("edss", "t25f_z", "hp_z", "pasat_z", "msfc", "msqol_physical_composite", "msqol_mental_composite")) {
        tryCatch({ 
                df_tx <- df_tx[c("baid", "visit", "gender", "age_base", "dmt_class", i, "time")] %>% filter(!is.na(!!as.symbol(i))) %>%
                        group_by(baid) %>% filter(n()>=2)
                formula <- as.formula(paste0(i, " ~ time"))
                model_tx <- lme(formula,
                                random = ~ time|baid,
                                data = df_tx,
                                na.action = na.omit,
                                control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(outcome = i, as.data.frame(summary(model_tx)$tTable)[2,], as.data.frame(unclass(intervals(model_tx, which = "fixed"))$fixed)[2,], tudca_n = length(levels(getGroups(model_tx))))
                clinical_results_tx <- rbind(clinical_results_tx, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for(i in c("edss")) {
        tryCatch({ 
                df_tx <- df_tx[c("baid", "visit", "gender", "age_base", "dmt_class", i, "time")] %>% filter(!is.na(!!as.symbol(i))) %>%
                        group_by(baid) %>% filter(n()>=2)
                formula <- as.formula(paste0(i, " ~ time"))
                model_tx <- lme(formula,
                                random = ~ 1|baid,
                                data = df_tx,
                                na.action = na.omit,
                                control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(outcome = i, as.data.frame(summary(model_tx)$tTable)[2,], as.data.frame(unclass(intervals(model_tx, which = "fixed"))$fixed)[2,], tudca_n = length(levels(getGroups(model_tx))))
                clinical_results_tx <- rbind(clinical_results_tx, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
clinical_results_tx <- clinical_results_tx %>% mutate(tudca_beta_ci = paste0(round(Value, 2), " (", round(lower, 2), " to ", round(upper, 2), ")")) %>% rename(tudca_pval = "p-value") %>% relocate(tudca_pval, .after = tudca_beta_ci)

clinical_results_p <- c()
for(i in c("edss", "t25f_z", "hp_z", "pasat_z", "msfc", "msqol_physical_composite", "msqol_mental_composite")) {
        tryCatch({ 
                df_p <- df_p[c("baid", "visit", "gender", "age_base", "dmt_class", i, "time")] %>% filter(!is.na(!!as.symbol(i))) %>%
                        group_by(baid) %>% filter(n()>=2)
                formula <- as.formula(paste0(i, " ~ time"))
                model_p <- lme(formula,
                               random = ~ time|baid,
                               data = df_p,
                               na.action = na.omit,
                               control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(outcome = i, as.data.frame(summary(model_p)$tTable)[2,], as.data.frame(unclass(intervals(model_p, which = "fixed"))$fixed)[2,], placebo_n = length(levels(getGroups(model_p))))
                clinical_results_p <- rbind(clinical_results_p, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

for(i in c("pasat_z")) {
        tryCatch({ 
                df_p <- df_p[c("baid", "visit", "gender", "age_base", "dmt_class", i, "time")] %>% filter(!is.na(!!as.symbol(i))) %>%
                        group_by(baid) %>% filter(n()>=2)
                formula <- as.formula(paste0(i, " ~ time"))
                model_p <- lme(formula,
                               random = ~ 1|baid,
                               data = df_p,
                               na.action = na.omit,
                               control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(outcome = i, as.data.frame(summary(model_p)$tTable)[2,], as.data.frame(unclass(intervals(model_p, which = "fixed"))$fixed)[2,], placebo_n = length(levels(getGroups(model_p))))
                clinical_results_p <- rbind(clinical_results_p, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

clinical_results_p <- clinical_results_p %>% mutate(placebo_beta_ci = paste0(round(Value, 2), " (", round(lower, 2), " to ", round(upper, 2), ")")) %>% rename(placebo_pval = "p-value") %>% relocate(placebo_pval, .after = placebo_beta_ci)

clinical_table <- full_join(clinical_results_p %>% select(outcome, placebo_n, placebo_beta_ci, placebo_pval), 
                            full_join(clinical_results_tx %>% select(outcome, tudca_n, tudca_beta_ci, tudca_pval), clinical_results %>% select(outcome, n, dif_beta_ci, dif_pval)))

# write_csv(clinical_table, "./out/clinical_outcomes.csv")

#AEs ----

redcap_aes <- redcap %>% filter(repeat_instrument == "Adverse Events") %>% 
        select(baid, event_name, repeat_instance, description_of_adverse_event:sae)

redcap_aes <- redcap_aes %>% left_join(ba_main[c("baid", "tx_group")]) %>%
        relocate(tx_group, .after = baid) %>%
        rename(ae = description_of_adverse_event,
               start_date = date_of_onset_if_day_not_available_put_1st_of_month,
               stop_date = when_did_it_stop_if_day_not_available_put_1st_of_month,
               ongoing = is_the_event_ongoing) %>%
        mutate(any_ae = ifelse(any(!is.na(ae)), 1, 0),
               gi_ae = if_else(any(str_detect(ae, "bdominal|Cholecystitis|Constipation|Diarrhea|GERD|flatul|digest|Nausea|Reflux")), 1, 0, missing = 0),
               other_ae = if_else(any(str_detect(ae, "Bronchitis|DVT|L facial|UTI")), 1, 0, missing = 0),
               hospitalization = if_else(any(str_detect(seriousness, "Hosp")), 1, 0, missing = 0),
               serious = if_else(any(sae == "Yes"), 1, 0, missing = 0))

redcap_aes <- redcap_aes %>% mutate(cum_aes = sum(!is.na(ae)))

redcap_aes_distinct <- redcap_aes %>% distinct(baid, .keep_all = T)
aes_fda <- redcap_aes %>% filter(!is.na(ae)) %>% arrange(baid) %>% left_join(df %>% filter(visit_date == first(visit_date)) %>% select(baid, visit_date))
aes_fda$start_date <- as.Date(aes_fda$start_date)
aes_fda <- aes_fda %>% mutate(days_post_base = time_length(difftime(start_date, visit_date), unit = "day"))
aes_fda_single <- aes_fda %>% distinct(baid, .keep_all = T) %>% select(baid, tx_group, any_ae:cum_aes) 
df_2visits <- df %>% group_by(baid) %>% filter(n() >=2)
aes_fda_single <- bind_rows(aes_fda_single, df_2visits %>% distinct(baid, .keep_all = T) %>% filter(!(baid %in% aes_fda_single$baid)) %>% select(baid, tx_group) %>%
                                    mutate(any_ae = 0,
                                           gi_ae = 0,
                                           other_ae = 0,
                                           hospitalization = 0,
                                           serious = 0,
                                           cum_aes = 0))

##Adverse events table----
# aes_fda_single %>% ungroup(baid) %>% select(-baid) %>% tbl_summary(by = tx_group) %>% add_p(all_categorical() ~ "fisher.test") %>%
#         bold_labels() %>%
#         as_flex_table() %>%
#         flextable::save_as_docx(path = "./out/adverse_events_fisher.docx")


aes_last_visit <- df %>% group_by(baid) %>% arrange(baid, visit_date) %>% mutate(years_from_baseline = time_length(difftime(visit_date, first(visit_date)), unit = "year")) %>% 
        filter(visit_date == last(visit_date))
aes_last_visit_sum <- aes_last_visit %>% 
        left_join(redcap_aes_distinct %>% select(any_ae:cum_aes)) %>% 
        ungroup() %>% 
        group_by(tx_group) %>%
        filter(!is.na(cum_aes)) %>% 
        mutate(cum_non_sae = cum_aes - serious) %>% 
        summarise(number_at_risk = n(),
                  person_yrs = sum(years_from_baseline),
                  gi_ae_sum = sum(gi_ae),
                  other_ae_sum = sum(other_ae),
                  hosp_sum = sum(hospitalization),
                  sae_sum = sum(serious),
                  cum_aes_sum = sum(cum_aes),
                  cum_non_sae_sum = sum(cum_non_sae)) %>% 
        rowwise() %>% 
        mutate(across(gi_ae_sum:cum_non_sae_sum, ~ .x*1000/person_yrs, .names = "{.col}_1000_person_yrs"))

aes_dummy <- data.frame(tx = c(rep(0,21), rep(1,26)), gi = c(rep(1,4), rep(0,17), rep(1,10), rep(0,16)), non_gi = c(rep(1,4), rep(0,17), 1,rep(0,25)), sae = c(rep(0,46), 1), any = c(rep(1,7), rep(0,14), rep(1,10), rep(0,16)))
aes_dummy <- aes_dummy %>% mutate(total = gi + non_gi)


#GFAP - nfl ----
##OLD----


plates_1 <- read_excel("nfl_gfap.xlsx", sheet = 3, na = "NaN") %>% clean_names %>% filter(!is.na(mean_conc))
plates_2 <- read_excel("nfl_gfap.xlsx", sheet = 4, na = "NaN") %>% clean_names %>% filter(!is.na(mean_conc))
plates_3 <- read_excel("nfl_gfap.xlsx", sheet = 5, na = "NaN") %>% clean_names %>% filter(!is.na(mean_conc))

plates <- bind_rows(plates_1, plates_2, plates_3)
plates_final <- plates %>% rename(visit_id = sample_barcode) %>% 
        dplyr::select(c(1:3, mean_conc, batch_name)) %>% 
        spread(key = plex, value = mean_conc) %>%
        mutate(plate = str_extract(location, "Plate [:digit:]"))

nfl_gfap <- read.csv("nfl_gfap_ba_cohort.csv") %>% clean_names() %>% left_join(plates_final[c("visit_id", "plate")]) %>% mutate(visit_id_original = visit_id) %>% mutate(visit_id = str_replace_all(visit_id, "_1", "")) %>% 
        mutate(visit_id = reduce2(c("BL", "W8", "W16"), c("V1", "V2", "V3"), .init = visit_id, str_replace))
nfl_gfap <- nfl_gfap %>% mutate(os_id = str_extract(visit_id, "^[:alnum:]{5}"),
                                visit = as.numeric(str_extract(visit_id, "[:alnum:]{1}$")))

nfl_gfap <- nfl_gfap %>% left_join(ba_main[c("os_id", "tx_group")]) %>%
        group_by(os_id) %>% filter(n() >= 2)

nfl_gfap <- nfl_gfap %>% filter(!(os_id == "BA021" & plate == "Plate 1"))  %>% left_join(ba_main[c("baid", "os_id")])


nfl_gfap <- nfl_gfap %>% arrange(os_id, visit) %>% 
        left_join((ba_levels %>% filter(analyte == "Tauroursodeoxycholic Acid"))[c("visit_id", "analyte", "result")])


nfl_gfap_same_plate <- nfl_gfap %>% ungroup() %>% group_by(os_id, plate) %>% filter(any(visit == 1)) %>% filter(n() >= 2)
nfl_gfap_dif_plate <- setdiff(nfl_gfap, nfl_gfap_same_plate) 

nfl_gfap_same_plate <- nfl_gfap_same_plate %>% filter(!os_id %in% unique(nfl_gfap_dif_plate$os_id))
nfl_gfap_dif_plate <- setdiff(nfl_gfap, nfl_gfap_same_plate) 




##Rerun samples ----
rerun <- read_excel("101122 _pavan_n2p_rerun.xlsx", sheet = 5) %>% clean_names()
rerun <- rerun %>% mutate(os_id = str_to_upper(str_extract(sample_barcode, "^[:alnum:]{5}")),
                          visit_id = str_to_upper(sample_barcode),
                          visit_id = str_to_upper(visit_id),
                          visit_id = str_replace(visit_id, " ", ""),
                          visit_id = str_replace(visit_id, "_1", ""),
                          visit_id = reduce2(c("BL", "W8", "W16"), c("V1", "V2", "V3"), .init = visit_id, str_replace),
                          visit = as.numeric(str_extract(visit_id, "[:alnum:]{1}$")),
                          plate = "Plate 4") %>%
        rename(visit_id_original = sample_barcode) %>% 
        left_join(ba_main[c("os_id", "baid", "tx_group")]) %>%
        group_by(os_id) %>% filter(n() >= 2) %>% 
        left_join((ba_levels %>% filter(analyte == "Tauroursodeoxycholic Acid"))[c("visit_id", "analyte", "result")]) %>% 
        dplyr::select(-c(mean_conc_2, sd_conc_4, cv_conc_5, mean_conc_6, sd_conc_8, cv_conc_9))

df_labs <- bind_rows(nfl_gfap_same_plate, rerun) %>% group_by(os_id) %>% arrange(os_id, visit) %>% mutate(nfl_base = first(nfl),
                                                                                                             gfap_base = first(gfap),
                                                                                                             tudca_change_2 = nth(result, 2) - first(result),
                                                                                                             tudca_change_3 = nth(result, 3) - first(result),
                                                                                                             visit_to_choose = if_else(tudca_change_2 > tudca_change_3, 2, 3, missing = 2),
                                                                                                             tudca_change = ifelse(visit_to_choose == 2, tudca_change_2, tudca_change_3),
                                                                                                             gfap_change = ifelse(visit_to_choose == 2, nth(gfap, 2) - first(gfap), nth(gfap, 3) - first(gfap)),
                                                                                                             nfl_change = ifelse(visit_to_choose == 2, nth(nfl, 2) - first(nfl), nth(nfl, 3) - first(nfl)))

df_labs <- df_labs %>% mutate(week_num = as.numeric(ifelse(visit == 1, 0,
                                                                 ifelse(visit == 2, 8, 16)))/16,
                                    visit_2 = as.numeric(ifelse(visit == 1, 0,
                                                                ifelse(visit == 2, 1, 2))))
df_labs <- df_labs %>% mutate(visit = ifelse(os_id == "BA025" & visit == 2, 3, visit)) %>% left_join(df %>% select(baid, visit, visit_date, time))

df_labs$tx_group <- as.factor(df_labs$tx_group)

##NFL ----
nfl_model <- lme(log(nfl) ~ time*tx_group + time*log(nfl_base), 
                 random = ~ time | os_id, 
                 data = df_labs,
                 control = lmeControl(opt = "optim", maxIter = 500),
                 method = "REML",
                 na.action = na.omit)

summary(nfl_model)$tTable; intervals(nfl_model, which = "fixed")

nfl <- cbind(variable = "nfl", N = length(levels(getGroups(nfl_model))),as.data.frame(summary(nfl_model)$tTable)[5,], as.data.frame(unclass(intervals(nfl_model, which = "fixed"))$fixed)[5,])

###NFL plot ----
nfl_df <- as.data.frame(Effect(c("time", "tx_group"), nfl_model, 
                               xlevels = list(time = c(0,0.5,1))))


nfl_plot <- ggplot(data = nfl_df, aes(x = time, y = fit, fill = factor(tx_group), color = factor(tx_group), shape = factor(tx_group))) + 
        geom_line(size = 1, show.legend = F, position = position_dodge(width = 0.2), alpha = 0.8)+
        geom_point(size = 4, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = fit - se, ymax = fit +se), size = 1.2, width = 0.1, alpha = 0.8, position = position_dodge(width = 0.2), show.legend = F) +
        theme_bw() + 
        scale_fill_manual(name = "Group",
                          limits = c("0", "1"),
                          labels = c("Placebo", "TUDCA"),
                          values = c("#3872b0", "#9c5338")) + 
        scale_color_manual(name = "Group",
                           limits = c("0", "1"),
                           labels = c("Placebo", "TUDCA"),
                           values = c("#3872b0", "#9c5338")) +
        scale_shape_manual(name = "Group",
                           limits = c("0", "1"),
                           labels = c("Placebo", "TUDCA"),
                           values = c(16, 15)) +
        scale_x_continuous(breaks = c(0,0.5,1),
                           labels = c("Baseline", 8, 16)) +
        ylab("log(nfl)") +
        xlab("Week") 
# annotate(geom = "label", x = 0.9, y = 12.5, label = "Beta interaction = 0.18 per 16 weeks\np=0.96" )

##GFAP ----
# df_labs <- df_labs %>% left_join(df_tudca %>% select(baid, tudca_change))

gfap_model <- lme(log(gfap) ~ time*tx_group + time*log(gfap_base), 
                  random = ~ time | os_id, 
                  data = df_labs,
                  control = lmeControl(opt = "optim", maxIter = 500),
                  method = "REML",
                  na.action = na.omit)

# gfap_change <- df_labs %>% group_by(baid) %>% arrange(baid, week_num) %>% mutate(gfap_roc = last(gfap) - first(gfap))
# gfap_change <- gfap_change %>% filter(tx_group == 1) %>% distinct(baid, .keep_all = T)
# cor.test(gfap_change$gfap_roc, gfap_change$tudca_change)
# 
# 
summary(gfap_model)$tTable; intervals(gfap_model, which = "fixed")

gfap <- cbind(variable = "gfap", N = length(levels(getGroups(gfap_model))), as.data.frame(summary(gfap_model)$tTable)[5,], as.data.frame(unclass(intervals(gfap_model, which = "fixed"))$fixed)[5,])

###GFAP plot ----
gfap_df <- as.data.frame(Effect(c("time", "tx_group"), gfap_model, 
                                xlevels = list(time = c(0,0.5,1))))



gfap_plot <- ggplot(data = gfap_df, aes(x = time, y = exp(fit), fill = factor(tx_group), color = factor(tx_group), shape = factor(tx_group))) + 
        geom_line(size = 1, show.legend = F, position = position_dodge(width = 0.2), alpha = 0.8)+
        geom_point(size = 4, position = position_dodge(width = 0.2)) +
        geom_errorbar(aes(ymin = exp(fit - se), ymax = exp(fit +se)), size = 1.2, width = 0.1, alpha = 0.8, position = position_dodge(width = 0.2), show.legend = F) +
        theme_bw() + 
        scale_fill_manual(name = "Group",
                          limits = c("0", "1"),
                          labels = c("Placebo", "TUDCA"),
                          values = c("#3872b0", "#9c5338")) + 
        scale_color_manual(name = "Group",
                           limits = c("0", "1"),
                           labels = c("Placebo", "TUDCA"),
                           values = c("#3872b0", "#9c5338")) +
        scale_shape_manual(name = "Group",
                           limits = c("0", "1"),
                           labels = c("Placebo", "TUDCA"),
                           values = c(16, 15)) +
        scale_x_continuous(breaks = c(0,0.5,1),
                           labels = c("Baseline", 8, 16)) +
        ylab("log(GFAP)") +
        xlab("Week") +
        coord_cartesian(ylim = c(4.5,5.2)) #+
# annotate(geom = "label", x = 0.85, y = 115, label = "Beta interaction = -22.0 per 16 weeks\np=0.24" )

nfl_plot + gfap_plot + plot_layout(guides = "collect") & theme(legend.position = "right")


lab_results <- c()
for(i in c("nfl", "gfap")) {
        formula <- as.formula(paste0("log(",i, ") ~ time*tx_group + time*",i,"_base"))
        model <- lme(formula,
                     random = ~ time | os_id, 
                     data = df_labs,
                     control = lmeControl(opt = "optim", maxIter = 500),
                     method = "REML",
                     na.action = na.omit)
        coef <- cbind(outcome = i, n = model$dims$ngrps[[1]], as.data.frame(summary(model)$tTable)[5,], as.data.frame(unclass(intervals(model, which = "fixed"))$fixed)[5,])
        lab_results <- rbind(lab_results, coef)
}

lab_results <- lab_results %>% mutate(dif_beta_ci = paste0(round(Value, 2), " (", round(lower, 2), " to ", round(upper, 2), ")")) %>% rename(dif_pval = "p-value") %>% relocate(dif_pval, .after = dif_beta_ci)
lab_results_percent <- lab_results %>% mutate(dif_beta_ci = paste0(round((exp(Value)-1)*100, 2), " (", round((exp(lower)-1)*100, 2), " to ", round((exp(upper)-1)*100, 2), ")")) 

lab_results_tx <- c()
for(i in c("nfl", "gfap")) {
        tryCatch({ 
                formula <- as.formula(paste0("log(",i, ") ~ week_num"))
                model_tx <- lme(formula,
                                random = ~ week_num | os_id,
                                data = df_labs %>% filter(tx_group == 1),
                                na.action = na.omit,
                                control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(outcome = i, as.data.frame(summary(model_tx)$tTable)[2,], as.data.frame(unclass(intervals(model_tx, which = "fixed"))$fixed)[2,], tudca_n = length(levels(getGroups(model_tx))))
                lab_results_tx <- rbind(lab_results_tx, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

lab_results_tx <- lab_results_tx %>% mutate(tudca_beta_ci = paste0(round(Value, 2), " (", round(lower, 2), " to ", round(upper, 2), ")")) %>% rename(tudca_pval = "p-value") %>% relocate(tudca_pval, .after = tudca_beta_ci)
lab_results_tx_percent <- lab_results_tx %>% mutate(dif_beta_ci = paste0(round((exp(Value)-1)*100, 2), " (", round((exp(lower)-1)*100, 2), " to ", round((exp(upper)-1)*100, 2), ")"))

lab_results_p <- c()
for(i in c("nfl", "gfap")) {
        tryCatch({ 
                formula <- as.formula(paste0("log(",i, ") ~ week_num"))
                model_p <- lme(formula,
                               random = ~ week_num | os_id,
                               data = df_labs %>% filter(tx_group == 0),
                               na.action = na.omit,
                               control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(outcome = i, as.data.frame(summary(model_p)$tTable)[2,], as.data.frame(unclass(intervals(model_p, which = "fixed"))$fixed)[2,], placebo_n = length(levels(getGroups(model_p))))
                lab_results_p <- rbind(lab_results_p, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
for(i in c("nfl")) {
        tryCatch({ 
                formula <- as.formula(paste0("log(",i, ") ~ week_num"))
                model_p <- lme(formula,
                               random = ~ 1 | os_id,
                               data = df_labs %>% filter(tx_group == 0),
                               na.action = na.omit,
                               control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(outcome = i, as.data.frame(summary(model_p)$tTable)[2,], as.data.frame(unclass(intervals(model_p, which = "fixed"))$fixed)[2,], placebo_n = length(levels(getGroups(model_p))))
                lab_results_p <- rbind(lab_results_p, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
lab_results_p <- lab_results_p %>% mutate(placebo_beta_ci = paste0(round(Value, 2), " (", round(lower, 2), " to ", round(upper, 2), ")")) %>% rename(placebo_pval = "p-value") %>% relocate(placebo_pval, .after = placebo_beta_ci)
lab_results_p_percent <- lab_results_p %>% mutate(dif_beta_ci = paste0(round((exp(Value)-1)*100, 2), " (", round((exp(lower)-1)*100, 2), " to ", round((exp(upper)-1)*100, 2), ")"))

lab_table <- full_join(lab_results_p %>% select(outcome, placebo_n, placebo_beta_ci, placebo_pval), 
                       full_join(lab_results_tx %>% select(outcome, tudca_n, tudca_beta_ci, tudca_pval), lab_results %>% select(outcome, n, dif_beta_ci, dif_pval)))

# write_csv(lab_table, "./out/lab_results.csv")

#Microbiota ----
stool_ids <- read_excel("P-00KA JH Bhargava_updated.xlsx", sheet = 2, na = "NA")
stool_ids$visit_date <- as.Date(stool_ids$visit_date)
stool_ids <-  stool_ids %>% 
        mutate(visit = ifelse(visit == "BL", 1,
                              ifelse(visit == "W8", 2, 3))) %>%
        mutate(visit = ifelse(visit_id == "BA06_W16", 2, visit)) %>% #mislabeled as visit 3 while date closer to visit 2
        mutate(baid = str_replace(baid, "BA", "BA0")) %>% 
        distinct(baid, visit, .keep_all = T) %>%
        left_join(df %>% select(baid, os_id, visit, visit_date, tx_group, tudca_levels), by = c("baid", "visit"), suffix = c("", "_binder"))  %>% 
        group_by(baid) %>%
        arrange(baid, visit_date) %>%
        mutate(time = (time_length(difftime(visit_date, first(visit_date)), unit = "week"))/16) %>%
        relocate(tx_group:time, .after = visit) %>%
        mutate(time = ifelse(is.na(time) & visit == 3, 1, time)) %>% 
        arrange(baid, visit)

stool_long <- stool_ids %>% group_by(baid) %>% filter(n() >=2) %>% rename(sample = stool_id)

metaphlan <- read.table("merged_abundance_table.txt") %>% data.frame() %>% row_to_names(row_number = 1)
metaphlan[2:ncol(metaphlan)] <- lapply(metaphlan[2:ncol(metaphlan)], as.numeric)
names(metaphlan) <- str_replace(names(metaphlan), "d_metagenome", "")
sum(names(metaphlan) %in% stool_long$sample)
metaphlan <- metaphlan %>% select(clade_name, all_of(stool_long$sample))
rownames(metaphlan) <- NULL
metaphlan <- metaphlan %>% column_to_rownames(var = "clade_name")
metaphlan <- metaphlan/100

str_excl = grep("\\|t__", row.names(metaphlan), invert = T)
tmp = metaphlan[str_excl,]
bacteria_metaphlan <- rownames(tmp)
sp_excl = grep("\\|s__", row.names(tmp), invert = T)
sp_ind = grep("\\|s__", row.names(tmp))
c_excl = grep("\\|c__", row.names(tmp), invert = T)

phyla = tmp[c_excl,]
p_ind = grep("\\|p__", row.names(phyla))
phyla = phyla[p_ind,]
row.names(phyla) = gsub(".*\\|", "", row.names(phyla))
colSums(phyla)[1:6]
phyla_filt = phyla[apply(phyla, 1, function(x) sum(x > 0.0001) > 0.05 * ncol(phyla)), ]
phyla_filt = as.data.frame(t(phyla_filt), check.names = F)
phyla = as.data.frame(t(phyla), check.names = F)

phyla_final <- phyla_filt %>% rownames_to_column(var = "stool_id") %>% left_join(stool_ids) %>% relocate(visit_id:visit_date_binder, .before = 1)
phyla_final <- phyla_final %>% mutate(b_f_ratio = p__Bacteroidetes/p__Firmicutes)
summary(lme(b_f_ratio ~ time*tx_group, random = ~ 1|baid, data = phyla_final))

# phyla_results <- c()
# for(i in colnames(phyla_final)[str_detect(colnames(phyla_final), "^p__")]){
#         tryCatch({
#                 formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time*tx_group"))
#                 nbzimm <- lme(formula, random = ~ 1|baid,  data=phyla_final)
#                 coef <- cbind(i, n = length(levels(getGroups(nbzimm))), as.data.frame(summary(nbzimm)$tTable)[4,], as.data.frame(unclass(intervals(nbzimm, which = "fixed"))$fixed)[4,])
#                 phyla_results <- rbind(phyla_results, coef)
#         }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }

genera = tmp[sp_excl,]
g_ind = grep("\\|g__", row.names(genera))
genera = genera[g_ind,]
row.names(genera) = gsub(".*\\|", "", row.names(genera))
colSums(genera)[1:6]
genera_filt = genera[apply(genera, 1, function(x) sum(x > 0.0001) > 0.05 * ncol(genera)), ]
genera_filt = as.data.frame(t(genera_filt), check.names = F)
genera = as.data.frame(t(genera), check.names = F)

genera_final <- genera_filt %>% rownames_to_column(var = "stool_id") %>% left_join(stool_ids) %>% relocate(visit_id:visit_date_binder, .before = 1)

# genera_results <- c()
# for(i in colnames(genera_final)[str_detect(colnames(genera_final), "^g__")]){
#         tryCatch({
#                 formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time*tx_group"))
#                 nbzimm <- lme.zig(formula, random = ~ 1|baid,  data=genera_final)
#                 coef <- cbind(i, n = length(levels(getGroups(nbzimm))), as.data.frame(summary(nbzimm)$tTable)[4,], as.data.frame(unclass(intervals(nbzimm, which = "fixed"))$fixed)[4,])
#                 genera_results <- rbind(genera_results, coef)
#         }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# 
# genera_results <- genera_results %>% 
#         dplyr::rename(p_value = `p-value`,
#                       beta = Value) %>%
#         mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
#         arrange(p_value) %>%
#         relocate(c(p_value, fdr), .after = beta)




species = tmp[sp_ind,]

#clean names
species_names <- data_frame(full_name = rownames(species), sp = gsub(".*\\|", "", row.names(species)))
species_names <- species_names %>% mutate(kingdom = str_extract(full_name, "(?<=k__)[:alnum:]+"),
                                          phylum = str_extract(full_name, "(?<=p__)[:alnum:]+"),
                                          class = str_extract(full_name, "(?<=c__)[:alnum:]+"),
                                          order = str_extract(full_name, "(?<=o__)[:alnum:]+"),
                                          family = str_extract(full_name, "(?<=f__)[:alnum:]+"),
                                          genus = str_extract(full_name, "(?<=g__)[:alnum:]+"),
                                          species = str_replace(sp, "s__", ""),
                                          species = str_replace(species, "sp_", "")) %>%
        mutate(clean_sp = ifelse(str_detect(species, "^GGB"), paste0(family, " (f) ", str_extract(species, "SGB[:graph:]+")), species),
               clean_sp = ifelse(str_detect(clean_sp, "^FGB"), paste0(order, " (o) ", str_extract(species, "SGB[:graph:]+")), clean_sp),
               clean_sp = ifelse(str_detect(clean_sp, "^OFG"), paste0(class, " (c) ", str_extract(species, "SGB[:graph:]+")), clean_sp),
               clean_sp = ifelse(str_detect(clean_sp, "^CFG"), paste0(phylum, " (p) ", str_extract(species, "SGB[:graph:]+")), clean_sp),
               clean_sp = ifelse(str_detect(clean_sp, "^Bacteria"), paste0(kingdom, " (k) ", str_extract(species, "SGB[:graph:]+")), clean_sp),
               clean_sp = str_replace_all(clean_sp, "_", " "))


row.names(species) = gsub(".*\\|", "", row.names(species))
colSums(species)[1:6]
species_filt = species[apply(species, 1, function(x) sum(x > 0.0001) > 0.05 * ncol(species)), ]
species_filt = as.data.frame(t(species_filt), check.names = F)
species = as.data.frame(t(species), check.names = F)

metadata <- stool_ids %>% left_join(food_q %>% select(baid, visit, hei)) %>% filter(stool_id %in% rownames(species)) %>% left_join(df %>% select(baid, visit, age, gender)) %>% column_to_rownames(var = "stool_id")
metadata <- metadata %>% mutate(visit_bin = ifelse(visit == 1, "base", "post"),
                                group = ifelse(tx_group == 1, "TUDCA", "Placebo"),
                                unique_group = paste0(group,"_",visit_bin))
species <- species[match(rownames(species), rownames(metadata)),]


##Shannon diversity----
alpha_div = data.frame(shannon = vegan::diversity(species, index = "shannon"))
alpha_meta <- merge(alpha_div, metadata, by = "row.names")

shannon_plot<-ggplot(alpha_meta, aes(x = factor(visit), y = shannon, group = interaction(factor(visit), factor(tx_group)))) +
        stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(width = 0.8, preserve = "total")) +
        geom_boxplot(position = position_dodge(width = 0.8, preserve = "total"), width = 0.6, outlier.shape = NA, coef = 0) +
        geom_point(position = position_dodge(width = 0.8, preserve = "total"),aes(color = factor(tx_group)), alpha = 0.7, size = 1.8) +
        theme_bw() +
        scale_x_discrete(breaks = c(1,2,3),
                         labels = c(0, 8, 16)) +
        scale_color_viridis_d(name = "Group",
                              labels = c("Placebo", "TUDCA"),
                              begin = 0.3,
                              end = 0.6,
                              option = "B") +
        theme(plot.title = element_text(hjust = 0.5, vjust = 0.1, size = 15),
              axis.title = element_text(color = "black", face = "bold", size = 15),
              axis.text = element_text(color = "black", size = 13),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 15, face = "bold"),
              legend.position = "bottom",
              panel.grid = element_blank()) +
        ylab("Shannon index") +
        xlab("Week")+ 
        ggtitle("n.s.") +
        coord_cartesian(ylim = c(2, 4.5))

pdf("shannon_boxplots_2.pdf", height = 5, width = 3)
shannon_plot
dev.off()

# write_csv(alpha_meta, "alpha_diversity_data_points.csv")

summary(lme(shannon ~ tx_group*time, 
            random = ~ time|baid,
            data = alpha_meta %>% group_by(baid) %>% filter(n()>=2),
            control = lmeControl(maxIter = 500, opt = "optim")))


##Bray-Curtis dissimilarity----
bray_curtis_d <- vegdist(species, method = "bray")
result <- adonis2(bray_curtis_d ~ unique_group + age + gender, data=metadata, method = "bray", permutations = 999)
result
pairwise.adonis2(bray_curtis_d ~ unique_group + age + gender, data = metadata)
#remove na of hei
metadata_no_na <- metadata %>% filter(!is.na(hei))
species_no_na <- species %>% rownames_to_column(var = "stool_id") %>% filter(stool_id %in% rownames(metadata_no_na)) %>% column_to_rownames(var = "stool_id")
species_no_na <- species_no_na[match(rownames(species_no_na), rownames(metadata_no_na)),]

bray_curtis_no_na <- vegdist(species_no_na, method = "bray")
result_no_na <- adonis2(bray_curtis_no_na ~ factor(tx_group)*time + hei, data=metadata_no_na, method = "bray", permutations = 999)
result_no_na


#bray curtis only base
metadata_no_na <- metadata %>% filter(unique_group %in% c("TUDCA_base", "TUDCA_post"))
species_no_na <- species %>% rownames_to_column(var = "stool_id") %>% filter(stool_id %in% rownames(metadata_no_na)) %>% column_to_rownames(var = "stool_id")
species_no_na <- species_no_na[match(rownames(species_no_na), rownames(metadata_no_na)),]

bray_curtis_no_na <- vegdist(species_no_na, method = "bray")
result_no_na <- adonis2(bray_curtis_no_na ~ factor(unique_group), data=metadata_no_na, method = "bray", permutations = 999)
result_no_na


MDS <- cmdscale(vegdist(species, method = "bray"), k = 2, eig = T, x.ret = T, list = T)
MDS$eig[MDS$eig < 0 ] = 0
# calculate the % of variance explained
percent = as.numeric(MDS$eig/sum(MDS$eig))[1:2]
pca_plot_df <- as.data.frame(MDS$points) %>% clean_names()
pca_plot_df <- base::merge(pca_plot_df, metadata, by = "row.names")
ggplot(pca_plot_df, aes(x = v1, y = v2, group = unique_group, fill = unique_group)) +
        geom_point(size = 3, alpha = 0.8, shape = 21, color = "black") +
        scale_fill_viridis_d(direction = -1,
                             option = "G",
                             name = "Visit",
                             limits = c("Placebo_base", "Placebo_post", "TUDCA_base", "TUDCA_post"),
                             labels = c("Placebo Baseline", "Placebo Follow-up (Weeks 8 & 16)", "TUDCA Baseline", "TUDCA Follow-up (Weeks 8 & 16)")) +
        theme_bw() +
        xlab(paste0("PC1 (", round(percent[1]*100, digits = 2), "%)")) +
        ylab(paste0("PC2 (", round(percent[2]*100, digits = 2), "%)")) +
        theme(axis.text = element_text(color = "black", size = 9),
              axis.title = element_text(size = 11),
              legend.title = element_text(face = "bold", size = 11),
              legend.text = element_text(size = 10),
              panel.grid = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line()) 
pdf("bray_curtis_2d.pdf", width = 6, height = 5)
ggplot(pca_plot_df, aes(x = v1, y = v2, group = unique_group, fill = unique_group)) +
        geom_point(size = 3, alpha = 0.7, shape = 21, color = "black") +
        scale_fill_viridis_d(direction = -1,
                             option = "G",
                             name = "Visit",
                             limits = c("Placebo_base", "Placebo_post", "TUDCA_base", "TUDCA_post"),
                             labels = c("Placebo Baseline", "Placebo Follow-up\n(Weeks 8 & 16)", "TUDCA Baseline", "TUDCA Follow-up\n(Weeks 8 & 16)")) +
        theme_bw() +
        xlab(paste0("PC1 (", round(percent[1]*100, digits = 2), "%)")) +
        ylab(paste0("PC2 (", round(percent[2]*100, digits = 2), "%)")) +
        theme(axis.text = element_text(color = "black", size = 13),
              axis.title = element_text(size = 15),
              legend.title = element_text(face = "bold", size = 15),
              legend.text = element_text(size = 13),
              panel.grid = element_blank(),
              panel.border = element_blank(), 
              axis.line = element_line()) 
dev.off()

# write_csv(pca_plot_df, "bray_curtis_data_points.csv")
# plot3d <- scatter3d(MDS$points[,1],-MDS$points[,2],MDS$points[,3],
#           groups=factor(metadata$unique_group),ellipsoid = F,
#           ellipsoid.alpha = 0.8,grid = FALSE, surface = FALSE,
#           axis.col = c("black", "black", "black"),
#           xlab = paste("PC1", "(", round(percent[1]*100, digits = 2), "%", ")",sep=""),
#           ylab = paste("PC2", " (",round(percent[2]*100, digits = 2), "%", ")",sep=""),
#           zlab = paste("PC3", " (",round(percent[3]*100, digits = 2), "%", ")",sep=""),
#           axis.scales = FALSE, surface.col = viridis(4, begin = 0, end = 1, direction = -1),
#           sphere.size = 0.7)

# rgl.postscript("braycurtis_3D_PCoA.pdf",fmt="pdf")
set.seed(1)
nmds <- vegan::metaMDS(bray_curtis_d, distance = "bray",add = T)
nmds_table <- nmds$points %>% 
        as.data.frame() %>% 
        rownames_to_column(var = "stool_id") %>% 
        left_join(stool_ids) %>% 
        mutate(nmds_group = ifelse(tx_group == 1 & visit == 1, "TUDCA_base",
                                   ifelse(tx_group == 1 & visit %in% c(2,3), "TUDCA_post", 
                                          ifelse(tx_group == 0 & visit == 1, "Placebo_base", "Placebo_post"))))

ggplot(nmds_table, aes(x = MDS1, y = MDS2, group = factor(nmds_group), color = factor(nmds_group), fill = factor(nmds_group))) + 
        geom_point(size = 3) + 
        stat_ellipse(geom = "polygon", alpha = 0.05) + 
        theme_bw()


##Species analysis ----
species_final <- species_filt %>% rownames_to_column(var = "stool_id") %>% left_join(stool_ids) %>% relocate(visit_id:visit_date_binder, .before = 1)

#Difference
metaphlan_results <- c()
for(i in colnames(species_final)[str_detect(colnames(species_final), "^s__")]){
        tryCatch({
                formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time*tx_group"))
                nbzimm <- lme.zig(formula, random = ~ 1|baid,  data=species_final)
                coef <- cbind(species = i, n = length(levels(getGroups(nbzimm))), as.data.frame(summary(nbzimm)$tTable)[4,], as.data.frame(unclass(intervals(nbzimm, which = "fixed"))$fixed)[4,])
                metaphlan_results <- rbind(metaphlan_results, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

metaphlan_results <- metaphlan_results %>% 
        dplyr::rename(p_value = `p-value`,
                      beta = Value) %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
        arrange(p_value) %>%
        relocate(c(p_value, fdr), .after = beta)

metaphlan_results <- left_join(metaphlan_results, species_names, by = c("species" = "sp"))

metaphlan_results <- metaphlan_results %>% 
        mutate(estimate_ci = paste0(round(beta, 3), " (", round(lower, 4), " to ", round(upper,4),")")) 

# metaphlan_results <- metaphlan_results %>% mutate(species_name_clean = str_extract())

#TUDCA
metaphlan_results_tudca_group <- c()
for(i in colnames(species_final)[str_detect(colnames(species_final), "^s__")]){
        tryCatch({
                formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time"))
                nbzimm <- lme.zig(formula, random = ~ 1|baid,  data=species_final %>% filter(tx_group==1))
                coef <- cbind(species = i, n = length(levels(getGroups(nbzimm))), as.data.frame(summary(nbzimm)$tTable)[2,], as.data.frame(unclass(intervals(nbzimm, which = "fixed"))$fixed)[2,])
                metaphlan_results_tudca_group <- rbind(metaphlan_results_tudca_group, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

metaphlan_results_tudca_group <- metaphlan_results_tudca_group %>% 
        dplyr::rename(p_value = `p-value`,
                      beta = Value) %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
        arrange(p_value) %>%
        relocate(c(p_value, fdr), .after = beta)
metaphlan_results_tudca_group <- left_join(metaphlan_results_tudca_group, species_names, by = c("species" = "sp"))

metaphlan_results_tudca_group <- metaphlan_results_tudca_group %>%
        mutate(estimate_ci_tudca = paste0(round(beta, 3), " (", round(lower, 4), " to ", round(upper,4),")")) %>%
        rename(p_value_tudca = p_value)

#placebo
metaphlan_results_placebo_group <- c()
for(i in colnames(species_final)[str_detect(colnames(species_final), "^s__")]){
        tryCatch({
                formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time"))
                nbzimm <- lme.zig(formula, random = ~ 1|baid,  data=species_final %>% filter(tx_group==0))
                coef <- cbind(species = i, n = length(levels(getGroups(nbzimm))), as.data.frame(summary(nbzimm)$tTable)[2,], as.data.frame(unclass(intervals(nbzimm, which = "fixed"))$fixed)[2,])
                metaphlan_results_placebo_group <- rbind(metaphlan_results_placebo_group, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

metaphlan_results_placebo_group <- metaphlan_results_placebo_group %>% 
        dplyr::rename(p_value = `p-value`,
                      beta = Value) %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
        arrange(p_value) %>%
        relocate(c(p_value, fdr), .after = beta) 
metaphlan_results_placebo_group <- left_join(metaphlan_results_placebo_group, species_names, by = c("species" = "sp")) 

metaphlan_results_placebo_group <- metaphlan_results_placebo_group %>% 
        mutate(estimate_ci_placebo = paste0(round(beta, 3), " (", round(lower, 4), " to ", round(upper,4),")")) %>%
        rename(p_value_placebo = p_value)

# write_csv(metaphlan_results, "./out/microbiome_species.csv")

metaphlan_results_table <- metaphlan_results %>% 
        left_join(metaphlan_results_tudca_group %>% select(clean_sp, estimate_ci_tudca, p_value_tudca)) %>%
        left_join(metaphlan_results_placebo_group %>% select(clean_sp, estimate_ci_placebo, p_value_placebo)) %>% 
        select(kingdom:clean_sp, estimate_ci_placebo, p_value_placebo, estimate_ci_tudca, p_value_tudca, estimate_ci, p_value) %>% 
        arrange(p_value)


metaphlan_results_table <- metaphlan_results_table %>% mutate(fdr = p.adjust(p_value, method = "BH"))
# write_csv(metaphlan_results_table, "species_results_supp_table.csv")

##Heatmap microbiota----
significant_sp <- (metaphlan_results %>% filter(p_value < 0.05))$species
significant_sp_all <- unique(c((metaphlan_results %>% filter(p_value < 0.05))$species,(metaphlan_results_tudca_group %>% filter(p_value < 0.05))$species))
significant_micro_dif <- metaphlan_results %>% filter(species %in% significant_sp_all) %>% select(clean_sp, beta) %>% rename(Difference = beta)
significant_micro_tudca <- metaphlan_results_tudca_group %>% filter(species %in% significant_sp_all) %>% select(clean_sp, beta) %>% rename(TUDCA = beta)
significant_micro_placebo <- metaphlan_results_placebo_group %>% filter(species %in% significant_sp_all) %>% select(clean_sp, beta) %>% rename(Placebo = beta)
significant_micro_matrix <- as.matrix(full_join(significant_micro_placebo, full_join(significant_micro_tudca, significant_micro_dif)) %>% column_to_rownames(var = "clean_sp"))
significant_micro_matrix[is.na(significant_micro_matrix)] <- 0

#extract p values and stars
significant_micro_dif_pval <- metaphlan_results %>% 
        filter(species %in% significant_sp_all) %>% 
        mutate(sign = ifelse(p_value < 0.001, "***",
                             ifelse(p_value < 0.01, "**", 
                                    ifelse(p_value < 0.05, "*", "")))) %>% 
        select(clean_sp, sign) %>% 
        rename(Difference = sign)

significant_micro_tudca_pval <- metaphlan_results_tudca_group %>% 
        filter(species %in% significant_sp_all) %>% 
        mutate(sign = ifelse(p_value < 0.001, "***",
                             ifelse(p_value < 0.01, "**", 
                                    ifelse(p_value < 0.05, "*", "")))) %>% 
        select(clean_sp, sign) %>% 
        rename(TUDCA = sign)

significant_micro_placebo_pval <- metaphlan_results_placebo_group %>%
        filter(species %in% significant_sp_all) %>% 
        mutate(sign = ifelse(p_value < 0.001, "***",
                             ifelse(p_value < 0.01, "**", 
                                    ifelse(p_value < 0.05, "*", "")))) %>% 
        select(clean_sp, sign) %>% 
        rename(Placebo = sign)

sign_pval <- full_join(significant_micro_placebo_pval, full_join(significant_micro_tudca_pval, significant_micro_dif_pval))
sign_pval[is.na(sign_pval)] <- "-"
sign_pval <- sign_pval %>% column_to_rownames("clean_sp")

beta_full_range <- c(as.vector(significant_micro_matrix[,1]), as.vector(significant_micro_matrix[,2]), as.vector(significant_micro_matrix[,3]))
min(beta_full_range[beta_full_range>0])
palette_color <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11)
palette_breaks <- c(seq(min(beta_full_range), max(beta_full_range[beta_full_range<0]), length.out = 6), 0,
                    seq(max(beta_full_range)/11, max(beta_full_range), length.out = 5))



micro_heatmap <- pheatmap::pheatmap(significant_micro_matrix, cellwidth = 23, cellheight = 15, gaps_col = 2, cluster_cols = F, color = palette_color, breaks = palette_breaks, display_numbers = as.matrix(sign_pval), fontsize_number = 15, fontsize = 12)

sign_color <- ifelse(significant_micro_matrix<= -0.0636571057 | significant_micro_matrix >= 0.0434999907, "white", "black")
sign_color <- sign_color[micro_heatmap$tree_row$order,]

micro_heamap <- pheatmap::pheatmap(significant_micro_matrix, cellwidth = 15, cellheight = 10, gaps_col = 2, cluster_cols = F, 
                                   breaks = palette_breaks, color = palette_color,  angle_col = 90,
                                   display_numbers = as.matrix(sign_pval), number_color = sign_color, fontsize_number = 10, fontsize = 10)



pdf("micro_heatmap.pdf", width = 7, height = 10)
micro_heamap
dev.off()

svg("micro_heatmap.svg", width = 7, height = 10)
micro_heamap
dev.off()


##Functional----

humann_paths <- read.csv("merged_pathabundance_relab.tsv", sep = "\t")
humann_paths <- humann_paths %>% column_to_rownames(var = "X..Pathway")
path_ind = grep("\\|", row.names(humann_paths), invert = T)
pathways = humann_paths[path_ind,]

pathways_filt = pathways[apply(pathways, 1, function(x) sum(x > 0) > 0.05 * ncol(pathways)), ]
pathways_filt = as.data.frame(t(pathways_filt), check.names = F)
pathways = as.data.frame(t(pathways), check.names = F)

pathways_final <- pathways_filt %>% rownames_to_column(var = "stool_id") %>% mutate(stool_id = str_replace(stool_id, "_R1_kneaddata_Abundance.RELAB", ""),
                                                                                    stool_id = str_replace(stool_id, "\\.", "-")) %>%
        filter(stool_id %in% stool_long$sample) %>%
        left_join(stool_ids) %>% relocate(visit_id:visit_date_binder, .before = 1)

names(pathways_final) <- make_clean_names(names(pathways_final))

pathway_names <- data.frame(path = names(pathways), path_clean = make_clean_names(names(pathways)))
pathway_names <- pathway_names %>% mutate(path = str_extract(path, "(?<=\\: ).+"))

#Difference
pathways_results <- c()
for(i in colnames(pathways_final)[15:ncol(pathways_final)]){
        tryCatch({
                formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time*tx_group"))
                lme_model <- lme(formula, random = ~ 1|baid,  data=pathways_final)
                coef <- cbind(path_clean = i, n = length(levels(getGroups(lme_model))), as.data.frame(summary(lme_model)$tTable)[4,], as.data.frame(unclass(intervals(lme_model, which = "fixed"))$fixed)[4,])
                pathways_results <- rbind(pathways_results, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

pathways_results <- pathways_results %>% 
        dplyr::rename(p_value = `p-value`,
                      beta = Value) %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
        arrange(p_value) %>%
        relocate(c(p_value, fdr), .after = beta) %>% 
        left_join(pathway_names) %>%
        relocate(path, .before = 1) %>%
        select(-path_clean)

pathways_results <- pathways_results %>% 
        mutate(estimate_ci = paste0(round(beta, 3), " (", round(lower, 4), " to ", round(upper,4),")")) 

#TUDCA
pathways_results_tudca_group <- c()
for(i in colnames(pathways_final)[15:ncol(pathways_final)]){
        tryCatch({
                formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time"))
                lme_model <- lme(formula, random = ~ 1|baid,  data=pathways_final %>% filter(tx_group == 1))
                coef <- cbind(path_clean = i, n = length(levels(getGroups(lme_model))), as.data.frame(summary(lme_model)$tTable)[2,], as.data.frame(unclass(intervals(lme_model, which = "fixed"))$fixed)[2,])
                pathways_results_tudca_group <- rbind(pathways_results_tudca_group, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

pathways_results_tudca_group <- pathways_results_tudca_group %>% 
        dplyr::rename(p_value = `p-value`,
                      beta = Value) %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
        arrange(p_value) %>%
        relocate(c(p_value, fdr), .after = beta) %>% 
        left_join(pathway_names) %>%
        relocate(path, .before = 1) %>%
        select(-path_clean)

pathways_results_tudca_group <- pathways_results_tudca_group %>%
        mutate(estimate_ci_tudca = paste0(round(beta, 3), " (", round(lower, 4), " to ", round(upper,4),")")) %>%
        rename(p_value_tudca = p_value)

#Placebo
pathways_results_placebo_group <- c()
for(i in colnames(pathways_final)[15:ncol(pathways_final)]){
        tryCatch({
                formula <- as.formula(paste0("asin(sqrt(", i, ")) ~ time"))
                lme_model <- lme(formula, random = ~ 1|baid,  data=pathways_final)
                coef <- cbind(path_clean = i, n = length(levels(getGroups(lme_model))), as.data.frame(summary(lme_model)$tTable)[2,], as.data.frame(unclass(intervals(lme_model, which = "fixed"))$fixed)[2,])
                pathways_results_placebo_group <- rbind(pathways_results_placebo_group, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

pathways_results_placebo_group <- pathways_results_placebo_group %>% 
        dplyr::rename(p_value = `p-value`,
                      beta = Value) %>%
        mutate(fdr = p.adjust(p_value, method = "BH")) %>% 
        arrange(p_value) %>%
        relocate(c(p_value, fdr), .after = beta) %>% 
        left_join(pathway_names) %>%
        relocate(path, .before = 1) %>%
        select(-path_clean)

pathways_results_placebo_group <- pathways_results_placebo_group %>% 
        mutate(estimate_ci_placebo = paste0(round(beta, 3), " (", round(lower, 4), " to ", round(upper,4),")")) %>%
        rename(p_value_placebo = p_value)

pathways_results_table <- pathways_results %>% 
        left_join(pathways_results_tudca_group %>% select(path, estimate_ci_tudca, p_value_tudca)) %>%
        left_join(pathways_results_placebo_group %>% select(path, estimate_ci_placebo, p_value_placebo)) %>% 
        select(path, estimate_ci_placebo, p_value_placebo, estimate_ci_tudca, p_value_tudca, estimate_ci, p_value) %>% 
        arrange(p_value)
pathways_results_table <- pathways_results_table %>% mutate(fdr = p.adjust(p_value, method = "BH"))
# write_csv(pathways_results_table, "pathway_results_supp_table_all.csv")

# write_csv(pathways_results_table, "./out/microbiome_pathways_all.csv")

###Heatmap functional----
significant_paths <- (pathways_results %>% filter(p_value < 0.05))$path
significant_paths_tudca <- (pathways_results_tudca_group %>% filter(p_value_tudca < 0.05))$path[!(pathways_results_tudca_group %>% filter(p_value_tudca < 0.05))$path %in% (pathways_results_placebo_group %>% filter(p_value_placebo < 0.05))$path]
significant_paths_all <- unique(c(significant_paths, significant_paths_tudca))
significant_path_dif <- pathways_results %>% filter(path %in% significant_paths_all) %>% select(path, beta) %>% rename(Difference = beta)
significant_path_tudca <- pathways_results_tudca_group %>% filter(path %in% significant_paths_all) %>% select(path, beta) %>% rename(TUDCA = beta)
significant_path_placebo <- pathways_results_placebo_group %>% filter(path %in% significant_paths_all) %>% select(path, beta) %>% rename(Placebo = beta)
significant_path_matrix <- as.matrix(full_join(significant_path_placebo, full_join(significant_path_tudca, significant_path_dif)) %>% column_to_rownames(var = "path"))
significant_path_matrix[is.na(significant_path_matrix)] <- 0

#extract p values and stars
significant_path_dif_pval <- pathways_results %>% 
        filter(path %in% significant_paths_all) %>% 
        mutate(sign = ifelse(p_value < 0.001, "***",
                             ifelse(p_value < 0.01, "**", 
                                    ifelse(p_value < 0.05, "*", "")))) %>% 
        select(path, sign) %>% 
        rename(Difference = sign)

significant_path_tudca_pval <- pathways_results_tudca_group %>% 
        filter(path %in% significant_paths_all) %>% 
        mutate(sign = ifelse(p_value_tudca < 0.001, "***",
                             ifelse(p_value_tudca < 0.01, "**", 
                                    ifelse(p_value_tudca < 0.05, "*", "")))) %>% 
        select(path, sign) %>% 
        rename(TUDCA = sign)

significant_path_placebo_pval <- pathways_results_placebo_group %>%
        filter(path %in% significant_paths_all) %>% 
        mutate(sign = ifelse(p_value_placebo < 0.001, "***",
                             ifelse(p_value_placebo < 0.01, "**", 
                                    ifelse(p_value_placebo < 0.05, "*", "")))) %>% 
        select(path, sign) %>% 
        rename(Placebo = sign)

sign_path_pval <- full_join(significant_path_placebo_pval, full_join(significant_path_tudca_pval, significant_path_dif_pval))
sign_path_pval[is.na(sign_path_pval)] <- "-"
sign_path_pval <- sign_path_pval %>% column_to_rownames("path")


path_beta_full_range <- c(as.vector(significant_path_matrix[,1]), as.vector(significant_path_matrix[,2]), as.vector(significant_path_matrix[,3]))

palette_color <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(11)
path_palette_breaks <- c(seq(min(path_beta_full_range), max(path_beta_full_range[path_beta_full_range<0]), length.out = 6), 0,
                         seq(max(path_beta_full_range)/11, max(path_beta_full_range), length.out = 5))

path_heamap <- pheatmap::pheatmap(significant_path_matrix, gaps_col = 2,cluster_cols = F, cellheight = 20, cellwidth = 20, color = palette_color, breaks = path_palette_breaks, display_numbers = as.matrix(sign_path_pval), number_color = "black", fontsize_number = 15, fontsize = 12)

path_sign_color <- ifelse(significant_path_matrix< -0.0019221512 | significant_path_matrix >= 0.0014304924, "white", "black")
path_sign_color <- path_sign_color[path_heamap$tree_row$order,]

path_heamap <- pheatmap::pheatmap(significant_path_matrix, gaps_col = 2, cluster_cols = F, cellheight = 17, cellwidth = 17, color = palette_color, breaks = path_palette_breaks,
                                  display_numbers = as.matrix(sign_path_pval), angle_col = 90, number_color = path_sign_color, fontsize_number = 13, fontsize = 10)

pdf("path_heatmap_all.pdf", height = 6, width = 7)
path_heamap
dev.off()

svg("path_heatmap_all.svg", height = 6, width = 6)
path_heamap
dev.off()

pdf("red_blue_palette.pdf", height = 2, width = 3)
display.brewer.pal(n = 11, name = "RdBu")
dev.off()

#Flow analysis ----
t_cells <- read_excel("20230327 PC151 TCell Quants 20230326 Unmixing.xlsx") %>% clean_names() #duplicates BA039
t_cells <- t_cells %>% dplyr::rename(ccr4_2_viable = ccr4_viable_2,
                                     ccr4_2_parent = ccr4_parent_2) %>% 
        mutate(os_id = paste0("BA", str_pad(ppid, pad = "0", width = 3, side = "left"))) %>%
        mutate(th_memory_parent = th_em_parent + th_cm_parent,
               tc_memory_parent = tc_em_parent + tc_cm_parent,
               th_memory_viable = th_em_viable + th_cm_viable,
               tc_memory_viable = tc_em_viable + tc_cm_viable)
t_cells_dups <- t_cells %>% select(os_id, visit) %>% mutate(unique_id = paste0(os_id, "_", visit)) %>% relocate(unique_id, .before = 1) %>% filter(duplicated(unique_id))
t_cells_dups_complete <- inner_join(t_cells, t_cells_dups) 
# write_csv(t_cells_dups_complete %>% arrange(os_id, visit), "t_cells_dups_complete.csv")

myeloid <- read_excel("20230212 PC151 Myeloid Data for Dimitrios.xlsx") %>% clean_names()
myeloid <- myeloid %>% dplyr::rename(os_id = ppid) %>% 
        filter(batch != 1) %>% #remove batch 1 per Matt
        mutate(visit = ifelse(visit == "BL", 1,
                              ifelse(visit == "W16", 3, 2)),
               os_id = str_to_upper(os_id)) 
myeloid_dups <- myeloid %>% select(os_id, visit) %>% mutate(unique_id = paste0(os_id, "_", visit)) %>% relocate(unique_id, .before = 1) %>% filter(duplicated(unique_id))
myeloid_dups_complete <- inner_join(myeloid, myeloid_dups) 
# write_csv(myeloid_dups_complete %>% arrange(os_id, visit), "myeloid_dups_complete.csv")



b_cells <- read_excel("20230416 PC151 BCell Data for Dimitrios.xlsx") %>% clean_names() %>% 
        mutate(os_id = paste0("BA", str_pad(ppid, pad = "0", width = 3, side = "left")),
               visit = ifelse(str_detect(sample, "BL"), 1,
                              ifelse(str_detect(sample, "W08"), 2 , 3))) %>% 
        relocate(os_id, visit, .before = 1)
b_cells_dups <- b_cells %>% select(os_id, visit) %>% mutate(unique_id = paste0(os_id, "_", visit)) %>% relocate(unique_id, .before = 1) %>% filter(duplicated(unique_id))
b_cells_dups_complete <- inner_join(b_cells, b_cells_dups) 
# write_csv(b_cells_dups_complete %>% arrange(os_id, visit), "b_cells_dups_complete.csv")

b_cells <- b_cells %>% mutate(b_memory = cd27_ig_d + ig_d,
                              am_b = am*memory_bcells_class_switch*ig_d/10000,
                              at_m_b = at_m*memory_bcells_class_switch*ig_d/10000,
                              im_b = im*memory_bcells_class_switch*ig_d/10000,
                              rm_b = rm*memory_bcells_class_switch*ig_d/10000,
                              asc_b = asc*ig_d/100)

t_cells_2 <- read.csv("20230606 PC151 TCell Deep Markers for Dimitrios.csv") %>% clean_names()
t_cells_2 <- t_cells_2 %>% mutate(across(.cols = str_which(names(t_cells_2), "glut1"), .fns = ~ ifelse(batch == 2, NA, .x)),
                                  os_id = paste0("BA", str_pad(ppid, pad = "0", width = 3, side = "left")))



#OS BA025 entered as 2nd visit but it is End of Study Visit
myeloid<- myeloid %>% mutate(visit = ifelse(os_id == "BA025" & visit == 2, 3, visit)) %>%
        distinct(os_id, visit, .keep_all = T) 
myeloid <- myeloid %>% distinct(os_id, visit, .keep_all = T) %>% 
        rename(myeloid_sample = sample,
               myeloid_batch = batch,
               myeloid_index = index,
               myeloid_index2 = index2) %>%
        group_by(os_id) %>% 
        arrange(os_id, visit) %>% 
        filter(n()>=2)
t_cells <- t_cells %>% mutate(visit = ifelse(os_id == "BA025" & visit == 2, 3, visit))
t_cells <- t_cells %>% distinct(os_id, visit, .keep_all = T) %>% 
        dplyr::rename(t_cell_batch = batch,
                      t_cell_index = index,
                      t_cell_sample = sample) %>%
        group_by(os_id) %>% 
        arrange(os_id, visit) %>% 
        filter(n()>=2)
b_cells <- b_cells %>% mutate(visit = ifelse(os_id == "BA025" & visit == 2, 3, visit)) %>%
        dplyr::rename(b_cell_inderx = index,
                      b_cell_sample = sample) %>%
        group_by(os_id) %>% 
        arrange(os_id, visit) %>% 
        filter(bcd == 0) %>% 
        filter(n()>=2)
t_cells_2 <- t_cells_2 %>% mutate(visit = ifelse(os_id == "BA025" & visit == 2, 3, visit)) %>% 
        distinct(os_id, visit, .keep_all = T) %>% 
        dplyr::rename(t_cell_2_batch = batch,
                      t_cell_2_index = index,
                      t_cell_2_sample = sample) %>%
        group_by(os_id) %>% 
        arrange(os_id, visit) %>% 
        filter(n()>=2)


# 
# b_cells_results <- c()
# for(i in names(b_cells %>% ungroup() %>% select(c(b_cells, naive)))) {
#         tryCatch({
#                 formula <- as.formula(paste0(i, " ~ time*factor(group)"))
#                 model <- lme(formula,
#                              random = ~ 1|os_id,
#                              data = b_cells,
#                              na.action = na.omit,
#                              correlation = corCAR1(form = ~ 1 | os_id),
#                              control = lmeControl(maxIter = 5000, opt = "optim"))
#                 coef <- cbind(cell_type = i, n = length(levels(getGroups(model))), as.data.frame(summary(model)$tTable)[4,], as.data.frame(unclass(intervals(model, which = "fixed"))$fixed)[4,])
#                 b_cells_results <- rbind(b_cells_results, coef)
#         }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }

table((b_cells %>% distinct(os_id, .keep_all = T))$group)
cells <- full_join(myeloid %>% select(-c(myeloid_index:myeloid_index2, replicate)), full_join(t_cells %>% select(-c(t_cell_sample:ppid)), full_join(t_cells_2 %>% select(-c(t_cell_2_sample:ppid)), b_cells %>% select(-c(b_cell_sample:bcd)))))

cells <- cells %>% left_join(df %>% select(baid, os_id, visit, tx_group, visit_date, gender, age_base, dmt_class, tudca_levels)) %>% relocate(baid:tudca_levels, .before = 1) 
cells <- cells %>% ungroup() %>% group_by(baid) %>% arrange(baid, visit_date)


viable_names <- names(cells)[str_detect(names(cells), "viable")]
cells <- cells %>% select(-c(all_of(viable_names)))
cells_base <- cells %>% group_by(baid) %>% select(-age_base) %>% filter(visit_date == first(visit_date)) %>% dplyr::rename_with(.cols = myeloid:asc_b, .fn = ~paste0(.x, "_base"))

cells <- cells %>% left_join(cells_base %>% select(baid, all_of(names(cells_base)[str_detect(names(cells_base), "_base")]))) %>%
        group_by(baid) %>% arrange(baid, visit_date) %>% mutate(time = (time_length(difftime(visit_date, first(visit_date)), unit = "week"))/16)

cells_results <- c()
for(i in names(cells %>% ungroup() %>% select(c(myeloid:asc_b)))) {
        tryCatch({
                formula <- as.formula(paste0(i, " ~ time*factor(tx_group) + time*", i, "_base"))
                model <- lme(formula,
                             random = ~ 1|baid,
                             data = cells,
                             na.action = na.omit,
                             control = lmeControl(maxIter = 5000, opt = "optim"))
                coef <- cbind(cell_type = i, n = length(levels(getGroups(model))), as.data.frame(summary(model)$tTable)[5,], as.data.frame(unclass(intervals(model, which = "fixed"))$fixed)[5,])
                cells_results <- rbind(cells_results, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

cells_results <- cells_results %>% dplyr::rename(dif_pval = "p-value") %>% mutate(dif_beta_ci = paste0(round(Value, 2), "% (", round(lower, 2), "% to ", round(upper, 2), "%)"),
                                                                                  dif_fdr = p.adjust(dif_pval, method = "BH")) %>% relocate(c(dif_pval,dif_fdr), .after = dif_beta_ci)


cells_results_tx <- c()
for(i in  names(cells %>% ungroup() %>% select(c(myeloid:asc_b)))) {
        tryCatch({ 
                formula <- as.formula(paste0(i, " ~ time"))
                model <- lme(formula,
                             random = ~ 1|baid,
                             data = cells %>% filter(tx_group == 1),
                             na.action = na.omit,
                             control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(cell_type = i, tudca_n = length(levels(getGroups(model))), as.data.frame(summary(model)$tTable)[2,], as.data.frame(unclass(intervals(model, which = "fixed"))$fixed)[2,])
                cells_results_tx <- rbind(cells_results_tx, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

cells_results_tx <- cells_results_tx %>% dplyr::rename(tudca_pval = "p-value") %>% mutate(tudca_beta_ci = paste0(round(Value, 2), "% (", round(lower, 2), "% to ", round(upper, 2), "%)"),
                                                                                          tudca_fdr = p.adjust(tudca_pval, method = "BH")) %>% relocate(c(tudca_pval,tudca_fdr), .after = tudca_beta_ci)

cells_results_p <- c()
for(i in  names(cells %>% ungroup() %>% select(c(myeloid:asc_b)))) {
        tryCatch({
                formula <- as.formula(paste0(i, " ~ time"))
                model <- lme(formula,
                             random = ~ 1|baid,
                             data = cells %>% filter(tx_group == 0),
                             na.action = na.omit,
                             control = lmeControl(maxIter = 500, opt = "optim"))
                coef <- cbind(cell_type = i, placebo_n = length(levels(getGroups(model))), as.data.frame(summary(model)$tTable)[2,], as.data.frame(unclass(intervals(model, which = "fixed"))$fixed)[2,])
                cells_results_p <- rbind(cells_results_p, coef)
        }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

cells_results_p <- cells_results_p %>% dplyr::rename(placebo_pval = "p-value") %>% mutate(placebo_beta_ci = paste0(round(Value, 2), "% (", round(lower, 2), "% to ", round(upper, 2), "%)"),
                                                                                          placebo_fdr = p.adjust(placebo_pval, method = "BH")) %>% relocate(c(placebo_pval,placebo_fdr), .after = placebo_beta_ci)


cells_table <- full_join(cells_results_p %>% select(cell_type, placebo_n, placebo_beta_ci, placebo_pval), 
                         full_join(cells_results_tx %>% select(cell_type, tudca_n, tudca_beta_ci, tudca_pval), cells_results %>% select(cell_type, n, dif_beta_ci, dif_pval)))

cells_table <- cells_table %>% mutate(type = ifelse(cell_type %in% names(myeloid), "myeloid",
                                                    ifelse(cell_type %in% names(t_cells), "t_cell_1",
                                                           ifelse(cell_type %in% names(t_cells_2), "t_cell_2",
                                                                  ifelse(cell_type %in% names(b_cells), "b_cell", NA)))), .after = cell_type) 
t_cells_table <- cells_table %>% filter(type == "t_cell_1") %>% 
        filter(!cell_type %in% c("ccr4_parent", "ccr4_2_parent","th_memory_parent", "tc_memory_parent")) %>%
        mutate(placebo_fdr = p.adjust(placebo_pval, method = "BH"),
               tudca_fdr = p.adjust(tudca_pval, method = "BH"),
               dif_fdr = p.adjust(dif_pval, method = "BH"))
# write_csv(t_cells_table, "./out/t_cells_outcomes.csv")

#plots for cells----
#Spaghetti plots t cells
significant_t_cells <- (t_cells_table %>% filter(dif_pval<0.05))$cell_type
t_cell_plot <- cells %>% select(baid:visit, significant_t_cells)
t_cell_plot_long <- t_cell_plot %>% pivot_longer(-c(baid:visit), names_to = "t_cell", values_to = "value")
# t_cell_plot_long <- t_cell_plot_long %>% mutate(visit = visit + 0.1*if_else(tx_group == 1, 1, -1))


ggplot(t_cell_plot_long, aes(x = visit, y = value, group = baid, color = factor(tx_group))) + 
        geom_point( ) + 
        geom_line() + theme(legend.position = "none") +
        scale_x_continuous(breaks = c(1,2,3), labels = c(1,2,3)) +
        facet_wrap(~ t_cell, scale = "free") +
        theme_bw() +
        ylab("T cell population (%)") +
        xlab("Visit") +
        scale_color_viridis_d(name = "Group",
                              labels = c("Placebo", "TUDCA"),
                              begin = 0.3,
                              end = 0.6,
                              option = "B") +
        theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
              axis.text = element_text(size = 12, color = "black"))


pdf("cell_plot.pdf", width = 12, height = 8)
ggplot(t_cell_plot_long, aes(x = factor(visit), y = value, group = interaction(factor(tx_group),factor(visit)))) + 
        # stat_boxplot(geom = "errorbar", width = 0.3, position = position_dodge(width = 0.7, preserve = "total")) +
        geom_boxplot(position = position_dodge(width = 0.7, preserve = "total"), width = 0.5, outlier.shape = NA, aes(fill = factor(tx_group)), alpha = 0.5) +
        geom_point(position = position_dodge(width = 0.7), aes(color = factor(tx_group)), size = 1.7) + 
        facet_wrap(~ t_cell, scale = "free", labeller = labeller(t_cell = 
                                                                         c("th1_th17_parent" = "Th1/17 (double positive)",
                                                                           "th_naive_parent" = "Th Naive", 
                                                                           "th_cm_parent" = "CD4+ CM",
                                                                           "tc_cm_parent" = "CD8+ CM",
                                                                           "th1_parent" = "Th1",
                                                                           "treg_parent" = "Treg"))) +
        theme_bw() +
        theme(axis.title = element_text(size = 14, face = "bold", color = "black"),
              axis.text = element_text(size = 12, color = "black"),
              strip.background = element_rect(fill="black"),
              strip.text=element_text(color="white",size = 14, face = "bold"),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14, face = "bold")) +
        
        ylab("T cell population (% of parent cell)") +
        xlab("Visit") +
        scale_color_viridis_d(name = "Group",
                              labels = c("Placebo", "TUDCA"),
                              begin = 0.3,
                              end = 0.6,
                              option = "B") +
        scale_fill_viridis_d(name = "Group",
                             labels = c("Placebo", "TUDCA"),
                             begin = 0.3,
                             end = 0.6,
                             option = "B") 
dev.off()
