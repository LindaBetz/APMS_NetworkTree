# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for    ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#              Disentangling heterogeneity in the psychosis spectrum: 
#          sex-specific moderation effects of environmental risk factors 
#                               on symptom networks  
#
#
#
#                   - Analysis reported in supplementary materials -
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# ---------------------------------- 0: Reproducibility  -----------------------------------

# for reproducibility, one can use the "checkpoint" package
# in a temporal directory, it will *install* those package versions used when the script was written
# these versions are then used to run the script
# to this end, a server with snapshot images of archived package versions needs to be contacted
# for more info visit: https://mran.microsoft.com/documents/rro/reproducibility

library(checkpoint)
checkpoint(
  snapshotDate = "2021-03-08",
  R.version = "4.0.4",
  checkpointLocation = tempdir()
)
# ---------------------------------- 1: Load packages & data -----------------------------------
# ------- 1.1 load packages -------
library(haven)
library(networktree)
library(bootnet)
library(qgraph)
library(BGGM)
library(coin)
library(tidyverse)

# ------- 1.2 define custom function -------
check_stability <- function(data, iters = 5000) {
  out <-  estimateNetwork(data, default = "pcor", alpha = 1)
  set.seed(1)
  boot_out <-
    bootnet(
      out,
      type = "case",
      statistics = "edge",
      iters = iters,
      verbose = FALSE,
      replacement = TRUE
    )
  res <- corStability(boot_out, verbose = FALSE)
  return(res)
}

# ------- 1.3 load 2007 AMPS data -------
apms07arch <-
  read_sav("UKDA-6379-spss/spss/spss19/apms07arch.sav")

# ------- 1.4 data preparation -------
# ------- 1.4 data preparation -------
data_apms <- apms07arch %>%
  transmute(
    # symptoms
    worry = case_when(I10 > 1 ~ 1, is.na(I1) ~ NA_real_, TRUE ~ 0),
    sleep_pr = case_when(D10 > 1 ~ 1,  is.na(D1) ~ NA_real_, TRUE ~ 0),
    anx = case_when(J11 > 1 ~ 1, is.na(J1) ~ NA_real_, TRUE ~ 0),
    depr = case_when(G10 > 1 ~ 1, is.na(G1) ~ NA_real_, TRUE ~ 0),
    # psychosis
    per = ifelse(PSQ3 == 1 & PSQ3a == 1, 1, 0),
    hal = ifelse(PSQ5 == 1 & PSQ5a == 1, 1, 0),
    # potential moderators/split variables
    age = as.numeric(ResAge),
    sex = as.factor(ifelse(ResSex == 2, "female", "male")),
    bullying =  as.factor(
      case_when(Trauma31 == 1 ~ "yes",
                is.na(Trauma31) ~ NA_character_,
                TRUE ~ "no")
    ),
    violence =  as.factor(
      case_when(Trauma33 == 1 ~ "yes",
                is.na(Trauma33) ~ NA_character_,
                TRUE ~ "no")
    ),
    separation = as.factor(
      case_when(
        LACare == 1 | ChldInst == 1 ~ "yes",
        is.na(LACare) & is.na(ChldInst) ~ NA_character_,
        TRUE ~ "no"
      )
    ),
    sexual_abuse = as.factor(
      case_when(
        VBa == 1 | VBb == 1 | VBc == 1 ~ "yes",
        is.na(VBa) &
          is.na(VBb) & is.na(VBc) & is.na(VBd) ~ NA_character_,
        TRUE ~ "no"
      )
    ),
    physical_abuse = as.factor(
      case_when(
        VBd == 1 ~ "yes",
        is.na(VBa) &
          is.na(VBb) & is.na(VBc) & is.na(VBd) ~ NA_character_,
        TRUE ~ "no"
      )
    ),
    cannabis = as.factor(ifelse(Cannyear == 1, "yes", "no")),
    alcohol = as.numeric(DVAudit1),
    ethnicity = as.factor(ETHNIC4),
    deprivation = as.numeric(qimd)
  )

data_missings_removed <- data_apms  %>% na.omit()


# ---------------------------------- 2: Supplementary Table 1: differences in included vs. excluded participants  -----------------------------------
# ------- 2.1 data preparation -------
data_included_excluded <- data_apms %>%
  mutate(included_excluded = if_else(rowSums(is.na(.)) > 0, "excluded", "included"))

# ------- 2.2 descriptive statistics -------
# interval data
data_included_excluded %>%
  select(., c(age, alcohol, deprivation, included_excluded)) %>%
  group_by(included_excluded) %>%
  summarise(across(everything(), c(
    median = ~ median(.x, na.rm = TRUE),
    IQR = ~ IQR(.x, na.rm = TRUE)
  )))

# frequency data
data_included_excluded %>%
  select(., -c(age, alcohol, deprivation, ethnicity)) %>%
  mutate(across(
    where(is.factor),
    ~ case_when(
      . %in% c("yes", "female") ~ 1,
      . %in% c("no", "male") ~ 0,
      TRUE ~ NA_real_
    )
  )) %>%
  group_by(included_excluded) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

excluded <- data_included_excluded %>%
  filter(included_excluded == "excluded")

round(table(excluded$ethnicity) / 96, 3)

# ------- 2.3 statistical comparison -------
# interval data
set.seed(1)
data_included_excluded %>%
  select(., c(age, alcohol, deprivation)) %>%
  map(
    ~ wilcox_test(
      as.numeric(.) ~ as.factor(data_included_excluded$included_excluded),
      data = data_included_excluded,
      distribution = "approximate"
    )
  )

# frequency data
set.seed(1)
data_included_excluded %>%
  select(.,-c(age, alcohol, deprivation, included_excluded)) %>%
  map(
    ~ chisq_test(
      as.factor(.) ~ as.factor(data_included_excluded$included_excluded),
      data = data_included_excluded,
      distribution = "approximate"
    )
  )

# ---------------------------------- 3: Supplementary Table 2: stability of results -----------------------------------
# ------ 3.1 full sample ------
full_sample <-
  check_stability(data_missings_removed %>% .[, 1:6]) # 0.7500344

# ------ 3.2 first split: sex differences ------

first_split_grp1 <-
  check_stability(data_missings_removed %>% filter(sex == "female") %>% .[, 1:6]) # 0.7500606


first_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "male") %>% .[, 1:6]) # 0.75

# ------ 3.3 second split: childhood sexual abuse in women ------

second_split_grp1 <-
  check_stability(data_missings_removed %>% filter(sex == "female" &
                                                     sexual_abuse == "yes") %>% .[, 1:6]) # 0.6728045

second_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "female" &
                                                     sexual_abuse == "no") %>% .[, 1:6]) # 0.7499269

# ------ 3.4 third split: in men: domestic violence ------

third_split_grp1 <-
  check_stability(data_missings_removed %>% filter(sex == "male" &
                                                     violence == "yes") %>% .[, 1:6]) # 0.4358974

third_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "male" &
                                                     violence == "no") %>% .[, 1:6])  # 0.75

# ------ 3.5 fourth split: in men & domestic violence==no: cannabis use ------

fourth_split_grp1 <-
  check_stability(
    data_missings_removed %>% filter(sex == "male" & violence == "no" &
                                       cannabis == "yes") %>% .[, 1:6]
  ) # 0.4368932

fourth_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "male" &
                                                     violence == "no" &
                                                     cannabis == "no") %>% .[, 1:6]) # 0.7501805

# ------ 3.6 fifth split: in men & domestic violence==no & cannabis==no: ethnicity ------
fifth_split_grp1 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "male" & violence == "no" &
        cannabis == "no" & ethnicity %in% c(1, 4)
    ) %>% .[, 1:6]
    
  ) # 0.7499041


fifth_split_grp2 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "male" & violence == "no" &
        cannabis == "no" & ethnicity %in% c(2, 3)
    ) %>% .[, 1:6]
    
  ) # 0.5153374

# ------ 3.7 sixth split: in men & domestic violence==no & cannabis==no & ethnicity == 1,4: lifetime bullying  ------
sixth_split_grp1 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "male" & violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) & bullying == "yes"
    ) %>% .[, 1:6]
  ) # 0.5943396

sixth_split_grp2 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "male" & violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) & bullying == "no"
    ) %>% .[, 1:6]
  ) # 0.7498855



### -- Supplementary figure 1

data_missings_removed %>%
  mutate(
    group = case_when(
      sex == "female" & sexual_abuse == "yes" ~ "women: sexual abuse",
      sex == "female" &
        sexual_abuse == "no" ~ "women: no sexual abuse",
      sex == "male" &
        violence == "yes" ~ "men: domestic violence",
      sex == "male" &
        violence == "no" &
        cannabis == "yes" ~ "men: no domestic violence, cannabis use",
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(2, 3) ~ "men: no domestic violence, no cannabis use, ethnic minority",
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) &
        bullying == "yes" ~ "men: no domestic violence, no cannabis use, ethnic majority, bullying",
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) &
        bullying == "no" ~ "men: no domestic violence, no cannabis use, ethnic majority, no bullying",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(group) %>%
  nest() %>%
  mutate(nets = map(
    data,
    . %>% select(c(
      "worry" , "sleep_pr", "anx", "depr", "per", "hal"
    )) %>% bootnet::estimateNetwork(
      .,
      default = "pcor",
      corMethod = "cor",
      alpha = 1
    ) %>% .$graph
  )) %>%
  mutate(global_strength = map_dbl(nets, ~ abs(sum(abs(
    .[upper.tri(.)]
  ))))) %>%
  ungroup() %>%
  mutate(group = reorder(group, global_strength)) %>%
  ggplot(., aes(x = group, y = global_strength)) +
  geom_bar(stat = "identity", fill = "#6884d5") +
  theme_classic() +
  scale_fill_viridis_d() +
  geom_hline(yintercept =  1.740785,
             linetype = "dashed",
             size = 1) +
  coord_flip() +
  ylab("\nGlobal Strength") +
  xlab("") +
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black")
  )
ggsave(
  "supplementary_figure_1.png",
  width = 12,
  height = 6,
  dpi = 400
)
