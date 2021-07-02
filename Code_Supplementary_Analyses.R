# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for    ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#              Disentangling heterogeneity of psychosis expression
#             in the general population: sex-specific moderation effects
#               of environmental risk factors on symptom networks
#
#
#
#                 - Analysis reported in supplementary material -
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ---------------------- 0: Reproducibility -----------------------------

# for reproducibility, we use the "checkpoint" package
# in a temporary directory, it will *install* package versions used when the script was written
# these versions are then used to run the script
# to this end, a server with snapshot images of archived package versions needs to be contacted
# for more info visit: https://mran.microsoft.com/documents/rro/reproducibility


library(checkpoint)
checkpoint(
  snapshotDate = "2021-03-08",
  R.version = "4.0.4",
  checkpointLocation = tempdir()
)
# ---------------------- 1: Load packages & data -----------------------------------
# ------ 1.1 load packages
library(haven)
library(networktree)
library(bootnet)
library(qgraph)
library(BGGM)
library(coin)
library(tidyverse)

# ------ 1.2 define custom function
# for checking stability of estimated networks
check_stability <- function(data, iters = 5000) {
  sub_network <- estimateNetwork(data, default = "pcor", alpha = 1)
  set.seed(1)
  boot_sub_network <-
    bootnet(
      sub_network,
      type = "case",
      statistics = "edge",
      iters = iters,
      verbose = FALSE,
      replacement = TRUE
    )
  cor_stability <- corStability(boot_sub_network, verbose = FALSE)
  return(cor_stability)
}


# for checking accuracy of edges in estimated networks
check_accuracy <- function(data, iters = 5000) {
  sub_network <- estimateNetwork(data, default = "pcor", alpha = 1)
  set.seed(1)
  boot_sub_network <-
    bootnet(
      sub_network,
      iters = iters,
      statistics = "edge",
      verbose = FALSE,
      replacement = TRUE
    )
  return(boot_sub_network)
}

# ------ 1.3 load 2007 AMPS data
apms07arch <-
  read_sav("UKDA-6379-spss/spss/spss19/apms07arch.sav")

# ------- 1.4 data preparation
data_apms <- apms07arch %>%
  transmute(
    # definition of symptoms following Moffa et al., 2017
    # affective symptoms
    worry = case_when(I10 > 1 ~ 1, is.na(I1) ~ NA_real_, TRUE ~ 0),
    sleep_pr = case_when(D10 > 1 ~ 1,  is.na(D1) ~ NA_real_, TRUE ~ 0),
    anx = case_when(J11 > 1 ~ 1, is.na(J1) ~ NA_real_, TRUE ~ 0),
    depr = case_when(G10 > 1 ~ 1, is.na(G1) ~ NA_real_, TRUE ~ 0),
    # psychosis symptoms
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
          is.na(VBb) & is.na(VBc) ~ NA_character_,
        TRUE ~ "no"
      )
    ),
    physical_abuse = as.factor(case_when(
      VBd == 1 ~ "yes",
      is.na(VBd) ~ NA_character_,
      TRUE ~ "no"
    )),
    cannabis = as.factor(ifelse(Cannyear == 1, "yes", "no")),
    alcohol = as.numeric(DVAudit1),
    ethnicity = as.factor(ETHNIC4),
    deprivation = as.numeric(qimd)
  )

data_missings_removed <- data_apms %>% na.omit()


# ---------------------- 2: Supplementary Table 1: differences in included vs. excluded participants ----------------------
# ------ 2.1 data preparation
data_included_excluded <- data_apms %>%
  mutate(included_excluded = if_else(rowSums(is.na(.)) > 0, "excluded", "included"))

# ------- 2.2 descriptive statistics
# n per group
data_included_excluded %>%
  group_by(included_excluded) %>%
  summarise(n = n())

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
  select(.,-c(age, alcohol, deprivation, ethnicity)) %>%
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

round(table(excluded$ethnicity) / sum(!is.na(excluded$ethnicity)), 3)

# ------- 2.3 statistical comparison
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
  select(., -c(age, alcohol, deprivation, included_excluded)) %>%
  map(
    ~ chisq_test(
      as.factor(.) ~ as.factor(data_included_excluded$included_excluded),
      data = data_included_excluded,
      distribution = "approximate"
    )
  )

# ---------------------- 3: Supplementary Table 2: stability of results ----------------------
# ------- 3.0 full sample
full_sample <-
  check_stability(data_missings_removed %>% .[, 1:6]) # 0.750069

# ------- 3.1 first split: sex differences

first_split_grp1 <-
  check_stability(data_missings_removed %>% filter(sex == "female") %>% .[, 1:6]) # 0.7499392


first_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "male") %>% .[, 1:6]) # 0.7499201

# ------- 3.2 second split: childhood sexual abuse in women

second_split_grp1 <-
  check_stability(data_missings_removed %>% filter(sex == "female" &
                                                     sexual_abuse == "yes") %>% .[, 1:6]) # 0.7503546

second_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "female" &
                                                     sexual_abuse == "no") %>% .[, 1:6]) # 0.7501466


# ------- 3.3 third split: childhood physical abuse in women

third_split_grp1 <-
  check_stability(
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "no" &
                                       physical_abuse == "yes") %>% .[, 1:6]
  ) # 0.2790698

third_split_grp2 <-
  check_stability(
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "no" &
                                       physical_abuse == "no") %>% .[, 1:6]
  ) # 0.75

# ------ 3.4 fourth split: in women: domestic violence

fourth_split_grp1 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "female" &
        sexual_abuse == "no" &
        physical_abuse == "no" &
        violence == "yes"
    ) %>% .[, 1:6]
  ) # 0.5181518

fourth_split_grp2 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "female" &
        sexual_abuse == "no" &
        physical_abuse == "no" &
        violence == "no"
    ) %>% .[, 1:6]
  )  # 0.75

# ------ 3.5 fifth split: in men: domestic violence

fifth_split_grp1 <-
  check_stability(data_missings_removed %>% filter(sex == "male" &
                                                     violence == "yes") %>% .[, 1:6]) # 0.4358974

fifth_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "male" &
                                                     violence == "no") %>% .[, 1:6])  # 0.7499159

# ------ 3.6 sixth split: in men: cannabis
sixth_split_grp1 <-
  check_stability(
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "no" &
                                       cannabis == "yes") %>% .[, 1:6]
  ) # 0.4368932


sixth_split_grp2 <-
  check_stability(data_missings_removed %>% filter(sex == "male" &
                                                     violence == "no" &
                                                     cannabis == "no") %>% .[, 1:6]) # 0.7500904

# ------ 3.7 seventh split: in men: ethnicity
seventh_split_grp1 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) # White, Other
    ) %>% .[, 1:6]
  ) # 0.75

seventh_split_grp2 <-
  check_stability(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(2, 3)  # Black, South Asian
    ) %>% .[, 1:6]
  ) # 0.515528

# ---------------------- 4: Supplementary Figures ----------------------
# ------ Supplementary Figure 1
# NetworkTree for women and men separately
# --- women
set.seed(234)
apms_networktree_women <- networktree(
  nodevars = data_missings_removed %>% filter(sex == "female") %>% .[, 1:6],
  # symptoms
  splitvars = data_missings_removed %>% filter(sex == "female") %>% .[, 7:ncol(data_missings_removed)],
  # risk factors
  method = "mob",
  # default method
  model = c("correlation"),
  transform = "pcor",
  minsize = ceiling(nrow(data_missings_removed) * 0.01) # minimal size terminal node
)

plot(
  apms_networktree_women,
  vsize = 10,
  edge.width = 2.5,
  label.cex = 2.5,
  cut = 0,
  borders = T,
  height = 1,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)

# --- men
set.seed(234)
apms_networktree_men <- networktree(
  nodevars = data_missings_removed %>% filter(sex == "male") %>% .[, 1:6],
  # symptoms
  splitvars = data_missings_removed %>% filter(sex == "male") %>% .[, 7:ncol(data_missings_removed)],
  # risk factors
  method = "mob",
  # default method
  model = c("correlation"),
  transform = "pcor",
  minsize = ceiling(nrow(data_missings_removed) * 0.01) # minimal size terminal node
)

plot(
  apms_networktree_men,
  vsize = 10,
  edge.width = 2.5,
  label.cex = 2.5,
  cut = 0,
  borders = T,
  height = 1,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)


# ------ Supplementary Figure 2

# ------- 3.0 full sample
full_sample_edge_accuracy <-
  check_accuracy(data_missings_removed %>% .[, 1:6])

plot(full_sample_edge_accuracy, order = "sample")
ggsave("full_sample_edge_accuracy.tiff")

# ------- 3.1 first split: sex differences

first_split_grp1_edge_accuracy <-
  check_accuracy(data_missings_removed %>% filter(sex == "female") %>% .[, 1:6])

plot(first_split_grp1_edge_accuracy, order = "sample")
ggsave("first_split_grp1_edge_accuracy.tiff")

first_split_grp2_edge_accuracy <-
  check_accuracy(data_missings_removed %>% filter(sex == "male") %>% .[, 1:6])

plot(first_split_grp2_edge_accuracy, order = "sample")
ggsave("first_split_grp2_edge_accuracy.tiff")

# ------- 3.2 second split: childhood sexual abuse in women

second_split_grp1_edge_accuracy <-
  check_accuracy(data_missings_removed %>% filter(sex == "female" &
                                                    sexual_abuse == "yes") %>% .[, 1:6])

plot(second_split_grp1_edge_accuracy, order = "sample")
ggsave("second_split_grp1_edge_accuracy.tiff")


second_split_grp2_edge_accuracy <-
  check_accuracy(data_missings_removed %>% filter(sex == "female" &
                                                    sexual_abuse == "no") %>% .[, 1:6])

plot(second_split_grp2_edge_accuracy, order = "sample")
ggsave("second_split_grp2_edge_accuracy.tiff")



# ------- 3.3 third split: childhood physical abuse in women

third_split_grp1_edge_accuracy <-
  check_accuracy(
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "no" &
                                       physical_abuse == "yes") %>% .[, 1:6]
  ) 


plot(third_split_grp1_edge_accuracy, order = "sample")
ggsave("third_split_grp1_edge_accuracy.tiff")


third_split_grp2_edge_accuracy <-
  check_accuracy(
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "no" &
                                       physical_abuse == "no") %>% .[, 1:6]
  )

plot(third_split_grp2_edge_accuracy, order = "sample")
ggsave("third_split_grp2_edge_accuracy.tiff")

# ------ 3.4 fourth split: in women: domestic violence

fourth_split_grp1_edge_accuracy <-
  check_accuracy(
    data_missings_removed %>% filter(
      sex == "female" &
        sexual_abuse == "no" &
        physical_abuse == "no" &
        violence == "yes"
    ) %>% .[, 1:6]
  )

plot(fourth_split_grp1_edge_accuracy, order = "sample")
ggsave("fourth_split_grp1_edge_accuracy.tiff")

fourth_split_grp2_edge_accuracy <-
  check_accuracy(
    data_missings_removed %>% filter(
      sex == "female" &
        sexual_abuse == "no" &
        physical_abuse == "no" &
        violence == "no"
    ) %>% .[, 1:6]
  )

plot(fourth_split_grp2_edge_accuracy, order = "sample")
ggsave("fourth_split_grp2_edge_accuracy.tiff")

# ------ 3.5 fifth split: in men: domestic violence

fifth_split_grp1_edge_accuracy <-
  check_accuracy(data_missings_removed %>% filter(sex == "male" &
                                                    violence == "yes") %>% .[, 1:6])

plot(fifth_split_grp1_edge_accuracy, order = "sample")
ggsave("fifth_split_grp1_edge_accuracy.tiff")

fifth_split_grp2_edge_accuracy <-
  check_accuracy(data_missings_removed %>% filter(sex == "male" &
                                                    violence == "no") %>% .[, 1:6])

plot(fifth_split_grp2_edge_accuracy, order = "sample")
ggsave("fifth_split_grp2_edge_accuracy.tiff")

# ------ 3.6 sixth split: in men: cannabis
sixth_split_grp1_edge_accuracy <-
  check_accuracy(
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "no" &
                                       cannabis == "yes") %>% .[, 1:6]
  )

plot(sixth_split_grp1_edge_accuracy, order = "sample")
ggsave("sixth_split_grp1_edge_accuracy.tiff")


sixth_split_grp2_edge_accuracy <-
  check_accuracy(data_missings_removed %>% filter(sex == "male" &
                                                    violence == "no" &
                                                    cannabis == "no") %>% .[, 1:6])


plot(sixth_split_grp2_edge_accuracy, order = "sample")
ggsave("sixth_split_grp2_edge_accuracy.tiff")

# ------ 3.7 seventh split: in men: ethnicity
seventh_split_grp1_edge_accuracy <-
  check_accuracy(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) # White, Other
    ) %>% .[, 1:6]
  )

plot(seventh_split_grp1_edge_accuracy, order = "sample")
ggsave("seventh_split_grp1_edge_accuracy.tiff")


seventh_split_grp2_edge_accuracy <-
  check_accuracy(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(2, 3)  # Black, South Asian
    ) %>% .[, 1:6]
  )

plot(seventh_split_grp2_edge_accuracy, order = "sample")
ggsave("seventh_split_grp2_edge_accuracy.tiff")

