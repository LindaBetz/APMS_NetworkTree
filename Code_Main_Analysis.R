# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     Code for    ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#
#              Disentangling heterogeneity of psychosis expression
#             in the general population: sex-specific moderation effects
#               of environmental risk factors on symptom networks
#
#
#
#                   - Analysis reported in main manuscript -
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ---------------------- 0: Reproducibility ----------------------

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

# ---------------------- 1: Load packages & data ----------------------
# ------ 1.1: load packages
library(haven)
library(networktree)
library(bootnet)
library(qgraph)
library(BGGM)
library(tidyverse)

# ------ 1.2: define custom functions 
global_strength <- function(data1, data2) {
  # global strength as defined in:
  # https://github.com/cvborkulo/NetworkComparisonTest/blob/master/R/NCT.R
  # estimate graphs
  fit1 <-  BGGM::estimate(data1, analytic = TRUE)
  fit2 <-  BGGM::estimate(data2, analytic = TRUE)
  
  # select graphs
  sel1 <-
    BGGM::select(fit1, cred = 0) # cred = 0 => no "pruning" of edges
  sel2 <- BGGM::select(fit2, cred = 0)
  
  abs(sum(abs(sel1$pcor_adj[upper.tri(sel1$pcor_adj)])) - sum(abs(sel2$pcor_adj[upper.tri(sel2$pcor_adj)])))
}

individual_global_strength <- function(data1, data2) {
  fit1 <-  BGGM::estimate(data1, analytic = TRUE)
  fit2 <-  BGGM::estimate(data2, analytic = TRUE)
  
  # select graphs
  sel1 <-
    BGGM::select(fit1, cred = 0) # cred = 0 => no "pruning" of edges
  sel2 <- BGGM::select(fit2, cred = 0)
  
  return(list(sum(abs(sel1$pcor_adj[upper.tri(sel1$pcor_adj)])),
              sum(abs(sel2$pcor_adj[upper.tri(sel2$pcor_adj)]))))
  
}

# custom function to compare networks in relevant aspects (global strength, individual connections)
# & store results
compare_networks <- function(data1,
                             data2,
                             iter = 5000,
                             node_names =  c("worry" , "sleep_pr", "anx", "depr", "per", "hal")) {
  individual_edges <- ggm_compare_estimate(data1,
                                           data2,
                                           progress = FALSE,
                                           iter = iter,
                                           seed = 1)
  
  individual_edges <-
    as.data.frame(BGGM::select(individual_edges, cred = 0.95)$pcor_adj)
  colnames(individual_edges) <-
    rownames(individual_edges) <- node_names
  
  
  obs_gs <-
    global_strength(data1,
                    data2)
  
  ind_gs <- individual_global_strength(data1,
                                       data2)
  
  res_gs <- ggm_compare_ppc(
    data1,
    data2,
    FUN = global_strength,
    iter = iter,
    custom_obs  = obs_gs,
    progress = FALSE
  )
  
  return(
    list(
      individual_edges = individual_edges,
      global_strength_diff = obs_gs,
      global_strength_group_1 = ind_gs[[1]],
      global_strength_group_2 = ind_gs[[2]],
      global_strength_p_val = res_gs$ppp_custom
    )
  )
}


# ------ 1.3: load 2007 AMPS data 
apms07arch <-
  read_sav("UKDA-6379-spss/spss/spss19/apms07arch.sav")

# ------ 1.4: data preparation
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

# ---------------------- 2: Sample Descriptives ----------------------
# n participants with missing data =
data_apms %>%
  mutate(any_NA = rowSums(is.na(.))) %>%
  filter(any_NA > 0) %>% nrow(.) # 146 participants have any variable missing

# % participants with missing data
146 / nrow(data_apms)  # 0.01972173


data_missings_removed <- data_apms  %>% na.omit()
nrow(data_missings_removed) # 7257 (= total sample size)

# Table 1: whole sample
# frequency data
data_missings_removed %>%
  select(., -c(age, alcohol, deprivation, ethnicity)) %>%
  mutate(across(
    where(is.factor),
    ~ case_when(
      . %in% c("yes", "female") ~ 1,
      . %in% c("no", "male") ~ 0,
      TRUE ~ NA_real_
    )
  )) %>%
  summarise(across(where(is.numeric), ~ mean(.x))) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

# ethnicity
round(table(data_missings_removed$ethnicity) / nrow(data_missings_removed),
      3)


# interval data
data_missings_removed %>%
  select(age, alcohol, deprivation) %>%
  summarise(across(everything(), c(median = median, IQR = IQR)))

# Table 1: men vs. women
# frequency data
data_missings_removed %>%
  group_by(sex) %>%
  select(., -c(age, alcohol, deprivation, ethnicity)) %>%
  mutate(across(
    where(is.factor),
    ~ case_when(
      . %in% c("yes", "female") ~ 1,
      . %in% c("no", "male") ~ 0,
      TRUE ~ NA_real_
    )
  )) %>%
  summarise(across(where(is.numeric), ~ mean(.x))) %>%
  mutate(across(where(is.numeric), ~ round(., 3)))

# ethnicity women
round(table(data_missings_removed[data_missings_removed$sex == "female", ]$ethnicity) /
        nrow(data_missings_removed[data_missings_removed$sex == "female", ]),
      3)

# ethnicity men
round(table(data_missings_removed[data_missings_removed$sex == "male", ]$ethnicity) /
        nrow(data_missings_removed[data_missings_removed$sex == "male", ]),
      3)

# interval data
data_missings_removed %>%
  group_by(sex) %>%
  select(age, alcohol, deprivation) %>%
  summarise(across(everything(), c(median = median, IQR = IQR)))

# ---------------------- 3: Network Estimation ----------------------
# estimate & plot network based on whole sample
node_names <-
  c(
    "worry",
    "sleep disturbance",
    "anxiety",
    "depression",
    "persecutory ideation",
    "hallucinatory experiences"
  )

png(filename = "whole_sample.png",
    width = 600,
    height = 600)
qgraph(
  cor(data_missings_removed %>% .[, 1:6],
      method = "pearson"),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()

# ---------------------- 4: Recursive Partitioning ----------------------
set.seed(234)
apms_networktree <- networktree(
  nodevars = data_missings_removed[, 1:6],
  # symptoms
  splitvars = data_missings_removed[, 7:ncol(data_missings_removed)],
  # risk factors
  method = "mob",
  # default method
  model = c("correlation"),
  transform = "pcor",
  minsize = ceiling(nrow(data_missings_removed) * 0.01) # minimal size terminal node
)


# ---------------------- 5: Subgroup Differences ----------------------
# ------ 5.1: first split: women vs. men 
# statistical comparison
set.seed(123)
first_split <-
  compare_networks(
    data_missings_removed %>% filter(sex == "female") %>% .[, 1:6],
    data_missings_removed %>% filter(sex == "male") %>% .[, 1:6]
  )

# ------ 5.2: second split: sexual abuse in women 
# statistical comparison
set.seed(123)
second_split <-
  compare_networks(
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "yes") %>% .[, 1:6],
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "no") %>% .[, 1:6]
  )

png(filename = "women_sexual_abuse_yes.png",
    width = 600,
    height = 600)
qgraph(
  cor(
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "yes") %>% .[, 1:6],
    method = "pearson"
  ),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()

png(filename = "women_sexual_abuse_no.png",
    width = 600,
    height = 600)
qgraph(
  cor(
    data_missings_removed %>% filter(sex == "female" &
                                       sexual_abuse == "no") %>% .[, 1:6],
    method = "pearson"
  ),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()


# ------ 5.3: third split: domestic violence in men
# statistical comparison
set.seed(123)
third_split <-
  compare_networks(
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "yes") %>% .[, 1:6],
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "no") %>% .[, 1:6]
  )

png(filename = "men_domestic_violence_yes.png",
    width = 600,
    height = 600)
qgraph(
  cor(
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "yes") %>% .[, 1:6],
    method = "pearson"
  ),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()

# ------ 5.4: fourth split: cannabis in men 
# statistical comparison
set.seed(123)
fourth_split <-
  compare_networks(
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "no" &
                                       cannabis == "yes") %>% .[, 1:6],
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "no" &
                                       cannabis == "no") %>% .[, 1:6]
  )

# plotting the respective subgroup network
png(filename = "men_cannabis_yes.png",
    width = 600,
    height = 600)
qgraph(
  cor(
    data_missings_removed %>% filter(sex == "male" &
                                       violence == "no" &
                                       cannabis == "yes") %>% .[, 1:6],
    method = "pearson"
  ),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()


# ------5.5: fifth split: ethnic background in men 
# statistical comparison
set.seed(123)
fifth_split <-
  compare_networks(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(2, 3)
    ) %>% .[, 1:6],
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4)
    ) %>% .[, 1:6]
  )

# plotting the respective subgroup network
png(filename = "men_ethnic_minority.png",
    width = 600,
    height = 600)
qgraph(
  cor(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(2, 3)
    ) %>% .[, 1:6],
    
    method = "pearson"
  ),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()

# ------ 5.6: sixth split: bullying in men
# statistical comparison
set.seed(123)
sixth_split <-
  compare_networks(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) & bullying == "yes"
    ) %>% .[, 1:6],
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" & ethnicity %in% c(1, 4) & bullying == "no"
    ) %>% .[, 1:6]
  )

png(filename = "men_ethnic_majority_bullying_yes.png",
    width = 600,
    height = 600)
qgraph(
  cor(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) & bullying == "yes"
    ) %>% .[, 1:6],
    
    method = "pearson"
  ),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()


png(filename = "men_ethnic_majority_bullying_no.png",
    width = 600,
    height = 600)
qgraph(
  cor(
    data_missings_removed %>% filter(
      sex == "male" &
        violence == "no" &
        cannabis == "no" &
        ethnicity %in% c(1, 4) & bullying == "no"
    ) %>% .[, 1:6],
    
    method = "pearson"
  ),
  graph = "pcor",
  vsize = 13.5,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 4,
  width = 1,
  minimum = 0.01,
  maximum = 0.4,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()

# ---------------------- 6: Figures ----------------------
# ------ Figure 2
# plotting for overview only; Figure 2 was edited in PowerPoint
plot(
  apms_networktree,
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
# ------ Figure 3: networks of women vs. men 
# plotting for overview only; Figure 3 was edited in PowerPoint

png(filename = "women.png",
    width = 900,
    height = 900)
qgraph(
  cor(data_missings_removed %>% filter(sex == "female") %>% .[, 1:6],
      method = "pearson"),
  graph = "pcor",
  vsize = 10,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 1,
  width = 1,
  minimum = 0,
  maximum = 0.4,
  edge.labels = TRUE,
  edge.label.color = "black",
  edge.label.position = 0.55,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()

png(filename = "men.png",
    width = 900,
    height = 900)
qgraph(
  cor(data_missings_removed %>% filter(sex == "male") %>% .[, 1:6],
      method = "pearson"),
  graph = "pcor",
  vsize = 10,
  edge.width = 1.75,
  label.cex = 2.25,
  cut = 0,
  borders = T,
  border.width  = 1,
  width = 1,
  minimum = 0,
  maximum = 0.4,
  edge.labels = TRUE,
  edge.label.color = "black",
  edge.label.position = 0.55,
  theme = "colorblind",
  layout = "circle",
  labels = c("1", "2", "3", "4", "5", "6"),
  color =  c(rep("#f7f5f2", 4),
             rep("#c5ceed", 2))
)
dev.off()
