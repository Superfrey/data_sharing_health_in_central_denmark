# ============================================================
# Public Attitudes Toward Health Data Sharing
# Full Analysis Script
# ============================================================
#
# Pipeline:
# 1. Load data
# 2. Define randomization groups and outcomes
# 3. Link registry-based covariates (age, sex, education, diabetes)
# 4. Apply inclusion/exclusion criteria
# 5. Estimate inverse probability weights
# 6. Descriptive analyses (weighted and unweighted)
# 7. Multinomial regression models
#
# NOTE:
# File paths are specific to Statistics Denmark (DST)
# ============================================================


# ---------------- Load necessary libraries ----------------
library(dplyr)
library(here)
library(readr)
library(data.table)
library(WeightIt)
library(finalfit)
library(nnet)


#---------------------- Functions defined--------------
# These functions are used to:
# - Link registry-based covariates
# - Format model outputs (odds ratios)


add_covariates <-
  function(study_population,
           index_year = 2024,
           education_year = 2021) {
    
    # Load data:
    # Bef (population register: age and sex)
    
    bef_data <-
      readRDS(
        paste0(
          "path_to/bef", #Specific to DST
          index_year - 1,
          "12.sas7bdat.rds"
        )
      )
    
    setDT(bef_data)
    names(bef_data) <- tolower(names(bef_data))
    
    bef_data <-
      bef_data[, .(
        pnr,
        age = alder,
        sex = as.factor(ifelse(koen ==
                                 '1', 'Men', 'Women'))
      )]
    
    # Education:
    # Function to derive years of education from registry codes
    
    edu_length_years <- function(education_year) {
      
      edu_data <-
        haven::read_sas(
          paste0(
            "path_to/udda", #Specific to DST
            education_year,
            "09.sas7bdat"
          )
        )
      
      formats <-
        read_sas(
          "path_to/n_audd_pria_l1lx_k.sas7bdat" #Specific to DST
        )
      
      setDT(edu_data)
      setDT(formats)
      
      edu_length <- merge(
        edu_data[, .(pnr = PNR, edu_code = HFAUDD)],
        formats[, .(
          edu_code = START,
          education_years = as.numeric(AUDD_PRIA_L1LX_K) / 12
        )],
        by = "edu_code",
        all.x = TRUE
      )
      
      return(edu_length[, .(pnr, education_years)])
    }
    
    education_data <- edu_length_years(education_year)
    
    # Diabetes
    # Define prevalent diabetes at index date
    
    diabetes_data <-
      readr::read_csv(
        "path_to/diabetes_onset_imputed_hicd_2023.csv",
        col_types = cols(
          diabetes_type = col_factor(levels = c("1",
                                                "2")),
          date_onset_diabetes_hicd = col_date(format = "%Y-%m-%d")
        )
      )
    
    setDT(diabetes_data)
    names(diabetes_data) <- tolower(names(diabetes_data))
    
    diabetes_data <- diabetes_data[, .(
      pnr,
      prevalent_diabetes = date_onset_diabetes_hicd < as.Date(paste0(index_year, "-01-01")),
      diabetes_type
    )]
    
    # standardize colnames:
    names(study_population) <- tolower(names(study_population))
    
    
    # Merge together:
    # Combine study population with registry covariates via pnr
    
    output <- Reduce(
      function(x, y)
        merge(x, y, by = 'pnr', all.x = TRUE),
      list(study_population,
           bef_data,
           education_data,
           diabetes_data)
    )
    
    # Replace missing diabetes status with FALSE
    output[, prevalent_diabetes := fifelse(is.na(prevalent_diabetes), FALSE, prevalent_diabetes)]
    
    return(output)
  }


estimates_OR <- function(model) {
  
  # Extract odds ratios with confidence intervals
  
  broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  ) %>%
    select(-c("std.error", "statistic")) %>%
    dplyr::filter(!term == "(Intercept)") %>%
    dplyr::rename(responds = y.level) %>%
    dplyr::mutate(across(c(estimate, conf.low, conf.high), ~ round_tidy(., 2)))
}


# ----------------------- Load Data---------------------------------
# Survey data and background population

hicd_2024 <- readRDS(here::here("path/to/hicd_cohort.rds"))
pop_2024 <- readRDS(here::here("path/to/background_population.rds"))


# ----------------------- Define randomisation groups ---------------------------------
# Combine randomization question into exposure groups (public vs private)

# Define public private randomization in HICD
# Combine randomization question in to single categories

hicd_2024 <- hicd_2024 %>%
  dplyr::mutate(
    
    aiq5 = coalesce(aiq5_a, aiq5_b),
    
    group_data_sharing_q6 = ifelse(
      (!is.na(aiq6a_1) | !is.na(aiq6a_2) | !is.na(aiq6a_3) | !is.na(aiq6a_4)) &
        (is.na(aiq6b_5) & is.na(aiq6b_6) & is.na(aiq6b_7) & is.na(aiq6b_8)), "Group public", "Group private"
    ),
    
    group_data_sharing_q6 = ifelse(
      (is.na(aiq6a_1) & is.na(aiq6a_2) & is.na(aiq6a_3) & is.na(aiq6a_4)) &
        (is.na(aiq6b_5) & is.na(aiq6b_6) & is.na(aiq6b_7) & is.na(aiq6b_8)), NA, group_data_sharing_q6
    ),
    
    aiq6_1 = coalesce(aiq6a_1, aiq6b_5),
    aiq6_2 = coalesce(aiq6a_2, aiq6b_6),
    aiq6_3 = coalesce(aiq6a_3, aiq6b_7),
    aiq6_4 = coalesce(aiq6a_4, aiq6b_8)
  ) %>%
  dplyr::select(-c(aiq6a_1, aiq6a_2, aiq6a_3, aiq6a_4, aiq6b_5, aiq6b_6, aiq6b_7, aiq6b_8))

# Current population: 56060


#-----------------------Label questionnarie data------------------
# Add descriptive labels for output tables

var_label(hicd_2024$aiq6_1) <- "Tekst og tal fra dine journaler"
var_label(hicd_2024$aiq6_2) <- "Scanningsbilleder fra hospitalet"
var_label(hicd_2024$aiq6_3) <- "Lydoptagelser af konsultationer"
var_label(hicd_2024$aiq6_4) <- "Genetisk information"


# ----------------------- Exclude participant without response ---------------------------------
# Restrict to complete cases on exposure and outcomes

hicd_2024 <- hicd_2024 %>%
  filter(across(c(group_data_sharing_q6, aiq6_1, aiq6_2, aiq6_3, aiq6_4), ~ !is.na(.)))

# Current population: 39048


# ----------------------- Add variables defined in the registries ---------------------------------
# Link registry-based covariates to study and background populations

responses_covariates <- add_covariates(
  study_population = as.data.table(hicd_2024),
  index_year = 2023,
  education_year = 2022
) %>%
  rename(diabetes_type = diabetes_type.y) %>%
  mutate(
    diabetes = factor(
      ifelse(
        is.na(prevalent_diabetes) |
          prevalent_diabetes == FALSE,
        "No diabetes",
        diabetes_type
      )
    ),
    age = ifelse(is.na(age), alder + 1, age)
  )

background_covariates <- add_covariates(
  study_population = as.data.table(pop_2024),
  index_year = 2023,
  education_year = 2022
) %>%
  mutate(
    diabetes = factor(ifelse(is.na(prevalent_diabetes) | prevalent_diabetes == FALSE, "No diabetes", diabetes_type))
  )


# ----------------------- Exclude participant without variable on education in HICD---------------------------------
# Ensure complete data for weighting variables

responses_covariates <- responses_covariates[!is.na(education_years)]
# 38740 left


background_covariates[, in_study_pop := pnr %in% responses_covariates$pnr]

background_covariates <- background_covariates[, .(pnr, age, sex, education_years, diabetes, in_study_pop)]

background_covariates_clean <- background_covariates[
  complete.cases(background_covariates) &
    age <= max(responses_covariates$age, na.rm = TRUE) &
    age >= min(responses_covariates$age, na.rm = TRUE)
]


# --------------Calculate weights based on background population on study population --------
# Estimate inverse probability weights using demographic and clinical variables

weightit_model_background <-
  WeightIt::weightit(
    in_study_pop ~ age + sex + diabetes + education_years,
    data = background_covariates_clean,
    method = "ps",
    stabilize = TRUE
  )


background_covariates_clean$weights_background <- weightit_model_background$weights


# Merge weight estimates to the study population and call the new df responses_weights
# This is the df we are going to use in our analysis

responses_weights <- Reduce(
  function(x, y)
    merge(x, y, by = 'pnr', all.x = TRUE),
  list(
    responses_covariates,
    background_covariates_clean[, .(pnr, weights_background)]
  )
)

responses_weights <- responses_weights[!is.na(weights_background)]


# Cap weights
responses_weights[, weights_background_capped := pmin(weights_background, quantile(weights_background, 0.95))]

# Rescale weights
responses_weights[, weights_background_rescaled := weights_background_capped / mean(weights_background_capped)]


# ------------------ Descriptive analysis ------
# Baseline characteristics and data sharing responses

## ------------------ Unweighted baseline characteristics  ------

basetable_unweighted <-
  responses_weights %>%
  finalfit::summary_factorlist(
    dependent = "group_data_sharing_q6",
    explanatory = c("sex",
                    "age",
                    "education_years",
                    "diabetes"),
    column = TRUE,
    total_col = TRUE,
    add_col_totals = TRUE,
    weights = NULL
  )

## ------------------ Weighted baseline characteristics  ------

basetable_weighted_background <-
  responses_weights %>%
  finalfit::summary_factorlist(
    dependent = "group_data_sharing_q6",
    explanatory = c("sex", "age", "education_years", "diabetes"),
    column = TRUE,
    total_col = TRUE,
    add_col_totals = TRUE,
    weights = "weights_background_rescaled"
  )


## ------------------ Unweighted response on data sharing ------
# Overall and by randomization to public or private institutions

data_sharing_unweighted <-
  responses_weights %>%
  finalfit::summary_factorlist(
    dependent = "group_data_sharing_q6",
    explanatory = c("aiq6_1", "aiq6_2", "aiq6_3", "aiq6_4"),
    column = TRUE,
    total_col = TRUE,
    add_col_totals = TRUE,
    weights = NULL,
    p = TRUE,
    p_cat = "chisq"
  )

## ------------------ Weighted response on data sharing ------
# Overall and by randomization to public or private institutions

data_sharing_weighted_background <-
  responses_weights %>%
  finalfit::summary_factorlist(
    dependent = "group_data_sharing_q6",
    explanatory = c("aiq6_1", "aiq6_2", "aiq6_3", "aiq6_4"),
    column = TRUE,
    total_col = TRUE,
    add_col_totals = TRUE,
    weights = "weights_background_rescaled",
    p = TRUE,
    p_cat = "chisq"
  )

# ------------------ Multinominal regression of randomiztion between private and public ------
# Model willingness to share different types of health data
# Exposure: group_data_sharing_q6 (public vs private)
# Reference category: "Ved ikke" (Don't know)

## ----------------- Weighted models-------

# Weighted by the entire population in Central Region Denmark age 18-92

### ---------------- Electronic Health records-----

model_text_weighted <- multinom(
  relevel(aiq6_1, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = weights_background_rescaled
)

text_weighted_estimates <- estimates_OR(model_text_weighted) %>%
  mutate(type = "text")


### ---------------- Images-----

model_images_weighted <- multinom(
  relevel(aiq6_2, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = weights_background_rescaled
)

images_weighted_estimates <- estimates_OR(model_images_weighted) %>%
  mutate(type = "images")


### -------Audio---------------

model_audio_weighted <- multinom(
  relevel(aiq6_3, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = weights_background_rescaled
)

audio_weighted_estimates <- estimates_OR(model_audio_weighted) %>%
  mutate(type = "audio")


### ---- Genetic data -----------

model_genetic_unweighted <- multinom(
  relevel(aiq6_4, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = NULL
)

genetic_unweighted_estimates <- estimates_OR(model_genetic_unweighted) %>%
  mutate(type = "genetic")


## ----------------- Unweighted models-------

# Estimates in the survey population in Health in Central Denmark

### ---------------- Electronic Health records-----

model_text_unweighted <- multinom(
  relevel(aiq6_1, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = NULL
)

text_unweighted_estimates <- estimates_OR(model_text_unweighted) %>%
  mutate(type = "text")


### ---------------- Images----

model_images_unweighted <- multinom(
  relevel(aiq6_2, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = NULL
)

images_unweighted_estimates <- estimates_OR(model_images_unweighted) %>%
  mutate(type = "images")


### -------Audio---------------

model_audio_unweighted <- multinom(
  relevel(aiq6_3, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = NULL
)

audio_unweighted_estimates <- estimates_OR(model_audio_unweighted) %>%
  mutate(type = "audio")


### ---- Genetic data ------------

model_genetic_weighted <- multinom(
  relevel(aiq6_4, ref = "Ved ikke") ~ group_data_sharing_q6,
  data = responses_weights,
  weights = weights_background_rescaled
)

genetic_weighted_estimates <- estimates_OR(model_genetic_weighted) %>%
  mutate(type = "genetic")