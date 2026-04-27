#This script will be for the ozone analysis


# Ozone mirror analysis (same structure as your PM2.5/PM10 workflow)

# =========================
# 0) SETTINGS
# =========================
fast_mode <- FALSE          # TRUE for quick test, FALSE for final
run_imputation <- TRUE      # FALSE to skip imputation section while testing

ozone_dir  <- "/Users/aakashparekh/Desktop/Epi internship/Ozone values"  # <-- edit
scd_file   <- "/Users/aakashparekh/Desktop/Epi internship/R_cognitive_decline_by_state_2015_to_2023.csv"
ozone_cache <- "/Users/aakashparekh/Desktop/Epi internship/ozone_all_states_years.csv"

m_imp <- if (fast_mode) 3 else 10
maxit_imp <- if (fast_mode) 5 else 10

# =========================
# 1) LIBRARIES
# =========================
library(tidyverse)
library(janitor)
library(zoo)
library(mice)
library(lmerTest)

# =========================
# 2) STATE CODE LOOKUP
# =========================
state_lookup <- tibble(
  state_code = c(1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,48,49,50,51,53,54,55,56),
  StateName = c("Alabama","Alaska","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","District Of Columbia",
                "Florida","Georgia","Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland",
                "Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire",
                "New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania",
                "Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virginia","Washington",
                "West Virginia","Wisconsin","Wyoming")
)

# =========================
# 3) BUILD OR LOAD OZONE STATE-YEAR DATA
# =========================
if (file.exists(ozone_cache)) {
  message("Loading cached ozone state-year file...")
  ozone_all_states_years <- read_csv(ozone_cache, show_col_types = FALSE)
} else {
  message("Building ozone state-year file from raw files...")
  ozone_files <- list.files(
    ozone_dir,
    pattern = "^8hour_44201_[0-9]{4}(\\.csv)?$",
    full.names = TRUE
  )
  if (length(ozone_files) == 0) stop("No ozone files found. Check ozone_dir and filenames.")
  
  ozone_state_year_list <- lapply(ozone_files, function(f) {
    message("Reading: ", basename(f))
    df <- read_csv(f, show_col_types = FALSE) %>% clean_names()
    
    # If state_name exists, use it. Otherwise map from state_code.
    if ("state_name" %in% names(df)) {
      out <- df %>%
        transmute(
          StateName = str_squish(str_to_title(state_name)),
          DateLocalRaw = as.character(date_local),
          MeanOzone = as.numeric(mean_including_all_data)
        )
    } else {
      out <- df %>%
        transmute(
          state_code = as.integer(state_code),
          DateLocalRaw = as.character(date_local),
          MeanOzone = as.numeric(mean_including_all_data)
        ) %>%
        left_join(state_lookup, by = "state_code") %>%
        select(StateName, DateLocalRaw, MeanOzone)
    }
    
    out %>%
      mutate(
        DateLocal = suppressWarnings(as.Date(DateLocalRaw, format = "%m/%d/%y")),
        DateLocal = coalesce(DateLocal, suppressWarnings(as.Date(DateLocalRaw, format = "%m/%d/%Y"))),
        DateLocal = coalesce(DateLocal, suppressWarnings(as.Date(DateLocalRaw, format = "%Y-%m-%d"))),
        Year = as.integer(format(DateLocal, "%Y"))
      ) %>%
      filter(Year >= 2010, Year <= 2023, !is.na(StateName)) %>%
      group_by(StateName, Year) %>%
      summarise(Mean_O3 = mean(MeanOzone, na.rm = TRUE), .groups = "drop")
  })
  
  ozone_all_states_years <- bind_rows(ozone_state_year_list) %>%
    group_by(StateName, Year) %>%
    summarise(Mean_O3 = mean(Mean_O3, na.rm = TRUE), .groups = "drop") %>%
    arrange(StateName, Year) %>%
    group_by(StateName) %>%
    mutate(
      O3_lag1 = lag(Mean_O3, 1),
      O3_lag2 = lag(Mean_O3, 2),
      O3_lag3 = lag(Mean_O3, 3)
    ) %>%
    ungroup()
  
  write_csv(ozone_all_states_years, ozone_cache)
}

# =========================
# 4) SCD PREP
# =========================
Total_Cog_Decline <- read_csv(scd_file, show_col_types = FALSE)

Total_Cog_Decline_long <- Total_Cog_Decline %>%
  clean_names() %>%
  rename(StateName = state_name) %>%
  pivot_longer(
    cols = starts_with("x"),
    names_to = "Year",
    names_prefix = "x",
    values_to = "MeanCognitiveDecline"
  ) %>%
  mutate(
    Year = as.integer(Year),
    StateName = str_squish(str_to_title(StateName))
  )

ozone_all_states_years <- ozone_all_states_years %>%
  mutate(StateName = str_squish(str_to_title(StateName)))

state_data_counts <- Total_Cog_Decline_long %>%
  group_by(StateName) %>%
  summarise(years_with_data = sum(!is.na(MeanCognitiveDecline)), .groups = "drop")

states_with_enough_data <- state_data_counts %>%
  filter(years_with_data > 2) %>%
  pull(StateName)

Total_Cog_Decline_long_filtered <- Total_Cog_Decline_long %>%
  filter(StateName %in% states_with_enough_data)

Total_Cog_Decline_long_interp <- Total_Cog_Decline_long %>%
  group_by(StateName) %>%
  arrange(Year) %>%
  mutate(InterpolatedLine = na.approx(MeanCognitiveDecline, Year, na.rm = FALSE)) %>%
  ungroup()

# =========================
# 5) MERGE + RAW MODELS
# =========================
merged_long_data_o3 <- Total_Cog_Decline_long_interp %>%
  left_join(ozone_all_states_years, by = c("StateName", "Year"))

filtered_merged_long_data_o3 <- Total_Cog_Decline_long_filtered %>%
  left_join(ozone_all_states_years, by = c("StateName", "Year"))

mod1_full_ozone <- merged_long_data_o3 %>% rename(O3_lag = O3_lag1) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + O3_lag + (1 | StateName), data = .)
mod2_full_ozone <- merged_long_data_o3 %>% rename(O3_lag = O3_lag2) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + O3_lag + (1 | StateName), data = .)
mod3_full_ozone <- merged_long_data_o3 %>% rename(O3_lag = O3_lag3) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + O3_lag + (1 | StateName), data = .)

summary(mod1_full_ozone); summary(mod2_full_ozone); summary(mod3_full_ozone)

mod1_filt_ozone <- filtered_merged_long_data_o3 %>% rename(O3_lag = O3_lag1) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + O3_lag + (1 | StateName), data = .)
mod2_filt_ozone <- filtered_merged_long_data_o3 %>% rename(O3_lag = O3_lag2) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + O3_lag + (1 | StateName), data = .)
mod3_filt_ozone <- filtered_merged_long_data_o3 %>% rename(O3_lag = O3_lag3) %>%
  lmer(MeanCognitiveDecline * 100 ~ Year + O3_lag + (1 | StateName), data = .)

summary(mod1_filt_ozone); summary(mod2_filt_ozone); summary(mod3_filt_ozone)

# =========================
# 6) IMPUTATION + POOLED MODELS
# =========================
if (run_imputation) {
  wide_merged_data_o3 <- merged_long_data_o3 %>%
    select(-InterpolatedLine) %>%
    pivot_wider(
      values_from = c("MeanCognitiveDecline", "Mean_O3", "O3_lag1", "O3_lag2", "O3_lag3"),
      names_from = "Year"
    )
  
  ini <- mice(wide_merged_data_o3, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  
  meth[] <- ""
  meth[grep("^MeanCognitiveDecline_", names(wide_merged_data_o3))] <- "pmm"
  
  imputed_data_o3 <- mice(
    wide_merged_data_o3,
    m = m_imp,
    maxit = maxit_imp,
    method = meth,
    predictorMatrix = pred,
    seed = 8675309,
    printFlag = TRUE
  )
  
  imputed_complete_o3 <- complete(imputed_data_o3, action = "long", include = FALSE)
  
  cog_long_o3 <- imputed_complete_o3 %>%
    pivot_longer(cols = starts_with("MeanCognitiveDecline_"),
                 names_to = "Year", names_prefix = "MeanCognitiveDecline_",
                 values_to = "MeanCognitiveDecline")
  
  o3_long <- imputed_complete_o3 %>%
    pivot_longer(cols = starts_with("Mean_O3_"),
                 names_to = "Year", names_prefix = "Mean_O3_",
                 values_to = "Mean_O3")
  
  o3_lag1_long <- imputed_complete_o3 %>%
    pivot_longer(cols = starts_with("O3_lag1_"),
                 names_to = "Year", names_prefix = "O3_lag1_",
                 values_to = "O3_lag1")
  
  o3_lag2_long <- imputed_complete_o3 %>%
    pivot_longer(cols = starts_with("O3_lag2_"),
                 names_to = "Year", names_prefix = "O3_lag2_",
                 values_to = "O3_lag2")
  
  o3_lag3_long <- imputed_complete_o3 %>%
    pivot_longer(cols = starts_with("O3_lag3_"),
                 names_to = "Year", names_prefix = "O3_lag3_",
                 values_to = "O3_lag3")
  
  imputed_long_o3 <- cog_long_o3 %>%
    left_join(o3_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    left_join(o3_lag1_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    left_join(o3_lag2_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    left_join(o3_lag3_long, by = c("StateName", "Year", ".imp", ".id")) %>%
    select(.imp, .id, StateName, Year, MeanCognitiveDecline, Mean_O3, O3_lag1, O3_lag2, O3_lag3) %>%
    mutate(Year = as.integer(Year))
  
  imp_list <- split(imputed_long_o3, imputed_long_o3$.imp)
  
  fit_lag1_ozone <- lapply(imp_list, function(d) {
    lmer(MeanCognitiveDecline * 100 ~ as.numeric(Year) + O3_lag1 + (1 | StateName), data = d)
  })
  fit_lag2_ozone <- lapply(imp_list, function(d) {
    lmer(MeanCognitiveDecline * 100 ~ as.numeric(Year) + O3_lag2 + (1 | StateName), data = d)
  })
  fit_lag3_ozone <- lapply(imp_list, function(d) {
    lmer(MeanCognitiveDecline * 100 ~ as.numeric(Year) + O3_lag3 + (1 | StateName), data = d)
  })
  
  print(summary(pool(as.mira(fit_lag1_ozone))))
  print(summary(pool(as.mira(fit_lag2_ozone))))
  print(summary(pool(as.mira(fit_lag3_ozone))))
}

summary(mod1_full_ozone)
summary(mod2_full_ozone)
summary(mod3_full_ozone)
summary(mod1_filt_ozone)
summary(mod2_filt_ozone)
summary(mod3_filt_ozone)
summary(pool(as.mira(fit_lag1_ozone)))
summary(pool(as.mira(fit_lag2_ozone)))
summary(pool(as.mira(fit_lag3_ozone)))
