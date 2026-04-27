#new try to get the confounders in, redoing the entire analysys

library(tidyverse)
library(haven)  # for reading .XPT files

# ── Set your data folder path ──────────────────────────────────────────────────
data_path <- "~/Desktop/CH internship/New versions BRFSS data"

# ── Updated variables to keep — includes both old and new variable names ──────
vars_to_keep <- c(
  "SEQNO",
  "_STATE",
  "IYEAR",
  "_AGE65YR",
  "ADDEPEV3",   # 2019 onwards
  "ADDEPEV2",   # 2015–2018
  "_SMOKER3",
  "CVDSTRK3",
  "_EDUCAG",
  "_INCOMG1",   # 2021 onwards
  "_INCOMG",    # 2015–2020
  "CHCCOPD3",   # 2021 onwards
  "CHCCOPD2",   # 2019–2020
  "CHCCOPD1",   # 2015–2018
  "_RACE",      # 2015–2021
  "_RACE1",     # 2022 onwards
  "_LLCPWT"
)

# ── Helper function to read one XPT file and select only needed columns ────────
read_brfss <- function(filepath, year) {
  read_xpt(filepath) %>%
    rename_with(toupper) %>%
    select(any_of(vars_to_keep)) %>%
    mutate(Year = year,
           source_file = basename(filepath),
           SEQNO = as.character(SEQNO))  # force SEQNO to character to avoid type conflicts
}

# ── Read all files by year, combining V files with main files ─────────────────

brfss_2015 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2015.XPT"),  2015),
  read_brfss(file.path(data_path, "LLCP15V1.XPT"),  2015),
  read_brfss(file.path(data_path, "LLCP15V2.XPT"),  2015)
)

brfss_2016 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2016.XPT"),  2016),
  read_brfss(file.path(data_path, "LLCP16V1.XPT"),  2016),
  read_brfss(file.path(data_path, "LLCP16V2.XPT"),  2016),
  read_brfss(file.path(data_path, "LLCP16V3.XPT"),  2016)
)

brfss_2017 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2017.XPT"),  2017),
  read_brfss(file.path(data_path, "LLCP17V2.XPT"),  2017),
  read_brfss(file.path(data_path, "LLCP17V3.XPT"),  2017)
)

brfss_2018 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2018.XPT"),  2018),
  read_brfss(file.path(data_path, "LLCP18V1.XPT"),  2018),
  read_brfss(file.path(data_path, "LLCP18V2.XPT"),  2018)
)

brfss_2019 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2019.XPT"),  2019),
  read_brfss(file.path(data_path, "LLCP19V1.XPT"),  2019),
  read_brfss(file.path(data_path, "LLCP19V2.XPT"),  2019),
  read_brfss(file.path(data_path, "LLCP19V3.XPT"),  2019)
)

brfss_2020 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2020.XPT"),  2020),
  read_brfss(file.path(data_path, "LLCP20V1.XPT"),  2020),
  read_brfss(file.path(data_path, "LLCP20V2.XPT"),  2020)
)

brfss_2021 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2021.XPT"),  2021),
  read_brfss(file.path(data_path, "LLCP21V1.XPT"),  2021),
  read_brfss(file.path(data_path, "LLCP21V2.XPT"),  2021)
)

brfss_2022 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2022.XPT"),  2022),
  read_brfss(file.path(data_path, "LLCP22V1.XPT"),  2022),
  read_brfss(file.path(data_path, "LLCP22V2.XPT"),  2022)
)

brfss_2023 <- bind_rows(
  read_brfss(file.path(data_path, "LLCP2023.XPT"),  2023),
  read_brfss(file.path(data_path, "LLCP23V2.XPT"),  2023),
  read_brfss(file.path(data_path, "LLCP23V3.XPT"),  2023)
)

# ── Stack all years into one dataset ──────────────────────────────────────────
brfss_all <- bind_rows(
  brfss_2015, brfss_2016, brfss_2017, brfss_2018, brfss_2019,
  brfss_2020, brfss_2021, brfss_2022, brfss_2023
)

# ── Quick check ───────────────────────────────────────────────────────────────
brfss_all %>%
  count(Year) %>%
  print()


# ── Add a state name lookup (FIPS code to state name) ─────────────────────────
fips_to_state <- tibble(
  `_STATE` = c(1,2,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,
               26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,44,45,46,
               47,48,49,50,51,53,54,55,56,66,72,78),
  StateName = c("Alabama","Alaska","Arizona","Arkansas","California","Colorado",
                "Connecticut","Delaware","District of Columbia","Florida","Georgia",
                "Hawaii","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky",
                "Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota",
                "Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire",
                "New Jersey","New Mexico","New York","North Carolina","North Dakota",
                "Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island","South Carolina",
                "South Dakota","Tennessee","Texas","Utah","Vermont","Virginia",
                "Washington","West Virginia","Wisconsin","Wyoming","Guam",
                "Puerto Rico","Virgin Islands")
)

# ── Compute state-level confounder averages by year ───────────────────────────
state_confounders <- brfss_all %>%
  rename(STATE = `_STATE`) %>%
  left_join(fips_to_state, by = c("STATE" = "_STATE")) %>%
  mutate(
    # depression
    depression_any = case_when(
      !is.na(ADDEPEV3) ~ as.numeric(ADDEPEV3 == 1),
      !is.na(ADDEPEV2) ~ as.numeric(ADDEPEV2 == 1),
      TRUE ~ NA_real_
    ),
    # COPD — three versions
    copd_any = case_when(
      !is.na(CHCCOPD3) ~ as.numeric(CHCCOPD3 == 1),
      !is.na(CHCCOPD2) ~ as.numeric(CHCCOPD2 == 1),
      !is.na(CHCCOPD1) ~ as.numeric(CHCCOPD1 == 1),
      TRUE ~ NA_real_
    ),
    # income
    income_low = case_when(
      !is.na(`_INCOMG1`) ~ as.numeric(`_INCOMG1` %in% c(1, 2)),
      !is.na(`_INCOMG`)  ~ as.numeric(`_INCOMG`  %in% c(1, 2)),
      TRUE ~ NA_real_
    ),
    # race — combine _RACE and _RACE1 into one column
    race_combined = case_when(
      !is.na(`_RACE`)  ~ `_RACE`,
      !is.na(`_RACE1`) ~ `_RACE1`,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(StateName, Year) %>%
  summarise(
    pct_age65plus  = mean(`_AGE65YR` == 1,         na.rm = TRUE) * 100,
    pct_depression = mean(depression_any,           na.rm = TRUE) * 100,
    pct_smoker     = mean(`_SMOKER3` %in% c(1, 2), na.rm = TRUE) * 100,
    pct_stroke     = mean(CVDSTRK3 == 1,            na.rm = TRUE) * 100,
    pct_lowedu     = mean(`_EDUCAG` == 1,           na.rm = TRUE) * 100,
    pct_lowincome  = mean(income_low,               na.rm = TRUE) * 100,
    pct_copd       = mean(copd_any,                 na.rm = TRUE) * 100,
    pct_nonwhite   = mean(race_combined != 1,       na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  filter(!is.na(StateName))

# check for any remaining NaN values
state_confounders %>%
  summarise(across(everything(), ~ sum(is.nan(.x))))

# save
write_csv(state_confounders,
          "~/Desktop/CH internship/New versions BRFSS data/state_confounders_2015_2023.csv")
