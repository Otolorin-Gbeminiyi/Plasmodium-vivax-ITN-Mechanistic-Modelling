# Chapter 6: Mechanistic modelling of P. vivax transmission and ITN effects ----

# 6.1 Model configuration and assumptions ----
# 6.1.1 Runtime controls and required packages ----
MODEL_USE_SPEEDUPS <- TRUE
MODEL_VERIFY_SPEEDUPS <- TRUE
MODEL_CACHE_ENTRIES <- 64L

required_packages <- c(
  "deSolve",
  "dplyr",
  "ggplot2",
  "lubridate"
)
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
# 6.1.2 Reproducibility settings and output helpers ----
MODEL_RANDOM_SEED <- 20260811
set.seed(MODEL_RANDOM_SEED)
print_block <- function(title) {
  cat("\n")
  cat("====================================================================\n")
  cat(title, "\n")
  cat("====================================================================\n")
}
# 6.1.3 Human population sizes ----
N_h_control      <- 5332
N_h_intervention <- 5462
# 6.1.4 Model burn-in period ----
burnin_days <- 365 * 15
# 6.1.5 Survey-round aggregation and surveillance assumptions ----
# 
ROUND_COLUMN <- NULL

# To maintain a consistent observation process in the primary
# calibration likelihood, all October records are excluded
#(because malaria record were more of passive detection) 
# 
TREAT_ALL_OCTOBER_AS_PASSIVE <- TRUE
PASSIVE_MONTH <- 10
PASSIVE_POSITIVITY_THRESHOLD <- 0.999999
#
# Default is I_only. Change only if your case/test definition supports it.
DETECTION_MODE <- "I_only"
# Fixed diagnostic characteristics.
diagnostic_sensitivity <- 1.0
diagnostic_specificity <- 1.0
# 6.1.7 Seasonality and background infection assumptions ----

seasonal_amp_fixed      <- 0.00
seasonal_peak_day_fixed <- 90
lambda_import_fixed     <- 0.00
# 6.1.8 Fixed biological parameters ----
# Hypnozoite parameters are fixed sensitivity assumptions, NOT fitted.
bio_parms <- c(
  tau_mh = 0.05,       # mosquito -> human transmission probability / infectious bite
  tau_hm = 0.023,      # human -> mosquito transmission probability / infectious bite
  mu_m   = 0.10,       # baseline adult mosquito mortality / day
  e_m    = 0.10,       # mosquito exposed -> infectious progression / day
  rho_h  = 1 / 12,     # human exposed -> infectious progression / day
  c_h    = 0.10,       # human infectious -> treated/recovered / day
  pi_h   = 1e-5,       # other human removal / day
  w_h    = 1 / 14,     # fixed T_h exit rate / day; mean residence 14 days

  # Relapse-eligible latent carriers already included in S_h
  p_ret     = 0.50,        # fixed proportion leaving T_h retaining hypnozoites
  alpha_hyp = 1 / 90,      # activation hazard per eligible carrier / day
  delta_hyp = 1 / 365      # loss of latent carriage while in S_h / day
)
# Mosquito : human population ratio used only for initialisation
m_ratio_fixed <- 5
# Initial prevalence used only to seed burn-in
prev_h0   <- 0.03
prev_hyp0 <- 0.03  # H_h/N_h; these carriers are already included in S_h
prev_m0   <- 0.08
# 6.1.9 Optimisation settings and parameter bounds ----
n_multistart <- 40
# Baseline b bounds
log_b_lower <- log(0.02)
log_b_upper <- log(5.00)
# q_ITN = 2 * plogis(eta_q), allowing 0 < q_ITN < 2.
# q=1 means no intervention effect.
# q=0.70 means 30% reduction in effective biting.
# q=1.30 means 30% increase in effective biting (possible harm).
# Preserve the original numerical lower limit; the upper limit approaches 2.
eta_q_lower <- stats::qlogis(stats::plogis(-8) / 2)
eta_q_upper <-  8
# Beta-binomial precision phi.
# Intraclass/overdispersion correlation rho = 1/(phi+1).
# Higher phi => closer to ordinary binomial.
log_phi_lower <- log(2)
log_phi_upper <- log(1e5)
RUN_ITN_PROFILE <- TRUE
# 6.1.10 q_ITN profile-likelihood grid ----

q_profile_grid <- sort(
  unique(
    c(
      2 * stats::plogis(eta_q_lower),
      0.001,
      0.005,
      seq(0.01, 1.99, by = 0.01),
      1,
      2 * stats::plogis(eta_q_upper)
    )
  )
)
RUN_HYPNOZOITE_SENSITIVITY <- FALSE
# 6.1.11 Secondary sensitivity-analysis settings ----

RUN_LIMITED_SENSITIVITY <- TRUE
# One-way sensitivity values around the fixed primary assumptions.
SENSITIVITY_P_RET_VALUES <- c(0.25, 0.75)
SENSITIVITY_T_RESIDENCE_DAYS <- c(7, 28, 180)
SENSITIVITY_RELAPSE_MEAN_DAYS <- c(60, 120)
SENSITIVITY_RESERVOIR_MEAN_DAYS <- c(180, 730)
SENSITIVITY_SEASONAL_AMPLITUDES <- c(0.20)
SENSITIVITY_IMPORT_RATES <- c(1e-6)
# Diagnostic values should only be added if assay-specific values are defensible.
DIAGNOSTIC_SENSITIVITY_VALUES <- numeric(0)
DIAGNOSTIC_SPECIFICITY_VALUES <- numeric(0)
# 6.2 Trial screening data preparation ----
# 6.2.1 Control and intervention data import ----
cat("\nChoose CONTROL CSV file...\n")
control_raw <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
)
cat("\nChoose INTERVENTION CSV file...\n")
intervention_raw <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
)
# 6.2.2 Survey date parsing ----
parse_model_date <- function(x) {

  if (inherits(x, "Date")) {
    return(x)
  }

  x_chr <- trimws(as.character(x))

  candidates <- list(
    iso = suppressWarnings(as.Date(x_chr)),
    ymd = suppressWarnings(lubridate::ymd(x_chr)),
    dmy = suppressWarnings(lubridate::dmy(x_chr)),
    mdy = suppressWarnings(lubridate::mdy(x_chr))
  )

  na_counts <- vapply(
    candidates,
    function(z) sum(is.na(z)),
    integer(1)
  )

  preference <- c("iso", "ymd", "dmy", "mdy")
  best_na <- min(na_counts)

  best_name <- preference[
    preference %in% names(na_counts)[na_counts == best_na]
  ][1]

  out <- candidates[[best_name]]

  if (any(is.na(out))) {
    stop("One or more dates could not be parsed.")
  }

  cat("Date format selected:", best_name, "\n")

  out
}
# 6.2.3 Screening input validation ----
standardise_input <- function(tab, arm_name) {

  required_cols <- c(
    "Date_of_the_survey",
    "number_tested",
    "number_positive"
  )

  missing_cols <- setdiff(
    required_cols,
    names(tab)
  )

  if (length(missing_cols) > 0) {
    stop(
      arm_name,
      " missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  tab$Date_of_the_survey <- parse_model_date(
    tab$Date_of_the_survey
  )

  tab$number_tested <- suppressWarnings(
    as.numeric(tab$number_tested)
  )

  tab$number_positive <- suppressWarnings(
    as.numeric(tab$number_positive)
  )

  if ("cum_positive" %in% names(tab)) {
    tab$cum_positive <- suppressWarnings(
      as.numeric(tab$cum_positive)
    )
  }

  tab <- tab[
    order(tab$Date_of_the_survey),
    ,
    drop = FALSE
  ]

  bad <- which(
    !is.finite(tab$number_tested) |
      !is.finite(tab$number_positive) |
      tab$number_tested < 0 |
      tab$number_positive < 0 |
      tab$number_positive > tab$number_tested
  )

  if (length(bad) > 0) {
    stop(
      arm_name,
      " has invalid screening counts at row(s): ",
      paste(head(bad, 20), collapse = ", ")
    )
  }

  tab
}
control_raw <- standardise_input(
  control_raw,
  "CONTROL"
)
intervention_raw <- standardise_input(
  intervention_raw,
  "INTERVENTION"
)
study_start_date <- min(
  min(control_raw$Date_of_the_survey, na.rm = TRUE),
  min(intervention_raw$Date_of_the_survey, na.rm = TRUE)
)
study_start_doy <- lubridate::yday(
  study_start_date
)
# 6.2.4 Surveillance aggregation and passive-case exclusions ----
prepare_screening_data <- function(
    tab,
    arm_name,
    round_column = ROUND_COLUMN
) {

  screen <- tab %>%
    dplyr::filter(
      is.finite(number_tested),
      number_tested > 0
    ) %>%
    dplyr::arrange(
      Date_of_the_survey
    )

  if (nrow(screen) == 0) {
    stop(
      arm_name,
      ": no screening rows with number_tested > 0."
    )
  }

  if (!is.null(round_column)) {

    if (!round_column %in% names(screen)) {
      stop(
        "ROUND_COLUMN = '",
        round_column,
        "' is not present in ",
        arm_name,
        " data."
      )
    }

    screen <- screen %>%
      dplyr::group_by(
        .data[[round_column]]
      ) %>%
      dplyr::summarise(
        survey_date =
          as.Date(
            round(
              mean(
                as.numeric(
                  Date_of_the_survey
                )
              )
            ),
            origin = "1970-01-01"
          ),

        number_tested =
          sum(
            number_tested,
            na.rm = TRUE
          ),

        number_positive =
          sum(
            number_positive,
            na.rm = TRUE
          ),

        .groups = "drop"
      ) %>%
      dplyr::arrange(
        survey_date
      )

    cat(
      arm_name,
      ": aggregated to ",
      nrow(screen),
      " verified survey rounds.\n",
      sep = ""
    )

  } else {

    screen <- screen %>%
      dplyr::group_by(
        Date_of_the_survey
      ) %>%
      dplyr::summarise(
        survey_date =
          max(
            Date_of_the_survey
          ),

        number_tested =
          sum(
            number_tested,
            na.rm = TRUE
          ),

        number_positive =
          sum(
            number_positive,
            na.rm = TRUE
          ),

        .groups = "drop"
      ) %>%
      dplyr::arrange(
        survey_date
      )

    cat(
      arm_name,
      ": using ",
      nrow(screen),
      " screening dates directly in the likelihood.\n",
      sep = ""
    )
  }

  screen <- screen %>%
    dplyr::mutate(
      arm = arm_name,

      followup_day =
        as.numeric(
          survey_date -
            study_start_date
        ),

      observed_positivity =
        number_positive /
        number_tested
    )

  if (TREAT_ALL_OCTOBER_AS_PASSIVE) {

    screen$passive_case_finding <-
      lubridate::month(
        screen$survey_date
      ) == PASSIVE_MONTH

  } else {

    screen$passive_case_finding <-
      lubridate::month(
        screen$survey_date
      ) == PASSIVE_MONTH &
      screen$observed_positivity >=
        PASSIVE_POSITIVITY_THRESHOLD
  }

  screen$include_in_fit <-
    !screen$passive_case_finding

  screen$surveillance_mode <-
    ifelse(
      screen$passive_case_finding,
      "Passive case finding",
      "Active surveillance"
    )

  print_block(
    paste(
      "SCREENING CHECK:",
      arm_name
    )
  )

  cat(
    "Screening observations:",
    nrow(screen),
    "\n"
  )

  cat(
    "Active observations used in fit:",
    sum(screen$include_in_fit),
    "\n"
  )

  cat(
    "Passive observations excluded:",
    sum(screen$passive_case_finding),
    "\n"
  )

  cat(
    "Active number tested:",
    sum(
      screen$number_tested[
        screen$include_in_fit
      ]
    ),
    "\n"
  )

  cat(
    "Active number positive:",
    sum(
      screen$number_positive[
        screen$include_in_fit
      ]
    ),
    "\n"
  )

  if (
    any(
      screen$passive_case_finding
    )
  ) {

    cat(
      "\nPassive observations retained descriptively:\n"
    )

    print(
      screen %>%
        dplyr::filter(
          passive_case_finding
        ) %>%
        dplyr::select(
          survey_date,
          number_tested,
          number_positive,
          observed_positivity
        )
    )
  }

  screen
}
# 6.2.5 Arm-specific screening datasets ----
control_obs <- prepare_screening_data(
  control_raw,
  "Control"
)
intervention_obs <- prepare_screening_data(
  intervention_raw,
  "Intervention"
)
# 6.2.6 Screening-versus-trial episode reconciliation ----
FULL_TRIAL_CONTROL_EPISODES <- 246
FULL_TRIAL_INTERVENTION_EPISODES <- 112
screening_reconciliation <-
  data.frame(
    arm =
      c(
        "Control",
        "Intervention"
      ),

    active_positive =
      c(
        sum(
          control_obs$number_positive[
            control_obs$include_in_fit
          ],
          na.rm = TRUE
        ),

        sum(
          intervention_obs$number_positive[
            intervention_obs$include_in_fit
          ],
          na.rm = TRUE
        )
      ),

    october_descriptive_positive =
      c(
        sum(
          control_obs$number_positive[
            control_obs$passive_case_finding
          ],
          na.rm = TRUE
        ),

        sum(
          intervention_obs$number_positive[
            intervention_obs$passive_case_finding
          ],
          na.rm = TRUE
        )
      ),

    all_screening_positive =
      c(
        sum(
          control_obs$number_positive,
          na.rm = TRUE
        ),

        sum(
          intervention_obs$number_positive,
          na.rm = TRUE
        )
      ),

    full_trial_episode_total =
      c(
        FULL_TRIAL_CONTROL_EPISODES,
        FULL_TRIAL_INTERVENTION_EPISODES
      ),

    stringsAsFactors =
      FALSE
  )
screening_reconciliation$difference_from_full_trial_episode_total <-
  screening_reconciliation$all_screening_positive -
  screening_reconciliation$full_trial_episode_total
print_block(
  "SCREENING VERSUS FULL-TRIAL POSITIVE RECONCILIATION"
)
print(
  screening_reconciliation
)
write.csv(
  screening_reconciliation,
  "screening_vs_full_trial_positive_reconciliation.csv",
  row.names = FALSE
)
# 6.3 Mechanistic P. vivax transmission model ----
# 6.3.1 Human and mosquito initial conditions ----
make_initial_state <- function(N_h0) {

  if (length(N_h0) != 1L || !is.finite(N_h0) || N_h0 <= 0) {
    stop("Initial human population must be positive and finite.")
  }
  seed_fractions <- c(prev_h0, prev_hyp0, prev_m0)
  if (any(!is.finite(seed_fractions)) || any(seed_fractions < 0) ||
      any(seed_fractions > 1) || prev_h0 + prev_hyp0 > 1) {
    stop("Initial fractions must be in [0,1], with prev_hyp0 <= 1 - prev_h0.")
  }

  infected_h_total <-
    N_h0 *
    prev_h0

  E_h0 <-
    infected_h_total *
    0.30

  T_h0 <-
    infected_h_total *
    0.10

  I_h0 <-
    infected_h_total -
    E_h0 -
    T_h0

  S_h0 <-
    N_h0 -
    E_h0 -
    I_h0 -
    T_h0

  # H_h is the latent-carrier subset of S_h.
  # Do not add it to N_h or subtract it again when constructing S_h.
  H_h0 <-
    N_h0 *
    prev_hyp0

  N_m0 <-
    N_h0 *
    m_ratio_fixed

  infected_m_total <-
    N_m0 *
    prev_m0

  E_m0 <-
    infected_m_total *
    0.40

  I_m0 <-
    infected_m_total -
    E_m0

  S_m0 <-
    N_m0 -
    E_m0 -
    I_m0

  c(
    S_h = S_h0,
    E_h = E_h0,
    I_h = I_h0,
    T_h = T_h0,
    H_h = H_h0,
    inc_h = 0,

    S_m = S_m0,
    E_m = E_m0,
    I_m = I_m0
  )
}
init_control_seed <- make_initial_state(
  N_h_control
)
# 6.3.2 Parameter assembly for trial arms ----
make_model_parms <- function(
    b_m,
    q_itn = 1,
    bio = bio_parms
) {

  if (!"p_ret" %in% names(bio) || length(bio[["p_ret"]]) != 1L ||
      !is.finite(bio[["p_ret"]]) || bio[["p_ret"]] < 0 || bio[["p_ret"]] > 1) {
    stop("p_ret must be the proportion leaving T_h retaining hypnozoites, in [0,1].")
  }
  if (!is.finite(bio[["w_h"]]) || bio[["w_h"]] <= 0 ||
      !is.finite(bio[["alpha_hyp"]]) || bio[["alpha_hyp"]] < 0 ||
      !is.finite(bio[["delta_hyp"]]) || bio[["delta_hyp"]] < 0) {
    stop("w_h must be positive; latent-carrier activation and loss rates must be non-negative.")
  }

  c(
    bio,

    b_m =
      as.numeric(
        b_m
      ),

    # The ONLY fitted intervention mechanism:
    # multiplier on effective biting/contact b.
    q_itn =
      as.numeric(
        q_itn
      ),

    lambda_import =
      as.numeric(
        lambda_import_fixed
      ),

    seasonal_amp =
      as.numeric(
        seasonal_amp_fixed
      ),

    seasonal_peak_day =
      as.numeric(
        seasonal_peak_day_fixed
      )
  )
}
# 6.3.3 Human SEIT and mosquito SEI differential equations ----
vivax_model <- function(
    t,
    state,
    parms
) {

  with(
    as.list(
      c(
        state,
        parms
      )
    ),
    {

      N_h <-
        S_h +
        E_h +
        I_h +
        T_h

      N_m <-
        S_m +
        E_m +
        I_m

      if (
        !is.finite(N_h) ||
        !is.finite(N_m) ||
        N_h <= 0 ||
        N_m <= 0
      ) {
        stop(
          "Invalid population size inside ODE."
        )
      }

      # Fixed/calendar-aligned seasonality.
      calendar_doy <-
        (
          (
            study_start_doy -
            1 +
            t
          ) %% 365
        ) +
        1

      seasonal_multiplier <-
        1 +
        seasonal_amp *
        cos(
          2 *
          pi *
          (
            calendar_doy -
            seasonal_peak_day
          ) /
          365
        )

      seasonal_multiplier <-
        max(
          seasonal_multiplier,
          0.05
        )

      # ONE intervention effect: q_itn acts on b.
      b_eff <-
        b_m *
        seasonal_multiplier *
        q_itn

      # Baseline mosquito mortality retained;
      # no separate ITN mortality effect.
      mu_m_eff <- mu_m

      lambda_h_vector <-
        tau_mh *
        b_eff *
        I_m /
        N_h

      lambda_h <-
        lambda_h_vector +
        lambda_import

      lambda_m <-
        tau_hm *
        b_eff *
        I_h /
        N_h

      primary_inf_h <-
        lambda_h *
        S_h

      # P. vivax hypnozoite relapse pathway
      # H_h already contains only eligible carriers within S_h.
      # Relapse leaves S_h and enters E_h; no S_h/N_h multiplier is needed.
      relapse_h <-
        alpha_hyp *
        H_h

      dS_h <-
        -primary_inf_h -
        relapse_h +
        w_h *
        T_h

      dE_h <-
        primary_inf_h +
        relapse_h -
        rho_h *
        E_h

      dI_h <-
        rho_h *
        E_h -
        (
          pi_h +
          c_h
        ) *
        I_h

      dT_h <-
        c_h *
        I_h -
        w_h *
        T_h

      # H_h tracks eligible carriers, not all liver-stage parasites.
      # Infection/relapse removes a carrier from this subset; it need not
      # clear all hypnozoites. Retention on return is represented by p_ret.
      dH_h <-
        p_ret *
        w_h *
        T_h -
        lambda_h *
        H_h -
        relapse_h -
        delta_hyp *
        H_h

      # Cumulative infection entries into E_h: primary + relapse.
      # This counter is distinct from cumulative survey-positive observations.
      dinc_h <-
        primary_inf_h +
        relapse_h

      # Constant adult mosquito population.
      Lambda_m <-
        mu_m_eff *
        N_m

      dS_m <-
        Lambda_m -
        lambda_m *
        S_m -
        mu_m_eff *
        S_m

      dE_m <-
        lambda_m *
        S_m -
        (
          e_m +
          mu_m_eff
        ) *
        E_m

      dI_m <-
        e_m *
        E_m -
        mu_m_eff *
        I_m

      list(
        c(
          dS_h,
          dE_h,
          dI_h,
          dT_h,
          dH_h,
          dinc_h,
          dS_m,
          dE_m,
          dI_m
        ),

        b_eff =
          b_eff,

        primary_inf_h =
          primary_inf_h,

        relapse_h =
          relapse_h,

        lambda_h_vector =
          lambda_h_vector,

        N_h =
          N_h,

        N_m =
          N_m
      )
    }
  )
}
# 6.3.4 Numerical ODE solution and state validity checks ----
run_model <- function(
    init_state,
    parms,
    days,
    arm_name = "ARM"
) {

  # The solver's existing 1e-5 count tolerance also applies to H_h <= S_h.
  if (!all(c("S_h", "H_h") %in% names(init_state)) ||
      any(!is.finite(init_state)) || init_state[["H_h"]] < -1e-5 ||
      init_state[["H_h"]] > init_state[["S_h"]] + 1e-5) {
    return(NULL)
  }

  output <- tryCatch(
    deSolve::ode(
      y = init_state,

      times = seq(
        0,
        ceiling(days),
        by = 1
      ),

      func = vivax_model,

      parms = parms,

      method = "lsoda",

      rtol = 1e-7,

      atol = 1e-8
    ),

    error = function(e) {

      message(
        arm_name,
        " ODE error: ",
        e$message
      )

      NULL
    }
  )

  if (is.null(output)) {
    return(NULL)
  }

  output <-
    as.data.frame(
      output
    )

  state_cols <- c(
    "S_h",
    "E_h",
    "I_h",
    "T_h",
    "H_h",
    "S_m",
    "E_m",
    "I_m"
  )

  if (
    any(
      !is.finite(
        as.matrix(
          output[
            ,
            state_cols,
            drop = FALSE
          ]
        )
      )
    )
  ) {
    return(NULL)
  }

  if (
    any(
      as.matrix(
        output[
          ,
          state_cols,
          drop = FALSE
        ]
      ) < -1e-5
    )
  ) {
    return(NULL)
  }

  # Reject invalid trajectories rather than clipping the carrier count.
  if (any(output$H_h > output$S_h + 1e-5)) {
    return(NULL)
  }

  output
}
# 6.3.5 Extraction of post-burn-in model states ----
extract_end_state <- function(
    model_output
) {

  state_names <- c(
    "S_h",
    "E_h",
    "I_h",
    "T_h",
    "H_h",
    "inc_h",
    "S_m",
    "E_m",
    "I_m"
  )

  state <-
    as.numeric(
      model_output[
        nrow(model_output),
        state_names
      ]
    )

  names(state) <-
    state_names

  state["inc_h"] <- 0

  state
}
# 6.3.6 Population rescaling between trial arms ----
rescale_state <- function(
    state,
    target_N_h
) {

  human_names <- c(
    "S_h",
    "E_h",
    "I_h",
    "T_h"
  )

  mosquito_names <- c(
    "S_m",
    "E_m",
    "I_m"
  )

  old_N_h <-
    sum(
      state[
        human_names
      ]
    )

  old_N_m <-
    sum(
      state[
        mosquito_names
      ]
    )

  mosquito_human_ratio <-
    old_N_m /
    old_N_h

  out <- state

  out[
    human_names
  ] <-
    state[
      human_names
    ] /
    old_N_h *
    target_N_h

  # Scale H_h by the same factor as S_h, preserving the subset relation.
  out["H_h"] <-
    state["H_h"] /
    old_N_h *
    target_N_h

  target_N_m <-
    target_N_h *
    mosquito_human_ratio

  out[
    mosquito_names
  ] <-
    state[
      mosquito_names
    ] /
    old_N_m *
    target_N_m

  out["inc_h"] <- 0

  out
}
# 6.3.7 Interpolation of model states at survey dates ----
interp_state <- function(
    model_output,
    t,
    variable
) {

  stats::approx(
    x =
      model_output$time,

    y =
      model_output[[variable]],

    xout =
      t,

    rule = 2
  )$y
}
# 6.4 Observation model and calibration likelihood ----
# 6.4.1 Model-based screening positivity ----
predict_screening_probability <- function(
    model_output,
    obs_table
) {

  rows <-
    lapply(
      seq_len(
        nrow(
          obs_table
        )
      ),
      function(i) {

        t_i <-
          max(
            obs_table$followup_day[i],
            0
          )

        S_h_i <-
          interp_state(
            model_output,
            t_i,
            "S_h"
          )

        E_h_i <-
          interp_state(
            model_output,
            t_i,
            "E_h"
          )

        I_h_i <-
          interp_state(
            model_output,
            t_i,
            "I_h"
          )

        T_h_i <-
          interp_state(
            model_output,
            t_i,
            "T_h"
          )

        N_h_i <-
          S_h_i +
          E_h_i +
          I_h_i +
          T_h_i

        if (
          DETECTION_MODE ==
          "I_only"
        ) {

          true_prevalence <-
            I_h_i /
            N_h_i

        } else if (
          DETECTION_MODE ==
          "E_plus_I"
        ) {

          true_prevalence <-
            (
              E_h_i +
              I_h_i
            ) /
            N_h_i

        } else {

          stop(
            "DETECTION_MODE must be 'I_only' or 'E_plus_I'."
          )
        }

        true_prevalence <-
          pmin(
            pmax(
              true_prevalence,
              1e-9
            ),
            1 - 1e-9
          )

        # Apply fixed diagnostic performance.
        observed_probability <-
          diagnostic_sensitivity *
          true_prevalence +
          (
            1 -
            diagnostic_specificity
          ) *
          (
            1 -
            true_prevalence
          )

        observed_probability <-
          pmin(
            pmax(
              observed_probability,
              1e-9
            ),
            1 - 1e-9
          )

        data.frame(
          survey_date =
            obs_table$survey_date[i],

          followup_day =
            obs_table$followup_day[i],

          model_true_prevalence =
            true_prevalence,

          model_probability =
            observed_probability,

          expected_positive =
            obs_table$number_tested[i] *
            observed_probability
        )
      }
    )

  dplyr::bind_rows(
    rows
  )
}
# 6.4.2 Beta-binomial likelihood ----
# Parameterisation:
#   alpha = p * phi
#   beta  = (1-p) * phi
#
# rho = 1/(phi + 1)
#
# phi -> infinity approximates ordinary binomial.
#
beta_binomial_logpmf <- function(
    y,
    n,
    p,
    phi
) {

  p <-
    pmin(
      pmax(
        p,
        1e-9
      ),
      1 - 1e-9
    )

  phi <-
    pmax(
      phi,
      1e-6
    )

  alpha <-
    p *
    phi

  beta <-
    (
      1 -
      p
    ) *
    phi

  lchoose(
    n,
    y
  ) +
    lbeta(
      y +
        alpha,

      n -
        y +
        beta
    ) -
    lbeta(
      alpha,
      beta
    )
}
beta_binomial_nll <- function(
    y,
    n,
    p,
    phi
) {

  -sum(
    beta_binomial_logpmf(
      y = y,
      n = n,
      p = p,
      phi = phi
    ),
    na.rm = TRUE
  )
}
# 6.4.3 Fitted-parameter transformations and bounds ----
# Primary fitted parameters:
#
# b_m      : shared baseline mosquito biting/contact rate
# q_itn    : ONLY intervention effect, multiplier on b
# phi_obs  : shared observation-dispersion nuisance parameter
#
unpack_fit_par <- function(par) {

  list(
    b_m =
      exp(
        as.numeric(
          par[["log_b_m"]]
        )
      ),

    q_itn =
      2 * stats::plogis(
        as.numeric(
          par[["eta_q_itn"]]
        )
      ),

    phi_obs =
      exp(
        as.numeric(
          par[["log_phi_obs"]]
        )
      )
  )
}
lower_par <- c(
  log_b_m =
    log_b_lower,

  eta_q_itn =
    eta_q_lower,

  log_phi_obs =
    log_phi_lower
)
upper_par <- c(
  log_b_m =
    log_b_upper,

  eta_q_itn =
    eta_q_upper,

  log_phi_obs =
    log_phi_upper
)
# 6.4.4 Joint trial-arm negative log-likelihood ----
objective_joint_fit <- function(
    par,
    return_details = FALSE,
    bio = bio_parms
) {

  p <-
    unpack_fit_par(
      par
    )

  control_parms <-
    make_model_parms(
      b_m =
        p$b_m,

      q_itn =
        1,

      bio =
        bio
    )

  intervention_parms <-
    make_model_parms(
      b_m =
        p$b_m,

      q_itn =
        p$q_itn,

      bio =
        bio
    )

  # Baseline/control burn-in
  burnin <-
    run_model(
      init_state =
        init_control_seed,

      parms =
        control_parms,

      days =
        burnin_days,

      arm_name =
        "CONTROL BURN-IN"
    )

  if (is.null(burnin)) {

    return(
      if (
        return_details
      ) {
        NULL
      } else {
        1e100
      }
    )
  }

  control_start <-
    extract_end_state(
      burnin
    )

  max_control_day <-
    max(
      control_obs$followup_day,
      na.rm = TRUE
    )

  control_model <-
    run_model(
      init_state =
        control_start,

      parms =
        control_parms,

      days =
        max_control_day,

      arm_name =
        "CONTROL FOLLOW-UP"
    )

  if (is.null(control_model)) {

    return(
      if (
        return_details
      ) {
        NULL
      } else {
        1e100
      }
    )
  }

  intervention_start <-
    rescale_state(
      control_start,
      N_h_intervention
    )

  max_intervention_day <-
    max(
      intervention_obs$followup_day,
      na.rm = TRUE
    )

  intervention_model <-
    run_model(
      init_state =
        intervention_start,

      parms =
        intervention_parms,

      days =
        max_intervention_day,

      arm_name =
        "INTERVENTION FOLLOW-UP"
    )

  if (is.null(intervention_model)) {

    return(
      if (
        return_details
      ) {
        NULL
      } else {
        1e100
      }
    )
  }

  control_prediction <-
    predict_screening_probability(
      control_model,
      control_obs
    )

  intervention_prediction <-
    predict_screening_probability(
      intervention_model,
      intervention_obs
    )

  control_fit <-
    control_obs %>%
    dplyr::left_join(
      control_prediction,
      by = c(
        "survey_date",
        "followup_day"
      )
    ) %>%
    dplyr::filter(
      include_in_fit
    )

  intervention_fit <-
    intervention_obs %>%
    dplyr::left_join(
      intervention_prediction,
      by = c(
        "survey_date",
        "followup_day"
      )
    ) %>%
    dplyr::filter(
      include_in_fit
    )

  if (
    nrow(control_fit) < 2 ||
    nrow(intervention_fit) < 2
  ) {

    return(
      if (
        return_details
      ) {
        NULL
      } else {
        1e100
      }
    )
  }

  nll_control <-
    beta_binomial_nll(
      y =
        control_fit$number_positive,

      n =
        control_fit$number_tested,

      p =
        control_fit$model_probability,

      phi =
        p$phi_obs
    )

  nll_intervention <-
    beta_binomial_nll(
      y =
        intervention_fit$number_positive,

      n =
        intervention_fit$number_tested,

      p =
        intervention_fit$model_probability,

      phi =
        p$phi_obs
    )

  objective <-
    nll_control +
    nll_intervention

  if (!is.finite(objective)) {
    objective <- 1e100
  }

  if (!return_details) {
    return(objective)
  }

  list(
    objective =
      objective,

    nll_control =
      nll_control,

    nll_intervention =
      nll_intervention,

    fitted =
      p,

    control_parms =
      control_parms,

    intervention_parms =
      intervention_parms,

    control_model =
      control_model,

    intervention_model =
      intervention_model,

    control_fit =
      control_fit,

    intervention_fit =
      intervention_fit,

    control_all =
      control_obs %>%
      dplyr::left_join(
        control_prediction,
        by = c(
          "survey_date",
          "followup_day"
        )
      ),

    intervention_all =
      intervention_obs %>%
      dplyr::left_join(
        intervention_prediction,
        by = c(
          "survey_date",
          "followup_day"
        )
      )
  )
}
# 6.5 Numerical optimisation and primary model calibration ----
# 6.5.1 Multistart parameter generation ----
make_random_start <- function() {

  q_start <-
    runif(
      1,
      0.10,
      1.90
    )

  phi_start <-
    exp(
      runif(
        1,
        log(5),
        log(500)
      )
    )

  c(
    log_b_m =
      runif(
        1,
        log_b_lower,
        log_b_upper
      ),

    eta_q_itn =
      qlogis(
        q_start / 2
      ),

    log_phi_obs =
      log(
        phi_start
      )
  )
}
# 6.5.2 Single-start optimisation procedure ----
fit_one_start <- function(
    start,
    start_id
) {

  fit <-
    tryCatch(
      nlminb(
        start =
          start,

        objective =
          objective_joint_fit,

        lower =
          lower_par,

        upper =
          upper_par,

        control = list(
          iter.max = 500,
          eval.max = 2000,
          rel.tol = 1e-9,
          x.tol = 1e-9
        )
      ),

      error = function(e) {
        NULL
      }
    )

  if (is.null(fit)) {

    return(
      data.frame(
        start_id =
          start_id,

        objective =
          Inf,

        convergence =
          999,

        message =
          "fit failed",

        log_b_m =
          NA_real_,

        eta_q_itn =
          NA_real_,

        log_phi_obs =
          NA_real_,

        b_m =
          NA_real_,

        q_itn =
          NA_real_,

        phi_obs =
          NA_real_
      )
    )
  }

  p <-
    unpack_fit_par(
      fit$par
    )

  data.frame(
    start_id =
      start_id,

    objective =
      fit$objective,

    convergence =
      fit$convergence,

    message =
      fit$message,

    log_b_m =
      as.numeric(
        fit$par[["log_b_m"]]
      ),

    eta_q_itn =
      as.numeric(
        fit$par[["eta_q_itn"]]
      ),

    log_phi_obs =
      as.numeric(
        fit$par[["log_phi_obs"]]
      ),

    b_m =
      p$b_m,

    q_itn =
      p$q_itn,

    phi_obs =
      p$phi_obs
  )
}
# 6.5.3 ODE caching and vectorised screening predictions ----
# Exact reuse of deterministic ODE solutions. No parameter rounding is used.
# Full daily trajectories and all original validity checks are retained.
.model_run_original <- run_model
.model_predict_original <- predict_screening_probability
.model_original_rhs <- vivax_model
.model_speed <- new.env(parent = emptyenv())
.model_speed$cache <- list()
.model_speed$hits <- 0L
.model_speed$misses <- 0L
.model_speed$checks <- 0L
.model_speed$started <- proc.time()[['elapsed']]

clear_model_cache <- function() {
  .model_speed$cache <- list()
  invisible(NULL)
}

.model_cached_run <- function(init_state, parms, days, arm_name = "ARM") {
  # The RHS is deterministic in these inputs and study_start_doy. If the
  # function is replaced interactively, bypass the cache for safety.
  if (!identical(vivax_model, .model_original_rhs)) {
    clear_model_cache()
    return(.model_run_original(init_state, parms, days, arm_name))
  }
  key <- serialize(list(init_state = init_state, parms = parms,
                        days = days, arm_name = arm_name,
                        study_start_doy = study_start_doy), NULL, version = 2)
  entries <- .model_speed$cache
  if (length(entries)) {
    for (i in seq_along(entries)) {
      if (identical(entries[[i]]$key, key)) {
        .model_speed$hits <- .model_speed$hits + 1L
        # Promote the entry; keeping recent finite-difference evaluations is useful.
        if (i != length(entries)) {
          .model_speed$cache <- c(entries[-i], entries[i])
        }
        return(entries[[i]]$value)
      }
    }
  }
  .model_speed$misses <- .model_speed$misses + 1L
  value <- .model_run_original(init_state, parms, days, arm_name)
  # Do not cache failures: error/warning behaviour remains with the original solver.
  if (!is.null(value)) {
    entries <- .model_speed$cache
    if (length(entries) >= MODEL_CACHE_ENTRIES) entries <- entries[-1L]
    .model_speed$cache <- c(entries, list(list(key = key, value = value)))
  }
  value
}

.model_predict_vectorised <- function(model_output, obs_table) {
  if (nrow(obs_table) == 0L) {
    return(.model_predict_original(model_output, obs_table))
  }
  t_i <- pmax(obs_table$followup_day, 0)
  # Keep the same approx() method, default tie handling, and rule = 2.
  S_h_i <- interp_state(model_output, t_i, "S_h")
  E_h_i <- interp_state(model_output, t_i, "E_h")
  I_h_i <- interp_state(model_output, t_i, "I_h")
  T_h_i <- interp_state(model_output, t_i, "T_h")
  # Preserve the original floating-point order; do not use rowSums().
  N_h_i <- S_h_i + E_h_i + I_h_i + T_h_i
  if (DETECTION_MODE == "I_only") {
    true_prevalence <- I_h_i / N_h_i
  } else if (DETECTION_MODE == "E_plus_I") {
    true_prevalence <- (E_h_i + I_h_i) / N_h_i
  } else {
    stop("DETECTION_MODE must be 'I_only' or 'E_plus_I'.")
  }
  true_prevalence <- pmin(pmax(true_prevalence, 1e-9), 1 - 1e-9)
  observed_probability <- diagnostic_sensitivity * true_prevalence +
    (1 - diagnostic_specificity) * (1 - true_prevalence)
  observed_probability <- pmin(pmax(observed_probability, 1e-9), 1 - 1e-9)
  data.frame(
    survey_date = obs_table$survey_date,
    followup_day = obs_table$followup_day,
    model_true_prevalence = true_prevalence,
    model_probability = observed_probability,
    expected_positive = obs_table$number_tested * observed_probability
  )
}

# Use the unchanged objective and all its dependencies in a reference environment.
# This avoids swapping global functions during the verification comparisons.
.model_reference_environment <- new.env(parent = environment(objective_joint_fit))
.model_reference_environment$run_model <- .model_run_original
.model_reference_environment$predict_screening_probability <- .model_predict_original
.model_objective_original <- objective_joint_fit
environment(.model_objective_original) <- .model_reference_environment

if (MODEL_USE_SPEEDUPS) {
  run_model <- .model_cached_run
  predict_screening_probability <- .model_predict_vectorised
}

# 6.5.4 Numerical equivalence checks for runtime optimisations ----
verify_model_speedups <- function() {
  # Fixed test values consume no random draws, preserving all multistart seeds.
  points <- list(
    c(log_b_m = log(1), eta_q_itn = qlogis(0.60 / 2), log_phi_obs = log(50)),
    c(log_b_m = log(1), eta_q_itn = qlogis(0.60 / 2), log_phi_obs = log(75)),
    c(log_b_m = log(1), eta_q_itn = qlogis(0.40 / 2), log_phi_obs = log(50)),
    c(log_b_m = log(1.4594), eta_q_itn = qlogis(0.4885 / 2), log_phi_obs = log(1162.63))
  )
  rows <- vector("list", length(points))
  for (i in seq_along(points)) {
    before <- proc.time()[['elapsed']]
    reference <- .model_objective_original(points[[i]], return_details = TRUE)
    original_seconds <- proc.time()[['elapsed']] - before
    before <- proc.time()[['elapsed']]
    accelerated <- objective_joint_fit(points[[i]], return_details = TRUE)
    accelerated_seconds <- proc.time()[['elapsed']] - before
    if (is.null(reference) || is.null(accelerated)) {
      stop("Runtime verification could not obtain a valid trajectory at test point ",
           i, ". Set MODEL_USE_SPEEDUPS <- FALSE to investigate with the original functions.")
    }
    # Compare every returned object, including both trajectories, predictions,
    # likelihoods, parameters, fit rows, and descriptive October rows.
    comparison <- all.equal(reference, accelerated, tolerance = 0,
                            check.attributes = TRUE)
    exact <- isTRUE(comparison)
    rows[[i]] <- data.frame(test_point = i, exact_match = exact,
                            original_seconds = original_seconds,
                            accelerated_seconds = accelerated_seconds)
    if (!exact) {
      write.csv(do.call(rbind, rows[seq_len(i)]),
                "runtime_equivalence_checks.csv", row.names = FALSE)
      stop("Runtime verification found a difference: ",
           paste(comparison, collapse = "; "),
           ". Set MODEL_USE_SPEEDUPS <- FALSE to use the original computations.")
    }
    .model_speed$checks <- .model_speed$checks + 1L
  }
  check_table <- do.call(rbind, rows)
  write.csv(check_table, "runtime_equivalence_checks.csv", row.names = FALSE)
  cat("\nRuntime check: original and faster results matched exactly at all ",
      length(points), " test points.\n", sep = "")
  # Keep validation work outside fit timing and cache hit statistics.
  clear_model_cache()
  .model_speed$hits <- 0L
  .model_speed$misses <- 0L
  .model_speed$started <- proc.time()[['elapsed']]
  invisible(check_table)
}

if (MODEL_USE_SPEEDUPS && MODEL_VERIFY_SPEEDUPS) verify_model_speedups()

# Same serial order and same random starts; progress does not change fitting.
.model_fit_one_start_original <- fit_one_start
fit_one_start <- function(start, start_id) {
  cat("\nPrimary start ", start_id, " of ", n_multistart, "\n", sep = "")
  flush.console()
  .model_fit_one_start_original(start, start_id)
}

# 6.5.5 Primary 40-start joint calibration ----
print_block(
  "MULTISTART JOINT FIT"
)
starts <-
  vector(
    "list",
    n_multistart
  )
starts[[1]] <-
  c(
    log_b_m =
      log(1.0),

    eta_q_itn =
      qlogis(0.60 / 2),

    log_phi_obs =
      log(50)
  )
if (n_multistart > 1) {

  for (
    i in 2:n_multistart
  ) {

    starts[[i]] <-
      make_random_start()
  }
}
multistart_results <-
  dplyr::bind_rows(
    lapply(
      seq_len(
        n_multistart
      ),
      function(i) {

        fit_one_start(
          starts[[i]],
          i
        )
      }
    )
  ) %>%
  dplyr::arrange(
    objective
  )
write.csv(
  multistart_results,
  "multistart_screening_fit.csv",
  row.names = FALSE
)
successful_fits <-
  multistart_results %>%
  dplyr::filter(
    convergence == 0,
    is.finite(objective)
  )
if (
  nrow(successful_fits) ==
  0
) {

  print(
    head(
      multistart_results,
      10
    )
  )

  stop(
    "No multistart optimisation converged. ",
    "Do not interpret parameters."
  )
}
# 6.5.6 Best converged parameter estimates ----
best_row <-
  successful_fits[
    1,
  ]
best_par <-
  c(
    log_b_m =
      best_row$log_b_m,

    eta_q_itn =
      best_row$eta_q_itn,

    log_phi_obs =
      best_row$log_phi_obs
  )
best_result <-
  objective_joint_fit(
    best_par,
    return_details = TRUE
  )
best <-
  best_result$fitted
print_block(
  "BEST CONVERGED FIT"
)
cat(
  "NLL =",
  best_result$objective,
  "\n"
)
cat(
  "Converged starts =",
  nrow(successful_fits),
  "of",
  n_multistart,
  "\n"
)
cat(
  "b_m =",
  best$b_m,
  "\n"
)
cat(
  "q_ITN =",
  best$q_itn,
  "\n"
)
cat(
  "Reduction in effective biting =",
  100 *
    (
      1 -
      best$q_itn
    ),
  "%\n"
)
cat(
  "Beta-binomial phi =",
  best$phi_obs,
  "\n"
)
cat(
  "Beta-binomial rho =",
  1 /
    (
      best$phi_obs +
      1
    ),
  "\n"
)
# 6.5.7 Multistart robustness diagnostics ----
best_nll <-
  min(
    successful_fits$objective
  )
near_best <-
  successful_fits %>%
  dplyr::filter(
    objective <=
      best_nll +
      2
  )
robustness_summary <-
  data.frame(
    attempted =
      n_multistart,

    converged =
      nrow(
        successful_fits
      ),

    within_2_nll =
      nrow(
        near_best
      ),

    best_nll =
      best_nll,

    q_min_near_best =
      min(
        near_best$q_itn
      ),

    q_median_near_best =
      median(
        near_best$q_itn
      ),

    q_max_near_best =
      max(
        near_best$q_itn
      ),

    b_min_near_best =
      min(
        near_best$b_m
      ),

    b_max_near_best =
      max(
        near_best$b_m
      )
  )
print_block(
  "MULTISTART ROBUSTNESS"
)
print(
  robustness_summary
)
write.csv(
  robustness_summary,
  "multistart_robustness_summary.csv",
  row.names = FALSE
)
# 6.5.8 Fitted-parameter boundary diagnostics ----
boundary_check <-
  data.frame(
    parameter = c(
      "b_m",
      "q_itn",
      "phi_obs"
    ),

    estimate = c(
      best$b_m,
      best$q_itn,
      best$phi_obs
    ),

    near_lower = c(
      best$b_m <=
        exp(log_b_lower) *
        1.05,

      best$q_itn <=
        0.02,

      best$phi_obs <=
        exp(log_phi_lower) *
        1.05
    ),

    near_upper = c(
      best$b_m >=
        exp(log_b_upper) /
        1.05,

      best$q_itn >=
        1.98,

      best$phi_obs >=
        exp(log_phi_upper) /
        1.05
    )
  )
print_block(
  "BOUNDARY CHECK"
)
print(
  boundary_check
)
write.csv(
  boundary_check,
  "parameter_boundary_check.csv",
  row.names = FALSE
)
# 6.6 Calibration diagnostics and empirical fit ----
# 6.6.1 Arm-level observed and expected screening positives ----
summarise_arm <- function(
    fit_df,
    arm_name
) {

  total_tested <-
    sum(
      fit_df$number_tested
    )

  total_observed_positive <-
    sum(
      fit_df$number_positive
    )

  total_expected_positive <-
    sum(
      fit_df$expected_positive
    )

  observed_overall_positivity <-
    total_observed_positive /
    total_tested

  model_overall_positivity <-
    total_expected_positive /
    total_tested

  data.frame(
    arm =
      arm_name,

    active_screening_observations =
      nrow(
        fit_df
      ),

    total_tested =
      total_tested,

    observed_positive =
      total_observed_positive,

    model_expected_positive =
      total_expected_positive,

    observed_overall_positivity =
      observed_overall_positivity,

    model_overall_positivity =
      model_overall_positivity,

    absolute_positivity_error =
      model_overall_positivity -
      observed_overall_positivity,

    percent_error_in_positive_count =
      100 *
      (
        total_expected_positive -
        total_observed_positive
      ) /
      max(
        total_observed_positive,
        1
      )
  )
}
fit_summary <-
  dplyr::bind_rows(
    summarise_arm(
      best_result$control_fit,
      "Control"
    ),

    summarise_arm(
      best_result$intervention_fit,
      "Intervention"
    )
  )
print_block(
  "ACTIVE-SURVEILLANCE FIT SUMMARY"
)
print(
  fit_summary
)
write.csv(
  fit_summary,
  "active_screening_fit_summary.csv",
  row.names = FALSE
)
# 6.6.2 Beta-binomial Pearson residuals ----
make_residuals <- function(
    fit_df,
    arm_name
) {

  rho <-
    1 /
    (
      best$phi_obs +
      1
    )

  fit_df %>%
    dplyr::mutate(
      arm =
        arm_name,

      variance_beta_binomial =
        number_tested *
        model_probability *
        (
          1 -
          model_probability
        ) *
        (
          1 +
          (
            number_tested -
            1
          ) *
          rho
        ),

      pearson_residual =
        (
          number_positive -
          expected_positive
        ) /
        sqrt(
          pmax(
            variance_beta_binomial,
            1e-9
          )
        )
    )
}
residual_data <-
  dplyr::bind_rows(
    make_residuals(
      best_result$control_fit,
      "Control"
    ),

    make_residuals(
      best_result$intervention_fit,
      "Intervention"
    )
  )
write.csv(
  residual_data,
  "active_screening_fit_and_residuals.csv",
  row.names = FALSE
)
# 6.6.3 Extreme date-level residual diagnostics ----
extreme_residuals <-
  residual_data %>%
  dplyr::mutate(
    absolute_pearson_residual =
      abs(
        pearson_residual
      ),

    sparse_denominator =
      number_tested < 20,

    expected_positive_below_one =
      expected_positive < 1
  ) %>%
  dplyr::arrange(
    dplyr::desc(
      absolute_pearson_residual
    )
  ) %>%
  dplyr::select(
    arm,
    survey_date,
    number_tested,
    number_positive,
    observed_positivity,
    model_probability,
    expected_positive,
    pearson_residual,
    absolute_pearson_residual,
    sparse_denominator,
    expected_positive_below_one
  )
print_block(
  "LARGEST ABSOLUTE PEARSON RESIDUALS"
)
print(
  head(
    extreme_residuals,
    20
  )
)
write.csv(
  extreme_residuals,
  "date_level_residual_diagnostics.csv",
  row.names = FALSE
)
# 6.6.4 Primary-period checks and Wilson intervals ----
# Primary empirical fitting plots must contain 2001 data only.
# This is a validation check, not a device for hiding erroneous dates.
primary_fit_dates_2001_ok <-
  all(
    lubridate::year(
      residual_data$survey_date
    ) ==
      2001
  )
if (!primary_fit_dates_2001_ok) {
  stop(
    "Primary fitting diagnostics contain survey dates outside 2001. ",
    "Check the input data/date parsing before plotting."
  )
}
fit_plot_start_date <-
  as.Date("2001-01-01")
fit_plot_end_date <-
  as.Date("2001-12-31")
wilson_interval <- function(
    y,
    n,
    z = 1.96
) {

  p <-
    y /
    n

  denominator <-
    1 +
    z^2 /
    n

  centre <-
    (
      p +
      z^2 /
      (
        2 *
        n
      )
    ) /
    denominator

  half <-
    z *
    sqrt(
      p *
      (
        1 -
        p
      ) /
      n +
      z^2 /
      (
        4 *
        n^2
      )
    ) /
    denominator

  data.frame(
    lower =
      pmax(
        centre -
        half,
        0
      ),

    upper =
      pmin(
        centre +
        half,
        1
      )
  )
}
wilson <-
  wilson_interval(
    residual_data$number_positive,
    residual_data$number_tested
  )
residual_data$observed_lower95 <-
  wilson$lower
residual_data$observed_upper95 <-
  wilson$upper
# 6.6.5 Observed-versus-modelled positivity figure ----
p_positivity <-
  ggplot2::ggplot(
    residual_data,
    ggplot2::aes(
      x =
        survey_date
    )
  ) +

  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin =
        observed_lower95,

      ymax =
        observed_upper95,

      group =
        interaction(
          arm,
          survey_date
        )
    ),
    alpha = 0.10
  ) +

  ggplot2::geom_point(
    ggplot2::aes(
      y =
        observed_positivity,

      size =
        number_tested,

      colour =
        "Observed"
    ),
    alpha = 0.8
  ) +

  ggplot2::geom_line(
    ggplot2::aes(
      y =
        model_probability,

      colour =
        "Model"
    ),
    linewidth = 0.9
  ) +

  ggplot2::facet_wrap(
    ~ arm
  ) +

  ggplot2::scale_x_date(
    limits =
      c(
        fit_plot_start_date,
        fit_plot_end_date
      ),
    date_breaks =
      "2 months",
    date_labels =
      "%b"
  ) +

  ggplot2::scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +

  ggplot2::labs(
    title =
      "Observed versus modelled active-surveillance positivity",

    subtitle =
      "Point size reflects number screened; October descriptive observations excluded from the primary fit",

    x =
      "Survey date (2001)",

    y =
      "Proportion positive",

    colour =
      "",

    size =
      "Number screened"
  ) +

  ggplot2::theme_minimal(
    base_size = 12
  )
print(
  p_positivity
)
ggplot2::ggsave(
  "fit_active_screening_positivity.png",
  p_positivity,
  width = 11,
  height = 6,
  dpi = 300
)
# 6.6.6 Observed-versus-expected positive-count figure ----
p_positive_counts <-
  ggplot2::ggplot(
    residual_data,
    ggplot2::aes(
      x =
        survey_date
    )
  ) +

  ggplot2::geom_point(
    ggplot2::aes(
      y =
        number_positive,

      colour =
        "Observed"
    ),
    size = 2
  ) +

  ggplot2::geom_line(
    ggplot2::aes(
      y =
        expected_positive,

      colour =
        "Model"
    ),
    linewidth = 0.9
  ) +

  ggplot2::facet_wrap(
    ~ arm,
    scales = "free_y"
  ) +

  ggplot2::scale_x_date(
    limits =
      c(
        fit_plot_start_date,
        fit_plot_end_date
      ),
    date_breaks =
      "2 months",
    date_labels =
      "%b"
  ) +

  ggplot2::labs(
    title =
      "Observed versus expected positives per active screening date",

    subtitle =
      "Expected positives = number screened × model-predicted positivity",

    x =
      "Survey date (2001)",

    y =
      "Number positive",

    colour =
      ""
  ) +

  ggplot2::theme_minimal(
    base_size = 12
  )
print(
  p_positive_counts
)
ggplot2::ggsave(
  "fit_active_positive_counts.png",
  p_positive_counts,
  width = 11,
  height = 6,
  dpi = 300
)
# 6.6.7 Residual diagnostic figure ----
p_residual <-
  ggplot2::ggplot(
    residual_data,
    ggplot2::aes(
      x =
        survey_date,

      y =
        pearson_residual
    )
  ) +

  ggplot2::geom_hline(
    yintercept = 0,
    linetype = 2
  ) +

  ggplot2::geom_hline(
    yintercept = c(
      -2,
      2
    ),
    linetype = 3
  ) +

  ggplot2::geom_point(
    size = 2
  ) +

  ggplot2::facet_wrap(
    ~ arm
  ) +

  ggplot2::scale_x_date(
    limits =
      c(
        fit_plot_start_date,
        fit_plot_end_date
      ),
    date_breaks =
      "2 months",
    date_labels =
      "%b"
  ) +

  ggplot2::labs(
    title =
      "Beta-binomial Pearson residuals",

    x =
      "Survey date (2001)",

    y =
      "Pearson residual"
  ) +

  ggplot2::theme_minimal(
    base_size = 12
  )
print(
  p_residual
)
ggplot2::ggsave(
  "fit_active_screening_residuals.png",
  p_residual,
  width = 11,
  height = 5.5,
  dpi = 300
)
# 6.6.8 Cumulative active-surveillance positives ----
cumulative_data <-
  residual_data %>%
  dplyr::group_by(
    arm
  ) %>%
  dplyr::arrange(
    survey_date,
    .by_group = TRUE
  ) %>%
  dplyr::mutate(
    observed_cumulative =
      cumsum(
        number_positive
      ),

    model_cumulative =
      cumsum(
        expected_positive
      )
  ) %>%
  dplyr::ungroup()
# Confirm that the fitted observations belong to the 2001 trial period.
if (
  any(
    lubridate::year(
      cumulative_data$survey_date
    ) !=
      2001
  )
) {
  stop(
    "Cumulative fitting diagnostic contains survey dates outside 2001. ",
    "Check the input data before plotting."
  )
}
plot_start_date <-
  as.Date("2001-01-01")
plot_end_date <-
  as.Date("2001-12-31")
p_cumulative <-
  ggplot2::ggplot(
    cumulative_data,
    ggplot2::aes(
      x =
        survey_date
    )
  ) +

  ggplot2::geom_step(
    ggplot2::aes(
      y =
        observed_cumulative,

      colour =
        "Observed cumulative"
    ),
    linewidth = 0.9
  ) +

  ggplot2::geom_line(
    ggplot2::aes(
      y =
        model_cumulative,

      colour =
        "Model cumulative"
    ),
    linewidth = 0.9
  ) +

  ggplot2::facet_wrap(
    ~ arm,
    scales = "free_y"
  ) +

  ggplot2::scale_x_date(
    limits =
      c(
        plot_start_date,
        plot_end_date
      ),
    date_breaks =
      "2 months",
    date_labels =
      "%b"
  ) +

  ggplot2::labs(
    title =
      "Cumulative active-surveillance positives",

    subtitle =
      paste0(
        "Observed and modelled cumulative positives during the ",
        "2001 fitting period; October descriptive observations excluded"
      ),

    x =
      "Survey date (2001)",

    y =
      "Cumulative positives",

    colour =
      ""
  ) +

  ggplot2::theme_minimal(
    base_size = 12
  ) +

  ggplot2::theme(
    legend.position =
      "bottom",

    axis.text.x =
      ggplot2::element_text(
        angle = 0,
        hjust = 0.5
      )
  )
print(
  p_cumulative
)
ggplot2::ggsave(
  "fit_active_cumulative_positives.png",
  p_cumulative,
  width = 10,
  height = 5.5,
  dpi = 300,
  bg = "white"
)
# 6.6.9 Descriptive reporting of October passive observations ----
passive_october <-
  dplyr::bind_rows(
    best_result$control_all %>%
      dplyr::filter(
        passive_case_finding
      ),

    best_result$intervention_all %>%
      dplyr::filter(
        passive_case_finding
      )
  ) %>%
  dplyr::select(
    arm,
    survey_date,
    number_tested,
    number_positive,
    observed_positivity,
    model_probability,
    expected_positive,
    surveillance_mode
  )
print_block(
  "OCTOBER OBSERVATIONS: DESCRIPTIVE ONLY"
)
print(
  passive_october
)
write.csv(
  passive_october,
  "passive_october_descriptive_only.csv",
  row.names = FALSE
)
# 6.7 Profile likelihood for the ITN contact multiplier ----
# 6.7.1 Profile-likelihood analysis setup ----
if (RUN_ITN_PROFILE) {
  
  print_block(
    "PROFILE LIKELIHOOD FOR q_ITN"
  )
  
  # 6.7.2 Profile-likelihood grid and fitting assumptions ----
  # The primary likelihood and model are unchanged. q_ITN is
  # fixed at each profile point while b_m and phi_obs are
  # re-optimised. 
  
  q_admissible_lower <- 2 * stats::plogis(eta_q_lower)
  q_admissible_upper <- 2 * stats::plogis(eta_q_upper)
  profile_threshold_95 <- 1.92
  primary_minimum_nll <- best_result$objective
  
 
  q_profile_grid_run <- sort(
    unique(
      c(
        q_profile_grid,
        best$q_itn
      )
    )
  )
  
  q_profile_grid_run <- q_profile_grid_run[
    q_profile_grid_run >= q_admissible_lower &
      q_profile_grid_run <= q_admissible_upper
  ]
  
  # 6.7.3 Conditional optimisation at fixed q_ITN ----
  profile_one_q <- function(q_value) {
    cat("\nProfile q_ITN = ", signif(q_value, 7), "\n", sep = "")
    flush.console()
    
    if (
      !is.finite(q_value) ||
      q_value < q_admissible_lower ||
      q_value > q_admissible_upper
    ) {
      return(
        data.frame(
          q_itn = q_value,
          nll = Inf,
          convergence = 999,
          message = "q outside admissible bounds",
          b_m = NA_real_,
          phi_obs = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    
    # q_value is held fixed. Only b_m and phi_obs are re-optimised.
    objective_profile <- function(shared_par) {
      
      full_par <- c(
        log_b_m =
          shared_par[["log_b_m"]],
        
        eta_q_itn =
          qlogis(q_value / 2),
        
        log_phi_obs =
          shared_par[["log_phi_obs"]]
      )
      
      objective_joint_fit(
        full_par
      )
    }
    
    # Start from the best full-model values for the shared parameters.
    start_shared <- c(
      log_b_m =
        best_par[["log_b_m"]],
      
      log_phi_obs =
        best_par[["log_phi_obs"]]
    )
    
    fit_q <- tryCatch(
      nlminb(
        start =
          start_shared,
        
        objective =
          objective_profile,
        
        lower = c(
          log_b_m =
            log_b_lower,
          
          log_phi_obs =
            log_phi_lower
        ),
        
        upper = c(
          log_b_m =
            log_b_upper,
          
          log_phi_obs =
            log_phi_upper
        ),
        
        control = list(
          iter.max = 500,
          eval.max = 2000,
          rel.tol = 1e-9,
          x.tol = 1e-9
        )
      ),
      
      error = function(e) {
        message(
          "Profile optimisation failed for q_ITN = ",
          signif(q_value, 7),
          ": ",
          e$message
        )
        NULL
      }
    )
    
    if (
      is.null(fit_q) ||
      !is.finite(fit_q$objective)
    ) {
      return(
        data.frame(
          q_itn = q_value,
          nll = Inf,
          convergence = 999,
          message = "profile fit failed",
          b_m = NA_real_,
          phi_obs = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    }
    
    data.frame(
      q_itn =
        q_value,
      
      nll =
        fit_q$objective,
      
      convergence =
        fit_q$convergence,
      
      message =
        if (is.null(fit_q$message)) "" else fit_q$message,
      
      b_m =
        exp(
          fit_q$par[["log_b_m"]]
        ),
      
      phi_obs =
        exp(
          fit_q$par[["log_phi_obs"]]
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  # 6.7.4 Profile evaluation and likelihood differences ----
  profile_itn <- dplyr::bind_rows(
    lapply(
      q_profile_grid_run,
      profile_one_q
    )
  )
  
  # Delta NLL must be relative to the full-model primary MLE.
  profile_itn <- profile_itn %>%
    dplyr::mutate(
      delta_nll =
        nll - primary_minimum_nll,
      
      usable_for_ci =
        is.finite(nll) &
        convergence == 0
    )
  
  finite_profile <- profile_itn %>%
    dplyr::filter(
      usable_for_ci
    ) %>%
    dplyr::arrange(
      q_itn
    )
  
  if (nrow(finite_profile) == 0) {
    stop(
      "All q_ITN profile likelihood fits failed or did not converge."
    )
  }
  

  if (
    min(
      finite_profile$delta_nll,
      na.rm = TRUE
    ) < -1e-4
  ) {
    warning(
      "A profile point has NLL materially below the stored primary MLE. ",
      "Re-check the primary multistart fit before interpreting the profile."
    )
  }
  
  # 6.7.5 Likelihood-threshold crossing interpolation ----
  interpolate_profile_crossing <- function(
      q1,
      d1,
      q2,
      d2,
      target = profile_threshold_95
  ) {
    
    if (
      any(!is.finite(c(q1, d1, q2, d2))) ||
      q1 == q2
    ) {
      return(NA_real_)
    }
    
    if (d1 == d2) {
      return(
        mean(
          c(q1, q2)
        )
      )
    }
    
    q_cross <-
      q1 +
      (target - d1) *
      (q2 - q1) /
      (d2 - d1)
    
    pmin(
      pmax(
        q_cross,
        min(q1, q2)
      ),
      max(q1, q2)
    )
  }
  
  q_hat <- best$q_itn
  
  left_profile <- finite_profile %>%
    dplyr::filter(
      q_itn <= q_hat
    ) %>%
    dplyr::arrange(
      q_itn
    )
  
  right_profile <- finite_profile %>%
    dplyr::filter(
      q_itn >= q_hat
    ) %>%
    dplyr::arrange(
      q_itn
    )
  
  lower_cross_indices <- integer(0)
  if (nrow(left_profile) >= 2) {
    lower_cross_indices <- which(
      head(
        left_profile$delta_nll,
        -1
      ) > profile_threshold_95 &
      tail(
        left_profile$delta_nll,
        -1
      ) <= profile_threshold_95
    )
  }
  
  if (length(lower_cross_indices) > 0) {
   
    lower_i <- max(lower_cross_indices)
    
    q_lower_95 <- interpolate_profile_crossing(
      left_profile$q_itn[lower_i],
      left_profile$delta_nll[lower_i],
      left_profile$q_itn[lower_i + 1],
      left_profile$delta_nll[lower_i + 1]
    )
    
    lower_bound_identified <- is.finite(q_lower_95)
    lower_bound_status <- "identified by profile-threshold crossing"
    
  } else {
    q_lower_95 <- NA_real_
    lower_bound_identified <- FALSE
    
    lower_reached_admissible_bound <-
      nrow(left_profile) > 0 &&
      abs(
        min(left_profile$q_itn) - q_admissible_lower
      ) <= 1e-10
    
    if (
      lower_reached_admissible_bound &&
      left_profile$delta_nll[1] <= profile_threshold_95
    ) {
      lower_bound_status <-
        "not identified within admissible q_ITN parameter range"
    } else {
      lower_bound_status <-
        "not identified within successfully profiled lower range"
    }
  }
  
  
  upper_cross_indices <- integer(0)
  if (nrow(right_profile) >= 2) {
    upper_cross_indices <- which(
      head(
        right_profile$delta_nll,
        -1
      ) <= profile_threshold_95 &
      tail(
        right_profile$delta_nll,
        -1
      ) > profile_threshold_95
    )
  }
  
  if (length(upper_cross_indices) > 0) {
    # Use the crossing nearest the MLE.
    upper_i <- min(upper_cross_indices)
    
    q_upper_95 <- interpolate_profile_crossing(
      right_profile$q_itn[upper_i],
      right_profile$delta_nll[upper_i],
      right_profile$q_itn[upper_i + 1],
      right_profile$delta_nll[upper_i + 1]
    )
    
    upper_bound_identified <- is.finite(q_upper_95)
    upper_bound_status <- "identified by profile-threshold crossing"
    
  } else {
    q_upper_95 <- NA_real_
    upper_bound_identified <- FALSE
    
    upper_reached_admissible_bound <-
      nrow(right_profile) > 0 &&
      abs(
        max(right_profile$q_itn) - q_admissible_upper
      ) <= 1e-10
    
    if (
      upper_reached_admissible_bound &&
      tail(
        right_profile$delta_nll,
        1
      ) <= profile_threshold_95
    ) {
      upper_bound_status <-
        "not identified within admissible q_ITN parameter range"
    } else {
      upper_bound_status <-
        "not identified within successfully profiled upper range"
    }
  }
  
  both_profile_bounds_identified <-
    lower_bound_identified &&
    upper_bound_identified
  
  reduction_estimate_percent <-
    100 *
    (
      1 - q_hat
    )
  
  # Because reduction = 1 - q, the CI direction reverses.
  reduction_lower_95_percent <-
    if (upper_bound_identified) {
      100 *
        (
          1 - q_upper_95
        )
    } else {
      NA_real_
    }
  
  reduction_upper_95_percent <-
    if (lower_bound_identified) {
      100 *
        (
          1 - q_lower_95
        )
    } else {
      NA_real_
    }
  
  # 6.7.6 Profile interval and reduction-in-biting summary ----
  profile_ci <- data.frame(
    q_itn_estimate =
      q_hat,
    
    q_itn_lower_95_profile =
      q_lower_95,
    
    q_itn_upper_95_profile =
      q_upper_95,
    
    lower_bound_identified =
      lower_bound_identified,
    
    upper_bound_identified =
      upper_bound_identified,
    
    both_bounds_identified =
      both_profile_bounds_identified,
    
    lower_bound_status =
      lower_bound_status,
    
    upper_bound_status =
      upper_bound_status,
    
    reduction_in_b_percent_estimate =
      reduction_estimate_percent,
    
    reduction_in_b_percent_lower =
      reduction_lower_95_percent,
    
    reduction_in_b_percent_upper =
      reduction_upper_95_percent,
    
    primary_minimum_nll =
      primary_minimum_nll,
    
    delta_nll_threshold =
      profile_threshold_95,
    
    stringsAsFactors = FALSE
  )
  
  # 6.7.7 Profile result tables ----
  print_block(
    "q_ITN PROFILE RESULTS"
  )
  
  print(
    profile_itn
  )
  
  cat(
    "\nApproximate 95% profile-likelihood interval:\n"
  )
  
  print(
    profile_ci
  )
  
  write.csv(
    profile_itn,
    "profile_likelihood_q_itn.csv",
    row.names = FALSE
  )
  
  write.csv(
    profile_ci,
    "profile_likelihood_q_itn_CI.csv",
    row.names = FALSE
  )
  
  # 6.7.8 Profile-likelihood figure ----
  p_profile <- ggplot2::ggplot(
    profile_itn,
    ggplot2::aes(
      x = q_itn,
      y = delta_nll
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 1.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept =
        profile_threshold_95,
      linetype = 2
    ) +
    ggplot2::geom_vline(
      xintercept =
        q_hat,
      linetype = 3
    ) +
    ggplot2::labs(
      title =
        "Profile likelihood for the single ITN effect on effective biting",
      
      subtitle =
        "Approximate 95% interval defined by Delta NLL = 1.92 threshold crossings",
      
      x =
        "q_ITN: intervention multiplier on effective biting b",
      
      y =
        "Delta negative log-likelihood"
    ) +
    ggplot2::theme_minimal(
      base_size = 12
    )
  
  print(
    p_profile
  )
  
  ggplot2::ggsave(
    filename =
      "profile_likelihood_q_itn.png",
    
    plot =
      p_profile,
    
    width =
      8,
    
    height =
      5.5,
    
    dpi =
      300
  )
  
  cat(
  # 6.7.9 Profile interval interpretation and convergence checks ----
    "\n------------------------------------------------------------\n"
  )
  cat(
    "PROFILE-LIKELIHOOD INTERPRETATION\n"
  )
  cat(
    "------------------------------------------------------------\n"
  )
  
  cat(
    "Primary minimum NLL = ",
    signif(
      primary_minimum_nll,
      8
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Best q_ITN estimate = ",
    round(
      q_hat,
      4
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Estimated reduction in effective biting = ",
    round(
      reduction_estimate_percent,
      1
    ),
    "%\n",
    sep = ""
  )
  
  if (lower_bound_identified) {
    cat(
      "Lower 95% q_ITN profile limit = ",
      signif(
        q_lower_95,
        6
      ),
      "\n",
      sep = ""
    )
  } else {
    cat(
      "Lower 95% q_ITN profile limit: NOT IDENTIFIED (",
      lower_bound_status,
      ")\n",
      sep = ""
    )
  }
  
  if (upper_bound_identified) {
    cat(
      "Upper 95% q_ITN profile limit = ",
      signif(
        q_upper_95,
        6
      ),
      "\n",
      sep = ""
    )
  } else {
    cat(
      "Upper 95% q_ITN profile limit: NOT IDENTIFIED (",
      upper_bound_status,
      ")\n",
      sep = ""
    )
  }
  
  if (both_profile_bounds_identified) {
    cat(
      "Approximate 95% profile interval for reduction in effective biting = ",
      round(
        reduction_lower_95_percent,
        1
      ),
      "% to ",
      round(
        reduction_upper_95_percent,
        1
      ),
      "%\n",
      sep = ""
    )
  } else {
    cat(
      "A complete two-sided reduction interval is not reportable until both q_ITN profile bounds are identified.\n"
    )
  }
  
  if (!both_profile_bounds_identified) {
    warning(
      "At least one q_ITN 95% profile-likelihood bound was not identified. ",
      "Do not substitute a grid endpoint for the missing confidence limit."
    )
  }
  
  nonconverged_profile_points <- profile_itn %>%
    dplyr::filter(
      is.finite(nll),
      convergence != 0
    )
  
  if (nrow(nonconverged_profile_points) > 0) {
    warning(
      nrow(nonconverged_profile_points),
      " profile point(s) returned finite NLL values but non-zero nlminb convergence codes; ",
      "these points were excluded from CI crossing calculations."
    )
  }
  
  cat(
    "\nProfile likelihood completed.\n"
  )
}
# 6.8 Additional hypnozoite-parameter sensitivity analysis ----
# 6.8.1 Optional parameter-variation grid ----
if (
  RUN_HYPNOZOITE_SENSITIVITY
) {

  print_block(
    "HYPNOZOITE SENSITIVITY"
  )

  sensitivity_grid <-
    expand.grid(
      p_ret =
        c(
          0.25,
          0.75
        ),

      relapse_mean_days =
        c(
          60,
          90,
          120
        ),

      reservoir_mean_days =
        c(
          180,
          365
        ),

      KEEP.OUT.ATTRS =
        FALSE
    )

  # 6.8.2 Scenario-specific model refitting ----
  sensitivity_rows <-
    vector(
      "list",
      nrow(
        sensitivity_grid
      )
    )

  for (
    i in seq_len(
      nrow(
        sensitivity_grid
      )
    )
  ) {

    g <-
      sensitivity_grid[
        i,
      ]

    bio_tmp <-
      bio_parms

    bio_tmp[
      "p_ret"
    ] <-
      g$p_ret

    bio_tmp[
      "alpha_hyp"
    ] <-
      1 /
      g$relapse_mean_days

    bio_tmp[
      "delta_hyp"
    ] <-
      1 /
      g$reservoir_mean_days

    objective_tmp <-
      function(par) {

        objective_joint_fit(
          par,
          bio =
            bio_tmp
        )
      }

    fit_tmp <-
      tryCatch(
        nlminb(
          start =
            best_par,

          objective =
            objective_tmp,

          lower =
            lower_par,

          upper =
            upper_par,

          control = list(
            iter.max = 300,
            eval.max = 1000
          )
        ),

        error = function(e) {
          NULL
        }
      )

    if (
      is.null(
        fit_tmp
      )
    ) {

      sensitivity_rows[[i]] <-
        data.frame(
          p_ret =
            g$p_ret,

          relapse_mean_days =
            g$relapse_mean_days,

          reservoir_mean_days =
            g$reservoir_mean_days,

          convergence =
            999,

          nll =
            Inf,

          b_m =
            NA_real_,

          q_itn =
            NA_real_,

          phi_obs =
            NA_real_
        )

    } else {

      p_tmp <-
        unpack_fit_par(
          fit_tmp$par
        )

      sensitivity_rows[[i]] <-
        data.frame(
          p_ret =
            g$p_ret,

          relapse_mean_days =
            g$relapse_mean_days,

          reservoir_mean_days =
            g$reservoir_mean_days,

          convergence =
            fit_tmp$convergence,

          nll =
            fit_tmp$objective,

          b_m =
            p_tmp$b_m,

          q_itn =
            p_tmp$q_itn,

          phi_obs =
            p_tmp$phi_obs
        )
    }
  }

  # 6.8.3 Sensitivity result export ----
  sensitivity_results <-
    dplyr::bind_rows(
      sensitivity_rows
    )

  print(
    sensitivity_results
  )

  write.csv(
    sensitivity_results,
    "hypnozoite_sensitivity_results.csv",
    row.names = FALSE
  )
}
# 6.9 Model parameter reporting and prerequisites ----
# 6.9.1 Primary parameter summary table ----
parameter_table <-
  data.frame(
    parameter = c(
      "b_m",
      "q_itn",
      "phi_obs",
      "rho_obs",
      "mu_m",
      "lambda_import",
      "seasonal_amp",
      "seasonal_peak_day",
      "w_h",
      "p_ret",
      "alpha_hyp",
      "delta_hyp",
      "diagnostic_sensitivity",
      "diagnostic_specificity"
    ),

    value = c(
      best$b_m,
      best$q_itn,
      best$phi_obs,
      1 /
        (
          best$phi_obs +
          1
        ),
      bio_parms[
        "mu_m"
      ],
      lambda_import_fixed,
      seasonal_amp_fixed,
      seasonal_peak_day_fixed,
      bio_parms[
        "w_h"
      ],
      bio_parms[
        "p_ret"
      ],
      bio_parms[
        "alpha_hyp"
      ],
      bio_parms[
        "delta_hyp"
      ],
      diagnostic_sensitivity,
      diagnostic_specificity
    ),

    status = c(
      "Fitted jointly",
      "Fitted: ONLY intervention effect on effective biting",
      "Fitted shared observation-dispersion nuisance parameter",
      "Derived from phi_obs",
      "Fixed baseline mosquito mortality; same in both arms",
      "Fixed primary-analysis assumption",
      "Fixed primary-analysis assumption",
      "Fixed primary-analysis assumption",
      "Fixed T_h exit rate; mean residence is 1/w_h days",
      "Fixed proportion leaving T_h retaining hypnozoites; not measured treatment failure",
      "Fixed activation hazard among eligible carriers in S_h",
      "Fixed latent-carriage loss rate while in S_h",
      "Fixed observation assumption",
      "Fixed observation assumption"
    ),

    stringsAsFactors =
      FALSE
  )
print_block(
  "FINAL PARAMETER TABLE"
)
print(
  parameter_table
)
write.csv(
  parameter_table,
  "final_parameter_table.csv",
  row.names = FALSE
)
# 6.9.2 Additional plotting dependencies ----
for (pkg in c("ggplot2", "dplyr", "scales", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
# 6.9.3 Prerequisite object validation ----
required_objects_after_step25 <- c(
  "best",
  "best_result",
  "control_obs",
  "intervention_obs",
  "N_h_intervention",
  "init_control_seed",
  "burnin_days",
  "m_ratio_fixed",
  "bio_parms",
  "DETECTION_MODE",
  "diagnostic_sensitivity",
  "diagnostic_specificity",
  "seasonal_amp_fixed",
  "seasonal_peak_day_fixed",
  "study_start_doy",
  "make_model_parms",
  "run_model",
  "extract_end_state",
  "rescale_state"
)
missing_objects_after_step25 <- required_objects_after_step25[
  !vapply(
    required_objects_after_step25,
    exists,
    logical(1),
    inherits = TRUE
  )
]
if (length(missing_objects_after_step25) > 0) {
  stop(
    "Run the original script through STEP 25 first. Missing object(s): ",
    paste(missing_objects_after_step25, collapse = ", ")
  )
}
# 6.10 Transmission potential and long-run projections ----
# 6.10.1 Mechanistic reproduction-number formulation ----
# is second order. One infectious episode returns to another with probability
#   r = [c_h/(c_h + pi_h)] * p_ret * alpha_hyp/(alpha_hyp + delta_hyp).
# The new I_h -> T_h -> H_h -> E_h -> I_h pathway can repeat, so the
# expected number of infectious episodes is 1/(1-r), not a one-relapse sum.
# Positive rho_h and w_h affect waiting times, but not this episode count:
# E_h and T_h have no competing exit in the specified equations.


calculate_R0_current_model <- function(
    b_eff,
    m_ratio = m_ratio_fixed,
    bio = bio_parms
) {
  
  tau_mh <- as.numeric(bio[["tau_mh"]])
  tau_hm <- as.numeric(bio[["tau_hm"]])
  mu_m   <- as.numeric(bio[["mu_m"]])
  e_m    <- as.numeric(bio[["e_m"]])
  
  gamma_h <- as.numeric(
    bio[["c_h"]] +
      bio[["pi_h"]]
  )
  
  p_ret <- as.numeric(
    bio[["p_ret"]]
  )
  
  alpha_hyp <- as.numeric(
    bio[["alpha_hyp"]]
  )
  
  delta_hyp <- as.numeric(
    bio[["delta_hyp"]]
  )
  
  relapse_probability <-
    alpha_hyp /
    (
      alpha_hyp +
        delta_hyp
    )
  
  # Recovery must occur before another relapse cycle can be entered.
  return_probability <-
    as.numeric(bio[["c_h"]]) / gamma_h *
    p_ret * relapse_probability

  if (!is.finite(return_probability) || return_probability < 0 ||
      return_probability >= 1 || p_ret < 0 || p_ret > 1 ||
      !is.finite(bio[["rho_h"]]) || bio[["rho_h"]] <= 0 ||
      !is.finite(bio[["w_h"]]) || bio[["w_h"]] <= 0) {
    stop("R0 requires positive progression/return rates and a relapse-cycle probability in [0,1).")
  }

  expected_blood_stage_episodes <-
    1 / (1 - return_probability)
  
  # Human infectious host -> newly infected mosquitoes
  R_human_to_mosquito <-
    tau_hm *
    b_eff *
    m_ratio /
    gamma_h
  
  # Probability an exposed mosquito survives the incubation period
  mosquito_survival_to_infectious <-
    e_m /
    (
      e_m +
        mu_m
    )
  
  # Infectious mosquito -> new human primary infections
  R_mosquito_to_human <-
    tau_mh *
    b_eff /
    mu_m
  
  # Spectral-radius convention for vector-borne transmission:
  R0 <-
    sqrt(
      R_human_to_mosquito *
        mosquito_survival_to_infectious *
        R_mosquito_to_human *
        expected_blood_stage_episodes
    )
  
  as.numeric(R0)
}
# 6.10.2 Arm-specific reproduction-number estimates ----
R0_control <-
  calculate_R0_current_model(
    b_eff =
      best$b_m
  )
R0_intervention <-
  calculate_R0_current_model(
    b_eff =
      best$b_m *
      best$q_itn
  )
R0_summary <- data.frame(
  arm = c(
    "Control",
    "Intervention"
  ),
  
  b_effective = c(
    best$b_m,
    best$b_m *
      best$q_itn
  ),
  
  q_itn = c(
    1,
    best$q_itn
  ),
  
  R0 = c(
    R0_control,
    R0_intervention
  )
)
cat("\n============================================================\n")
cat("STEP 26: MODEL-SPECIFIC R0\n")
cat("============================================================\n")
print(
  R0_summary
)
cat(
  "\nControl R0 = ",
  round(R0_control, 3),
  "\n",
  sep = ""
)
cat(
  "Intervention R0 = ",
  round(R0_intervention, 3),
  "\n",
  sep = ""
)
write.csv(
  R0_summary,
  "R0_control_intervention.csv",
  row.names = FALSE
)
# 6.10.3 Intervention threshold for R0 below one ----
# Approximate ITN effect required to move R0 below 1.
# Because R0 is proportional to b in this formulation:
if (R0_control > 1) {
  
  q_R0_threshold <-
    1 /
    R0_control
  
  reduction_b_threshold_percent <-
    100 *
    (
      1 -
        q_R0_threshold
    )
  
  cat(
    "Approximate q_ITN required for R0 < 1: ",
    round(
      q_R0_threshold,
      3
    ),
    "\n",
    sep = ""
  )
  
  cat(
    "Approximate minimum reduction in effective biting required for R0 < 1: ",
    round(
      reduction_b_threshold_percent,
      1
    ),
    "%\n",
    sep = ""
  )
  
} else {
  
  q_R0_threshold <- NA_real_
  reduction_b_threshold_percent <- 0
  
  cat(
    "Control R0 is already below 1 under the fitted assumptions.\n"
  )
}
# Change this if you want a different projection period.
LONG_RUN_YEARS <- 5
long_run_days <-
  ceiling(
    LONG_RUN_YEARS *
      365.25
  )
# 6.10.4 Reconstruct fitted baseline states ----
longrun_burnin <-
  run_model(
    init_state =
      init_control_seed,
    
    parms =
      best_result$control_parms,
    
    days =
      burnin_days,
    
    arm_name =
      "LONG-RUN BURN-IN"
  )
if (is.null(longrun_burnin)) {
  stop(
    "Long-run burn-in failed."
  )
}
longrun_control_start <-
  extract_end_state(
    longrun_burnin
  )
longrun_intervention_start <-
  rescale_state(
    longrun_control_start,
    N_h_intervention
  )
# 6.10.5 Long-run simulations by trial arm ----
longrun_control_model <-
  run_model(
    init_state =
      longrun_control_start,
    
    parms =
      best_result$control_parms,
    
    days =
      long_run_days,
    
    arm_name =
      "CONTROL LONG RUN"
  )
longrun_intervention_model <-
  run_model(
    init_state =
      longrun_intervention_start,
    
    parms =
      best_result$intervention_parms,
    
    days =
      long_run_days,
    
    arm_name =
      "INTERVENTION LONG RUN"
  )
if (
  is.null(longrun_control_model) ||
  is.null(longrun_intervention_model)
) {
  stop(
    "Long-run simulation failed."
  )
}
# 6.10.6 Positivity, cumulative episodes and reproduction number ----
make_longrun_output <- function(
    model_output,
    arm_name,
    q_value
) {
  
  N_h_series <-
    model_output$S_h +
    model_output$E_h +
    model_output$I_h +
    model_output$T_h
  
  
  if (DETECTION_MODE == "I_only") {
    
    true_prevalence <-
      model_output$I_h /
      N_h_series
    
  } else if (
    DETECTION_MODE ==
    "E_plus_I"
  ) {
    
    true_prevalence <-
      (
        model_output$E_h +
          model_output$I_h
      ) /
      N_h_series
    
  } else {
    
    stop(
      "DETECTION_MODE must be 'I_only' or 'E_plus_I'."
    )
  }
  
  
  model_probability <-
    diagnostic_sensitivity *
    true_prevalence +
    (
      1 -
        diagnostic_specificity
    ) *
    (
      1 -
        true_prevalence
    )
  
  
  calendar_doy <-
    (
      (
        study_start_doy -
          1 +
          model_output$time
      ) %%
        365
    ) +
    1
  
  
  seasonal_multiplier <-
    1 +
    seasonal_amp_fixed *
    cos(
      2 *
        pi *
        (
          calendar_doy -
            seasonal_peak_day_fixed
        ) /
        365
    )
  
  
  seasonal_multiplier <-
    pmax(
      seasonal_multiplier,
      0.05
    )
  
  
  b_eff_series <-
    best$b_m *
    seasonal_multiplier *
    q_value
  
  
  R0_series <-
    vapply(
      b_eff_series,
      calculate_R0_current_model,
      numeric(1)
    )
  
  
  data.frame(
    arm =
      arm_name,
    
    day =
      model_output$time,
    
    year =
      model_output$time /
      365.25,
    
    true_prevalence =
      true_prevalence,
    
    model_probability =
      model_probability,
    
    b_effective =
      b_eff_series,
    
    R0 =
      R0_series,
    
    cumulative_blood_stage_episodes =
      model_output$inc_h
  )
}
longrun_control <-
  make_longrun_output(
    longrun_control_model,
    "Control",
    1
  )
longrun_intervention <-
  make_longrun_output(
    longrun_intervention_model,
    "Intervention",
    best$q_itn
  )
longrun_data <-
  dplyr::bind_rows(
    longrun_control,
    longrun_intervention
  )
write.csv(
  longrun_data,
  "long_run_model_R0.csv",
  row.names = FALSE
)
# 6.10.7 Observed and fitted surveillance data ----
observed_fitted_data <-
  dplyr::bind_rows(
    best_result$control_fit %>%
      dplyr::mutate(
        arm =
          "Control"
      ),
    
    best_result$intervention_fit %>%
      dplyr::mutate(
        arm =
          "Intervention"
      )
  )
active_end_year <-
  max(
    observed_fitted_data$followup_day,
    na.rm = TRUE
  ) /
  365.25
# 6.10.8 Observed-versus-fitted model panel ----
p_observed_fitted <-
  ggplot2::ggplot(
    observed_fitted_data,
    ggplot2::aes(
      x =
        survey_date
    )
  ) +
  
  ggplot2::geom_point(
    ggplot2::aes(
      y =
        observed_positivity,
      
      size =
        number_tested,
      
      colour =
        "Observed"
    ),
    alpha =
      0.80
  ) +
  
  ggplot2::geom_line(
    ggplot2::aes(
      y =
        model_probability,
      
      colour =
        "Model"
    ),
    linewidth =
      0.9
  ) +
  
  ggplot2::facet_wrap(
    ~ arm,
    scales =
      "free_x"
  ) +
  
  ggplot2::scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy =
          0.1
      )
  ) +
  
  ggplot2::labs(
    title =
      "Observed versus fitted active-surveillance positivity",
    
    subtitle =
      "Point size reflects number screened; passive October observations excluded",
    
    x =
      "Survey date",
    
    y =
      "Proportion positive",
    
    colour =
      "",
    
    size =
      "Number screened"
  ) +
  
  ggplot2::theme_minimal(
    base_size =
      12
  )
# 6.10.9 Five-year projection panel ----
p_longrun <-
  ggplot2::ggplot(
    longrun_data,
    ggplot2::aes(
      x =
        year,
      
      y =
        model_probability,
      
      colour =
        arm
    )
  ) +
  
  ggplot2::geom_line(
    linewidth =
      1
  ) +
  
  ggplot2::geom_vline(
    xintercept =
      active_end_year,
    
    linetype =
      3
  ) +
  
  ggplot2::scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy =
          0.1
      )
  ) +
  
  ggplot2::labs(
    title =
      paste0(
        "Long-run counterfactual projection (",
        LONG_RUN_YEARS,
        " years)"
      ),
    
    subtitle =
      "Model-based counterfactual under fixed assumptions; vertical dashed line marks end of active surveillance",
    
    x =
      "Years from study start",
    
    y =
      "Model-predicted positivity",
    
    colour =
      "Arm"
  ) +
  
  ggplot2::theme_minimal(
    base_size =
      12
  )
# 6.10.10 Reproduction-number projection panel ----
p_R0_longrun <-
  ggplot2::ggplot(
    longrun_data,
    ggplot2::aes(
      x =
        year,
      
      y =
        R0,
      
      colour =
        arm
    )
  ) +
  
  ggplot2::geom_hline(
    yintercept =
      1,
    
    linetype =
      2
  ) +
  
  ggplot2::geom_vline(
    xintercept =
      active_end_year,
    
    linetype =
      3
  ) +
  
  ggplot2::geom_line(
    linewidth =
      1
  ) +
  
  ggplot2::labs(
    title =
      "Model-specific R0",
    
    subtitle =
      "Horizontal dashed line marks the R0 = 1 threshold",
    
    x =
      "Years from study start",
    
    y =
      "R0",
    
    colour =
      "Arm"
  ) +
  
  ggplot2::theme_minimal(
    base_size =
      12
  )
# 6.10.11 Combined observed-fit and projection figure ----
combined_observed_longrun_R0 <-
  p_observed_fitted /
  (
    p_longrun |
      p_R0_longrun
  ) +
  patchwork::plot_annotation(
    title =
      "P. vivax ITN model: observed fit, counterfactual projection and model-specific R0"
  )
print(
  combined_observed_longrun_R0
)
ggplot2::ggsave(
  filename =
    "observed_fitted_longrun_R0_combined.png",
  
  plot =
    combined_observed_longrun_R0,
  
  width =
    14,
  
  height =
    10,
  
  dpi =
    300
)
# 6.11 Alternative ITN contact-reduction scenarios ----
# This is used only for the final integrity check.
primary_integrity_snapshot <-
  list(
    b_m =
      best$b_m,

    q_itn =
      best$q_itn,

    phi_obs =
      best$phi_obs,

    nll =
      best_result$objective,

    control_observed_positive =
      best_result$control_fit$number_positive,

    intervention_observed_positive =
      best_result$intervention_fit$number_positive,

    control_expected_positive =
      best_result$control_fit$expected_positive,

    intervention_expected_positive =
      best_result$intervention_fit$expected_positive
  )
fitted_effect_label <-
  paste0(
    "Fitted ITN effect (",
    round(
      100 *
        (
          1 -
            best$q_itn
        ),
      1
    ),
    "% reduction)"
  )
# 6.11.2 Effective-biting scenario definitions ----
scenario_definitions <-
  data.frame(
    scenario = c(
      "No ITN effect",
      "20% reduction in effective biting",
      fitted_effect_label,
      "50% reduction in effective biting",
      "70% reduction in effective biting"
    ),

    q_itn = c(
      1.00,
      0.80,
      best$q_itn,
      0.50,
      0.30
    ),

    stringsAsFactors =
      FALSE
  )
scenario_definitions$reduction_in_b_percent <-
  100 *
  (
    1 -
      scenario_definitions$q_itn
  )
scenario_definitions$scenario <-
  factor(
    scenario_definitions$scenario,
    levels =
      scenario_definitions$scenario
  )
# 6.11.3 Five-year scenario model simulations ----
run_one_scenario <- function(
    scenario_name,
    q_value
) {

  scenario_parms <-
    make_model_parms(
      b_m =
        best$b_m,

      q_itn =
        q_value,

      bio =
        bio_parms
    )

  scenario_model <-
    run_model(
      init_state =
        longrun_intervention_start,

      parms =
        scenario_parms,

      days =
        long_run_days,

      arm_name =
        paste0(
          "SCENARIO: ",
          scenario_name
        )
    )

  if (is.null(scenario_model)) {
    return(NULL)
  }

  scenario_output <-
    make_longrun_output(
      scenario_model,
      as.character(
        scenario_name
      ),
      q_value
    )

  scenario_output$scenario <-
    as.character(
      scenario_name
    )

  scenario_output$q_itn <-
    q_value

  scenario_output$reduction_in_b_percent <-
    100 *
    (
      1 -
        q_value
    )

  scenario_output
}
scenario_longrun <-
  dplyr::bind_rows(
    lapply(
      seq_len(
        nrow(
          scenario_definitions
        )
      ),
      function(i) {

        run_one_scenario(
          scenario_name =
            as.character(
              scenario_definitions$scenario[i]
            ),

          q_value =
            scenario_definitions$q_itn[i]
        )
      }
    )
  )
if (nrow(scenario_longrun) == 0) {
  stop(
    "All effective-biting scenario simulations failed."
  )
}
scenario_longrun$scenario <-
  factor(
    scenario_longrun$scenario,
    levels =
      levels(
        scenario_definitions$scenario
      )
  )
write.csv(
  scenario_longrun,
  "scenario_longrun_trajectories.csv",
  row.names = FALSE
)
# 6.11.4 Extract scenario outcomes at selected years ----
value_at_year <- function(
    data,
    variable,
    selected_year
) {

  stats::approx(
    x =
      data$year,

    y =
      data[[variable]],

    xout =
      selected_year,

    rule =
      2
  )$y
}
# 6.11.5 Alternative intervention-strength scenario summaries ----
scenario_summary <-
  dplyr::bind_rows(
    lapply(
      seq_len(
        nrow(
          scenario_definitions
        )
      ),
      function(i) {

        scenario_name <-
          as.character(
            scenario_definitions$scenario[i]
          )

        q_value <-
          scenario_definitions$q_itn[i]

        scenario_data <-
          scenario_longrun %>%
          dplyr::filter(
            scenario ==
              scenario_name
          )

        data.frame(
          scenario =
            scenario_name,

          q_itn =
            q_value,

          reduction_in_effective_biting_percent =
            100 *
            (
              1 -
                q_value
            ),

          R0 =
            calculate_R0_current_model(
              best$b_m *
                q_value
            ),

          predicted_positivity_year1 =
            value_at_year(
              scenario_data,
              "model_probability",
              1
            ),

          predicted_positivity_year3 =
            value_at_year(
              scenario_data,
              "model_probability",
              3
            ),

          predicted_positivity_year5 =
            value_at_year(
              scenario_data,
              "model_probability",
              LONG_RUN_YEARS
            ),

          cumulative_blood_stage_episodes_per_1000_year5 =
            value_at_year(
              scenario_data,
              "cumulative_blood_stage_episodes",
              LONG_RUN_YEARS
            ) /
            N_h_intervention *
            1000
        )
      }
    )
  )
scenario_summary$scenario <-
  factor(
    scenario_summary$scenario,
    levels =
      levels(
        scenario_definitions$scenario
      )
  )
write.csv(
  scenario_summary,
  "scenario_analysis_summary.csv",
  row.names = FALSE
)
# 6.11.6 Effective-biting scenario plots ----
p_scenario_longrun <-
  ggplot2::ggplot(
    scenario_longrun,
    ggplot2::aes(
      x = year,
      y = model_probability * 1000,
      group = scenario,
      colour = scenario,
      linetype = scenario
    )
  ) +
  ggplot2::geom_line(
    linewidth = 1
  ) +
  ggplot2::labs(
    x = "Years after intervention",
    y = "Predicted malaria positivity per 1,000",
    colour = "Intervention-strength scenario",
    linetype = "Intervention-strength scenario",
    title = "Predicted long-run malaria positivity"
  ) +
  ggplot2::guides(
    colour =
      ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      ),
    linetype =
      ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title =
      ggplot2::element_text(
        face = "bold"
      ),
    plot.margin =
      ggplot2::margin(
        8,
        12,
        8,
        8
      )
  )
p_scenario_R0 <-
  ggplot2::ggplot(
    scenario_summary,
    ggplot2::aes(
      x = scenario,
      y = R0
    )
  ) +
  ggplot2::geom_col(
    width = 0.7
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::labs(
    x = NULL,
    y = expression(R[0]),
    title = "Model-specific reproduction number"
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_classic() +
  ggplot2::theme(
    plot.margin =
      ggplot2::margin(
        8,
        18,
        8,
        8
      )
  )
if (!requireNamespace(
  "patchwork",
  quietly = TRUE
)) {
  install.packages("patchwork")
}
combined_scenario_plot <-
  patchwork::wrap_plots(
    p_scenario_longrun,
    p_scenario_R0,
    ncol = 2,
    widths = c(1.55, 1.15)
  ) +
  patchwork::plot_annotation(
    title =
      "Long-run P. vivax transmission under alternative ITN-effect scenarios",

    subtitle =
      paste0(
        "Fitted ITN effect: ",
        round(
          100 *
            (
              1 -
                best$q_itn
            ),
          1
        ),
        "% reduction in effective biting"
      )
  )
print(
  combined_scenario_plot
)
ggplot2::ggsave(
  filename =
    "scenario_analysis_longrun_R0.png",

  plot =
    combined_scenario_plot,

  width =
    16.5,

  height =
    7.5,

  dpi =
    300,

  bg =
    "white"
)
# 6.12 Coverage-dependent counterfactual scenarios ----
# 6.12.1 Coverage–contact multiplier calibration ----
OBSERVED_ITN_COVERAGE <- 0.35
coverage_gamma <-
  -log(
    best$q_itn
  ) /
  OBSERVED_ITN_COVERAGE
q_from_coverage <- function(
    coverage
) {

  exp(
    -coverage_gamma *
      coverage
  )
}
coverage_q0_ok <-
  isTRUE(
    all.equal(
      q_from_coverage(0),
      1,
      tolerance = 1e-12
    )
  )
coverage_q35_ok <-
  isTRUE(
    all.equal(
      q_from_coverage(
        OBSERVED_ITN_COVERAGE
      ),
      best$q_itn,
      tolerance = 1e-12
    )
  )
if (!coverage_q0_ok) {
  stop(
    "Coverage mapping failed: 0% coverage does not give q_ITN = 1."
  )
}
if (!coverage_q35_ok) {
  stop(
    "Coverage mapping failed: 35% coverage does not reproduce fitted q_ITN."
  )
}
# 6.12.2 ITN coverage levels and parameter mapping ----
coverage_levels <- c(
  0.00,
  0.20,
  0.35,
  0.40,
  0.60,
  0.80,
  1.00
)
coverage_scenario_definitions <-
  data.frame(
    coverage =
      coverage_levels,

    coverage_percent =
      coverage_levels *
      100,

    q_itn =
      q_from_coverage(
        coverage_levels
      ),

    stringsAsFactors =
      FALSE
  )
coverage_scenario_definitions$reduction_in_effective_biting_percent <-
  100 *
  (
    1 -
      coverage_scenario_definitions$q_itn
  )
coverage_scenario_definitions$scenario <-
  ifelse(
    abs(
      coverage_scenario_definitions$coverage -
        OBSERVED_ITN_COVERAGE
    ) <
      1e-12,

    paste0(
      round(
        100 *
          OBSERVED_ITN_COVERAGE
      ),
      "% ITN coverage (observed/fitted)"
    ),

    paste0(
      round(
        coverage_scenario_definitions$coverage_percent
      ),
      "% ITN coverage"
    )
  )
coverage_scenario_definitions$scenario <-
  factor(
    coverage_scenario_definitions$scenario,
    levels =
      coverage_scenario_definitions$scenario
  )
print_block(
  "ITN COVERAGE SCENARIO DEFINITIONS"
)
print(
  coverage_scenario_definitions
)
write.csv(
  coverage_scenario_definitions,
  "ITN_coverage_scenario_definitions.csv",
  row.names = FALSE
)
# 6.12.3 Coverage-specific model trajectories ----
coverage_scenario_longrun <-
  dplyr::bind_rows(
    lapply(
      seq_len(
        nrow(
          coverage_scenario_definitions
        )
      ),
      function(i) {

        run_one_scenario(
          scenario_name =
            as.character(
              coverage_scenario_definitions$scenario[i]
            ),

          q_value =
            coverage_scenario_definitions$q_itn[i]
        )
      }
    )
  )
if (nrow(coverage_scenario_longrun) == 0) {
  stop(
    "All ITN coverage scenario simulations failed."
  )
}
coverage_scenario_longrun <-
  coverage_scenario_longrun %>%
  dplyr::left_join(
    coverage_scenario_definitions %>%
      dplyr::select(
        scenario,
        coverage,
        coverage_percent
      ) %>%
      dplyr::mutate(
        scenario =
          as.character(
            scenario
          )
      ),
    by = "scenario"
  )
coverage_scenario_longrun$scenario <-
  factor(
    coverage_scenario_longrun$scenario,
    levels =
      levels(
        coverage_scenario_definitions$scenario
      )
  )
write.csv(
  coverage_scenario_longrun,
  "ITN_coverage_scenario_trajectories.csv",
  row.names = FALSE
)
# 6.12.4 Coverage scenario outcome tables ----
coverage_scenario_summary <-
  dplyr::bind_rows(
    lapply(
      seq_len(
        nrow(
          coverage_scenario_definitions
        )
      ),
      function(i) {

        scenario_name <-
          as.character(
            coverage_scenario_definitions$scenario[i]
          )

        q_value <-
          coverage_scenario_definitions$q_itn[i]

        coverage_value <-
          coverage_scenario_definitions$coverage[i]

        scenario_data <-
          coverage_scenario_longrun %>%
          dplyr::filter(
            scenario ==
              scenario_name
          )

        data.frame(
          scenario =
            scenario_name,

          ITN_coverage_percent =
            100 *
            coverage_value,

          q_itn =
            q_value,

          reduction_in_effective_biting_percent =
            100 *
            (
              1 -
                q_value
            ),

          R0 =
            calculate_R0_current_model(
              best$b_m *
                q_value
            ),

          predicted_positivity_year1 =
            value_at_year(
              scenario_data,
              "model_probability",
              1
            ),

          predicted_positivity_year3 =
            value_at_year(
              scenario_data,
              "model_probability",
              3
            ),

          predicted_positivity_year5 =
            value_at_year(
              scenario_data,
              "model_probability",
              LONG_RUN_YEARS
            ),

          cumulative_blood_stage_episodes_per_1000_year5 =
            value_at_year(
              scenario_data,
              "cumulative_blood_stage_episodes",
              LONG_RUN_YEARS
            ) /
            N_h_intervention *
            1000
        )
      }
    )
  )
coverage_scenario_summary$scenario <-
  factor(
    coverage_scenario_summary$scenario,
    levels =
      levels(
        coverage_scenario_definitions$scenario
      )
  )
print_block(
  "ITN COVERAGE SCENARIO SUMMARY"
)
print(
  coverage_scenario_summary
)
write.csv(
  coverage_scenario_summary,
  "ITN_coverage_scenario_summary.csv",
  row.names = FALSE
)
# 6.12.5 Coverage-dependent projection figures ----
p_coverage_longrun <-
  ggplot2::ggplot(
    coverage_scenario_longrun,
    ggplot2::aes(
      x = year,
      y = model_probability * 1000,
      group = scenario,
      colour = scenario,
      linetype = scenario
    )
  ) +
  ggplot2::geom_line(
    linewidth = 1
  ) +
  ggplot2::geom_vline(
    xintercept =
      active_end_year,
    linetype =
      "dashed",
    linewidth =
      0.7
  ) +
  ggplot2::labs(
    x = "Years after intervention",
    y = "Predicted malaria positivity per 1,000",
    colour = "ITN coverage",
    linetype = "ITN coverage",
    title = "Five-year counterfactual malaria positivity",
    caption = "Dashed vertical line marks the end of observed 2001 active-surveillance follow-up."
  ) +
  ggplot2::guides(
    colour =
      ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      ),
    linetype =
      ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE
      )
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position =
      "bottom",
    plot.margin =
      ggplot2::margin(
        8,
        12,
        8,
        8
      )
  )
p_coverage_R0 <-
  ggplot2::ggplot(
    coverage_scenario_summary,
    ggplot2::aes(
      x = scenario,
      y = R0
    )
  ) +
  ggplot2::geom_col(
    width = 0.7
  ) +
  ggplot2::geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::labs(
    x = NULL,
    y = expression(R[0]),
    title = "Model-specific reproduction number"
  ) +
  ggplot2::coord_flip() +
  ggplot2::theme_classic() +
  ggplot2::theme(
    plot.margin =
      ggplot2::margin(
        8,
        18,
        8,
        8
      )
  )
# 6.12.6 Composite ITN coverage scenario figure ----
combined_coverage_scenario_plot <-
  patchwork::wrap_plots(
    p_coverage_longrun,
    p_coverage_R0,
    ncol = 2,
    widths = c(1.55, 1.15)
  ) +
  patchwork::plot_annotation(
    title =
      "Five-year model-based P. vivax projections under alternative ITN coverage scenarios",

    subtitle =
      "Coverage-effect relationship anchored to the fitted intervention effect at 35% observed ITN coverage"
  )
print(
  combined_coverage_scenario_plot
)
ggplot2::ggsave(
  filename =
    "ITN_coverage_scenario_longrun_R0.png",

  plot =
    combined_coverage_scenario_plot,

  width =
    16.5,

  height =
    7.5,

  dpi =
    300,

  bg =
    "white"
)
# 6.13 One-way sensitivity analysis and parameter refitting ----
# 6.13.1 Multistart refitting of sensitivity scenarios ----
# Each one-way sensitivity scenario uses the same 40-start nlminb
# strategy and the same parameter bounds as the primary fit. Only
# b_m, q_ITN and phi_obs are fitted; the sensitivity parameter is
# held fixed through the appropriate biological/observation override.
refit_secondary_model <- function(
    control_obs_override = NULL,
    intervention_obs_override = NULL,
    bio_override = bio_parms,
    seasonal_amp_override = seasonal_amp_fixed,
    seasonal_peak_day_override = seasonal_peak_day_fixed,
    lambda_import_override = lambda_import_fixed,
    diagnostic_sensitivity_override = diagnostic_sensitivity,
    diagnostic_specificity_override = diagnostic_specificity,
    n_starts = n_multistart
) {
  
  original_control_obs <-
    control_obs
  
  original_intervention_obs <-
    intervention_obs
  
  original_seasonal_amp <-
    seasonal_amp_fixed
  
  original_seasonal_peak_day <-
    seasonal_peak_day_fixed
  
  original_lambda_import <-
    lambda_import_fixed
  
  original_diagnostic_sensitivity <-
    diagnostic_sensitivity
  
  original_diagnostic_specificity <-
    diagnostic_specificity
  
  on.exit(
    {
      control_obs <<-
        original_control_obs
      
      intervention_obs <<-
        original_intervention_obs
      
      seasonal_amp_fixed <<-
        original_seasonal_amp
      
      seasonal_peak_day_fixed <<-
        original_seasonal_peak_day
      
      lambda_import_fixed <<-
        original_lambda_import
      
      diagnostic_sensitivity <<-
        original_diagnostic_sensitivity
      
      diagnostic_specificity <<-
        original_diagnostic_specificity
    },
    add = TRUE
  )
  
  if (!is.null(control_obs_override)) {
    control_obs <<-
      control_obs_override
  }
  
  if (!is.null(intervention_obs_override)) {
    intervention_obs <<-
      intervention_obs_override
  }
  
  seasonal_amp_fixed <<-
    seasonal_amp_override
  
  seasonal_peak_day_fixed <<-
    seasonal_peak_day_override
  
  lambda_import_fixed <<-
    lambda_import_override
  
  diagnostic_sensitivity <<-
    diagnostic_sensitivity_override
  
  diagnostic_specificity <<-
    diagnostic_specificity_override
  
  objective_secondary <- function(par) {
    objective_joint_fit(
      par,
      bio =
        bio_override
    )
  }
  
  # Use the same starting-value strategy as the primary multistart fit:
  # one prespecified start followed by random starts generated by
  # make_random_start() within the same fitted-parameter bounds.
  secondary_starts <- vector(
    "list",
    n_starts
  )
  
  secondary_starts[[1]] <- c(
    log_b_m =
      log(1.0),
    
    eta_q_itn =
      qlogis(0.60 / 2),
    
    log_phi_obs =
      log(50)
  )
  
  if (n_starts > 1) {
    for (i in 2:n_starts) {
      secondary_starts[[i]] <-
        make_random_start()
    }
  }
  
  fit_secondary_one_start <- function(
      start,
      start_id
  ) {
    cat("  Sensitivity start ", start_id, " of ", n_starts, "\n", sep = "")
    flush.console()
    
    fit_tmp <- tryCatch(
      nlminb(
        start =
          start,
        
        objective =
          objective_secondary,
        
        lower =
          lower_par,
        
        upper =
          upper_par,
        
        control = list(
          iter.max = 500,
          eval.max = 2000,
          rel.tol = 1e-9,
          x.tol = 1e-9
        )
      ),
      
      error = function(e) {
        NULL
      }
    )
    
    if (
      is.null(fit_tmp) ||
      !is.finite(fit_tmp$objective)
    ) {
      return(
        data.frame(
          start_id =
            start_id,
          
          objective =
            Inf,
          
          convergence =
            999,
          
          message =
            "fit failed",
          
          log_b_m =
            NA_real_,
          
          eta_q_itn =
            NA_real_,
          
          log_phi_obs =
            NA_real_,
          
          b_m =
            NA_real_,
          
          q_itn =
            NA_real_,
          
          phi_obs =
            NA_real_,
          
          stringsAsFactors = FALSE
        )
      )
    }
    
    p_tmp <- unpack_fit_par(
      fit_tmp$par
    )
    
    data.frame(
      start_id =
        start_id,
      
      objective =
        fit_tmp$objective,
      
      convergence =
        fit_tmp$convergence,
      
      message =
        if (is.null(fit_tmp$message)) "" else fit_tmp$message,
      
      log_b_m =
        as.numeric(
          fit_tmp$par[["log_b_m"]]
        ),
      
      eta_q_itn =
        as.numeric(
          fit_tmp$par[["eta_q_itn"]]
        ),
      
      log_phi_obs =
        as.numeric(
          fit_tmp$par[["log_phi_obs"]]
        ),
      
      b_m =
        p_tmp$b_m,
      
      q_itn =
        p_tmp$q_itn,
      
      phi_obs =
        p_tmp$phi_obs,
      
      stringsAsFactors = FALSE
    )
  }
  
  secondary_multistart_results <- dplyr::bind_rows(
    lapply(
      seq_len(n_starts),
      function(i) {
        fit_secondary_one_start(
          secondary_starts[[i]],
          i
        )
      }
    )
  ) %>%
    dplyr::arrange(
      objective
    )
  
  secondary_successful <- secondary_multistart_results %>%
    dplyr::filter(
      convergence == 0,
      is.finite(objective)
    )
  
  convergence_count <-
    nrow(secondary_successful)
  
  if (convergence_count == 0) {
    return(
      list(
        fit = NULL,
        fitted = NULL,
        multistart_results =
          secondary_multistart_results,
        convergence_count =
          0L,
        near_best_count =
          0L,
        status =
          "NO VALID CONVERGED SOLUTION"
      )
    )
  }
  
  secondary_best_row <-
    secondary_successful[1, ]
  
  secondary_best_par <- c(
    log_b_m =
      secondary_best_row$log_b_m,
    
    eta_q_itn =
      secondary_best_row$eta_q_itn,
    
    log_phi_obs =
      secondary_best_row$log_phi_obs
  )
  
  # Retain the best valid converged solution across all starts.
  # Do not replace it with a single additional optimiser run.
  fit_secondary_best <- list(
    par =
      secondary_best_par,
    objective =
      secondary_best_row$objective,
    convergence =
      secondary_best_row$convergence,
    message =
      secondary_best_row$message
  )
  
  fitted_secondary <- unpack_fit_par(
    fit_secondary_best$par
  )
  
  best_secondary_nll <-
    fit_secondary_best$objective
  
  near_best_count <- sum(
    secondary_successful$objective <=
      best_secondary_nll + 2
  )
  
  list(
    fit =
      fit_secondary_best,
    
    fitted =
      fitted_secondary,
    
    multistart_results =
      secondary_multistart_results,
    
    convergence_count =
      convergence_count,
    
    near_best_count =
      near_best_count,
    
    status =
      "MULTISTART COMPLETED"
  )
}
# 6.13.2 One-way sensitivity scenario specifications ----
sensitivity_status_text <- "NOT REQUESTED"
primary_model_sensitivity_summary <- NULL
if (RUN_LIMITED_SENSITIVITY) {
  
  print_block(
    "LIMITED ONE-WAY SENSITIVITY ANALYSIS"
  )
  
  sensitivity_specs <- list()
  
  add_sensitivity_spec <- function(
      name,
      assumption,
      value,
      bio_value = bio_parms,
      seasonal_amp_value = seasonal_amp_fixed,
      seasonal_peak_day_value = seasonal_peak_day_fixed,
      lambda_import_value = lambda_import_fixed,
      diagnostic_sensitivity_value = diagnostic_sensitivity,
      diagnostic_specificity_value = diagnostic_specificity
  ) {
    
    sensitivity_specs[[
      length(sensitivity_specs) + 1
    ]] <<-
      list(
        name =
          name,
        
        assumption =
          assumption,
        
        value =
          as.character(value),
        
        bio =
          bio_value,
        
        seasonal_amp =
          seasonal_amp_value,
        
        seasonal_peak_day =
          seasonal_peak_day_value,
        
        lambda_import =
          lambda_import_value,
        
        diagnostic_sensitivity =
          diagnostic_sensitivity_value,
        
        diagnostic_specificity =
          diagnostic_specificity_value
      )
  }
  
  # p_ret sensitivity.
  for (p_ret_value in SENSITIVITY_P_RET_VALUES) {
    
    bio_tmp <-
      bio_parms
    
    bio_tmp["p_ret"] <-
      p_ret_value
    
    add_sensitivity_spec(
      name =
        paste0(
          "p_ret = ",
          p_ret_value
        ),
      
      assumption =
        "Proportion leaving T_h retaining hypnozoites",
      
      value =
        p_ret_value,
      
      bio_value =
        bio_tmp
    )
  }
  
  # T_h residence-time sensitivity.
  for (t_residence_days in SENSITIVITY_T_RESIDENCE_DAYS) {
    
    bio_tmp <-
      bio_parms
    
    bio_tmp["w_h"] <-
      1 / t_residence_days
    
    add_sensitivity_spec(
      name =
        paste0(
          "T_h mean residence = ",
          t_residence_days,
          " days"
        ),
      
      assumption =
        "Mean treated/recovered residence time (1/w_h)",
      
      value =
        paste0(
          t_residence_days,
          " days"
        ),
      
      bio_value =
        bio_tmp
    )
  }
  
  # Relapse/eligible-carrier activation sensitivity.
  for (relapse_mean_days in SENSITIVITY_RELAPSE_MEAN_DAYS) {
    
    bio_tmp <-
      bio_parms
    
    bio_tmp["alpha_hyp"] <-
      1 / relapse_mean_days
    
    add_sensitivity_spec(
      name =
        paste0(
          "Eligible-carrier activation timescale = ",
          relapse_mean_days,
          " days"
        ),
      
      assumption =
        "Inverse eligible-carrier activation hazard (1/alpha_hyp); not time since treatment",
      
      value =
        paste0(
          relapse_mean_days,
          " days"
        ),
      
      bio_value =
        bio_tmp
    )
  }
  
  # Latent-carriage duration sensitivity.
  for (reservoir_mean_days in SENSITIVITY_RESERVOIR_MEAN_DAYS) {
    
    bio_tmp <-
      bio_parms
    
    bio_tmp["delta_hyp"] <-
      1 / reservoir_mean_days
    
    add_sensitivity_spec(
      name =
        paste0(
          "Latent-carriage loss timescale = ",
          reservoir_mean_days,
          " days"
        ),
      
      assumption =
        "Inverse latent-carriage loss rate while in S_h (1/delta_hyp)",
      
      value =
        paste0(
          reservoir_mean_days,
          " days"
        ),
      
      bio_value =
        bio_tmp
    )
  }
  
  # Seasonality sensitivity.
  for (seasonal_amp_value in SENSITIVITY_SEASONAL_AMPLITUDES) {
    
    add_sensitivity_spec(
      name =
        paste0(
          "Seasonal amplitude = ",
          seasonal_amp_value
        ),
      
      assumption =
        "Seasonal amplitude",
      
      value =
        seasonal_amp_value,
      
      seasonal_amp_value =
        seasonal_amp_value
    )
  }
  
  # Importation sensitivity.
  for (import_value in SENSITIVITY_IMPORT_RATES) {
    
    add_sensitivity_spec(
      name =
        paste0(
          "Importation force = ",
          format(
            import_value,
            scientific = TRUE
          )
        ),
      
      assumption =
        "Imported/background force of infection per susceptible per day",
      
      value =
        format(
          import_value,
          scientific = TRUE
        ),
      
      lambda_import_value =
        import_value
    )
  }
  
  # Retain the existing optional diagnostic sensitivity hooks unchanged.
  if (length(DIAGNOSTIC_SENSITIVITY_VALUES) > 0) {
    for (
      diagnostic_sensitivity_value in
        DIAGNOSTIC_SENSITIVITY_VALUES
    ) {
      add_sensitivity_spec(
        name =
          paste0(
            "Diagnostic sensitivity = ",
            diagnostic_sensitivity_value
          ),
        
        assumption =
          "Diagnostic sensitivity",
        
        value =
          diagnostic_sensitivity_value,
        
        diagnostic_sensitivity_value =
          diagnostic_sensitivity_value
      )
    }
  }
  
  if (length(DIAGNOSTIC_SPECIFICITY_VALUES) > 0) {
    for (
      diagnostic_specificity_value in
        DIAGNOSTIC_SPECIFICITY_VALUES
    ) {
      add_sensitivity_spec(
        name =
          paste0(
            "Diagnostic specificity = ",
            diagnostic_specificity_value
          ),
        
        assumption =
          "Diagnostic specificity",
        
        value =
          diagnostic_specificity_value,
        
        diagnostic_specificity_value =
          diagnostic_specificity_value
      )
    }
  }
  
  # 6.13.3 Primary-model reference for sensitivity comparisons ----
  sensitivity_rows <- list(
    data.frame(
      scenario =
        "Primary fitted model",
      
      assumption_varied =
        "None",
      
      parameter_varied =
        "None",
      
      value =
        "Primary assumptions",
      
      value_tested =
        "Primary assumptions",
      
      convergence =
        0,
      
      convergence_count =
        nrow(successful_fits),
      
      nll =
        best_result$objective,
      
      best_nll =
        best_result$objective,
      
      delta_nll =
        0,
      
      b_m =
        best$b_m,
      
      q_itn =
        best$q_itn,
      
      reduction_in_effective_biting_percent =
        100 *
        (
          1 - best$q_itn
        ),
      
      phi_obs =
        best$phi_obs,
      
      R0_control =
        R0_control,
      
      R0_intervention =
        R0_intervention,
      
      multistart_status =
        "PRIMARY FIT",
      
      stringsAsFactors = FALSE
    )
  )
  
  # 6.13.4 Refit and evaluate all sensitivity scenarios ----
  if (length(sensitivity_specs) > 0) {
    
    for (i in seq_along(sensitivity_specs)) {
      
      spec <-
        sensitivity_specs[[i]]
      
      cat(
        "\nSensitivity scenario ",
        i,
        " of ",
        length(sensitivity_specs),
        ": ",
        spec$name,
        "\n",
        sep = ""
      )
      
      fit_sensitivity <-
        refit_secondary_model(
          bio_override =
            spec$bio,
          
          seasonal_amp_override =
            spec$seasonal_amp,
          
          seasonal_peak_day_override =
            spec$seasonal_peak_day,
          
          lambda_import_override =
            spec$lambda_import,
          
          diagnostic_sensitivity_override =
            spec$diagnostic_sensitivity,
          
          diagnostic_specificity_override =
            spec$diagnostic_specificity,
          
          n_starts =
            n_multistart
        )
      
      scenario_multistart_file <- paste0(
        "sensitivity_multistart_",
        sprintf("%02d", i),
        ".csv"
      )
      
      write.csv(
        fit_sensitivity$multistart_results,
        scenario_multistart_file,
        row.names = FALSE
      )
      
      if (
        is.null(fit_sensitivity$fit) ||
        is.null(fit_sensitivity$fitted)
      ) {
        
        sensitivity_rows[[
          length(sensitivity_rows) + 1
        ]] <-
          data.frame(
            scenario =
              spec$name,
            
            assumption_varied =
              spec$assumption,
            
            parameter_varied =
              spec$assumption,
            
            value =
              spec$value,
            
            value_tested =
              spec$value,
            
            convergence =
              999,
            
            convergence_count =
              fit_sensitivity$convergence_count,
            
            nll =
              NA_real_,
            
            best_nll =
              NA_real_,
            
            delta_nll =
              NA_real_,
            
            b_m =
              NA_real_,
            
            q_itn =
              NA_real_,
            
            reduction_in_effective_biting_percent =
              NA_real_,
            
            phi_obs =
              NA_real_,
            
            R0_control =
              NA_real_,
            
            R0_intervention =
              NA_real_,
            
            multistart_status =
              fit_sensitivity$status,
            
            stringsAsFactors = FALSE
          )
        
      } else {
        
        p_sensitivity <-
          fit_sensitivity$fitted
        
        scenario_best_nll <-
          fit_sensitivity$fit$objective
        
        sensitivity_rows[[
          length(sensitivity_rows) + 1
        ]] <-
          data.frame(
            scenario =
              spec$name,
            
            assumption_varied =
              spec$assumption,
            
            parameter_varied =
              spec$assumption,
            
            value =
              spec$value,
            
            value_tested =
              spec$value,
            
            convergence =
              fit_sensitivity$fit$convergence,
            
            convergence_count =
              fit_sensitivity$convergence_count,
            
            nll =
              scenario_best_nll,
            
            best_nll =
              scenario_best_nll,
            
            delta_nll =
              scenario_best_nll -
              best_result$objective,
            
            b_m =
              p_sensitivity$b_m,
            
            q_itn =
              p_sensitivity$q_itn,
            
            reduction_in_effective_biting_percent =
              100 *
              (
                1 - p_sensitivity$q_itn
              ),
            
            phi_obs =
              p_sensitivity$phi_obs,
            
            R0_control =
              calculate_R0_current_model(
                p_sensitivity$b_m,
                bio =
                  spec$bio
              ),
            
            R0_intervention =
              calculate_R0_current_model(
                p_sensitivity$b_m *
                p_sensitivity$q_itn,
                bio =
                  spec$bio
              ),
            
            multistart_status =
              fit_sensitivity$status,
            
            stringsAsFactors = FALSE
          )
      }
    }
  }
  
  # 6.13.5 Sensitivity results, convergence and export ----
  primary_model_sensitivity_summary <- dplyr::bind_rows(
    sensitivity_rows
  )
  
  print(
    primary_model_sensitivity_summary
  )
  
  write.csv(
    primary_model_sensitivity_summary,
    "primary_model_sensitivity_summary.csv",
    row.names = FALSE
  )
  
  sensitivity_successes <- primary_model_sensitivity_summary %>%
    dplyr::filter(
      scenario != "Primary fitted model",
      convergence == 0,
      convergence_count > 0,
      is.finite(q_itn),
      is.finite(best_nll)
    )
  
  if (
    nrow(sensitivity_successes) ==
    (
      nrow(primary_model_sensitivity_summary) - 1
    )
  ) {
    sensitivity_status_text <-
      "COMPLETED"
    
  } else if (nrow(sensitivity_successes) > 0) {
    sensitivity_status_text <-
      "PARTIAL: one or more multistart sensitivity refits had no valid converged solution"
    
  } else {
    sensitivity_status_text <-
      "FAILED: no alternative multistart sensitivity refits converged"
  }
}
# 6.14 Final model validation and thesis reporting ----
# 6.14.1 Integrity checks for primary fitted outputs ----
print_block(
  "FINAL THESIS-READINESS CHECKLIST"
)
# The primary fitted objects must remain unchanged after all secondary analyses.
primary_fit_integrity_ok <-
  isTRUE(
    all.equal(
      best$b_m,
      best_result$fitted$b_m,
      tolerance = 1e-12
    )
  ) &&
  isTRUE(
    all.equal(
      best$q_itn,
      best_result$fitted$q_itn,
      tolerance = 1e-12
    )
  ) &&
  isTRUE(
    all.equal(
      best$phi_obs,
      best_result$fitted$phi_obs,
      tolerance = 1e-12
    )
  )
observed_positive_integrity_ok <-
  identical(
    as.numeric(
      best_result$control_fit$number_positive
    ),
    as.numeric(
      primary_integrity_snapshot$control_observed_positive
    )
  ) &&
  identical(
    as.numeric(
      best_result$intervention_fit$number_positive
    ),
    as.numeric(
      primary_integrity_snapshot$intervention_observed_positive
    )
  ) &&
  sum(
    best_result$control_fit$number_positive
  ) ==
    sum(
      control_obs$number_positive[
        control_obs$include_in_fit
      ]
    ) &&
  sum(
    best_result$intervention_fit$number_positive
  ) ==
    sum(
      intervention_obs$number_positive[
        intervention_obs$include_in_fit
      ]
    )
model_expected_positive_integrity_ok <-
  isTRUE(
    all.equal(
      best_result$control_fit$expected_positive,
      primary_integrity_snapshot$control_expected_positive,
      tolerance = 1e-12
    )
  ) &&
  isTRUE(
    all.equal(
      best_result$intervention_fit$expected_positive,
      primary_integrity_snapshot$intervention_expected_positive,
      tolerance = 1e-12
    )
  )
profile_status <-
  if (
    RUN_ITN_PROFILE &&
    exists(
      "profile_ci",
      inherits = TRUE
    )
  ) {
    "COMPLETED"
  } else if (
    RUN_ITN_PROFILE
  ) {
    "NOT COMPLETED"
  } else {
    "NOT REQUESTED"
  }
boundary_status <-
  if (
    any(
      boundary_check$near_lower |
        boundary_check$near_upper
    )
  ) {
    "CAUTION"
  } else {
    "PASS"
  }
sensitivity_check_status <-
  if (!RUN_LIMITED_SENSITIVITY) {
    "NOT REQUESTED"
  } else if (
    identical(
      sensitivity_status_text,
      "COMPLETED"
    )
  ) {
    "PASS"
  } else {
    "CHECK"
  }
reconciliation_status <-
  if (
    exists(
      "screening_reconciliation",
      inherits = TRUE
    )
  ) {
    "PASS"
  } else {
    "FAIL"
  }
residual_diagnostic_status <-
  if (
    exists(
      "extreme_residuals",
      inherits = TRUE
    ) &&
    nrow(
      extreme_residuals
    ) > 0
  ) {
    "PASS"
  } else {
    "FAIL"
  }
plot_2001_status <-
  if (
    isTRUE(
      primary_fit_dates_2001_ok
    ) &&
    identical(
      fit_plot_start_date,
      as.Date("2001-01-01")
    ) &&
    identical(
      fit_plot_end_date,
      as.Date("2001-12-31")
    )
  ) {
    "PASS"
  } else {
    "FAIL"
  }
# 6.14.2 Thesis-readiness checklist and verification ----
thesis_checklist <-
  data.frame(
    domain = c(
      "PRIMARY MODEL",
      "PRIMARY MODEL",
      "PRIMARY MODEL",
      "PRIMARY MODEL",
      "DATA HANDLING",
      "DATA HANDLING",
      "DATA HANDLING",
      "UNCERTAINTY",
      "MODEL VALIDATION",
      "MODEL VALIDATION",
      "SENSITIVITY",
      "TRANSMISSION METRICS",
      "TRANSMISSION METRICS",
      "SCENARIO ANALYSIS",
      "SCENARIO ANALYSIS",
      "SCENARIO ANALYSIS",
      "SCENARIO ANALYSIS",
      "SCENARIO ANALYSIS",
      "PLOTTING",
      "PLOTTING",
      "CORE OUTPUT INTEGRITY",
      "CORE OUTPUT INTEGRITY",
      "CORE OUTPUT INTEGRITY"
    ),

    item = c(
      "Primary calibration completed",
      "Multistart optimisation converged",
      "Primary fitted parameters retained",
      "No fitted parameter at practical boundary",
      "Active-surveillance observations used in likelihood",
      "October observations excluded from primary likelihood and retained descriptively",
      "Screening-versus-full-trial reconciliation produced",
      "Profile likelihood completed",
      "Aggregate/residual diagnostics produced",
      "Extreme date-level residuals identified",
      "Limited one-way sensitivity analysis completed",
      "Model-specific control R0 calculated",
      "Model-specific intervention R0 calculated",
      "Effective-biting scenario analysis completed",
      "ITN percentage-coverage scenario analysis completed",
      "0% coverage gives q_ITN approximately 1",
      "35% coverage reproduces fitted q_ITN",
      "Five-year counterfactual projections completed",
      "Primary fitting plots restricted to 2001",
      "Five-year scenario plots retained as counterfactual projections",
      "Original observed positives unchanged",
      "Primary model-predicted positives retained",
      "Primary fitted model unchanged by secondary analyses"
    ),

    status = c(
      if (!is.null(best_result)) "PASS" else "FAIL",

      if (
        nrow(
          successful_fits
        ) > 0
      ) "PASS" else "FAIL",

      if (
        primary_fit_integrity_ok
      ) "PASS" else "FAIL",

      boundary_status,

      if (
        all(
          c(
            "number_positive",
            "number_tested",
            "include_in_fit"
          ) %in%
            names(
              control_obs
            )
        ) &&
        all(
          c(
            "number_positive",
            "number_tested",
            "include_in_fit"
          ) %in%
            names(
              intervention_obs
            )
        )
      ) "PASS" else "FAIL",

      if (
        any(
          control_obs$passive_case_finding
        ) ||
        any(
          intervention_obs$passive_case_finding
        )
      ) "PASS" else "CHECK",

      reconciliation_status,

      if (
        profile_status == "COMPLETED"
      ) "PASS" else if (
        profile_status == "NOT REQUESTED"
      ) "NOT REQUESTED" else "CHECK",

      if (
        exists(
          "fit_summary",
          inherits = TRUE
        ) &&
        exists(
          "residual_data",
          inherits = TRUE
        )
      ) "PASS" else "FAIL",

      residual_diagnostic_status,

      sensitivity_check_status,

      if (
        is.finite(
          R0_control
        )
      ) "PASS" else "FAIL",

      if (
        is.finite(
          R0_intervention
        )
      ) "PASS" else "FAIL",

      if (
        exists(
          "scenario_summary",
          inherits = TRUE
        )
      ) "PASS" else "FAIL",

      if (
        exists(
          "coverage_scenario_summary",
          inherits = TRUE
        )
      ) "PASS" else "FAIL",

      if (
        coverage_q0_ok
      ) "PASS" else "FAIL",

      if (
        coverage_q35_ok
      ) "PASS" else "FAIL",

      if (
        LONG_RUN_YEARS == 5 &&
        exists(
          "longrun_data",
          inherits = TRUE
        )
      ) "PASS" else "CHECK",

      plot_2001_status,

      if (
        LONG_RUN_YEARS == 5 &&
        exists(
          "combined_coverage_scenario_plot",
          inherits = TRUE
        )
      ) "PASS" else "CHECK",

      if (
        observed_positive_integrity_ok
      ) "PASS" else "FAIL",

      if (
        model_expected_positive_integrity_ok
      ) "PASS" else "FAIL",

      if (
        primary_fit_integrity_ok
      ) "PASS" else "FAIL"
    ),

    details = c(
      paste0(
        "Primary NLL = ",
        signif(
          best_result$objective,
          7
        )
      ),

      paste0(
        nrow(
          successful_fits
        ),
        "/",
        n_multistart,
        " starts converged"
      ),

      paste0(
        "b_m=",
        signif(
          best$b_m,
          7
        ),
        "; q_ITN=",
        signif(
          best$q_itn,
          7
        ),
        "; phi=",
        signif(
          best$phi_obs,
          7
        )
      ),

      if (
        boundary_status == "PASS"
      ) {
        "No primary fitted parameter near a practical boundary"
      } else {
        "At least one fitted parameter is near a practical boundary"
      },

      "number_positive / number_tested used directly in beta-binomial likelihood",

      "October excluded because ascertainment differed from the primary active-screening process; retained descriptively",

      paste0(
        "Reconciliation table rows = ",
        if (
          exists(
            "screening_reconciliation",
            inherits = TRUE
          )
        ) {
          nrow(
            screening_reconciliation
          )
        } else {
          0
        }
      ),

      profile_status,

      "Observed/modelled totals, positivity and residual diagnostics retained",

      paste0(
        "Residual diagnostic rows = ",
        if (
          exists(
            "extreme_residuals",
            inherits = TRUE
          )
        ) {
          nrow(
            extreme_residuals
          )
        } else {
          0
        }
      ),

      sensitivity_status_text,

      paste0(
        "Control model-specific R0 = ",
        signif(
          R0_control,
          6
        )
      ),

      paste0(
        "Intervention model-specific R0 = ",
        signif(
          R0_intervention,
          6
        )
      ),

      "Mechanistic intervention-strength scenarios change q_ITN only",

      "Coverage scenarios: 0, 20, 35, 40, 60, 80 and 100%",

      paste0(
        "q(0) = ",
        signif(
          q_from_coverage(0),
          7
        )
      ),

      paste0(
        "q(0.35) = ",
        signif(
          q_from_coverage(
            OBSERVED_ITN_COVERAGE
          ),
          7
        ),
        "; fitted q_ITN = ",
        signif(
          best$q_itn,
          7
        )
      ),

      paste0(
        LONG_RUN_YEARS,
        "-year model-based counterfactual projection under fixed assumptions"
      ),

      "Observed/fitted positivity, positive-count, residual and cumulative diagnostic displays restricted to calendar year 2001",

      "Coverage scenario remains a five-year model-based counterfactual projection; observed follow-up end marked on the trajectory panel",

      paste0(
        "Control active positives = ",
        sum(
          best_result$control_fit$number_positive
        ),
        "; intervention active positives = ",
        sum(
          best_result$intervention_fit$number_positive
        )
      ),

      paste0(
        "Control expected positives = ",
        signif(
          sum(
            best_result$control_fit$expected_positive
          ),
          7
        ),
        "; intervention expected positives = ",
        signif(
          sum(
            best_result$intervention_fit$expected_positive
          ),
          7
        )
      ),

      "best and best_result retained as the primary fitted model"
    ),

    stringsAsFactors = FALSE
  )
print(
  thesis_checklist,
  row.names = FALSE
)
write.csv(
  thesis_checklist,
  "thesis_modelling_final_checklist.csv",
  row.names = FALSE
)
acceptable_statuses <- c(
  "PASS",
  "COMPLETED",
  "NOT REQUESTED",
  "NOT APPLICABLE"
)
thesis_ready_for_final_interpretation <-
  all(
    thesis_checklist$status %in%
      acceptable_statuses
  )
cat(
  "\nTHESIS_READY_FOR_FINAL_INTERPRETATION = ",
  thesis_ready_for_final_interpretation,
  "\n",
  sep = ""
)
cat(
  "MODEL_RANDOM_SEED = ",
  MODEL_RANDOM_SEED,
  "\n",
  sep = ""
)
# 6.14.3 Final profile-likelihood and sensitivity diagnostics ----
print_block(
  "FINAL PROFILE AND SENSITIVITY DIAGNOSTIC SUMMARY"
)
profile_lower_final <-
  if (
    RUN_ITN_PROFILE &&
    exists("profile_ci", inherits = TRUE)
  ) {
    profile_ci$q_itn_lower_95_profile[1]
  } else {
    NA_real_
  }
profile_upper_final <-
  if (
    RUN_ITN_PROFILE &&
    exists("profile_ci", inherits = TRUE)
  ) {
    profile_ci$q_itn_upper_95_profile[1]
  } else {
    NA_real_
  }
profile_both_identified_final <-
  if (
    RUN_ITN_PROFILE &&
    exists("profile_ci", inherits = TRUE)
  ) {
    isTRUE(profile_ci$both_bounds_identified[1])
  } else {
    FALSE
  }
reduction_lower_final <-
  if (
    RUN_ITN_PROFILE &&
    exists("profile_ci", inherits = TRUE)
  ) {
    profile_ci$reduction_in_b_percent_lower[1]
  } else {
    NA_real_
  }
reduction_upper_final <-
  if (
    RUN_ITN_PROFILE &&
    exists("profile_ci", inherits = TRUE)
  ) {
    profile_ci$reduction_in_b_percent_upper[1]
  } else {
    NA_real_
  }
primary_fit_unchanged_after_secondary <-
  isTRUE(
    all.equal(
      best$b_m,
      primary_integrity_snapshot$b_m,
      tolerance = 1e-12
    )
  ) &&
  isTRUE(
    all.equal(
      best$q_itn,
      primary_integrity_snapshot$q_itn,
      tolerance = 1e-12
    )
  ) &&
  isTRUE(
    all.equal(
      best$phi_obs,
      primary_integrity_snapshot$phi_obs,
      tolerance = 1e-12
    )
  ) &&
  isTRUE(
    all.equal(
      best_result$objective,
      primary_integrity_snapshot$nll,
      tolerance = 1e-12
    )
  )
final_primary_diagnostic_summary <- data.frame(
  metric = c(
    "Primary minimum NLL",
    "Primary b_m",
    "Primary q_ITN",
    "Reduction in effective biting (%)",
    "Primary phi",
    "q_ITN lower 95% profile limit",
    "q_ITN upper 95% profile limit",
    "Reduction lower 95% profile limit (%)",
    "Reduction upper 95% profile limit (%)",
    "Both q_ITN profile limits identified",
    "Primary fit unchanged after secondary analyses"
  ),
  
  value = c(
    as.character(best_result$objective),
    as.character(best$b_m),
    as.character(best$q_itn),
    as.character(100 * (1 - best$q_itn)),
    as.character(best$phi_obs),
    as.character(profile_lower_final),
    as.character(profile_upper_final),
    as.character(reduction_lower_final),
    as.character(reduction_upper_final),
    as.character(profile_both_identified_final),
    as.character(primary_fit_unchanged_after_secondary)
  ),
  
  stringsAsFactors = FALSE
)
print(
  final_primary_diagnostic_summary,
  row.names = FALSE
)
write.csv(
  final_primary_diagnostic_summary,
  "final_primary_profile_diagnostic_summary.csv",
  row.names = FALSE
)
if (
  RUN_LIMITED_SENSITIVITY &&
  !is.null(primary_model_sensitivity_summary)
) {
  
  final_sensitivity_diagnostic_table <-
    primary_model_sensitivity_summary %>%
    dplyr::select(
      parameter_varied,
      value_tested,
      best_nll,
      delta_nll,
      convergence_count,
      b_m,
      q_itn,
      reduction_in_effective_biting_percent,
      phi_obs
    )
  
  cat(
    "\nSensitivity diagnostic table:\n"
  )
  
  print(
    final_sensitivity_diagnostic_table
  )
  
  write.csv(
    final_sensitivity_diagnostic_table,
    "final_sensitivity_diagnostic_table.csv",
    row.names = FALSE
  )
}
if (!primary_fit_unchanged_after_secondary) {
  warning(
    "Primary fitted values changed after secondary analyses. ",
    "Do not interpret profile or sensitivity results until the mutation is identified."
  )
}
cat(
  "\nFINAL MODEL RUN COMPLETE.\n"
)

# 6.15 Runtime diagnostics and output reproducibility ----
# 6.15.1 Exact runtime equivalence and performance summary ----
# Additional runtime checks and reporting. Scientific output tables above retain
# their original names and calculations.
if (MODEL_USE_SPEEDUPS) {
  .runtime_reference_best <- .model_objective_original(best_par, return_details = TRUE)
  .runtime_best_comparison <- all.equal(
    .runtime_reference_best, best_result, tolerance = 0, check.attributes = TRUE
  )
  .runtime_best_exact <- isTRUE(.runtime_best_comparison)
  write.csv(data.frame(exact_match = .runtime_best_exact),
            "runtime_best_fit_check.csv", row.names = FALSE)
  if (!.runtime_best_exact) {
    stop("The original computation did not exactly reproduce the final fitted result: ",
         paste(.runtime_best_comparison, collapse = "; "))
  }
}
.runtime_report <- data.frame(
  speedups_enabled = MODEL_USE_SPEEDUPS,
  exact_preflight_checks_passed = .model_speed$checks,
  cached_ODE_calls_reused = .model_speed$hits,
  ODE_calls_computed = .model_speed$misses,
  elapsed_seconds_after_preflight = proc.time()[['elapsed']] - .model_speed$started,
  cache_entries_retained = length(.model_speed$cache),
  cache_size_MiB = as.numeric(object.size(.model_speed$cache)) / 1024^2
)
write.csv(.runtime_report, "runtime_performance_summary.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "runtime_session_info.txt")
cat("\nRuntime summary (timing is specific to this computer and run):\n")
print(.runtime_report, row.names = FALSE)

# 6.15.2 Comparison of model output files between runs ----
# This checks CSV results, not image encodings or elapsed times.
compare_model_output_folders <- function(original_dir, faster_dir, tolerance = 0) {
  original_dir <- normalizePath(original_dir, mustWork = TRUE)
  faster_dir <- normalizePath(faster_dir, mustWork = TRUE)
  if (identical(original_dir, faster_dir)) stop("Choose two different output folders.")
  original_files <- list.files(original_dir, pattern = "\\.csv$", full.names = FALSE)
  faster_files <- list.files(faster_dir, pattern = "\\.csv$", full.names = FALSE)
  # Runtime diagnostic files are intentionally additional and are not model outputs.
  files <- sort(unique(c(original_files, faster_files)))
  files <- files[!grepl("^runtime_", files)]
  if (!length(files)) stop("No model CSV files were found to compare.")
  checks <- lapply(files, function(nm) {
    a_path <- file.path(original_dir, nm)
    b_path <- file.path(faster_dir, nm)
    if (!file.exists(a_path) || !file.exists(b_path)) {
      return(data.frame(file = nm, match = FALSE, status = "MISSING FILE"))
    }
    if (identical(unname(tools::md5sum(a_path)), unname(tools::md5sum(b_path)))) {
      return(data.frame(file = nm, match = TRUE, status = "BYTE IDENTICAL"))
    }
    a <- read.csv(a_path, stringsAsFactors = FALSE, check.names = FALSE)
    b <- read.csv(b_path, stringsAsFactors = FALSE, check.names = FALSE)
    difference <- all.equal(a, b, tolerance = tolerance, check.attributes = TRUE)
    data.frame(file = nm, match = isTRUE(difference),
               status = if (isTRUE(difference)) "VALUES MATCH" else
                 paste(difference, collapse = "; "))
  })
  report <- do.call(rbind, checks)
  print(report, row.names = FALSE)
  if (!all(report$match)) warning("At least one output is missing or different.")
  invisible(report)
}


# 6.16 Biological-parameter uncertainty: 200 refits ----
# Independent uniform draws: 11 fixed biological parameters +/-10%.
# For each draw, refit b_m, q_itn and phi_obs using the existing
# beta-binomial likelihood, ODE equations and current fit bounds.
# Results are a conditional-MLE sensitivity ensemble, not a
# conventional confidence interval or Bayesian posterior.

# 6.16.1 Uncertainty-refit function and configuration ----
chapter6_run_200_refits <- function(
    n_draws = 200L,
    n_starts = 10L,
    seed = 20261024L,
    checkpoint_every = 10L,
    output_dir = file.path(getwd(), "Chapter6_200_refit_uncertainty_qITN_0_to_2")
) {
  required_objects <- c(
    "bio_parms", "best", "best_par", "best_result",
    "lower_par", "upper_par", "control_obs", "intervention_obs",
    "objective_joint_fit", "unpack_fit_par"
  )
  missing_objects <- required_objects[
    !vapply(required_objects, exists, logical(1), envir = .GlobalEnv,
            inherits = FALSE)
  ]
  if (length(missing_objects)) {
    stop("Run the original Chapter 6 model first. Missing: ",
         paste(missing_objects, collapse = ", "))
  }
  if (length(n_draws) != 1L || !is.finite(n_draws) ||
      n_draws < 1L || n_draws != as.integer(n_draws) ||
      length(n_starts) != 1L || !is.finite(n_starts) ||
      n_starts < 1L || n_starts != as.integer(n_starts) ||
      length(checkpoint_every) != 1L || !is.finite(checkpoint_every) ||
      checkpoint_every < 1L || checkpoint_every != as.integer(checkpoint_every)) {
    stop("n_draws, n_starts and checkpoint_every must be positive integers.")
  }
  
  # 6.16.2 Biological-parameter uncertainty bounds ----
  bio_base <- bio_parms
  if (length(bio_base) != 11L || is.null(names(bio_base)) ||
      anyDuplicated(names(bio_base)) ||
      any(!is.finite(bio_base)) || any(bio_base <= 0)) {
    stop("Expected the model's 11 named, positive biological parameters.")
  }
  probability_names <- intersect(c("tau_mh", "tau_hm", "p_ret"),
                                 names(bio_base))
  if (any(bio_base[probability_names] * 1.1 > 1)) {
    stop("A +/-10% biological-parameter draw would exceed probability 1.")
  }
  
  # 6.16.3 Primary-fit snapshot and run configuration ----
  primary_snapshot <- list(
    best = best, best_par = best_par, best_result = best_result,
    bio_parms = bio_parms
  )
  configuration <- list(
    n_draws = as.integer(n_draws),
    n_starts = as.integer(n_starts),
    seed = as.integer(seed),
    base_biology = bio_base,
    original_fit = best_par,
    fit_bounds_lower = lower_par,
    fit_bounds_upper = upper_par,
    control_observations = control_obs,
    intervention_observations = intervention_obs,
    detection_mode = DETECTION_MODE,
    diagnostic_sensitivity = diagnostic_sensitivity,
    diagnostic_specificity = diagnostic_specificity,
    seasonal_amp = seasonal_amp_fixed,
    lambda_import = lambda_import_fixed
  )
  
  # 6.16.4 Output directories and parameter-range tables ----
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  refit_parameter_ranges <- data.frame(
    parameter = names(bio_base),
    baseline_value = as.numeric(bio_base),
    lower_bound = 0.90 * as.numeric(bio_base),
    upper_bound = 1.10 * as.numeric(bio_base),
    sampling = "Independent uniform within +/-10%",
    stringsAsFactors = FALSE
  )
  write.csv(refit_parameter_ranges,
            file.path(output_dir, "refit_biological_parameter_ranges.csv"),
            row.names = FALSE)
  checkpoint_file <- file.path(output_dir, "refit_checkpoint.rds")
  draw_rows <- vector("list", as.integer(n_draws))
  start_rows <- vector("list", as.integer(n_draws))
  
  # 6.16.5 Random seed preservation and checkpoint restoration ----
  # Restore the user's existing RNG state after these secondary refits.
  had_rng <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng) previous_rng <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_rng) {
      assign(".Random.seed", previous_rng, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  
  completed <- 0L
  if (file.exists(checkpoint_file)) {
    previous <- readRDS(checkpoint_file)
    if (!identical(previous$configuration, configuration)) {
      stop("Existing refit checkpoint belongs to a different model or ",
           "configuration. Choose another uncertainty output directory.")
    }
    completed <- as.integer(previous$completed)
    if (length(completed) != 1L || is.na(completed) ||
        completed < 0L || completed > n_draws ||
        length(previous$draw_rows) != completed ||
        length(previous$start_rows) != completed) {
      stop("Refit checkpoint is incomplete or invalid.")
    }
    if (completed > 0L) {
      draw_rows[seq_len(completed)] <- previous$draw_rows
      start_rows[seq_len(completed)] <- previous$start_rows
    }
    assign(".Random.seed", previous$rng_state, envir = .GlobalEnv)
    message("Resuming 200-refit uncertainty analysis after draw ", completed)
  } else {
    set.seed(seed)
  }
  
  # 6.16.6 Joint biological-parameter sampling and optimisation ----
  if (completed < n_draws) {
    for (d in seq.int(completed + 1L, as.integer(n_draws))) {
      bio_draw <- setNames(
        vapply(seq_along(bio_base), function(j) {
          runif(1L, min = 0.90 * bio_base[[j]],
                max = 1.10 * bio_base[[j]])
        }, numeric(1)),
        names(bio_base)
      )
      objective_draw <- function(par) {
        objective_joint_fit(par, bio = bio_draw)
      }
      one_draw_starts <- vector("list", as.integer(n_starts))
      
      for (s in seq_len(as.integer(n_starts))) {
        initial <- if (s == 1L) {
          best_par
        } else {
          setNames(runif(length(lower_par), min = lower_par,
                         max = upper_par), names(lower_par))
        }
        fit <- tryCatch(
          stats::nlminb(
            start = initial,
            objective = objective_draw,
            lower = lower_par,
            upper = upper_par,
            control = list(iter.max = 500L, eval.max = 2000L,
                           rel.tol = 1e-9)
          ),
          error = function(e) {
            list(par = initial, objective = Inf, convergence = 999L,
                 message = conditionMessage(e))
          }
        )
        one_draw_starts[[s]] <- data.frame(
          draw = d,
          start = s,
          nll = as.numeric(fit$objective),
          convergence = as.integer(fit$convergence),
          message = if (is.null(fit$message)) "" else as.character(fit$message),
          log_b_m = as.numeric(fit$par[["log_b_m"]]),
          eta_q_itn = as.numeric(fit$par[["eta_q_itn"]]),
          log_phi_obs = as.numeric(fit$par[["log_phi_obs"]]),
          stringsAsFactors = FALSE
        )
      }
      
      starts_d <- do.call(rbind, one_draw_starts)
      start_rows[[d]] <- starts_d
      good <- which(
        starts_d$convergence == 0L &
          is.finite(starts_d$nll) & starts_d$nll < 1e90 &
          is.finite(starts_d$log_b_m) &
          is.finite(starts_d$eta_q_itn) &
          is.finite(starts_d$log_phi_obs)
      )
      winner <- if (length(good)) good[which.min(starts_d$nll[good])] else NA_integer_
      
      if (is.na(winner)) {
        fitted_b <- NA_real_
        fitted_q <- NA_real_
        fitted_phi <- NA_real_
        best_nll_draw <- NA_real_
        near_boundary <- NA
      } else {
        par_draw <- c(
          log_b_m = starts_d$log_b_m[winner],
          eta_q_itn = starts_d$eta_q_itn[winner],
          log_phi_obs = starts_d$log_phi_obs[winner]
        )
        fitted <- unpack_fit_par(par_draw)
        fitted_b <- fitted$b_m
        fitted_q <- fitted$q_itn
        fitted_phi <- fitted$phi_obs
        best_nll_draw <- starts_d$nll[winner]
        scaled_distance <- (par_draw - lower_par) / (upper_par - lower_par)
        near_boundary <- any(scaled_distance < 0.01 |
                               scaled_distance > 0.99)
      }
      
      draw_rows[[d]] <- data.frame(
        draw = d,
        status = if (is.na(winner)) "FAILED" else "CONVERGED",
        nll = best_nll_draw,
        converged_starts = length(good),
        near_boundary = near_boundary,
        b_m = fitted_b,
        q_itn = fitted_q,
        phi_obs = fitted_phi,
        contact_reduction_percent = 100 * (1 - fitted_q),
        as.list(bio_draw),
        check.names = FALSE
      )
      
      if (d %% checkpoint_every == 0L || d == n_draws) {
        saveRDS(list(
          configuration = configuration,
          completed = d,
          draw_rows = draw_rows[seq_len(d)],
          start_rows = start_rows[seq_len(d)],
          rng_state = get(".Random.seed", envir = .GlobalEnv)
        ), checkpoint_file)
        message("Biological-parameter refits completed: ", d, "/", n_draws)
      }
    }
  }
  
  # 6.16.7 Refit results and empirical uncertainty percentiles ----
  refit_ensemble <- do.call(rbind, draw_rows)
  refit_optimisation_log <- do.call(rbind, start_rows)
  write.csv(refit_ensemble,
            file.path(output_dir, "refit_ensemble.csv"), row.names = FALSE)
  write.csv(refit_optimisation_log,
            file.path(output_dir, "refit_optimisation_log.csv"),
            row.names = FALSE)
  
  successful <- refit_ensemble[
    refit_ensemble$status == "CONVERGED" &
      is.finite(refit_ensemble$q_itn) &
      is.finite(refit_ensemble$b_m) &
      is.finite(refit_ensemble$phi_obs),
    , drop = FALSE
  ]
  if (nrow(successful) == 0L) {
    stop("No uncertainty refits converged. Inspect refit_optimisation_log.csv.")
  }
  
  metrics <- c("b_m", "q_itn", "phi_obs", "contact_reduction_percent")
  primary_values <- c(
    b_m = best$b_m,
    q_itn = best$q_itn,
    phi_obs = best$phi_obs,
    contact_reduction_percent = 100 * (1 - best$q_itn)
  )
  refit_uncertainty_summary <- do.call(rbind, lapply(metrics, function(metric) {
    values <- successful[[metric]]
    data.frame(
      parameter = metric,
      primary_estimate = unname(primary_values[[metric]]),
      refit_median = stats::median(values),
      lower_2.5_percentile = unname(stats::quantile(values, 0.025)),
      upper_97.5_percentile = unname(stats::quantile(values, 0.975)),
      minimum = min(values),
      maximum = max(values),
      converged_draws = nrow(successful),
      attempted_draws = n_draws,
      interval_type = "Empirical sensitivity ensemble percentiles",
      stringsAsFactors = FALSE
    )
  }))
  write.csv(refit_uncertainty_summary,
            file.path(output_dir, "refit_uncertainty_summary.csv"),
            row.names = FALSE)
  
  # 6.16.8 ITN contact-reduction distribution and integrity checks ----
  p_refit <- ggplot2::ggplot(
    successful,
    ggplot2::aes(x = contact_reduction_percent)
  ) +
    ggplot2::geom_histogram(bins = 25, fill = "#4477AA",
                            colour = "white") +
    ggplot2::geom_vline(xintercept = 100 * (1 - best$q_itn),
                        linetype = "solid", linewidth = 0.8) +
    ggplot2::geom_vline(
      xintercept = as.numeric(stats::quantile(
        successful$contact_reduction_percent, c(0.025, 0.975))),
      linetype = "dashed", linewidth = 0.7
    ) +
    ggplot2::labs(
      title = "ITN contact reduction across biological-parameter refits",
      subtitle = "200 intended draws; dashed lines show empirical 2.5th and 97.5th percentiles",
      x = "Reduction in effective biting (%)",
      y = "Number of converged refits"
    ) +
    ggplot2::theme_classic(base_size = 12)
  ggplot2::ggsave(
    file.path(output_dir, "refit_contact_reduction_distribution.png"),
    plot = p_refit, width = 9, height = 5.5, dpi = 300
  )
  
  unchanged <- identical(best, primary_snapshot$best) &&
    identical(best_par, primary_snapshot$best_par) &&
    identical(best_result, primary_snapshot$best_result) &&
    identical(bio_parms, primary_snapshot$bio_parms)
  if (!unchanged) {
    stop("Primary model objects changed during the uncertainty refits.")
  }
  
  # 6.16.9 Uncertainty-result object and saved outputs ----
  result <- list(
    draws = refit_ensemble,
    starts = refit_optimisation_log,
    summary = refit_uncertainty_summary,
    primary_model_unchanged = unchanged,
    seed = seed,
    output_dir = output_dir,
    interpretation = paste(
      "Empirical sensitivity ensemble conditional on independent +/-10%",
      "uniform variation of 11 fixed biological parameters.",
      "Not a statistical confidence interval or a Bayesian posterior."
    )
  )
  saveRDS(result, file.path(output_dir, "refit_ensemble.rds"))
  print(refit_uncertainty_summary, row.names = FALSE)
  cat("\nPrimary fitted model unchanged:", unchanged, "\n")
  cat("Converged refits:", nrow(successful), "of", n_draws, "\n")
  cat("Refit output folder:", normalizePath(output_dir), "\n")
  if (nrow(successful) < n_draws) {
    warning("Some refits did not converge. Inspect the optimisation log; ",
            "do not present the percentiles as a complete 200-refit range.")
  }
  invisible(result)
}

# 6.16.10 Execute 200 uncertainty refits ----
cat("\nCHAPTER 6: 200-REFIT BIOLOGICAL-PARAMETER UNCERTAINTY ANALYSIS\n")
chapter6_uncertainty_200 <- chapter6_run_200_refits(
  n_draws = 200L,
  n_starts = 10L,
  seed = 20261024L,
  checkpoint_every = 10L
)

