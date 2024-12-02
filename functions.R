library(survival)
library(foreach)

# Reformatting data ------------------------------------------------------- 
# Put data into required format, along with truncation (trunc) and minimum
# allowable surgery time (i.e. what do we set surgery at t=0 as for analysis?)
# For cloning/splitting, need to tell it the name of: surgery, time to
# surgery, death, and max follow-up. 
# The default arguments are names as they appear in ETTAA.
formatData <- function(data, trunc = 365.25 * 7, 
                       bump = 0.5,
                       surgery.id = 'intervention',
                       tts.id = 'final_proc_time',
                       fup.id = 'final_surv_total',
                       death.id = 'final_status',
                       id.id = 'ssid'){
  # copy
  .copy <- data
  data[,c(surgery.id, tts.id, fup.id, death.id, id.id)] <- NULL # remove from outputted set
  # Make these specific colnames so rest of script works
  data$id <- .copy[, id.id]
  data$surgery <- .copy[, surgery.id]
  data$timetosurgery <- .copy[, tts.id]
  data$death <- .copy[, death.id]
  data$fup_obs <- .copy[, fup.id]
  
  # Inflate procedure time of those who have it at t = 0 for analysis.
  data$timetosurgery[data$timetosurgery == 0] <- bump
  # Restrict outcome at truncation time (i.e. rmst)
  data$death[data$fup_obs > trunc] <- 0
  data$fup_obs[data$fup_obs > trunc] <- trunc
  
  # Neaten up the outputted data frame
  nn <- names(data)
  new.names <- nn[is.na(match(nn, names(.copy)))]
  subset(data, select = c(new.names, nn[-which(nn%in%new.names)]))
}


# Reducing data based on formulae -----------------------------------------
# One source of computational expense is in splitting the data (`survSplit`),
# this function aims to keep _only_ the variables specified in control.formula
# and surgery.formula i.e. the censoring models. Possible for these to be
# different.
reduceData <- function(data, surgery.formula, control.formula){
  data$.dummy <- rnorm(nrow(data))
  if(gsub("\\s", "", surgery.formula) == gsub("\\s", "", control.formula)){
    M <- as.data.frame(
      lm(as.formula(paste0(".dummy~", surgery.formula)), data, method = 'model.frame')
    )[,-1]
  }else{
    Ms <- as.data.frame(
      lm(as.formula(paste0(".dummy~", surgery.formula)), data, method = 'model.frame')
    )[,-1]
    Mc <- as.data.frame(
      lm(as.formula(paste0(".dummy~", control.formula)), data, method = 'model.frame')
    )[,-1]
    Mc$..id <- Ms$..id <- 1:nrow(Ms)
    nc <- colnames(Mc); ns <- colnames(Ms)
    c.in.s <- match(nc, ns)
    s.in.c <- match(ns, nc)
    
    # Add columns to control matrix which are in surgery, but not control
    allvars <- dplyr::select(dplyr::inner_join(Mc,Ms,'..id'), -`..id`)
    # if variable exists and is suffixed with '.y' and '.x' then it exists in
    # both Mc, Ms
    nv <- names(allvars)
    in.both.x <- nv[grepl("\\.x$", nv)]
    in.both.y <- nv[grepl("\\.y$", nv)]
    rest <- nv[-c(match(in.both.x, nv), match(in.both.y, nv))]
    M <- allvars[, c(in.both.x, rest)]
    names(M) <- gsub("\\.x$", "", names(M))
  }
  # Now that we have all of the variables, get the stuff from `data` which
  # allows us to clone data
  cbind(data.frame(
    id = data$id,
    surgery = data$surgery,
    timetosurgery = data$timetosurgery,
    death = data$death,
    fup_obs = data$fup_obs
  ), M)

}

# Cloning data ------------------------------------------------------------
# Takes a data.frame, grace period (default: 6 months), returns stacked
# data.frame
# Will need updating to handle other data set structures
# All credit to Maringe et al.
cloneData <- function(data, GP = 365.25/2){
  
  # Copy data to tab.control
  tab.control <- data
  # New variable names to be populated -->
  tab.control$outcome <- tab.control$fup <- tab.control$censoring <- 
    tab.control$fup.uncens <- NA
  tab.control$arm <- 'control'
  # And populate these for surgery while we're here.
  tab.surgery <- tab.control
  tab.surgery$arm <- 'surgery'
  
  # Control x {outcome, fup} ----
  
  # Outcome --> This is the the outcome in the emulated trial (which can be
  # different from the observed outcome).
  # fup --> the follow-up time in the emulated trial (which can be different from
  # the observed follow-up time).
  # These are the outcome and follow-up time in the outcome model!
  
  # CASE 1: They receive surgery in GP, survival up until surgery.
  # Set outcome = 0 since they don't die
  tab.control$outcome <- ifelse(tab.control$surgery == 1 & tab.control$timetosurgery <= GP, 
                                0, tab.control$outcome)
  # and set fup to surgery time
  tab.control$fup <- ifelse(tab.control$surgery == 1 & tab.control$timetosurgery <= GP, 
                            tab.control$timetosurgery, tab.control$fup)
  
  # CASE 2: Don't have surgery in GP, keep reported outcomes
  tab.control$outcome <- ifelse(tab.control$surgery == 0 | (tab.control$surgery == 1 & tab.control$timetosurgery > GP),
                                tab.control$death, tab.control$outcome)
  tab.control$fup <- ifelse(tab.control$surgery == 0 | (tab.control$surgery == 1 & tab.control$timetosurgery > GP),
                            tab.control$fup_obs, tab.control$fup)
  
  if(any(is.na(tab.control$outcome) | is.na(tab.control$fup))) stop("Some missing for {outcome, fup} control.\n")
  
  # Surgery x {outcome, fup} ----
  # Define the same new outcome variables {outcome, fup} for surgery arm.
  
  # CASE 1: They DID receive surgery in GP --> keep observed outcomes and tte.
  tab.surgery$outcome <- ifelse(tab.surgery$surgery == 1 & tab.control$timetosurgery <= GP, 
                                tab.surgery$death, tab.surgery$outcome)
  tab.surgery$fup <- ifelse(tab.surgery$surgery == 1 & tab.control$timetosurgery <= GP, 
                            tab.surgery$fup_obs, tab.surgery$fup)
  
  # CASE 2: They die or are censored in GP without surgery --> keep observed outcomes and tte.
  tab.surgery$outcome <- ifelse(tab.surgery$surgery == 0 & tab.control$fup_obs <= GP, 
                                tab.surgery$death, tab.surgery$outcome)
  tab.surgery$fup <- ifelse(tab.surgery$surgery == 0 & tab.control$fup_obs <= GP, 
                            tab.surgery$fup_obs, tab.surgery$fup)
  
  # CASE 3: They survive GP without surgery --> Alive with tte set as GP.
  tab.surgery$outcome <- ifelse((tab.surgery$surgery == 1 & tab.surgery$timetosurgery > GP) |
                                  (tab.surgery$surgery == 0 & tab.surgery$fup_obs > GP), 0, tab.surgery$outcome)
  tab.surgery$fup <- ifelse((tab.surgery$surgery == 1 & tab.surgery$timetosurgery > GP) |
                              (tab.surgery$surgery == 0 & tab.surgery$fup_obs > GP), GP, tab.surgery$fup)
  
  if(any(is.na(tab.surgery$outcome) | is.na(tab.surgery$fup))) stop("Some missing for {outcome, fup} surgery\n")
  
  # Defining censoring ----
  # Define new variables: {censoring, fup.uncensored} ----
  # censoring --> 0/1 to indicate if clone was censored in a given arm (either
  # because they receive surgery in the control arm or they didn't receive surgery
  # in the surgery arm).
  # fup.uncensored --> the follow-up time uncensored in the trial arm (can be
  # shorter than the follow-up time in the outcome model).
  # These are the outcome and follow-up time in the _weighting_ model!
  
  # Control x {censoring, fup.uncensored} ----
  # CASE 1: They receive surgery in GP, survival up until surgery.
  # censoring is 1 in control group at time of surgery
  tab.control$censoring <- ifelse(tab.control$surgery == 1 & tab.control$timetosurgery <= GP, 1, tab.control$censoring)
  tab.control$fup.uncens <-ifelse(tab.control$surgery == 1 & tab.control$timetosurgery <= GP, 
                                  tab.control$timetosurgery, tab.control$fup.uncens)
  
  # CASE 2: Die or LTFU in GP; keep true event time and leave uncensored
  tab.control$censoring <- ifelse(tab.control$surgery == 0 & tab.control$fup_obs <= GP, 0, tab.control$censoring)
  tab.control$fup.uncens <- ifelse(tab.control$surgery ==0 & tab.control$fup_obs <= GP, tab.control$fup_obs, tab.control$fup.uncens)
  
  # CASE 3: No surgery in GP and are still alive at GP; uncensored and their follow-up time is GP.
  tab.control$censoring <- ifelse((tab.control$surgery == 1 & tab.control$timetosurgery > GP) |
                                    (tab.control$surgery == 0 & tab.control$fup_obs > GP), 0, tab.control$censoring)
  tab.control$fup.uncens <- ifelse((tab.control$surgery == 1 & tab.control$timetosurgery > GP) |
                                     (tab.control$surgery == 0 & tab.control$fup_obs > GP), GP, tab.control$fup.uncens)
  
  if(any(is.na(tab.control$censoring) | is.na(tab.control$fup.uncens))) stop("Some missing for {censoring, fup.uncens} control.\n") # Check
  
  # Surgery x {censoring, fup.uncensored} ----
  # CASE 1: Patients receive surgery in GP; uncensored and have uncensored
  # follow-up until time of surgery.
  tab.surgery$censoring <- ifelse(tab.surgery$surgery == 1 & tab.surgery$timetosurgery <= GP, 0, tab.surgery$censoring)
  tab.surgery$fup.uncens <- ifelse(tab.surgery$surgery == 1 & tab.surgery$timetosurgery <= GP, 
                                   tab.surgery$timetosurgery, tab.surgery$fup.uncens)
  
  # CASE 2: Patients die or are LFTU in GP; uncensored and keep their follow-up
  # times
  tab.surgery$censoring <- ifelse(tab.surgery$surgery == 0 & tab.surgery$fup_obs <= GP, 0, tab.surgery$censoring)
  tab.surgery$fup.uncens <- ifelse(tab.surgery$surgery == 0 & tab.surgery$fup_obs <= GP, 
                                   tab.surgery$fup_obs, tab.surgery$fup.uncens)
  
  # CASE 3: No surgery in GP and are still alive at GP; censored and their follow-up time is GP.
  tab.surgery$censoring <- ifelse((tab.surgery$surgery == 0 & tab.surgery$fup_obs > GP) | 
                                    (tab.surgery$surgery == 1 & tab.surgery$timetosurgery > GP), 1, tab.surgery$censoring)
  tab.surgery$fup.uncens <- ifelse((tab.surgery$surgery == 0 & tab.surgery$fup_obs > GP) | 
                                     (tab.surgery$surgery == 1 & tab.surgery$timetosurgery > GP), GP, tab.surgery$fup.uncens)
  
  if(any(is.na(tab.surgery$censoring) | is.na(tab.surgery$fup.uncens))) stop("Some missing for {censoring, fup.uncens} surgery.\n") # Check
  
  out <- rbind(tab.control, tab.surgery)
  
  stopifnot(all(table(out$id)==2))
  out
}

# Splitting data by named arms
splitData <- function(data, arms = c('surgery', 'control'), t.events = NULL){
  # Can use this on e.g. bootstrapped data (hence finding the unique event times
  # and returning a df).
  t.events <- if(is.null(t.events)) sort(unique(data$fup)) else t.events
  t.events <- data.frame("Tevent" = t.events, "Tid" = 1:length(t.events))
  
  out <- do.call(rbind, splits <- lapply(arms, function(a){
    tab <- data[data$arm == a,]
    tab$Tstart <- 0
    
    # Splitting the dataset at each time of event until the event happens
    data.long <- survSplit(tab, cut = t.events$Tevent, end = "fup", 
                           start = "Tstart", event = "outcome", id = "ID")
    # Arranging by ID and event time.
    data.long <- dplyr::arrange(data.long, ID, Tstart)
    
    # Splitting the  dataset at each time of event until censoring happens. 
    # (this is to have the censoring status at each time of event).
    data.long.cens <- survSplit(tab, cut = t.events$Tevent, end = "fup", 
                                start = "Tstart", event = "censoring", id = "ID") 
    data.long.cens <- dplyr::arrange(data.long.cens, ID, Tstart)
    
    data.long.surv.cens <- survSplit(tab, cut = t.events$Tevent, end = "fup", 
                                     start = "Tstart", event = "eventCens", id = "ID") 
    
    # Replacing the censoring variable in data.long by the censoring variable obtained
    # in data.long.cens
    data.long$censoring <- data.long.cens$censoring
    data.long$eventCens <- data.long.surv.cens$eventCens
    
    # Creating Tstop (end of the interval) (is this needed??)
    data.long$Tstop <- data.long$fup
    
    # Merge, sort, fix identifiers for unique Tstarts
    data.long <- dplyr::left_join(data.long, t.events, 
                                  by = c("Tstart" = "Tevent"))
    data.long <- dplyr::arrange(data.long, ID, fup)
    data.long$Tid <- ifelse(is.na(data.long$Tid), 0, data.long$Tid)
    data.long
  }))
  
  list(tab = out,
       t.events = t.events)
}


# Estimating weights ------------------------------------------------------

# Fit Cox model, takes confounders in right hand side of formula. Can elect
# stabilised weights as well as whether IPTW/IPCW should be returned.
fitCox <- function(data, t.events,
                   rhs.form = "female + age + maxdiam + arch + copd1 + nyha + cad1 + ctd1 + diabetes1",
                   stabilised = FALSE, weights.only = TRUE, LTFU.mod = FALSE){
  if(!LTFU.mod){
    coxform <- as.formula(
      paste0("Surv(Tstart, Tstop, censoring) ~ ", rhs.form),
      env = parent.frame()
    )
  }else{
    coxform <- as.formula(
      paste0("Surv(Tstart, Tstop, eventCens) ~ ", rhs.form),
      env = parent.frame()
    )
  }
  
  
  mod <- coxph(coxform, data, x = T)
  rhs <- mod$x %*% coef(mod)
  
  BH <- as.data.frame(basehaz(mod, centered = FALSE))
  BH <- unique(dplyr::left_join(BH, t.events, by = c("time" = "Tevent")))
  dat2 <- dplyr::left_join(data, BH, "Tid")
  dat2$hazard <- ifelse(is.na(dat2$hazard), 0, dat2$hazard)
  dat2$pr.unc <- exp(-dat2$hazard*exp(rhs))
  
  if(stabilised){
    mod0 <- if(LTFU.mod) 
      coxph(Surv(Tstart, Tstop, eventCens) ~ 1, data) 
    else 
      coxph(Surv(Tstart, Tstop, censoring) ~ 1, data)
    BH0 <- as.data.frame(basehaz(mod0, centered = FALSE))
    names(BH0) <- c("haz0", "time")
    BH0 <- unique(dplyr::left_join(BH0, t.events, by = c("time" = "Tevent")))
    dat2 <- dplyr::left_join(dat2, BH0, "Tid")
    
    dat2$haz0 <-  ifelse(is.na(dat2$haz0), 0, dat2$haz0)
    dat2$pr.unc0 <- exp(-(dat2$haz0))
    dat2$w <- dat2$pr.unc0/dat2$pr.unc
  }else{
    dat2$w <- 1/dat2$pr.unc
  }
  
  # data$w[data$fup <= GP & data$arm == 'surgery'] <- 1
  
  # Return full data set as well as just the weights (for ease of merging)
  # or just the weights (the default), with full data then for debugging.
  if(weights.only)
    return(c(dat2$w))
  else
    return(list(
      data = dat2,
      mod = mod, mod0 = mod0,
      w = dat2$w
    ))  
  
}

# Fit Logistic regression model, put in right hand side of the formula for the censoring model
# and specify if these are the stabilised weights.
fitLogReg <- function(data,
                      rhs.form = "female + age + maxdiam + arch + copd1 + nyha + cad1 + ctd1 + diabetes1",
                      stabilised = FALSE, weights.only = TRUE, GP){
  glm.form <- as.formula(
    paste0("censoring ~ ", rhs.form),
    env = parent.frame()
  )
  
  # Surgery -> 1 up until the GP and then it's fixed based on their
  # characteristics thereafter. So, distinct version of the data dropping
  # Tstart, Tstop (occurs at Tid = 0).
  # Are they _ever_ censored
  to.merge <- data.frame(
    id = unique(data$id), 
    # Creating a new variable which is 1 if each id is _ever_ censored.
    is.censored = with(data, tapply(censoring, id, max))
  )
  # Merging on
  redx <- dplyr::left_join(dplyr::filter(data, Tid == 0),
                           to.merge, by = 'id')
  redx$censoring <- redx$is.censored # replace `censoring` (same as swapping-out formula)
  
  # Fit the model
  mod <- glm(glm.form, binomial, redx)
  redx$pr.unc <- 1 - plogis(predict(mod))
  
  if(stabilised){
    mod0 <- glm(censoring ~ 1, binomial, redx)
    redx$pr.unc0 <- 1 - plogis(predict(mod0))
    redx$w <- redx$pr.unc0/redx$pr.unc
  }else{
    redx$w <- 1/data$pr.unc
  }
  
  # Merge subject characteristic - specific weights back onto parent data
  data <- dplyr::left_join(data, dplyr::select(redx, id, w), by = 'id')
  
  # If follow-up is in GP then weight as one (arm: surgery consideration),
  # otherwise their weight is based on characteristics
  data$w[data$fup <= GP] <- 1
  
  # Return full data set as wellor just the weights (the default), with full
  # data then for debugging.
  if(weights.only)
    return(unname(data$w))
  else
    return(list(
      data = data,
      mod = mod, mod0 = mod0,
      w = data$w
    ))  
}


# Composite outcomes ------------------------------------------------------
estSurv <- function(data, rmst.time = 365.25){
  S.surgery <- survfit(Surv(Tstart, Tstop, outcome) ~ 1,
                   data = data[data$arm == "surgery",],
                   weights = w)
  # and for the control arm -->
  S.control <- survfit(Surv(Tstart, Tstop, outcome) ~ 1,
                   data = data[data$arm == "control",],
                   weights = w)
  
  # one-year survival
  S1 <- min(S.surgery$surv) # Survival to end of follow-up (i.e. last time point)
  S0 <- min(S.control$surv) # -"-
  
  # S(t_max) and HR don't change across different rmsts
  full.cox <- coxph(Surv(Tstart, Tstop, outcome) ~ arm, 
                    data = data, weights = w)
  
  # RMST --> Allow multiple RMSTs and get them all in one go!
  rmsts <- structure(sapply(rmst.time, function(r){
    RMST1 <- summary(S.surgery, rmean = r)$table['rmean']
    RMST0 <- summary(S.control, rmean = r)$table['rmean']
    c(RMST0, RMST1, RMST1 - RMST0, RMST1 / RMST0)
  }),
  dimnames = list(
    c("control", 
      "surgery",
      "surgery - control",
      "surgery / control"),
    paste0("t = ", round(rmst.time, 3))
  ))
  
  # Survival time taken at the RMSTs
  Sts.sur <- summary(S.surgery, times = rmst.time)
  Sts.con <- summary(S.control, times = rmst.time)
  S.rmsts <- structure(rbind(
    Sts.con$surv,
    Sts.sur$surv,
    Sts.sur$surv - Sts.con$surv,
    Sts.sur$surv / Sts.con$surv
  ), dimnames = list(
    c("S(t|control)", "S(t|surgery)", "surgery - control", "surgery / control"),
    paste0('t = ', rmst.time)
  ))
  
  out <- list(
    S0 = S0, S1 = S1, 
    S.diff = S1 - S0,
    RMSTs = rmsts,
    S.rmsts = S.rmsts,
    logHR = coefficients(full.cox),
    S.surgery = S.surgery,
    S.control = S.control
  )
  out
}


# Estimating competing risks ----------------------------------------------
estSurvCR <- function(data, rmst.time = 365.25, rmst.scale = 365.25){
  
  # Separate-out control and surgery sets to change ids, otherwise survfit fails
  # due to its `id` argument seeing two clusters.
  con <- data[data$arm == 'control',]; sur <- data[data$arm == 'surgery',]
  uids <- unique(con$id)
  con$id <- match(con$id, uids)
  sur$id <- match(sur$id, uids) + length(uids)
  data2 <- rbind(con, sur) # Rejoin
  
  sCR <- survfit(Surv(Tstart, Tstop, outcomeCR.fac) ~ arm,
                 data = data2, id = id, weights = w)
  
  # Pr(state shown) matrix
  pstate <- sCR$pstate
  colnames(pstate) <- sCR$states
  # sur.starts <- which(pstate==1,arr.ind = T)[2,1] # This splits after control arm
  sur.starts <- sCR$strata['arm=control'] + 1
  pstate.con <- pstate[1:(sur.starts-1),]
  pstate.sur <- pstate[sur.starts:nrow(pstate),]
  
  # S_0(t) and S_1(t) (i.e. probability of surviving to chosen horizon time).
  aa <- try(pstate.con[, 1], silent = T)
  if(inherits(aa, 'try-error')){
    cat("erroring: Here is pstate::\n")
    print(pstate)
    dbg <<- pstate
  }
  S0 <- min(pstate.con[,1]); S1 <- min(pstate.sur[,1])
  S.diff <- S1 - S0
  
  # RMSTs
  RMSTs <- do.call(cbind, 
                   lapply(seq_along(rmst.time), function(i){
                     r <- rmst.time[i]
                     temp <- survival:::survmean2(sCR, scale = rmst.scale, rmean = r)$matrix
                     colnames(temp)[3] <- sprintf("rmean: %.3f", r/rmst.scale)
                     if(i > 1) temp <- temp[,3,drop=F]
                     temp
                   })
  )
  
  # S(t) at each r \in RMSTs
  Sts <- summary(sCR, times = rmst.time)
  ps <- Sts$pstate
  row.names(ps) <- Sts$strata
  colnames(ps) <- Sts$states
  Sts <- cbind(ps, t = rmst.time)
  
  # Return list 
  list(
    sCR = sCR,
    S0 = S0, S1 = S1,
    S.diff = S.diff,
    RMSTs = RMSTs,
    rmst.scale = rmst.scale,
    St = Sts
  )
  
}

# Function to create data for Cox models above ----------------------------
createCoxdata <- function(control.data, surgery.data, 
                          c.weights, s.weights){
  rbind(
    cbind(dplyr::select(control.data, arm, id, Tstart, Tstop, outcome, censoring), w = c.weights), # control
    cbind(dplyr::select(surgery.data, arm, id, Tstart, Tstop, outcome, censoring), w = s.weights)  # surgery
  )
}

# Function to fit one set -------------------------------------------------
fixITB <- function(data, 
                   GP = 365.25, truncation.time = 365.25 * 7, bump = 0.5,
                   stabilised.weights = TRUE, 
                   rmsts = seq(1, 7, 2) * 365.25, rmst.scale = 1,
                   surgery.formula = "female + age + maxdiam + arch + copd1 + nyha + cad1 + ctd1 + diabetes1", 
                   control.formula = NULL, max.weight = NULL,
                   verbose = F, outputs = c("composite", "CR")){
  outputs <- unname(sapply(outputs, match.arg, outputs))
  # Format the data for use with cloneData function
  formatted.data <- formatData(data, trunc = truncation.time, bump = bump)
  # Assume control/surgery censoring model formulae are the same if it's not
  # provideds
  if(is.null(control.formula)) control.formula <- surgery.formula
  # Reduce data down to only necessary elements (largely for debugging)
  formatted.data <- reduceData(formatted.data, surgery.formula, control.formula)
  
  # Clone data
  clone.data <- cloneData(formatted.data, GP)
  clone.data$eventCens <- as.integer(clone.data$outcome == 0)
  if(verbose) cat("Data cloned.\n")
  # Long format data
  long.data <- splitData(clone.data)
  if(verbose) cat("Data in (Tstart, Tstop] format.\n")
  t.events <- long.data$t.events
  long.data <- long.data$tab
  
  if(verbose) cat("Obtaining weights...\n")
  # IPTWs
  con.weights <- fitCox(data = long.data[long.data$arm == 'control', ],
                        t.events = t.events, rhs.form = control.formula,
                        stabilised = stabilised.weights)
  
  sur.weights <- fitLogReg(data = long.data[long.data$arm == 'surgery', ], GP = GP,
                           rhs.form = surgery.formula, stabilised = stabilised.weights)
  
  if(verbose){ 
    cat("IPTWs obtained.\n --> Summary: Surgery arm:\n")
    print(summary(sur.weights))
    cat('\n--> Summary: Control arm:\n')
    print(summary(con.weights)); cat('\n')
  }
  
  # IPCWs
  con.c.weights <- fitCox(data = long.data[long.data$arm == 'control', ],
                          t.events = t.events, rhs.form = control.formula,
                          weights.only = T, LTFU.mod = T,
                          stabilised = stabilised.weights)
  sur.c.weights <- fitCox(data = long.data[long.data$arm == 'surgery', ],
                          t.events = t.events, rhs.form = surgery.formula,
                          weights.only = T, LTFU.mod = T,
                          stabilised = stabilised.weights)
  
  if(verbose){ 
    cat("IPCWs obtained.\n --> Summary: Surgery arm:\n")
    print(summary(sur.c.weights))
    cat('\n--> Summary: Control arm:\n')
    print(summary(con.c.weights)); cat('\n')
  }
  
  # Summary of the weights -->
  weight.summary <- rbind(IPTW.con = summary(con.weights),
                         IPCW.con = summary(con.c.weights),
                         comb.con = summary(con.weights * con.c.weights),
                         IPTW.sur = summary(sur.weights),
                         IPCW.sur = summary(sur.c.weights),
                         comb.sur = summary(sur.c.weights * sur.weights))
  
  # Create the final data set
  final.data <- createCoxdata(control.data = long.data[long.data$arm == 'control', ],
                              surgery.data = long.data[long.data$arm == 'surgery', ],
                              c.weights = cbind(iptw = con.weights, 
                                                ipcw = con.c.weights,
                                                combined = con.weights * con.c.weights), 
                              s.weights = cbind(iptw = sur.weights, 
                                                ipcw = sur.c.weights,
                                                combined = sur.weights * sur.c.weights))
  final.data$w <- final.data$w.combined
  
  # Setting max weight as needed
  max.weight <- if(is.null(max.weight)) Inf else max.weight
  final.data$w <- pmin(final.data$w, max.weight)
  if(!is.infinite(max.weight)) 
    ids.exceeding <- unique(final.data[final.data$w == max.weight, c('id', 'arm')]) 
  else 
    ids.exceeding <- data.frame()
  
  # Survival stuff --->
  if(verbose) cat("\nCalculating survival information of interest...\n")
  
  # Output 1 --> Stock survival analysis death vs no-death
  if("composite" %in% outputs) S.noCR <- estSurv(final.data, rmst.time = rmsts)
  
  # Output 2 --> Survival analysis no-death vs aneurysm death vs all-cause (a competing risk)
  if("CR" %in% outputs){
    # CR information --->
    df <- data.frame(
      id = formatted.data$id,
      AR = data$final_ar_status
    )
    
    final.data <- dplyr::left_join(final.data, df, 'id')
    final.data$outcomeCR <- 0
    final.data$outcomeCR <- ifelse(final.data$outcome == 1 & final.data$AR == 0, 1, final.data$outcomeCR)
    final.data$outcomeCR <- ifelse(final.data$outcome == 1 & final.data$AR == 1, 2, final.data$outcomeCR)
    # NB: final.data$outcomeCR is not directly comparable with e.g. table(data$death, data$AR)!
    
    final.data$outcomeCR.fac <- factor(final.data$outcomeCR, 0:2,
                                       c("censored", "other death", "aneurysm death"))
    
    # Survival items of interest
    S.CR <- estSurvCR(final.data, rmst.time = rmsts, rmst.scale = rmst.scale)
  }
  
  
  if(verbose) cat("\nDone!\n")
  out <- structure(
    list(weights = weight.summary, GP = GP, rmsts = rmsts, max.weight = max.weight,
         truncation.time = truncation.time, stabilised = stabilised.weights,
         surgery.formula = surgery.formula, ids.exceeding = ids.exceeding,
         control.formula = control.formula, data = final.data, outputs = outputs),
    class = 'fixITB'
  )
  
  if("composite" %in% outputs) out$noCR <- S.noCR
  if("CR" %in% outputs) out$CR <- S.CR
  
  out
}

# S3 print
print.fixITB <- function(x, show.RMSTs = TRUE, ...){
  cat("\n============================\n")
  cat("Trial emulation in the presence of immortal-time bias.\n")
  cat(sprintf("Grace-period used: %.2f. Original data truncated at %.2f\n", x$GP, x$truncation.time))
  cat(sprintf("%s weights were used.\n", ifelse(x$stabilised, "Stabilised", "Unstabilised")))
  if(length(x$outputs) == 2L){
    CRflag <- noCRflag <- T
    cat("Survival analysis was conducted on both the composite endpoint AND competing risks.")
  }else{
    cat("Survival analysis was carried out on", 
        ifelse(x$outputs == 'CR', "The competing risk of death due to aneurysm and other causes",
               "the composite endpoint"))
    CRflag <- if(x$outputs == 'CR') T else F
    noCRflag <- if(x$outputs == 'composite') T else F
  }
  cat("\n============================\n\n")
  cat("--> Summary of weights: Control clone arm:\n")
  print(x$weights['comb.con',])
  cat("--> Summary of weights: Surgery clone arm:\n")
  print(x$weights['comb.sur',])
  
  if(noCRflag){
    cat('\n', cli::rule(width = 10), '\n', sep = '')
    cat(cli::style_underline("Composite endpoints"), ":", sep = '')
    cat(sprintf("\nDifference in end-of-period survival: %.3f (RR: %.3f)\n", 
                x$noCR$S.diff,
                x$noCR$S1/x$noCR$S0))
    
    cat("\nSurvival probabilities at chosen RMSTs:\n")
    print(x$noCR$S.rmsts)
    
    if(show.RMSTs){
      cat("\nRMSTs:\n")
      print(x$noCR$RMSTs)
    }
    
    cat(sprintf("\nlog hazard ratio (hazard ratio): %.3f (%.3f)\n",
                x$noCR$logHR, exp(x$noCR$logHR)))
  }
  
  if(CRflag){
    cat('\n', cli::rule(width = 10), '\n', sep = '')
    cat(cli::style_underline("Competing Risks"), ":", sep = '')
    cat(sprintf("\nDifference in end-of-period survival: %.3f (RR: %.3f)\n", 
                x$CR$S.diff,
                x$CR$S1/x$CR$S0))
    
    cat("\nSurvival probabilities at chosen RMSTs:\n")
    tab <- x$CR$St
    con <- tab[1:length(x$rmsts),]; sur <- tab[(1 + length(x$rmsts)):nrow(tab),]
    invisible(lapply(x$CR$sCR$states, function(s){
      con <- con[, c('t', s)]
      sur <- sur[, c('t', s)]
      S1 <- sur[,s]; S0 <- con[,s]
      cat(cli::col_blue(paste0("Pr(", ifelse(s == '(s0)', "censored", s), ")",'\n')))
      print(structure(rbind(S0, S1, S1 - S0, S1 / S0),
                      dimnames = list(c("control", "surgery", "surgery - control", "surgery / control"),
                                      paste0("t = ", x$rmsts))))
      cat('\n')
    }))
    
    if(show.RMSTs){
      cat("Restricted mean time in state:\n")
      tab <- x$CR$RMSTs
      row.names(tab) <- gsub("arm\\=", '', row.names(tab))
      row.names(tab) <- gsub("\\, ", ' x ', row.names(tab))
      row.names(tab) <- gsub("(s0)", 'censoring', row.names(tab))
      colnames(tab) <- gsub("rmean\\:\\s", 't = ', colnames(tab))
      print(tab)
    }
  }
  
}

# S3 KM and incidence plots
plot.fixITB <- function(x, what = 'CR', xscale = 365.25,
                        # Legend arguments
                        show.legend = TRUE, legend.pos = 'topleft',
                        legend.size = 0.75, legend.border = FALSE,
                        # Arguments for controlling control/surgery
                        # line colours
                        control.col = 'black', surgery.col = 'red',
                        # Line types of competing risks
                        otherdeath.lty = 1, aneurysm.lty = 2,
                        # Show vertical line for GP?
                        show.GP = T, GP.lty = 2, GP.col = 'lightgrey',
                        # Show vertical line for truncation?
                        show.trunc = T, trunc.lty = 2, trunc.col = 'lightgrey',
                        # xlab, ylab
                        xlab = 'Time (years) from enrollment', ylab = NULL,
                        # Other arguments to pass to plot(...)
                        ...){
  if(what == 'CR' && "CR" %in% x$outputs){
    if(is.null(ylab)) ylab <- 'Cumulative incidence'
    leg <- paste0(rep(names(x$CR$sCR$strata), 2), ' x ', rep(x$CR$sCR$states[-1], each = 2))
    leg <- gsub("arm=", "", leg)
    plot(x$CR$sCR, col = rep(c(control.col, surgery.col), 2), 
         lty = c(otherdeath.lty, otherdeath.lty, 
                 aneurysm.lty, aneurysm.lty), 
         xscale = xscale, xlab = xlab, ylab = ylab, ...)
    # Vertical lines for GP and truncation if requested.
    if(show.GP) abline(v = x$GP, lty = GP.lty, col = GP.col)
    if(show.trunc) abline(v = x$truncation.time, lty = trunc.lty, col = trunc.col)
    # Legend
    if(show.legend){
      bty <- if(legend.border) 'o' else 'n'
      legend(legend.pos, bty = bty, 
             lty = c(otherdeath.lty, otherdeath.lty, 
                     aneurysm.lty, aneurysm.lty), 
             lwd = rep(1.25, 4),
             col = c(control.col, surgery.col, control.col, surgery.col), 
             cex = legend.size,
             legend = leg)
    }
  }
  
  if(what == 'KM' && 'CR' %in% x$outputs){
    if(is.null(ylab)) ylab <- expression(S(t))
    leg <- c('control', 'surgery')
    plot(x$CR$sCR, noplot = x$CR$sCR$states[-1],
         xscale = xscale, xlab = xlab, ylab = ylab,
         ylim = 0:1,
         col = c(control.col, surgery.col))
    # Vertical lines for GP and truncation if requested.
    if(show.GP) abline(v = x$GP, lty = GP.lty, col = GP.col)
    if(show.trunc) abline(v = x$truncation.time, lty = trunc.lty, col = trunc.col)
    # Legend
    if(show.legend){
      bty <- if(legend.border) 'o' else 'n'
      legend(legend.pos, bty = bty, 
             lty = c(1, 1), lwd = rep(1.25, 4),
             col = c(control.col, surgery.col),
             cex = legend.size, legend = leg)
    }
  }else if(what == 'KM' && x$outputs == 'composite'){
    if(is.null(ylab)) ylab <- expression(S(t))
    leg <- c('control', 'surgery')
    S <- survfit(Surv(Tstart, Tstop, outcome) ~ arm, x$data, weights = w) 
    plot(S, ylab = ylab, xlab = xlab,
         col = c(control.col, surgery.col), ylim = 0:1,
         xscale = xscale, ...)
    # Vertical lines for GP and truncation if requested.
    if(show.GP) abline(v = x$GP, lty = GP.lty, col = GP.col)
    if(show.trunc) abline(v = x$truncation.time, lty = trunc.lty, col = trunc.col)
    # Legend
    if(show.legend){
      bty <- if(legend.border) 'o' else 'n'
      legend(legend.pos, bty = bty, 
             lty = c(1, 1), lwd = rep(1.25, 4),
             col = c(control.col, surgery.col),
             cex = legend.size, legend = leg)
    }
  }
  
  if(what == 'CI'){
   if(is.null(ylab)) ylab <- ("Overall cumulative incidence") 
   if("CR" %in% x$outputs){
     plot(x$CR$sCR, noplot = x$CR$sCR$states[-1], fun = function(y) 1 - y,
          xscale = xscale, col = c(control.col, surgery.col),
          ylab = ylab, xlab = xlab, ...)
   }else if(x$outputs == 'composite'){
     S <- survfit(Surv(Tstart, Tstop, outcome) ~ arm, x$data, weights = w) 
     plot(S, ylab = ylab, xlab = xlab, fun = function(y) 1 - y,
          col = c(control.col, surgery.col), 
          xscale = xscale, ...)
   }
   
   # Vertical lines for GP and truncation if requested.
   if(show.GP) abline(v = x$GP, lty = GP.lty, col = GP.col)
   if(show.trunc) abline(v = x$truncation.time, lty = trunc.lty, col = trunc.col)
   # Legend
   if(show.legend){
     bty <- if(legend.border) 'o' else 'n'
     legend(legend.pos, bty = bty, 
            lty = c(1, 1), lwd = rep(1.25, 4),
            col = c(control.col, surgery.col),
            cex = legend.size, legend = c('control', 'surgery'))
   }
  }
}


# Two random plotter functions for cum. inc/haz + difference --------------
plot.cuminc <- function(x, lpos = 'topleft', ...){
  par(mfrow=c(1, 2))
  leg <- paste0(rep(names(x$CR$sCR$strata), 2), ' x ', rep(x$CR$sCR$states[-1], each = 2))
  leg <- gsub('arm=', '', leg)
  plot(x$CR$sCR, col = rep(c('black', 'red'), 2), lty = c(1, 1, 2, 2), 
       xscale = 365.25, xlab = 'Time (years) from enrollment',
       ylab = 'Cumulative incidence')
  legend(lpos, lty = c(1,1,2,2), col = c('black', 'red', 'black', 'red'),
         legend = leg, bty = 'n', cex = .75)
  # Working out difference in cumulative incidence
  pstates <- summary(x$CR$sCR, times = 0:x$truncation.time)
  ps <- pstates$pstate
  colnames(ps) <- pstates$states
  cc <- table(pstates$strata)[1]
  con <- ps[1:cc,]; sur <- ps[(cc + 1):nrow(ps),]
  diff <- sur - con
  ylims <- c(min(diff[, -1]), max(diff[, -1]))
  plot(0:x$truncation.time, diff[, pstates$states[2]],
       type = 'l', ylab = 'Difference in cumulative incidence', ylim = ylims,
       xlab = 'Time (years) from enrollment')
  abline(h = 0, lty = 'dotted')
  lines(0:x$truncation.time, diff[, pstates$states[3]], main = pstates$states[3],
        lty = 2)
  par(mfrow=c(1,1))
}

plot.cumhaz <- function(x, lpos = 'topleft', ...){
  par(mfrow=c(1, 2))
  leg <- paste0(rep(names(x$CR$sCR$strata), 2), ' x ', rep(x$CR$sCR$states[-1], each = 2))
  leg <- gsub('arm=', '', leg)
  plot(x$CR$sCR, col = rep(c('black', 'red'), 2), lty = c(1, 1, 2, 2), 
       xscale = 365.25, xlab = 'Time (years) from enrollment',
       ylab = 'Cumulative hazard', cumhaz = T)
  legend(lpos, lty = c(1,1,2,2), col = c('black', 'red', 'black', 'red'),
         legend = leg, bty = 'n', cex = .75)
  # Working out difference in cumulative incidence
  pstates <- summary(x$CR$sCR, times = 0:x$truncation.time)
  ps <- pstates$cumhaz
  colnames(ps) <- pstates$states[-1]
  cc <- table(pstates$strata)[1]
  con <- ps[1:cc,]; sur <- ps[(cc + 1):nrow(ps),]
  diff <- sur - con
  ylims <- c(min(diff), max(diff))
  plot(0:x$truncation.time, diff[, pstates$states[2]],
       type = 'l', ylab = 'Difference in cumulative hazard', ylim = ylims,
       xlab = 'Time (years) from enrollment')
  abline(h = 0, lty = 'dotted')
  lines(0:x$truncation.time, diff[, pstates$states[3]], main = pstates$states[3],
        lty = 2)
  par(mfrow=c(1,1))
}


# ////////////////////////////////////////////////////////////////////////////
# BOOTSTRAPPING FUNCTIONS                                                    /
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- /
# Defining function `bfit` which takes in a set of data from `cloneData` and /
# then performs a fixITB-style fit on a bootstrap resample. The wrapper      /
# function, bootfixITB is parallelised for macOS only. S3 methods are defined/
# for printing and plotting.                                                 /
# ////////////////////////////////////////////////////////////////////////////

# Bootstrapping functions -------------------------------------------------
# Function to obtain ONE bootstrap fit. This takes control arguments matching
# fixITB.CR and are meant to be supplied by boot.fixITB.CR.
bfit <- function(cloned.data, GP, truncation.time,
                 stabilised.weights, rmsts, 
                 surgery.formula, control.formula,
                 outputs, rmst.scale, AR.df, max.weight){
  con <- cloned.data[cloned.data$arm == 'control',]
  sur <- cloned.data[cloned.data$arm == 'surgery',]
  inds <- sample(1:nrow(con), nrow(con), replace = T)
  new.id <- rep(1:length(inds), 2)
  data <- rbind(con[inds,], sur[inds,])
  # Save old ids somewhere (a separate df in this case) for lookup later
  id.lookup <- data.frame(
    orig.id = con$id[inds],
    new.id = new.id[1:length(inds)]
  )
  data$id <- new.id
  # Split out the data
  long.data <- splitData(data)
  t.events <- long.data$t.events
  long.data <- long.data$tab
  if(is.null(control.formula)) control.formula <- surgery.formula
  
  # IPTWs
  con.weights <- fitCox(data = long.data[long.data$arm == 'control', ],
                        t.events = t.events, rhs.form = control.formula,
                        stabilised = stabilised.weights)
  
  sur.weights <- fitLogReg(data = long.data[long.data$arm == 'surgery', ], GP = GP,
                           rhs.form = surgery.formula, stabilised = stabilised.weights)
  
  # IPCWs
  con.c.weights <- fitCox(data = long.data[long.data$arm == 'control', ],
                          t.events = t.events, rhs.form = control.formula,
                          weights.only = T, LTFU.mod = T,
                          stabilised = stabilised.weights)
  sur.c.weights <- fitCox(data = long.data[long.data$arm == 'surgery', ],
                          t.events = t.events, rhs.form = surgery.formula,
                          weights.only = T, LTFU.mod = T,
                          stabilised = stabilised.weights)
  
  # Summary of the weights -->
  weight.summary <- rbind(IPTW.con = summary(con.weights),
                          IPCW.con = summary(con.c.weights),
                          comb.con = summary(con.weights * con.c.weights),
                          IPTW.sur = summary(sur.weights),
                          IPCW.sur = summary(sur.c.weights),
                          comb.sur = summary(sur.c.weights * sur.weights))
  
  # Create the final data set
  final.data <- createCoxdata(control.data = long.data[long.data$arm == 'control', ],
                              surgery.data = long.data[long.data$arm == 'surgery', ],
                              c.weights = cbind(iptw = con.weights, 
                                                ipcw = con.c.weights,
                                                combined = con.weights * con.c.weights), 
                              s.weights = cbind(iptw = sur.weights, 
                                                ipcw = sur.c.weights,
                                                combined = sur.weights * sur.c.weights))
  final.data$w <- final.data$w.combined
  
  # Setting max weight as needed
  max.weight <- if(is.null(max.weight)) Inf else max.weight
  final.data$w <- pmin(final.data$w, max.weight)
  
  # Output 1 --> Stock survival analysis death vs no-death
  if("composite" %in% outputs) S.noCR <- estSurv(final.data, rmst.time = rmsts)
  
  # Output 2 --> Survival analysis no-death vs aneurysm death vs all-cause (a competing risk)
  if("CR" %in% outputs){
    # CR information --->
    id.lookup <- dplyr::left_join(id.lookup, AR.df,
                                  by = c("orig.id" = 'ssid'))
    
    final.data <- dplyr::left_join(final.data, id.lookup, by = c('id' = 'new.id'))
    final.data$outcomeCR <- 0
    final.data$outcomeCR <- ifelse(final.data$outcome == 1 & final.data$final_ar_status == 0, 1, final.data$outcomeCR)
    final.data$outcomeCR <- ifelse(final.data$outcome == 1 & final.data$final_ar_status == 1, 2, final.data$outcomeCR)
    
    final.data$outcomeCR.fac <- factor(final.data$outcomeCR, 0:2,
                                       c("censored", "other death", "aneurysm death"))
    
    S.CR <- estSurvCR(final.data, rmst.time = rmsts, rmst.scale = rmst.scale)
  }
  
  out <- list(weights = weight.summary, GP = GP,
              truncation.time = truncation.time, stabilised = stabilised.weights,
              surgery.formula = surgery.formula, max.weight = max.weight,
              control.formula = control.formula, outputs = outputs, inds = inds)
  
  if("composite" %in% outputs) out$noCR <- S.noCR
  if("CR" %in% outputs) out$CR <- S.CR
  
  out
}

boot.fixITB <- function(data, GP = 365.25, truncation.time = 365.25 * 7, bump = 0.5,
                        stabilised.weights = TRUE, rmsts = seq(1, 7, 2) * 365.25, rmst.scale = 1,
                        surgery.formula = "female + age + maxdiam + arch + copd1 + nyha + cad1 + ctd1 + diabetes1",
                        control.formula = NULL, max.weight = NULL,
                        verbose = F, outputs = c("composite", "CR"),
                        B = 100L, ncores = floor(parallel::detectCores()/2)){
  outputs <- unname(sapply(outputs, match.arg, outputs))
  # Format the data
  formatted.data <- formatData(data, trunc = truncation.time, bump = bump)
  # Assume control/surgery censoring model formulae are the same if it's not
  # provided
  if(is.null(control.formula)) control.formula <- surgery.formula
  # Reduce data down to only necessary elements 
  formatted.data <- reduceData(formatted.data, surgery.formula, control.formula)
  
  # Clone data
  clone.data <- cloneData(formatted.data, GP)
  clone.data$eventCens <- as.integer(clone.data$outcome == 0)
  
  # Bootstrap in parallel 
  p <- proc.time()[3]
  if(ncores > 1){
    if(grepl('macOS', sessionInfo()$running)) 
      doMC::registerDoMC(cores = ncores)
    # else: <Register appropriate cluster for Windows use.>
    out <- foreach(i = 1:B, .verbose = T) %dopar% { # This is about 50% faster than doing sequentially
      bfit(clone.data, GP, truncation.time,
           stabilised.weights, rmsts, surgery.formula, control.formula,
           outputs, rmst.scale, AR.df, max.weight)
    }
  }else{ # If we want to do it sequentially instead...
    out <- lapply(1:B, function(x){
      bfit(clone.data, GP, truncation.time,
           stabilised.weights, rmsts, surgery.formula, control.formula,
           outputs, rmst.scale, AR.df, max.weight)
    })
  }
  
  class(out) <- "bootITB"
  attr(out, 'time.taken') <- proc.time()[3] - p
  out
}


# S3 Methods for these ----------------------------------------------------
print.bootITB <- function(x, ...){
  # Get information from the first bfit since they're all the same
  temp <- x[[1]]
  
  cat("\n============================\n")
  cat("Trial emulation in the presence of immortal-time bias.\n")
  cat(sprintf("This is based on %d bootstrapped fits (%.3f seconds).\n", length(x), attr(x, 'time.taken')))
  cat(sprintf("Grace-period used: %.2f. Original data truncated at %.2f\n", temp$GP, temp$truncation.time))
  cat(sprintf("%s weights were used.\n", ifelse(temp$stabilised, "Stabilised", "Unstabilised")))
  if(length(temp$outputs) == 2L){
    CRflag <- noCRflag <- T
    cat("Survival analysis was conducted on both the composite endpoint AND competing risks.")
  }else{
    cat("Survival analysis was carried out on", 
        ifelse(temp$outputs == 'CR', "The competing risk of deaths due to aneurysm and other causes",
               "the composite endpoint"))
    CRflag <- if(temp$outputs == 'CR') T else F
    noCRflag <- if(temp$outputs == 'composite') T else F
  }
  cat("\n============================\n\n")
  
  # Weights summary ->
  w <- do.call(rbind, lapply(x, function(y) rbind(con = y$weights['comb.con',], sur = y$weights['comb.sur',])))
  cat(cli::style_underline("Weight summary"), ":\n", sep = '')
  # min-of-min and max-of-max
  invisible(lapply(c("con", 'sur'), function(y){
    tab <- w[which(row.names(w) == y),]
    min.min <- min(tab[,'Min.']); max.max <- max(tab[,'Max.'])
    cat(sprintf("%s arm had minimum and maximum observed weights: %.3f, %.3f\n", y,
                min.min, max.max))
  }))
  cat('(Across all bootstraps)\n\n')
  
  # mean (SD) [95% CI] reporting
  .f <- function(y) format(round(y, 3), nsmall = 3)
  .report2 <- function(y){
    # m <- mean(y); s <- sd(y); qz <- qnorm(.975) * s
    # paste0(.f(m), " [", .f(m - qz), ', ', .f(m + qz), ']')
    m <- mean(y); qq <- quantile(y, c(.025, .975))
    paste0(.f(m), ' [', .f(qq[1]), ', ', .f(qq[2]), ']')
  }
  
  # Composite endpoints -->
  if(noCRflag){
    cat('\n', cli::rule(width = 10), '\n', sep = '')
    cat(cli::style_underline("Composite endpoints"), ":\n", sep = '')
    noCRs <- lapply(x, '[[', 'noCR')
    lHRs <- sapply(noCRs, '[[', 'logHR')
    cat(sprintf("Mean [95%% CI] log hazard ratio: %s\n\n", .report2(lHRs)))
    
    # End-of-period survival
    cat(sprintf("Mean [95%% CI] difference in end-of-period survival: %s\n", 
                .report2(sapply(noCRs, '[[', 'S.diff'))))
    S1 <- sapply(noCRs, '[[', 'S1'); S0 <- sapply(noCRs, '[[', 'S0')
    cat(sprintf("\t\tCorresponding percentage change: %s%%\n\n", 
                .report2((S1-S0)/S0 * 100)))
    
    # Survival probabilties at chosen RMSTs
    cat("Mean [95% CI] survival probabilities at chosen RMSTs:\n")
    tabs <- lapply(noCRs, '[[', 'S.rmsts')
    tabs.df <- do.call(rbind, lapply(seq_along(tabs), function(i){
      tab.i <- tabs[[i]]
      data.frame(
        b = i,
        quant = rep(row.names(tab.i), times = ncol(tab.i)),
        time = rep(colnames(tab.i), each = nrow(tab.i)),
        val = as.vector(tab.i)
      )
    }))
    med.Sts <- aggregate(val ~ time + quant, data = tabs.df, .report2)
    to.print <- do.call(cbind, lapply(unique(tabs.df$time), function(t){
      df <- med.Sts[med.Sts$time == t, -1]
      row.names(df) <- df$quant
      names(df)[2] <- t
      df[,-1,drop=F]
    }))
    print(to.print)
    
    # Now the actual RMSTs
    cat("\nMean [95% CI] RMSTs:\n")
    tabs <- lapply(noCRs, '[[', 'RMSTs')
    tabs.df <- do.call(rbind, lapply(seq_along(tabs), function(i){
      tab.i <- tabs[[i]]
      data.frame(
        b = i,
        quant = rep(row.names(tab.i), times = ncol(tab.i)),
        time = rep(colnames(tab.i), each = nrow(tab.i)),
        val = as.vector(tab.i)
      )
    }))
    med.Sts <- aggregate(val ~ time + quant, data = tabs.df, .report2)
    to.print <- do.call(cbind, lapply(unique(tabs.df$time), function(t){
      df <- med.Sts[med.Sts$time == t, -1]
      row.names(df) <- df$quant
      names(df)[2] <- t
      df[,-1,drop=F]
    }))
    print(to.print)
  }
  
  # Competing risks -->
  if(CRflag){
    cat('\n', cli::rule(width = 10), '\n', sep = '')
    cat(cli::style_underline("Competing risks"), ":\n", sep = '')
    CRs <- lapply(x, '[[', 'CR')
    
    # End-of-period survival
    cat(sprintf("Mean [95%% CI] difference in end-of-period survival: %s\n", 
                .report2(sapply(CRs, '[[', 'S.diff'))))
    S1 <- sapply(CRs, '[[', 'S1'); S0 <- sapply(CRs, '[[', 'S0')
    cat(sprintf("\t\tCorresponding percentage change: %s%%\n\n", 
                .report2((S1-S0)/S0 * 100)))
    
    cat("Mean [95% CI] survival probabilities at chosen RMSTs:\n")
    tabs <- lapply(CRs, '[[', 'St')
    tabs.df <- do.call(rbind, lapply(seq_along(tabs), function(i){
      tab.i <- tabs[[i]]
      rn <- row.names(tab.i); tab.i <- as.data.frame(tab.i)
      tab.i$arm <- gsub('arm=', '', rn); row.names(tab.i) <- NULL
      con <- tab.i[tab.i$arm == 'control',]; sur <- tab.i[tab.i$arm == 'surgery',]
      temp <- sur[, 1:3] - con[,1:3]
      temp$arm <- 'surgery - control'
      temp$t <- con$t
      tab.i <- rbind(tab.i, temp); tab.i$b <- i
      tab.i
    }))
    
    # now go over individual states and aggregate
    invisible(lapply(CRs[[1]]$sCR$states, function(s){
      tab.s <- tabs.df[, c(s,'t','arm','b')]
      names(tab.s)[1] <- 'val'
      agg <- aggregate(val ~ arm + t, tab.s, .report2)
      cat(cli::col_blue(paste0("Pr(", ifelse(s == '(s0)', "censored", s), ")",'\n')))
      out <- structure(sapply(unique(agg$t), function(t) agg[agg$t==t, 'val']),
                       dimnames = list(
                         unique(agg$arm), paste0('t = ', unique(agg$t))
                       ))
      print(as.data.frame(out))
      cat("\n")
    }))
    
  }
  
  invisible(x)
}

plot.bootITB <- function(x, what = 'CR',
                         # Legend arguments
                         show.legend = TRUE, legend.pos = 'topleft',
                         legend.size = 0.75, legend.border = FALSE,
                         # Arguments for controlling control/surgery
                         # line colours
                         control.col = 'black', surgery.col = 'red',
                         # thickness of mean line, thickness of individual lines
                         mean.lwd = 1, indiv.lwd = 0.15,
                         indiv.alpha = 0.25,
                         # Line types of competing risks
                         otherdeath.lty = 1, aneurysm.lty = 2,
                         # Show vertical line for GP?
                         show.GP = T, GP.lty = 2, GP.col = 'lightgrey',
                         # Show vertical line for truncation?
                         show.trunc = T, trunc.lty = 2, trunc.col = 'lightgrey',
                         # xlab, ylab
                         xlab = 'Time (years) from enrollment', ylab = NULL,
                         # Other arguments to pass to plot(...)
                         ...){
  CRs <- lapply(x, '[[', 'CR')
  sCRs <- lapply(CRs, '[[', 'sCR')
  states <- lapply(seq_along(sCRs), function(i){
    s <- sCRs[[i]]
    pstates <- summary(s, times = 0:x[[1]]$truncation.time)
    ps <- pstates$pstate
    colnames(ps) <- pstates$states
    cc <- table(pstates$strata)[1]
    con <- ps[1:cc,]; sur <- ps[(cc + 1):nrow(ps),]
    diff <- sur - con
    con <- as.data.frame(con); sur <- as.data.frame(sur); diff <- as.data.frame(diff)
    con$arm <- 'control'; sur$arm <- 'surgery'; diff$arm <- 'difference'
    out <- rbind(con, sur, diff)
    out$b <- i; out$t <- 0:x[[1]]$truncation.time
    out
  })
  states2 <- do.call(rbind, states)
  
  cens <- aggregate(`(s0)` ~ arm + t, states2, mean)
  ODs <- aggregate(`other death` ~ arm + t, states2, mean)
  ARs <- aggregate(`aneurysm death` ~ arm + t, states2, mean)
  
  t <- unique(cens$t)
  
  surgery.col.indiv <- grDevices::adjustcolor(surgery.col, alpha.f = indiv.alpha)
  control.col.indiv <- grDevices::adjustcolor(control.col, alpha.f = indiv.alpha)
  
  # KM plot
  if(what == 'KM'){
    if(is.null(ylab)) ylab <- expression(S(t))
    to.plot <- cens[cens$arm != 'difference', ]
    all.lines <- states2[states2$arm != 'difference', c('b', 't', 'arm', '(s0)')]
    ubs <- unique(all.lines$b)
    t <- unique(to.plot$t)
    # Surgery
    plot(pmin(1, to.plot[to.plot$arm == 'surgery', ]$`(s0)`) ~ t, type = 's',
         col = surgery.col, ylim = 0:1, xaxt = 'n', xlab = xlab,
         ylab = ylab, lwd = mean.lwd)
    axis(side = 1, at = seq(0, x[[1]]$truncation.time, by = 365.25), 
         labels = seq(0, x[[1]]$truncation.time%/%365.25, by = 1))
    for(b in ubs) lines(all.lines[all.lines$b == b & all.lines$arm == 'surgery', '(s0)'] ~ t, 
                        col = surgery.col.indiv, lwd = indiv.lwd)
    # Vertical lines for GP and truncation if requested.
    if(show.GP) abline(v = x[[1]]$GP, lty = GP.lty, col = GP.col)
    if(show.trunc) abline(v = x[[1]]$truncation.time, lty = trunc.lty, col = trunc.col)
    
    # Control
    lines(pmin(1, to.plot[to.plot$arm == 'control', ]$`(s0)`) ~ t, type = 's',
          col = control.col, lwd = mean.lwd)
    for(b in ubs) lines(all.lines[all.lines$b == b & all.lines$arm == 'control', '(s0)'] ~ t, 
                        col = control.col.indiv, lwd = indiv.lwd)
    
    # Legend
    if(show.legend){
      bty <- if(legend.border) 'o' else 'n'
      legend(legend.pos, bty = bty, 
             lty = c(1, 1), lwd = c(1, 1),
             col = c(control.col, surgery.col),
             cex = legend.size, legend = c('control', 'surgery'))
    }
  }
  
  if(what == 'CR'){
    if(is.null(ylab)) ylab <- 'Cumulative incidence'
    # Cumulative incidence plot
    to.plotOD <- ODs[ODs$arm != 'difference', ]
    to.plotAR <- ARs[ARs$arm != 'difference', ]
    all.lines <- states2[states2$arm != 'difference', c('b', 't', 'arm', 'aneurysm death', "other death")]
    ubs <- unique(all.lines$b)
    ylims <- c(0, max(c(all.lines$`aneurysm death`, all.lines$`other death`)))
    
    # Let's do Aneurysms first
    plot(to.plotAR[to.plotAR$arm=='surgery',]$`aneurysm death`~t, type = 's',
         xaxt = 'n', xlab = xlab, ylab = ylab,
         ylim = ylims, col = surgery.col, lwd = mean.lwd, lty = aneurysm.lty)
    axis(side = 1, at = seq(0, x[[1]]$truncation.time, by = 365.25), 
         labels = seq(0, x[[1]]$truncation.time%/%365.25, by = 1))
    # Vertical lines for GP and truncation if requested.
    if(show.GP) abline(v = x[[1]]$GP, lty = GP.lty, col = GP.col)
    if(show.trunc) abline(v = x[[1]]$truncation.time, lty = trunc.lty, col = trunc.col)
    for(b in ubs) lines(all.lines[all.lines$arm == 'surgery' & all.lines$b == b, 'aneurysm death'] ~ t,
                        col = surgery.col.indiv, lwd = indiv.lwd, lty = aneurysm.lty)
    
    lines(to.plotAR[to.plotAR$arm=='control',]$`aneurysm death`~t, col = control.col,
          lwd = mean.lwd, lty = aneurysm.lty)
    for(b in ubs) lines(all.lines[all.lines$arm == 'control' & all.lines$b == b, 'aneurysm death'] ~ t,
                        col = control.col.indiv, lwd = indiv.lwd, lty = aneurysm.lty)
    
    # Now other deaths
    lines(to.plotOD[to.plotOD$arm=='surgery',]$`other death` ~ t, 
          col = surgery.col, lwd = mean.lwd, lty = otherdeath.lty)
    for(b in ubs) lines(all.lines[all.lines$arm == 'surgery' & all.lines$b == b, 'other death'] ~ t,
                        col = surgery.col.indiv, lwd = indiv.lwd, lty = otherdeath.lty)
    
    lines(to.plotOD[to.plotOD$arm=='control',]$`other death` ~ t, col = control.col,
          lwd = mean.lwd, lty = otherdeath.lty)
    for(b in ubs) lines(all.lines[all.lines$arm == 'control' & all.lines$b == b, 'other death'] ~ t,
                        col = control.col.indiv, lwd = indiv.lwd, lty = otherdeath.lty)
    
  }
  
  if(what == 'CI'){ # This could be within the KM bracket, but im feeling lazy...
    if(is.null(ylab)) ylab <- 'Overall cumulative incidence'
    to.plot <- cens[cens$arm != 'difference', ]
    all.lines <- states2[states2$arm != 'difference', c('b', 't', 'arm', '(s0)')]
    ubs <- unique(all.lines$b)
    t <- unique(to.plot$t)
    ylims <- c(0, max(1 - all.lines$`(s0)`))
    # Surgery
    plot((1 - to.plot[to.plot$arm == 'surgery', ]$`(s0)`) ~ t, type = 's',
         col = surgery.col, ylim = ylims, xaxt = 'n', xlab = xlab,
         ylab = ylab, lwd = mean.lwd)
    axis(side = 1, at = seq(0, x[[1]]$truncation.time, by = 365.25), 
         labels = seq(0, x[[1]]$truncation.time%/%365.25, by = 1))
    for(b in ubs) lines(1 - all.lines[all.lines$b == b & all.lines$arm == 'surgery', '(s0)'] ~ t, 
                        col = surgery.col.indiv, lwd = indiv.lwd)
    # Vertical lines for GP and truncation if requested.
    if(show.GP) abline(v = x[[1]]$GP, lty = GP.lty, col = GP.col)
    if(show.trunc) abline(v = x[[1]]$truncation.time, lty = trunc.lty, col = trunc.col)
    
    # Control
    lines((1 - to.plot[to.plot$arm == 'control', ]$`(s0)`) ~ t, type = 's',
          col = control.col, lwd = mean.lwd)
    for(b in ubs) lines(1 - all.lines[all.lines$b == b & all.lines$arm == 'control', '(s0)'] ~ t, 
                        col = control.col.indiv, lwd = indiv.lwd)
    
    # Legend
    if(show.legend){
      bty <- if(legend.border) 'o' else 'n'
      legend(legend.pos, bty = bty, 
             lty = c(1, 1), lwd = c(1, 1),
             col = c(control.col, surgery.col),
             cex = legend.size, legend = c('control', 'surgery'))
    }
    
  }

}

# manual version
boot.fixITB2 <- function(data, GP = 365.25, truncation.time = 365.25 * 7, bump = 0.5,
                        stabilised.weights = TRUE, rmsts = seq(1, 7, 2) * 365.25,
                        surgery.formula = "female + age + maxdiam + arch + copd1 + nyha + cad1 + ctd1 + diabetes1", 
                        control.formula = NULL, max.weight = NULL,
                        verbose = F, outputs = c("composite", "CR"), rmst.scale = 1,
                        B = 100L, ncores = floor(parallel::detectCores()/2)){
  outputs <- unname(sapply(outputs, match.arg, outputs))
  # Format the data
  formatted.data <- formatData(data, trunc = truncation.time, bump = bump)
  # Assume control/surgery censoring model formulae are the same if it's not
  # provided
  if(is.null(control.formula)) control.formula <- surgery.formula
  # Reduce data down to only necessary elements
  formatted.data <- reduceData(formatted.data, surgery.formula, control.formula)
  
  # Clone data
  clone.data <- cloneData(formatted.data, GP)
  clone.data$eventCens <- as.integer(clone.data$outcome == 0)
  
  # Bootstrap in parallel 
  p <- proc.time()[3]
  pb <- utils::txtProgressBar(max = B, style = 3)
  out <- vector("list", B)
  for(i in 1:B){
    out[[i]] <- bfit(clone.data, GP, truncation.time,
                     stabilised.weights, rmsts, surgery.formula, control.formula,
                     outputs, rmst.scale, AR.df, max.weight)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  class(out) <- "bootITB"
  attr(out, 'time.taken') <- proc.time()[3] - p
  out
}



