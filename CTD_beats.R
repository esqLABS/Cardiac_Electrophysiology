CTD_beats <- function(time, X, xx, flag) {
  # CTDxx: Ca2+ transient duration at xx% of the initial base value (xx, e.g., 90%)
  # cai_syst: systolic calcium concentration (uM)
  # cai_diast: diastolic calcium concentration (uM)
  # RAsCa: repolarization abnormalities in the calcium transient
  # alternansCa: alternans in the calcium transient
  # t:time vector
  # X: output from modelRunner_drugs
  # xx: % of the initial base value to calculate the transient durantion
  # (eg. 90)
  # flag indicates if an EAD has been detected previously

  CTDxx <- numeric(length(X))
  cai_syst <- numeric(length(X))
  cai_diast <- numeric(length(X))
  RAsCa <- numeric(length(X))
  alternansCa <- numeric(length(X))

  for (i in seq_along(X)) {
    t <- time[[i]]
    cai <- X[[i]][, 6] * 1e3  # Assuming column 6 is calcium concentration

    result <- CTD(t, cai, xx, flag[i])
    CTDxx[i] <- result$CTDxx
    cai_syst[i] <- result$cai_syst
    cai_diast[i] <- result$cai_diast
    RAsCa[i] <- result$EAD
  }

  if (sum(diff(na.omit(CTDxx)) >= 5, na.rm = TRUE) >= 1 || sum(flag, na.rm = TRUE) >= 1) {
    alternansCa <- CTDxx
  }

  return(list(CTDxx = mean(CTDxx, na.rm = TRUE),
              cai_syst = mean(cai_syst, na.rm = TRUE),
              cai_diast = mean(cai_diast, na.rm = TRUE),
              RAsCa = RAsCa,
              alternansCa = alternansCa))
}

CTD <- function(t, cai, xx, flag) {
  cai_syst <- max(cai)
  cai_sys_index <- which.max(cai)

  t_ind_rep <- which(t > 300)[1]
  cai_diast <- min(cai[t_ind_rep:length(cai)])

  dCaidt <- diff(cai) / diff(t)

  if (sum(dCaidt[t_ind_rep:length(dCaidt)] > 0.0001) > 1 || flag == 1) {
    return(list(CTDxx = NA, cai_syst = NA, cai_diast = NA, EAD = 1))
  } else {
    Acai <- cai_syst - cai_diast
    cat_xxrec <- cai_syst - (xx / 100) * Acai
    t_catxx <- approx(cai[cai_sys_index:length(cai)], t[cai_sys_index:length(t)], xout = cat_xxrec)$y

    ind <- which.max(dCaidt)
    return(list(CTDxx = t_catxx - t[ind], cai_syst = cai_syst, cai_diast = cai_diast, EAD = 0))
  }
}
