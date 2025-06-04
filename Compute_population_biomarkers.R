Compute_population_biomarkers <- function(Simresults){
  #'simresults' contains several Times and X lists, one for each simulated individual --> extract them:
  Times <- Simresults[grep("time", names(Simresults))]
  X <-Simresults[grep("X", names(Simresults))]
  IsJs <-Simresults[grep("IsJs", names(Simresults))]

  # full_results <- mapply(Compute_biomarkers(),X, Times, IsJs, SIMPLIFY = FALSE)
  full_results <- Map(Compute_biomarkers,X, Times, IsJs)

  flattened_list <- purrr::map(full_results, ~ .x) %>% purrr::flatten()
}

Compute_biomarkers = function(X, time, IsJs){
  # Compute biomarkers for a given individual

  # Calculate biomarkers
  APD90 <- APD_beats(time, X, 90)
  APD50 <- APD_beats(time, X, 50)
  APD40 <- APD_beats(time, X, 40)
  APD30 <- APD_beats(time, X, 30)

  Tri90_40 <- APD90$APDxx - APD40$APDxx
  Tri90_30 <- APD90$APDxx - APD30$APDxx

  CTD50 <- CTD_beats(time, X, 50, APD40$RAs)
  CTD90 <- CTD_beats(time, X, 90, APD40$RAs)

  qnet <- getQnet(time, IsJs)
  qnet_apd <- qnet$qnet_mean / APD90$APDxx

  EMw <- CTD90$CTDxx - APD90$APDxx

  list(X= X,
       Time = time,
       IsJs= IsJs,
       APD40 = APD40$APDxx,
       APD30 = APD30$APDxx,
       APD50 = APD50$APDxx,
       APD90 = APD90$APDxx,
       Tri90_40 = Tri90_40,
       Tri90_30 = Tri90_30,
       Vpp = APD40$Vpp,
       dVdt = APD40$dVdt,
       Vr = APD40$Vr,
       Vp = APD40$Vp,
       RAs = APD40$RAs,
       alternans = APD40$alternans,
       catd90 = CTD90$CTDxx,
       catd50 = CTD50$CTDxx,
       cai_syst = CTD50$cai_syst,
       cai_diast = CTD50$cai_diast,
       RAsCa = CTD50$RAsCa,
       qnet = qnet$qnet_mean,
       qnet_apd = qnet_apd,
       EMw = EMw)
}
