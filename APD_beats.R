APD_beats <- function(time, X, percent) {
  APDxx <- numeric(length(X))
  Vpp <- numeric(length(X))
  dVdt <- numeric(length(X))
  Vr <- numeric(length(X))
  Vp <- numeric(length(X))
  RAs <- logical(length(X))
  alternans <- list() # unknown functionality, not used later, just returned

  for (i in seq_along(X)) {
    t <- time[[i]]
    V <- X[[i]][, 1]  # The first column represents voltage

    V_max <- max(V)
    V_min <- min(V)
    Vr[i] <- V_min
    Vp[i] <- V_max
    Vpp[i] <- V_max - V_min
    dVdt[i] <- max(diff(V) / diff(t))

    V_threshold <- V_max - (percent / 100) * (V_max - V_min)

    # Find indices where voltage crosses the threshold
    above_threshold <- which(V >= V_threshold)

    if (length(above_threshold) == 0) {
      return(NA)  # Return NA if no crossings found
    }

    t_depol <- t[min(above_threshold)]  # Time at depolarization
    t_repol <- t[max(above_threshold)]  # Time at repolarization

    APDxx[i] <- t_repol - t_depol # Compute APD duration

    # Detect EADs: if there are multiple positive slopes in repolarization phase
    dVdt_rep <- diff(V[t > min(t) + 100]) / diff(t[t > min(t) + 100]) #Calculate rate change of voltage during repolarization phase (assuming it occurs later than 100ms)
    RAs[i] <- sum(dVdt_rep > 0.0001) > 1 # Counts the number of instances where the slope is positive (threshold 0.0001 to filter out small fluctuations), single occurrence considered artifact.
  }

  return(list(APDxx = APDxx, Vpp = Vpp, dVdt = dVdt, Vr = Vr, Vp = Vp, RAs = RAs, alternans = alternans))
}
