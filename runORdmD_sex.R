runORdmD_sex <- function(BCL, scaling_factors, drug_effect, gender, BeatsSaved, model) {
  beats <- 300 # 500
  CL <- BCL
  phase <- 1

  # Load the initial condition (assuming it is saved as 'X0.rds' in R)
  X0 <- unname(as.matrix(read.csv("X0.csv", header = FALSE)))

  time <- vector("list", BeatsSaved)
  X <- vector("list", BeatsSaved)
  IsJs <- vector("list", BeatsSaved)

  options <- list()  # Placeholder for options, can be used for solver configuration

  pos <- 1
  for (n in 1:beats) {
    print(n)
    if (n <= beats - BeatsSaved) {
      # Solve the ODE (adjusted function call)
      out <- ode(y = X0, times = seq(0, CL, by = 0.5), func = ORdmD_sex_c, parms = list(
        flag_ode = 1,
        scalars = scaling_factors,
        drug_effect = drug_effect,
        cell = model,
        sex = gender,
        phase = phase
      ))
      X0 <- tail(out[, -1], 1)  # Update the initial condition for the next iteration

    } else if (n > beats - BeatsSaved) {
      # Solve and save the time and state variables
      out <- ode(y = X0, times = seq(0, CL, by = 0.5), func = ORdmD_sex_c, parms = list(
        flag_ode = 1,
        scalars = scaling_factors,
        drug_effect = drug_effect,
        cell = model,
        sex = gender,
        phase = phase
      ))
      time[[pos]] <- out[, 1]  # Time points
      X[[pos]] <- out[, -1]    # State variables
      X0 <- tail(out[, -1], 1)  # Update X0 for the next beat
      pos <- pos + 1
    }
  }

  # Calculate IsJs for each saved beat
  for (b in 1:BeatsSaved) {
    IsJs[[b]] <- matrix(NA, nrow = nrow(X[[b]]), ncol = 25)  # Initialize as matrix

    # Calculate fluxes and current for each time points
    for (i in seq_len(nrow(X[[b]]))){

      IsJs[[b]][i, ] <- ORdmD_sex_c(time[[b]][i], X[[b]][i,], list(
        flag_ode = 0,
        scalars = scaling_factors,
        drug_effect = drug_effect,
        cell = model,
        sex = gender,
        phase = phase
      ))
    }
  }


  # Return the results
  return(list(time = time, X = X, IsJs = IsJs))
}
