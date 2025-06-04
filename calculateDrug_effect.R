calculateDrug_effect <- function(D, IC50) {
  # Returns a vector with the non-block fraction for each channel
  # D is the drug concentration
  # IC50 is a vector containing IC50 values and the Hill coefficients

  ICNa <- IC50[3]
  hNa <- IC50[4]
  ICNaL <- IC50[5]
  hNaL <- IC50[6]
  ICCaL <- IC50[7]
  hCaL <- IC50[8]
  ICKr <- IC50[1]
  hKr <- IC50[2]

  drug_effect <- numeric(4)  # Initialize vector of length 4

  if (is.na(ICNa) && is.na(hNa)) {
    drug_effect[1] <- 1
  } else {
    drug_effect[1] <- 1 / (1 + (D / ICNa)^hNa)
  }

  if (is.na(ICNaL) && is.na(hNaL)) {
    drug_effect[2] <- 1
  } else {
    drug_effect[2] <- 1 / (1 + (D / ICNaL)^hNaL)
  }

  if (is.na(ICCaL) && is.na(hCaL)) {
    drug_effect[3] <- 1
  } else {
    drug_effect[3] <- 1 / (1 + (D / ICCaL)^hCaL)
  }

  if (is.na(ICKr) && is.na(hKr)) {
    drug_effect[4] <- 1
  } else {
    drug_effect[4] <- 1 / (1 + (D / ICKr)^hKr)
  }

  return(drug_effect)
}
