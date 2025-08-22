# Load required libraries
library(deSolve)   # For solving ODEs
library(readxl)    # For reading Excel files
library(doParallel) # For parallel execution
library(foreach)    # For parallel looping
library(compiler)

source("calculateDrug_effect.R")
source("ORdmD_sex.R")
source("runORdmD_sex.R")


#Pre-compile functions for computational efficiency
calculateDrug_effect_c <- cmpfun(calculateDrug_effect)
runORdmD_sex_c <- cmpfun(runORdmD_sex)
ORdmD_sex_c <- cmpfun(ORdmD_sex)

# Define parameters
drugs_to_simulate <- c('Control','Quinidine',"Clarithromycin",'Dofetilide')
# 'Azimilide', 'Dofetilide', 'Vandetanib', 'Quinidine', 'Disopyramide',
#'Sotalol', 'Chlorpromazine', 'Cisapride', 'Clarithromycin', 'Clozapine',
#'Droperidol', 'Domperidone', 'Pimozide', 'Risperidone', 'Ondansetron'

path_name <- "Drug_scenarios"  # Folder containing pharmacological data

# Conditions to simulate
BCLs <- 1000  # Basic cycle length (ms) (1000)
celltypes <- c("endo")  # Options: "endo", "mid", "epi"
BeatsSaved <- 5

# Load population scaling factors
scalars <- unname(as.matrix(read.csv("factors_population.csv", header = FALSE)))
# Ensure it's a matrix format (600x12)

# Set up parallel computing cluster
num_cores <- detectCores() - 1  # Use all cores except one
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Loop through each drug
for (drug_to_simulate in drugs_to_simulate) {

  # Read Excel data
  drug_data <- read_excel(paste0(path_name, "/", drug_to_simulate, "_maxC_IC50.xlsx"))
  # Set the first column as rownames and remove it from the data
  #drug_data <- drug_data[, -1]

  genders <- drug_data[[2]]
  FPCs <- drug_data[[3]]
  IC50s <- as.matrix(drug_data[,-(1:3)])
  scenarios <- drug_data[[1]]

  # Loop through each scenario
  for (cases in seq_along(scenarios)) {

    cat("Simulating:", drug_to_simulate, scenarios[cases], "\n")

    D <- FPCs[cases]
    IC50 <- IC50s[cases,]
    drug_effect <- calculateDrug_effect_c(D, IC50)
    gender <- genders[cases]

    if (gender == 1) {
      cat("Male population\n")
      scaling <- scalars[1:300, ] # Do 20 patients (1 and 2 already done)
      sex <- "male"
    } else if (gender == 2) {
      cat("Female population\n")
      scaling <- scalars[301:600, ] # Do 20 patients (1 and 2 already done)
      sex <- "female"
    }

    # Logging errors during parallel computing
    log_info <- function(message, index) {
      write(paste(Sys.time(), "Task", index, "Info:", message, "\n"),
            file = "debug_log.txt", append = TRUE)
    }


    # Loop through BCL and cell types
    for (BCL in BCLs) {
      for (celltype in celltypes) {

        My_function<- function(i){

          print(paste0("individual n°", i, "\n"))
          ind_scaling <- scaling[i, ]

          cat("Simulation starts..")
          tryCatch({
            log_info("Starting computation", i)
              # Run simulation
              sim_result <- runORdmD_sex_c(BCL, ind_scaling, drug_effect, gender, BeatsSaved, celltype)
            log_info(paste("Computation complete"), i)
            return(sim_result)
          }, error=function(e){
            log_error(conditionMessage(e),i)
            return(list(time= NA, X=NA, IsJs = NA))
          })

          #Time <- sim_result$time
          #X <- sim_result$X
          #IsJs <- sim_result$IsJs

          #sims_i <- X[[length(X)]][nrow(X[[length(X)]]),]

          #save_filename_i <- paste0(drug_to_simulate, "_", scenarios[cases], "_", sex, "_BCL", BCL, "_", celltype, "_ind", i, ".RDS")
          #saveRDS(sim_result, file = save_filename_i)  # Save each individual's results
        }

        tic <- Sys.time()  # Start timing

        # Run parallel computation for each individual
        results <- foreach(i = 1:nrow(scaling),
                            .packages = c("deSolve"),
                            .errorhandling = "pass",
                            .combine = "c",
                            .export= c("scaling", "BCL", "drug_effect", "gender", "BeatsSaved", "celltype","sex","cases" ,"drug_to_simulate", "runORdmD_sex_c")) %dopar% My_function(i)


        # Save results
        save_filename <- paste0("Results/", drug_to_simulate, "_", scenarios[cases], "_", sex, "_population_BCL", BCL, "_", celltype, ".RData")
        save(results, file = save_filename)

        toc <- Sys.time()  # End timing
        print(toc - tic)  # Print execution time
      }
    }
  }
}

# Stop parallel cluster
stopCluster(cl)
