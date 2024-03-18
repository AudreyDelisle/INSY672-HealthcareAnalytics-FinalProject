################# Health Care Assignment - Group 7 (Yvette and Audrey)

# Please note that to run this code, you must change the data file path below to match your file path to the COVIDdata file.

#############Step 0: Edit file path where covid data is saved
# File path to the Excel data
data_file_path <- "C:/Users/aadel/OneDrive/Documents/Mcgill Masters/3. Winter Term/Winter 1/HealthCare/Final Project/COVIDdata.xlsm"

########### Step 1: Install libraries
# Load required libraries
library(readxl)
library(shiny)
library(deSolve)

############################## Question 3 ####################################

plotSEIR <- function(beta = 0.00006, gamma = 0.15,theta = 0.005, S0 = 2040000, I0 = 80, R0 = 0,E0 = 10000, 
                    tmax = 100) {
  times <- seq(0, tmax, 1)
  
  parameters = c(beta = beta, gamma = gamma, theta = theta)
  
  iniState <- c(S = S0, E = E0, I = I0, R = R0)
  
  solution <- ode(iniState, times, SEIRmodel, parameters)
  
  par(mar = c(5, 5, 0, 0), oma = c(0, 0, 1, 1), mgp = c(2.5, 1, 0), xpd = T)
  plot(NA, xlab = "Time", ylab = "Number of individuals", ylim = c(0, 1000), 
       xlim = c(0, 60), cex.lab = 2)
  lines(x = solution[, "time"], y = solution[, "S"], col = "black", lwd = 3)
  lines(x = solution[, "time"], y = solution[, "E"], col = "green", lwd = 3)
  lines(x = solution[, "time"], y = solution[, "I"], col = "red", lwd = 3)
  lines(x = solution[, "time"], y = solution[, "R"], col = "blue", lwd = 3)
  
  legend(x = 'top', y = 1100, legend = c("S","E", "I", "R"), 
         col = c("black","green", "red", "blue"), lwd = 3, horiz = T, bty = "n", cex = 2) 
  
}

SEIRmodel <- function(times, state, parameters) {
  with (as.list(c(state, parameters)), {
    dS <- - beta * I * S 
    dE <- beta * I * S - theta * E
    dI <- theta * E - gamma * I
    dR <- gamma * I
    return(list(c(dS,dE, dI, dR)))
  })
}

app <- shinyApp(
  ui = fluidPage(
    titlePanel("SEIR Model"),
    
    sidebarLayout(
      sidebarPanel(      
        sliderInput("transmission", "Transmission Rate:",
                    min = 0, max = 0.02,
                    value = 0.000006, step = 0.000001),
        
        sliderInput("recovery", "Recovery Rate:",
                    min = 0, max = 1,
                    value = 0.15, step = 0.001),
        sliderInput("incubation", "Incubation Rate:",
                    min = 0, max = 1, 
                    value = 0.005, step = 0.001),
        
        sliderInput("S0", "Initial Susceptible Individuals:",
                    min = 100, max = 2040000,
                    value = 2040000, step = 10000),
        
        sliderInput("E0", "Initial Exposure Individuals",
                    min = 0, max = 10000, 
                    value = 10000, step = 100),
        
        sliderInput("I0", "Initial Infected Individuals:",
                    min = 1, max = 100,
                    value = 80, step = 1),
        
        sliderInput("tmax", "Time:",
                    min = 1, max = 100,
                    value = 21, step = 10)
        
      ),
      
      mainPanel(
        plotOutput("plot")
      ),
      
      "left"
    )
  ),
  
  server = function(input, output) {
    output$plot <- renderPlot({
      plotSEIR(beta = input$transmission, gamma = input$recovery, S0 = input$S0, E0 = input$E0, 
              I0 = input$I0, R0 = 0, tmax = input$tmax)
    })
  }
)

runApp(app)



#################################### Question 4 ##############################

############### Please note that to run this code, you must change the data file path at the beggining of the code to match your file path to the COVIDdata file.

# Read the COVID data from the specified Excel sheet
covid_statistics <- read_excel(data_file_path, sheet = "Sheet1")

# Define the SEIR differential equations model
SEIR_dynamics <- function(time_points, state_variables, model_params) {
  with(as.list(c(state_variables, model_params)), {
    rate_of_S_loss <- -beta * I * S 
    rate_of_E_gain_loss <- beta * I * S - theta * E
    rate_of_I_gain_loss <- theta * E - gamma * I
    rate_of_R_gain <- gamma * I
    
    return(list(c(rate_of_S_loss, rate_of_E_gain_loss, rate_of_I_gain_loss, rate_of_R_gain)))
  })
}

# Chosen parameters for SEIR model simulation
trans_rate <- 6e-5  # Transmission rate
exp_rate <- 0.005   # Rate from exposed to infectious
recov_rate <- 0.15  # Recovery rate
init_S <- 2.04e6    # Initial susceptible population
init_E <- 10000     # Initial exposed individuals
init_I <- 80        # Initial infected individuals
init_R <- 0         # Initial recovered individuals
simulation_length <- 21  # Number of time points to simulate

# Forecast function based on SEIR model dynamics
forecast_new_infections <- function(trans_rate, recov_rate, exp_rate, init_S, init_E, init_I, init_R, simulation_length) {
  simulation_times <- seq(0, simulation_length, by = 1)
  model_params <- c(beta = trans_rate, gamma = recov_rate, theta = exp_rate)
  initial_state <- c(S = init_S, E = init_E, I = init_I, R = init_R)
  
  # Solve the SEIR model
  model_output <- ode(y = initial_state, times = simulation_times, func = SEIR_dynamics, parms = model_params)
  
  # Calculate new positive cases
  new_positive_cases <- (model_output[,"I"][2:(simulation_length + 1)] - model_output[,"I"][1:simulation_length]) +
    (model_output[,"R"][2:(simulation_length + 1)] - model_output[,"R"][1:simulation_length])
  return(new_positive_cases)
}

# Execute the forecast function to get predicted cases
predicted_cases <- forecast_new_infections(trans_rate, recov_rate, exp_rate, init_S, init_E, init_I, init_R, simulation_length)

# Retrieve the actual newly positive cases from the data
actual_cases <- na.omit(covid_statistics$`Newly Positive (Confirmed Infections)`)

# Calculate RMSE between forecasted and actual cases
rmse_value <- sqrt(mean((predicted_cases - actual_cases) ^ 2))

# Display forecasted newly positive cases
cat("Forecasted Newly Positive Cases:", format(predicted_cases, scientific = FALSE, trim = TRUE), sep = "\n")

# Print RMSE
rmse_output <- paste("\nRoot Mean Square Error (RMSE):", rmse_value)
cat(rmse_output)
