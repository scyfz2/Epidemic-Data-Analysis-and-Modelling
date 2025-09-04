# Load required libraries
library(ggplot2)
library(dplyr)
library(deSolve)
library(lubridate)
library(tidyr)
# Load and prepare the data
measles_data <- data.frame(
  month = seq(as.Date("2023-01-01"), as.Date("2024-12-01"), by = "month"),
  cases = c(13,10,14,34,35,23,11,5,5,17,41,154,278,279,334,389,396,314,341,170,101,75,129,97)
)

# Plot case counts over time
ggplot(measles_data, aes(x = month, y = cases)) +
  geom_line(color = "red", linewidth = 1) +
  geom_point(color = "darkred") +
  labs(title = "Monthly Confirmed Measles Cases in England (2023–2024)",
       x = "Month", y = "Cases") +
  theme_minimal() +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = as.Date("2023-10-01"), linetype = "dashed", color = "blue") + 
  geom_vline(xintercept = as.Date("2024-05-01"), linetype = "dashed", color = "green")   

# Estimate exponential growth rate (log-linear model)
exp_data <- measles_data %>%
  filter(month >= as.Date("2023-11-01") & month <= as.Date("2024-05-01")) %>%
  mutate(t = 1:n(), log_cases = log(cases))

growth_model <- lm(log_cases ~ t, data = exp_data)
summary(growth_model)
r <- coef(growth_model)["t"]  # daily growth rate
D <- 12  # mean generation time in days
R <- 1 + r * D  # approximate R

# Estimate Vaccination Protection
p = 1 - 1 / R

# Use the growth rate to infer 𝛽
sigma <- 1/12
gamma <- 1/8
# r <- 0.315725

f <- function(beta) {
  left <- r
  right <- (- (sigma + gamma) + sqrt((sigma + gamma)^2 + 4 * sigma * beta)) / 2
  return(left - right)
}

# Solve the equation using uniroot
result <- uniroot(f, c(0, 10))
beta_est <- result$root
R0_est <- beta_est / gamma

cat("Estimated beta:", beta_est, "\n")
cat("Estimated R0:", R0_est, "\n")


#------------------------------------
#------------------------------------
# Parameter settings
N <- 56e6                  # Total population
I0 <- 13                   # Initial number of infected individuals (January 2023)
E0 <- 10                   # Initial number of exposed individuals (estimated from Jan–Feb 2023)
R0_initial <- 0.819 * N    # Initially immune population (81.9% coverage)
S0_initial <- N - I0 - E0 - R0_initial  # Initial susceptible population

# Model parameters
parameters <- c(
  beta = beta_est,     # Calibrated transmission rate
  sigma = 1/12,        # Transition rate from exposed to infectious (12-day incubation)
  gamma = 1/8,         # Recovery rate (8-day infectious period)
  v = 0.005            # Vaccination rate (0.5% per day)
)

# Time vector from 2023-01-01 to 2024-12-31 (730 days)
times <- seq(0, 730, by = 1)
start_vaccination <- 30  # Day when catch-up vaccination starts (2023-12-01)

# SEIRV model definition
seirv_model <- function(time, state, params) {
  with(as.list(c(state, params)), {
    # Start vaccination only after day 30
    vacc_rate <- ifelse(time >= start_vaccination, v, 0)
    
    dS <- -beta * S * I / N - vacc_rate * S
    dE <- beta * S * I / N - sigma * E
    dI <- sigma * E - gamma * I
    dR <- gamma * I + vacc_rate * S
    
    return(list(c(dS, dE, dI, dR)))
  })
}

# Initial states for the model compartments
initial_state <- c(S = S0_initial, E = E0, I = I0, R = R0_initial)

# Run the SEIRV model
output <- ode(
  y = initial_state,
  times = times,
  func = seirv_model,
  parms = parameters
)
output <- as.data.frame(output)

# Convert time to actual calendar date
output$date <- as.Date("2023-01-01") + output$time

# Plot infected population over time
ggplot(output, aes(x = date, y = I)) +
  geom_line(color = "red", linewidth = 1) +
  labs(title = "Adjusted Measles SEIRV Model Simulation (2023–2024)",
       x = "Date", y = "Number of Infected Individuals") +
  geom_vline(xintercept = as.Date("2023-12-01"), linetype = "dashed", color = "blue") +
  scale_x_date(date_breaks = "2 months", date_labels = "%Y-%m") +
  theme_minimal()

# Run scenario with no vaccination
parameters_no_vacc <- parameters
parameters_no_vacc["v"] <- 0
output_no_vacc <- ode(initial_state, times, seirv_model, parameters_no_vacc)
outTotal <- as.data.frame(output_no_vacc)

# Add total population column
outTotal$Total <- rowSums(outTotal[, c("S", "E", "I", "R")])

# Convert to long format for plotting

out_long <- pivot_longer(outTotal, cols = c("S", "E", "I", "R", "Total"),
                         names_to = "Compartment", values_to = "Count")

# Plot SEIR compartment counts
ggplot(out_long, aes(x = time, y = Count, color = Compartment)) +
  geom_line(size = 1.1) +
  labs(title = "SEIR Simulation of Measles in England (2023–2024)",
       x = "Time (Days)",
       y = "Population Count") +
  theme_minimal() +
  scale_color_manual(values = c(S = "#1f77b4", E = "#ff7f0e", I = "#d62728", R = "#2ca02c", Total = "black"))

# Filter for infected compartment only
out_I <- out_long[out_long$Compartment == "I", ]

# Plot infected population only
ggplot(out_I, aes(x = time, y = Count, color = Compartment)) +
  geom_line(size = 1.1) +
  labs(title = "Infection Count in SEIR Simulation of Measles in England (2023–2024)",
       x = "Time (Days)",
       y = "Infection Count") +
  theme_minimal() +
  scale_color_manual(values = c(I = "#d62728"))

# Function to calculate cumulative cases from model output
cumulative_cases <- function(output) {
  sum(output[, "I"] * parameters["gamma"])  # γI approximates daily new recoveries (i.e., incident cases)
}

# Calculate cumulative cases for both scenarios
cases_base <- cumulative_cases(output_no_vacc)  # ~944,000 cases
cases_vacc <- cumulative_cases(output)          # ~91,000 cases
cases_prevented <- cases_base - cases_vacc      # ~852,000 cases

# Key results
cat("Cumulative cases without vaccination:", round(cases_base / 1e5, 1), "x10^5 cases\n")
cat("Cumulative cases with catch-up vaccination:", round(cases_vacc / 1e5, 1), "x10^5 cases\n")
cat("Cases prevented:", round(cases_prevented / 1e5, 1), "x10^5 cases\n")

# Calculate monthly case counts
monthly_cases <- output %>%
  mutate(month = floor_date(date, "month")) %>%  # Group by month
  group_by(month) %>%
  summarise(
    total_I = sum(I),  # Total infected population in the month
    new_cases = sum(I * parameters["gamma"])  # Estimated new cases in the month
  ) %>%
  mutate(
    growth_rate = (total_I - lag(total_I)) / lag(total_I)  # Monthly growth rate of infections
  )

# View monthly results
print(monthly_cases)

#------------------------------------
# Calculate monthly new infections for both scenarios

# With vaccination
monthly_vacc <- output %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(new_cases = sum(I * parameters["gamma"])) %>%
  mutate(vaccination = "With Vaccination")

# Without vaccination
monthly_no_vacc <- output_no_vacc %>%
  as.data.frame() %>%
  mutate(date = as.Date("2023-01-01") + time) %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(new_cases = sum(I * parameters["gamma"])) %>%
  mutate(vaccination = "Without Vaccination")

# Combine both datasets
monthly_compare <- bind_rows(monthly_vacc, monthly_no_vacc)

# Plot comparison of monthly new cases
ggplot(monthly_compare, aes(x = month, y = new_cases, color = vaccination)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Monthly New Measles Cases: With vs. Without Catch-up Vaccination",
    x = "Month",
    y = "New Cases",
    color = "Scenario"
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%m-%Y", limits = c(as.Date("2023-01-01"), as.Date("2024-02-01"))) + # 修改为适合的日期范围
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

