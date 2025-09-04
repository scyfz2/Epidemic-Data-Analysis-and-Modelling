# **Epidemic Data Analysis and Modelling**

## Project Overview
This project aims to analyze the measles outbreak in England during 2023-2024 using national surveillance data. It estimates key epidemic parameters such as growth rate and reproduction number, and develops a mathematical SEIR model to evaluate the impact of increased vaccination coverage. Despite the existence of an effective universal vaccination program, measles outbreaks have persisted, making this research significant for public health.

## Table of Contents
1. [Background](#background)
2. [Data Analysis](#data-analysis)
3. [Model Methodology](#model-methodology)
4. [Results](#results)
5. [Conclusion](#conclusion)
6. [References](#references)

## Background
Measles is a highly infectious viral disease caused by the Morbillivirus hominis, primarily spread through respiratory droplets. Despite an effective vaccination program, outbreaks continue to occur. The Basic Reproduction Number (R₀) for measles is estimated to be between 12 and 18, indicating its high transmissibility.

## Data Analysis
### Surveillance Data Biases
National surveillance data may be subject to underreporting, particularly in mild or asymptomatic cases. Reporting delays and misclassification errors may also occur, along with regional disparities in healthcare access, leading to an underestimation of the actual disease burden.

### Exploring the Dataset
- **Case Trend Figure**: Shows the temporal trend of confirmed measles cases from 2023 to 2024.
- **Epidemic Stages**:
  ![Monthly Confirmed Measles Cases in England (2023–2024)](/image/Monthly_Confirmed_Measles_Cases_England_2023-2024.png)
  - Low-level period (Jan-Oct 2023): Average monthly cases were 16.7.
  - Exponential growth period (Nov 2023 - May 2024): Cases increased from 41 to 396.
  - Decline period (Jun-Dec 2024): Average monthly decrease of approximately 50 cases.

### Persistence Despite Vaccination
Despite high vaccination coverage, outbreaks can occur due to:
1. Two-dose coverage not meeting the herd immunity threshold.
2. Waning immunity in certain populations.
3. Imported cases due to international travel.

## Model Methodology
### SEIR Model
- **Model Structure**:
  - S (Susceptible): Individuals who are susceptible to infection.
  - E (Exposed): Individuals who have been exposed to the virus but are not yet infectious.
  - I (Infectious): Individuals who are infected and can transmit the virus to others.
  - R (Recovered): Individuals who have recovered from infection or were vaccinated and are now immune.

### Model Equations

![SEIR Simulation of Measles in England (2023–2024)](/image/SEIR_Simulation.png)

The SEIR model is described by the following differential equations:

dS/dt = -β⋅S⋅I/N   
dE/dt = β⋅S⋅I/N − σE   
dI/dt = σE − γI   
dR/dt = γI 

Where:  
• β = 1.9855: transmission rate  
• σ = 1/12: incubation rate  
• γ = 1/8: recovery rate  
• R₀ = β / γ = 15.884  

### Model Implementation
The model is implemented in R, with the code provided in the file `EDAM_code.R`. Detailed comments explain the code functionality.

## Results

![Monthly New Measles Cases: With vs. Without Catch-up Vaccination](/image/Monthly_New_Measles_Cases.png)

- **Without Vaccination**: Approximately 944,000 cumulative cases.
- **With Vaccination**: Approximately 91,000 cumulative cases.
- **Cases Prevented**: Approximately 852,000 cases.

### Monthly New Cases Comparison
Visualizations demonstrate the comparison of monthly new cases with and without vaccination.

## Conclusion
This report analyzed the 2023-2024 measles outbreak in England using national surveillance data and a mathematical SEIR model. Despite high baseline vaccination coverage, an epidemic occurred, highlighting the impact of suboptimal herd immunity, waning protection, and imported cases. The implemented vaccination intervention significantly reduced cumulative cases, emphasizing the importance of timely enhanced vaccination strategies during outbreaks.

## References
1. Centers for Disease Control and Prevention. (n.d.). Measles vaccines. Retrieved April 23, 2025, from [CDC](https://www.cdc.gov/measles/vaccines/index.html)
2. Rota, P., Moss, W., Takeda, M. et al. (2016). Measles. Nat Rev Dis Primers, 2, 16049. Retrieved from [Nature](https://doi.org/10.1038/nrdp.2016.49)
3. UK Government. (2025). Population estimates for the UK, England, Wales, Scotland, and Northern Ireland: mid-2023. Retrieved from [UK Government](https://assets.publishing.service.gov.uk/media/67d875029dc953ac3bfe937b/1_Population_18_03_2025.pdf)

---

For further information, please refer to the code file `EDAM_code.R` and the report file `EDAM_report.pdf`.
