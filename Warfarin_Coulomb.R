# -----------------------------------------------------------------------------
# Integrated Shiny App: Warfarin-Albumin Molecular Dynamics + System Dynamics
# Author: Your Name
# Packages: shiny, bio3d, rgl, tidyverse, plotly
# -----------------------------------------------------------------------------

# Install and load required packages
if (!require("shiny")) install.packages("shiny")
if (!require("bio3d")) install.packages("bio3d")
if (!require("rgl")) install.packages("rgl")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("plotly")) install.packages("plotly")

library(shiny)
library(bio3d)
library(rgl)
library(tidyverse)
library(plotly)

# ---------------------------
# 1. Molecular Dynamics Data
# ---------------------------

# Load albumin structure (PDB: 1AO6)
albumin_pdb <- read.pdb("1ao6")  # Human Serum Albumin (HSA)

# Generate simplified warfarin structure (manual coordinates)
warfarin <- data.frame(
  x = c(0, 1, 2, 3, 4, 5, 0.5, 1.5),  # X coordinates
  y = c(0, 1, 0, -1, 0, 1, 0.5, -0.5), # Y coordinates
  z = c(0, 0, 0, 0, 0, 0, 1, -1),     # Z coordinates
  atom = c("C", "C", "O", "O", "C", "C", "H", "H")  # Atom types
)

# Sudlow Site I (primary warfarin binding site)
sudlow_site <- albumin_pdb$atom %>%
  filter(chain == "A", resid >= 195 & resid <= 198)  # Adjust based on PDB data

# Generate warfarin docking coordinates
set.seed(123)
warfarin_docked <- warfarin %>%
  mutate(
    x = x + runif(n(), min = 15, max = 18),  # Align to Sudlow Site I
    y = y + runif(n(), min = 15, max = 18),
    z = z + runif(n(), min = 15, max = 18)
  )

# Generate electron cloud data
generate_electron_cloud <- function(atom_coords, sigma = 1) {
  grid <- expand.grid(
    x = seq(min(atom_coords$x) - 2, max(atom_coords$x) + 2, length.out = 50),
    y = seq(min(atom_coords$y) - 2, max(atom_coords$y) + 2, length.out = 50),
    z = seq(min(atom_coords$z) - 2, max(atom_coords$z) + 2, length.out = 50)
  )
  grid$density <- 0
  for (i in 1:nrow(atom_coords)) {
    dist <- sqrt((grid$x - atom_coords$x[i])^2 +
                   (grid$y - atom_coords$y[i])^2 +
                   (grid$z - atom_coords$z[i])^2)
    grid$density <- grid$density + exp(-dist^2 / (2 * sigma^2))
  }
  return(grid)
}

warfarin_cloud <- generate_electron_cloud(warfarin_docked)

# ---------------------------
# 2. System Dynamics Simulation
# ---------------------------

# Define simulation parameters
simulate_system <- function(steps = 100) {
  # Initialize variables
  time <- 1:steps
  free_warfarin <- numeric(steps)
  antibody_concentration <- numeric(steps)
  rbc_count <- numeric(steps)
  
  # Initial conditions
  free_warfarin[1] <- 100  # Initial free warfarin (µM)
  antibody_concentration[1] <- 50  # Initial antibody concentration (nM)
  rbc_count[1] <- 1000  # Initial RBC count (cells/µL)
  
  # Simulate dynamics
  for (i in 2:steps) {
    # Warfarin binding to albumin
    free_warfarin[i] <- free_warfarin[i - 1] * 0.95  # Decay rate
    
    # Antibody dynamics (increase if foreign RBCs present)
    if (i > 20) {  # Introduce foreign RBCs at step 20
      antibody_concentration[i] <- antibody_concentration[i - 1] * 1.05
    } else {
      antibody_concentration[i] <- antibody_concentration[i - 1] * 0.98
    }
    
    # RBC lysis by antibodies
    if (i > 20) {
      rbc_count[i] <- rbc_count[i - 1] - 10  # Lysis rate
    } else {
      rbc_count[i] <- rbc_count[i - 1]  # No lysis
    }
  }
  
  # Return results as a data frame
  data.frame(
    Time = time,
    Free_Warfarin = free_warfarin,
    Antibody_Concentration = antibody_concentration,
    RBC_Count = rbc_count
  )
}

# Run simulation
sim_data <- simulate_system(steps = 100)

# ---------------------------
# 3. Shiny App UI
# ---------------------------

ui <- fluidPage(
  titlePanel("Warfarin-Albumin Interaction & System Dynamics"),
  sidebarLayout(
    sidebarPanel(
      h4("Visualization Options"),
      checkboxInput("show_albumin", "Show Albumin", TRUE),
      checkboxInput("show_warfarin", "Show Warfarin", TRUE),
      checkboxInput("show_electron_cloud", "Show Electron Cloud", FALSE),
      checkboxInput("show_sudlow_site", "Show Sudlow Site I", TRUE),
      hr(),
      h4("System Dynamics Metrics"),
      checkboxGroupInput(
        "metrics", "Select Metrics to Display:",
        choices = c(
          "Free Warfarin (µM)" = "Free_Warfarin",
          "Antibody Concentration (nM)" = "Antibody_Concentration",
          "RBC Count (cells/µL)" = "RBC_Count"
        ),
        selected = c("Free_Warfarin", "Antibody_Concentration", "RBC_Count")
      ),
      hr(),
      h4("Descriptive Statistics"),
      verbatimTextOutput("stats_output")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Molecular Visualization", rglwidgetOutput("rgl_plot", width = "800px", height = "800px")),
        tabPanel("System Dynamics", plotlyOutput("system_plot"))
      )
    )
  )
)

# ---------------------------
# 4. Shiny App Server Logic
# ---------------------------

server <- function(input, output) {
  # Descriptive statistics
  output$stats_output <- renderPrint({
    albumin_stats <- albumin_pdb$atom %>%
      summarise(
        mean_x = mean(x, na.rm = TRUE),
        mean_y = mean(y, na.rm = TRUE),
        mean_z = mean(z, na.rm = TRUE),
        sd_x = sd(x, na.rm = TRUE),
        sd_y = sd(y, na.rm = TRUE),
        sd_z = sd(z, na.rm = TRUE)
      )
    
    warfarin_stats <- warfarin_docked %>%
      summarise(
        mean_x = mean(x, na.rm = TRUE),
        mean_y = mean(y, na.rm = TRUE),
        mean_z = mean(z, na.rm = TRUE),
        sd_x = sd(x, na.rm = TRUE),
        sd_y = sd(y, na.rm = TRUE),
        sd_z = sd(z, na.rm = TRUE)
      )
    
    cat("Albumin Descriptive Statistics (Å):\n")
    print(albumin_stats)
    cat("\nWarfarin Descriptive Statistics (Å):\n")
    print(warfarin_stats)
  })
  
  # 3D Molecular Visualization
  output$rgl_plot <- renderRglwidget({
    open3d()
    bg3d(color = "white")
    
    # Show albumin
    if (input$show_albumin) {
      plot3d(
        albumin_pdb$atom[, c("x", "y", "z")],
        col = "gray",
        type = "s",
        radius = 0.5,
        aspect = FALSE
      )
    }
    
    # Show warfarin
    if (input$show_warfarin) {
      points3d(
        x = warfarin_docked$x,
        y = warfarin_docked$y,
        z = warfarin_docked$z,
        col = ifelse(warfarin_docked$atom == "C", "black",
                     ifelse(warfarin_docked$atom == "O", "red", "blue")),
        size = 10
      )
    }
    
    # Show electron cloud
    if (input$show_electron_cloud) {
      scatter3D(
        warfarin_cloud$x, warfarin_cloud$y, warfarin_cloud$z,
        colvar = warfarin_cloud$density,
        col = ramp.col(c("blue", "red")),
        alpha = 0.5,
        pch = 16,
        cex = 0.5,
        add = TRUE
      )
    }
    
    # Show Sudlow Site I
    if (input$show_sudlow_site) {
      spheres3d(
        x = sudlow_site$x,
        y = sudlow_site$y,
        z = sudlow_site$z,
        col = "gold",
        radius = 0.8
      )
    }
    
    rglwidget()
  })
  
  # System Dynamics Plot
  output$system_plot <- renderPlotly({
    plot_data <- sim_data %>%
      select(Time, all_of(input$metrics)) %>%
      pivot_longer(-Time, names_to = "Metric", values_to = "Value")
    
    plot_ly(plot_data, x = ~Time, y = ~Value, color = ~Metric, type = 'scatter', mode = 'lines') %>%
      layout(
        title = "System Dynamics Over Time",
        xaxis = list(title = "Time (Steps)"),
        yaxis = list(title = "Value"),
        hovermode = "x unified"
      )
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)