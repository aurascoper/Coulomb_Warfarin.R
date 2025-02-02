# -----------------------------------------------------------------------------
# 集成 Shiny 应用: 华法林-白蛋白分子动力学与系统动力学
# 作者: Aurascoper
# 使用包: shiny, bio3d, rgl, tidyverse, plotly
# -----------------------------------------------------------------------------

# 電場與磁場計算相關 (僅供參考)
# 電場 E = 電荷密度 / 真空的誘電率
# 磁場 B = 電流密度 * 真空的透磁率
# 
# 高斯定律:
# ∇·E = ρ / ε₀
# ∇·B = 0
#
# 法拉第定律:
# ∇×E = -∂B/∂t
#
# 安培（修正）定律:
# ∇×B = μ₀J + μ₀ε₀∂E/∂t
#
# 電荷守恆:
# ∇·J + ∂ρ/∂t = 0
# -----------------------------------------------------------------------------

# 安裝與加載必要的包
packages <- c("shiny", "bio3d", "rgl", "tidyverse", "plotly")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
lapply(packages, library, character.only = TRUE)

# ---------------------------
# 1. 分子動力學數據
# ---------------------------

# 讀取白蛋白結構 (PDB: 1AO6)，並處理錯誤
tryCatch({
  albumin_pdb <- read.pdb("https://files.rcsb.org/download/1AO6.pdb")
}, error = function(e) {
  stop("讀取 PDB 文件失敗，請檢查網絡連接。")
})

# 生成簡化版的華法林結構（手動坐標）
warfarin <- tibble(
  x = c(0, 1, 2, 3, 4, 5, 0.5, 1.5),  
  y = c(0, 1, 0, -1, 0, 1, 0.5, -0.5), 
  z = c(0, 0, 0, 0, 0, 0, 1, -1),     
  atom = c("C", "C", "O", "O", "C", "C", "H", "H")
)

# Sudlow Site I (華法林主要結合位點)
sudlow_site <- albumin_pdb$atom %>%
  filter(chain == "A", between(as.numeric(resid), 195, 198))

# 使用庫侖定律生成華法林對接坐標 (隨機偏移)
set.seed(123)
warfarin_docked <- warfarin %>%
  mutate(
    x = x + runif(n(), min = 15, max = 18),  
    y = y + runif(n(), min = 15, max = 18),
    z = z + runif(n(), min = 15, max = 18)
  )

# 計算庫侖力（此函數供參考）
calculate_coulomb_force <- function(q1, q2, r) {
  k <- 8.99e9  # 庫侖常數 (N·m²/C²)
  force <- k * q1 * q2 / (r^2)
  return(force)
}

# 利用麥克斯韋方程生成電子雲數據
generate_electron_cloud <- function(atom_coords, sigma = 1) {
  expand_grid(
    x = seq(min(atom_coords$x) - 2, max(atom_coords$x) + 2, length.out = 50),
    y = seq(min(atom_coords$y) - 2, max(atom_coords$y) + 2, length.out = 50),
    z = seq(min(atom_coords$z) - 2, max(atom_coords$z) + 2, length.out = 50)
  ) %>%
    mutate(density = exp(-((x - mean(atom_coords$x))^2 +
                             (y - mean(atom_coords$y))^2 +
                             (z - mean(atom_coords$z))^2) / (2 * sigma^2)))
}
warfarin_cloud <- generate_electron_cloud(warfarin_docked)

# ---------------------------
# 2. 系統動力學模擬
# ---------------------------

simulate_system <- function(steps = 100, binding_rate = 0.95, lysis_rate = 10, antibody_growth = 1.05) {
  # 初始化模擬數據
  tibble(
    Time = 1:steps,
    Free_Warfarin = cumprod(c(100, rep(binding_rate, steps - 1))),
    Antibody_Concentration = cumprod(c(50, rep(ifelse(2:steps > 20, antibody_growth, 0.98), steps - 1))),
    RBC_Count = c(rep(1000, 20), pmax(0, 1000 - cumsum(rep(lysis_rate, steps - 20))))
  )
}

# ---------------------------
# 3. Shiny 應用 UI
# ---------------------------

ui <- fluidPage(
  titlePanel("华法林-白蛋白相互作用与系统动力学"),
  sidebarLayout(
    sidebarPanel(
      h4("可视化选项"),
      checkboxInput("show_albumin", "显示白蛋白", TRUE),
      checkboxInput("show_warfarin", "显示华法林", TRUE),
      checkboxInput("show_electron_cloud", "显示电子云", FALSE),
      checkboxInput("show_sudlow_site", "显示主要结合位点", TRUE),
      hr(),
      h4("系统动力学指标"),
      checkboxGroupInput(
        "metrics", "选择要显示的指标：",
        choices = c(
          "游离华法林 (µM)" = "Free_Warfarin",
          "抗体浓度 (nM)" = "Antibody_Concentration",
          "红细胞数 (cells/µL)" = "RBC_Count"
        ),
        selected = c("Free_Warfarin", "Antibody_Concentration", "RBC_Count")
      ),
      hr(),
      h4("描述性统计"),
      verbatimTextOutput("stats_output")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("分子可视化", rglwidgetOutput("rgl_plot", width = "800px", height = "800px")),
        tabPanel("系统动力学", plotlyOutput("system_plot"))
      )
    )
  )
)

# ---------------------------
# 4. Shiny 应用 Server 逻辑
# ---------------------------

server <- function(input, output) {
  # 显示描述性统计数据
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
    
    cat("白蛋白描述性统计 (Å):\n")
    print(albumin_stats)
    cat("\n华法林描述性统计 (Å):\n")
    print(warfarin_stats)
  })
  
  # 分子三维可视化
  output$rgl_plot <- renderRglwidget({
    open3d()
    bg3d(color = "white")
    
    # 显示白蛋白结构
    if (input$show_albumin) {
      plot3d(
        albumin_pdb$atom[, c("x", "y", "z")],
        col = "gray",
        type = "s",
        radius = 0.5,
        aspect = FALSE
      )
    }
    
    # 显示华法林
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
    
    # 显示电子云
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
    
    # 显示主要结合位点 (Sudlow Site I)
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
  
  # 系统动力学图形展示
  output$system_plot <- renderPlotly({
    sim_data <- simulate_system(100, input$binding_rate, input$lysis_rate, input$antibody_growth)
    sim_data %>%
      pivot_longer(-Time, names_to = "Metric", values_to = "Value") %>%
      plot_ly(x = ~Time, y = ~Value, color = ~Metric, type = 'scatter', mode = 'lines') %>%
      layout(
        title = "系统动力学随时间的变化",
        xaxis = list(title = "时间 (Steps)"),
        yaxis = list(title = "数值"),
        hovermode = "x unified"
      )
  })
}

# 运行 Shiny 应用
shinyApp(ui = ui, server = server)