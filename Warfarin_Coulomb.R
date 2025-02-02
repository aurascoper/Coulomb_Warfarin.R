# -----------------------------------------------------------------------------
# 統合Shinyアプリ：ワルファリン-アルブミン分子動力学 + システムダイナミクス
# 作者：あなたの名前
# パッケージ：shiny, bio3d, rgl, tidyverse, plotly, parallel
# -----------------------------------------------------------------------------

# 必要なパッケージを効率的にインストールおよびロード
packages <- c("shiny", "bio3d", "rgl", "tidyverse", "plotly", "parallel")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
lapply(packages, library, character.only = TRUE)

# ---------------------------
# 1. 分子動力学データ
# ---------------------------

# アルブミン構造をロード（PDB: 1AO6）、エラーハンドリングを追加
tryCatch({
  albumin_pdb <- read.pdb("https://files.rcsb.org/download/1AO6.pdb")
}, error = function(e) {
  stop("PDBファイルのロード中にエラーが発生しました。インターネット接続を確認してください。")
})

# ワルファリンの簡易構造を生成（手動座標）
warfarin <- data.frame(
  x = c(0, 1, 2, 3, 4, 5, 0.5, 1.5),  
  y = c(0, 1, 0, -1, 0, 1, 0.5, -0.5), 
  z = c(0, 0, 0, 0, 0, 0, 1, -1),     
  atom = c("C", "C", "O", "O", "C", "C", "H", "H")  
)

# Sudlow Site I（ワルファリンの主要結合部位）
sudlow_site <- albumin_pdb$atom %>%
  filter(chain == "A", as.numeric(resid) >= 195 & as.numeric(resid) <= 198)

# ワルファリンのドッキング座標を生成
set.seed(123)
warfarin_docked <- warfarin %>%
  mutate(
    x = x + runif(n(), min = 15, max = 18),  
    y = y + runif(n(), min = 15, max = 18),
    z = z + runif(n(), min = 15, max = 18)
  )

# 並列処理を使用して電子雲データを効率的に生成
generate_electron_cloud <- function(atom_coords, sigma = 1) {
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("atom_coords", "sigma"))
  cloud_data <- parLapply(cl, 1:50, function(i) {
    expand.grid(
      x = seq(min(atom_coords$x) - 2, max(atom_coords$x) + 2, length.out = 50),
      y = seq(min(atom_coords$y) - 2, max(atom_coords$y) + 2, length.out = 50),
      z = seq(min(atom_coords$z) - 2, max(atom_coords$z) + 2, length.out = 50)
    ) %>%
      mutate(density = exp(-((x - mean(atom_coords$x))^2 +
                               (y - mean(atom_coords$y))^2 +
                               (z - mean(atom_coords$z))^2) / (2 * sigma^2)))
  })
  stopCluster(cl)
  do.call(rbind, cloud_data)
}
warfarin_cloud <- generate_electron_cloud(warfarin_docked)

# ---------------------------
# 2. システムダイナミクスシミュレーション
# ---------------------------

simulate_system <- function(steps = 100, binding_rate = 0.95, lysis_rate = 10, antibody_growth = 1.05, dose = 100) {
  time <- 1:steps
  free_warfarin <- numeric(steps)
  antibody_concentration <- numeric(steps)
  rbc_count <- numeric(steps)
  plasma_concentration <- numeric(steps)
  
  free_warfarin[1] <- dose  
  antibody_concentration[1] <- 50  
  rbc_count[1] <- 1000  
  plasma_concentration[1] <- dose * 0.1  # 初期血漿濃度
  
  for (i in 2:steps) {
    free_warfarin[i] <- free_warfarin[i - 1] * binding_rate  
    
    antibody_concentration[i] <- if (i > 20) antibody_concentration[i - 1] * antibody_growth else antibody_concentration[i - 1] * 0.98
    
    rbc_count[i] <- if (i > 20) max(0, rbc_count[i - 1] - lysis_rate) else rbc_count[i - 1]
    
    plasma_concentration[i] <- plasma_concentration[i - 1] * 0.95  # 薬物クリアランスをシミュレート
  }
  
  data.frame(Time = time, Free_Warfarin = free_warfarin, Antibody_Concentration = antibody_concentration, RBC_Count = rbc_count, Plasma_Concentration = plasma_concentration)
}

# ---------------------------
# 3. Shinyアプリのユーザーインターフェース
# ---------------------------

ui <- fluidPage(
  titlePanel("ワルファリン-アルブミン相互作用とシステムダイナミクス"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("binding_rate", "ワルファリン結合率", min = 0.9, max = 1.0, value = 0.95, step = 0.005),
      sliderInput("lysis_rate", "赤血球溶血率", min = 5, max = 20, value = 10),
      sliderInput("antibody_growth", "抗体成長率", min = 1.0, max = 1.1, value = 1.05, step = 0.01),
      sliderInput("dose", "ワルファリン投与量", min = 50, max = 200, value = 100, step = 10),
      actionButton("update", "シミュレーションを更新")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("分子可視化", rglwidgetOutput("rgl_plot", width = "100%", height = "600px")),
        tabPanel("システムダイナミクス", plotlyOutput("system_plot"))
      )
    )
  )
)

# ---------------------------
# 4. Shinyアプリのサーバーロジック
# ---------------------------

server <- function(input, output) {
  output$rgl_plot <- renderRglwidget({
    open3d()
    bg3d(color = "white")
    points3d(warfarin_docked$x, warfarin_docked$y, warfarin_docked$z, col = "red", size = 10)
    spheres3d(sudlow_site$x, sudlow_site$y, sudlow_site$z, col = "gold", radius = 0.8)
    surface3d(warfarin_cloud$x, warfarin_cloud$y, warfarin_cloud$z, col = "blue", alpha = 0.1)
    rglwidget()
  })
  
  output$system_plot <- renderPlotly({
    sim_data <- simulate_system(100, input$binding_rate, input$lysis_rate, input$antibody_growth, input$dose)
    sim_data %>%
      pivot_longer(-Time, names_to = "Metric", values_to = "Value") %>%
      plot_ly(x = ~Time, y = ~Value, color = ~Metric, type = 'scatter', mode = 'lines') %>%
      layout(title = "時間経過によるシステムダイナミクス", xaxis = list(title = "時間"), yaxis = list(title = "値"))
  })
}

shinyApp(ui = ui, server = server)