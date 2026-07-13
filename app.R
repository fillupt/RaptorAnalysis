library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(plotly)
library(tidyr)
library(nlme)
library(shinycssloaders)
library(shinyjs)
library(mgcv)
library(stringr)

# Apply polynomial contrasts to non-ordered factors
options(contrasts = c(
   unordered = "contr.treatment",
   ordered = "contr.poly"
))

ui <- fluidPage(
   
   shinyjs::useShinyjs(),
   
   titlePanel("2026 Raptor Analysis"),
   div(id = "sidebarDiv",  
       sidebarPanel(
          fileInput("file", "Upload CSV File", accept = ".csv"),
          checkboxGroupInput("predictors", "Select Predictors for Regression:",
                             choices = c("TimeRel", "RoadSlopeAngle", "RoadHeadingChange", "LanePosition", "CamHeightm", "RoadAttribute"),
                             selected = c("TimeRel", "RoadSlopeAngle", "RoadHeadingChange", "LanePosition", "CamHeightm", "RoadAttribute")),
          hr(),
          strong("Participant Info:"),
          # Wrap the output in a div and hide it initially
          div(id = "fileDetailsContainer", style = "visibility: hidden;",
              verbatimTextOutput("fileDetails")
          )
       )
   ),
   mainPanel(
      # Fixed the pipe syntax bugs by properly nesting within withSpinner
      withSpinner(plotlyOutput("speedPlot"), color = "#0dc5c1"),
      withSpinner(verbatimTextOutput("modelSummary"), color = "#0dc5c1"),
      plotOutput("designMatrixPlot")
   )
)


server <- function(input, output, session) {
   
   # Observe the file input
   observeEvent(input$file, {
      # Show the participant info container when a file is uploaded
      shinyjs::show("fileDetailsContainer")
   })
   
   data <- reactive({
      req(input$file)
      
      read_csv(
         input$file$datapath,
         progress = FALSE,
         show_col_types = FALSE,
         col_types = cols(.default = col_character())
      )
   })
   
   
   processed <- reactive({
      df <- data()
      
      time_num <- suppressWarnings(as.numeric(df$Time))
      speed_num <- suppressWarnings(as.numeric(df$Speed))
      
      df <- df[!is.na(time_num) & !is.na(speed_num), , drop = FALSE]
      
      df %>%
         mutate(
            Time = as.numeric(Time),
            Speed = as.numeric(Speed),
            across(
               c(RoadSlopeAngle, RoadHeadingChange, LanePosition, CamHeightm),
               as.numeric
            ),
            SpeedKPH = Speed * 3.6,
            TimeRel = Time - first(Time),
            SampleInt = row_number(),
            RoadAttribute = factor(
               RoadAttribute,
               levels = c("Normal", "TransSharp", "TransSpike")
            ),
         )
   })
   
   regression <- reactive({
      df <- processed()
      predictors <- input$predictors
      
      validate(
         need(
            "TimeRel" %in% predictors,
            "TimeRel must be included in the model."
         )
      )
      
      continuous_vars <- c(
         "TimeRel",
         "RoadSlopeAngle",
         "RoadHeadingChange",
         "LanePosition",
         "CamHeightm"
      )
      
      categorical_vars <- "RoadAttribute"
      
      selected_continuous <- intersect(
         predictors,
         continuous_vars
      )
      
      selected_categorical <- intersect(
         predictors,
         categorical_vars
      )
      
      required_vars <- unique(c(
         "SpeedKPH",
         "SampleInt",
         selected_continuous,
         selected_categorical
      ))
      
      df <- df %>%
         filter(if_all(
            all_of(required_vars),
            ~ !is.na(.x)
         )) %>%
         mutate(
            RoadAttribute = droplevels(factor(RoadAttribute))
         )
      
      validate(
         need(nrow(df) > 20, "Insufficient complete observations.")
      )
      
      smooth_terms <- paste0(
         "s(",
         selected_continuous,
         ", k = 10)",
         collapse = " + "
      )
      
      formula_parts <- c(
         if (length(selected_continuous) > 0) smooth_terms,
         selected_categorical
      )
      
      validate(
         need(
            length(formula_parts) > 0,
            "No valid predictors selected."
         )
      )
      
      gam_formula <- as.formula(
         paste(
            "SpeedKPH ~",
            paste(formula_parts, collapse = " + ")
         )
      )
      
      fit <- mgcv::gamm(
         formula = gam_formula,
         data = df,
         method = "REML",
         correlation = nlme::corAR1(form = ~ SampleInt)
      )
      
      predicted <- predict(
         fit$gam,
         newdata = df,
         type = "response",
         se.fit = TRUE
      )
      
      df <- df %>%
         mutate(
            ModelSpeed = as.numeric(predicted$fit),
            LowerCI = ModelSpeed -
               1.96 * as.numeric(predicted$se.fit),
            UpperCI = ModelSpeed +
               1.96 * as.numeric(predicted$se.fit)
         )
      
      list(
         df = df,
         gam = fit$gam,
         lme = fit$lme,
         formula = gam_formula,
         selected_continuous = selected_continuous # Passed down for the design matrix plot
      )
   })
   
   
   fileInfo <- reactive({
      req(input$file)
      
      fname <- tools::file_path_sans_ext(input$file$name)
      
      matched <- stringr::str_match(
         fname,
         "^(.*)_(\\d{4}-\\d{2}-\\d{2})_(\\d{2}-\\d{2}-\\d{2})$"
      )
      
      validate(
         need(!is.na(matched[1, 1]), "Filename does not match the expected format.")
      )
      
      list(
         ID = matched[1, 2],
         Date = matched[1, 3],
         Time = gsub("-", ":", matched[1, 4])
      )
   })
   
   output$fileDetails <- renderPrint({
      info <- fileInfo()
      cat(info$ID, "on", info$Date, info$Time, "\n")
   })
   
   output$speedPlot <- renderPlotly({
      reg <- regression()
      if (is.null(reg) || is.null(reg$df)) {
         return(NULL)
      }
      df <- reg$df
      
      p <- plot_ly(df, x = ~TimeRel) %>%
         add_trace(y = ~SpeedKPH, name = "Actual", type = 'scatter', mode = 'lines',
                   line = list(color = '#1f77b4')) %>%
         add_trace(y = ~ModelSpeed, name = "Model", mode = 'lines', line = list(color = '#ff7f0e')) %>%
         add_ribbons(y = ~ModelSpeed, ymin = ~LowerCI, ymax = ~UpperCI, name = "Approx 95% CI",
                     line = list(color = 'rgba(255,127,14,0.2)'),
                     fillcolor = 'rgba(255,127,14,0.2)', hoverinfo = "none") %>%
         layout(title = "Speed Observed vs Model",
                xaxis = list(title = "Time (s)"),
                yaxis = list(title = "Speed (km/h)"),
                hovermode = "x unified")
      
      return(p)
   })
   
   output$modelSummary <- renderPrint({
      reg <- regression()
      
      validate(
         need(!is.null(reg), "TimeRel must be included.")
      )
      
      cat("GAM summary:\n")
      print(summary(reg$gam))
      
      cat("\nAIC for the full GAMM/AR(1) model:",
          AIC(reg$lme), "\n")
   })
   
   # FIXED: Wrapped the floating script block in a proper renderPlot call
   output$designMatrixPlot <- renderPlot({
      reg <- regression()
      df <- reg$df
      vars <- reg$selected_continuous
      
      validate(
         need(length(vars) > 0, "Select at least one continuous predictor to visualize Z-scores.")
      )
      
      Z <- df %>% select(all_of(vars)) %>% scale() %>% as.data.frame()
      Z <- Z + matrix(rep((1:length(vars)) * 5, each = nrow(Z)), ncol = length(vars))
      Z$TimeRel <- df$TimeRel
      Z_long <- tidyr::pivot_longer(Z, -TimeRel, names_to = "Variable", values_to = "Zscore")
      
      ggplot(Z_long, aes(x = TimeRel, y = Zscore, group = Variable)) +
         geom_line(color = "black", linewidth = 1) +
         scale_y_continuous(breaks = (1:length(vars)) * 5, labels = vars) +
         labs(x = "Time (s)", y = "Z-scored Variables") +
         theme_minimal()
   })
}

shinyApp(ui, server)