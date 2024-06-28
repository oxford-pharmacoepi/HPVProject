# Libraries
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(shinydashboard)
library(DT)
library(gt)
library(ggforestplot)
library(ggplot2)
library(plotly)
library(here)
library(dplyr)
library(stringr)
library(tidyr)
library(fresh)
library(rclipboard)
library(bit)
library(bit64)
library(rsconnect)
library(packrat)
library(visOmopResults)
library(DiagrammeR)

# Variables
dataFolder <- "data"

# Functions
source(here("functions.R"))
print("Functions added")
# Data
source(here("prepare_data.R"))
print("Data prepared")

# App
# theme -----------------------
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#002147" 
  ),
  adminlte_sidebar(
    # width = "400px",
    dark_bg = "#7eb3eb",
    dark_hover_bg = "#226584",
    dark_color = "black"
  ), 
  adminlte_global(
    content_bg = "#eaebea" 
  ),
  adminlte_vars(
    border_color = "#002147",
    active_link_hover_bg = "#FFF",
    active_link_hover_color = "#002147",
    active_link_hover_border_color = "#002147",
    link_hover_border_color = "#002147"
  )
)

shinyApp(ui = ui, server = server)
