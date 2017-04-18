#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(MoBspatial)
require(ggplot2) # for plotting
require(gridExtra) # for multiple plots using grid.arrange()
require(raster) # creating a grid polygon
require(broom) # converting spatial polygon to dataframe for plotting with ggplot

# Define UI for slider demo application
  
  fluidPage(theme="simplex.min.css",
            tags$style(type="text/css",
                       "label {font-size: 12px;}"),
            tags$h2("Visualization of biodiversity pattern"),
    fluidRow(
    column(4, 
           # Simple integer interval
           sliderInput("S", "Species Richness",
                       min=10, max=500, value=50, step=10),
           
           sliderInput("N", "Number of individuals",
                       min=500, max=10000, value=1000, step=100),
           
           sliderInput("cv.abund", "CV(SAD)",
                       min = 0.5, max = 3, value = 1, step= 0.5),
           
           sliderInput("spatagg", "Spatial Agregation",
                       min = 0.01, max = 1, value = 0.01, step= 0.01),
           
           sliderInput("Resolution", "Resolution",
                       min = 2, max = 6, value = 2),

           sliderInput("cell", "Cell",
                       min = 1, max = 36, value = 1)
           
    ),
    column(8, plotOutput("InteractivePlot", height="700px"))
  )  
)