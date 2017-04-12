shinyUI(fluidPage(
  
  # Application title
  titlePanel("Censored Statistics Simulator"),
  
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        
        selectInput("distribution", "Distribution Type:", c("Log-Normal", "Gamma", "Normal")),
        # This UI limits the reactivity to only this button and changing the distribution
        
      
        # IF LN this is the input values
        conditionalPanel(
          condition="input.distribution == 'Log-Normal'",
          numericInput("sdlog",
                       "SD of Logs:",
                       min=0.1,
                       max=2.5,
                       value=1,
                       step=0.5
          ),
          numericInput("meanlog",
                       "Mean of Logs:",
                       min=0.0,
                       max=10,
                       value=0.5,
                       step=0.5
          )),
        
        # IF GAMMA this is the input values
        conditionalPanel(
          condition="input.distribution == 'Gamma'",
          numericInput("shape",
                       "Shape:",
                       min=0.1,
                       max=4,
                       value=1,
                       step=0.5
          ),
          numericInput("scale",
                       "Scale:",
                       min=1,
                       max=200,
                       value=100,
                       step=10
          )),
        # IF NORMAL this is the input values
        conditionalPanel(
          condition="input.distribution == 'Normal'",
          numericInput("sd",
                       "Std.Dev:",
                       min=0.1,
                       max=20,
                       value=1,
                       step=0.5
          ),
          
          numericInput("mu",
                       "Mu:",
                       min=-10,
                       max=10,
                       value=10,
                       step=1
          )),
        
        numericInput("simNum",
                     "Number of Simulated Samples:",
                     min=1,
                     max=1000,
                     value=100,
                     step=10
        ),
        numericInput("ROPE",
                     "Region of Practical Equivalence:",
                     min=0.01,
                     max=1,
                     value=0.9,
                     step=0.05
        )
      ),
      
      wellPanel(
        
       h4("Number of Samples and Censoring Information"),
       
       wellPanel(
         h5("Sample subgroup 1"),
         numericInput("obs1", "Sample Size:", min=0, max = 200, value=10),
         numericInput("censQ1", "Censoring Quantile:", min=0, max = 1, value=0.25, step=0.01)
         ),
        
       wellPanel(
         h5("Sample subgroup 2"),
         numericInput("obs2", "Sample Size:", min=0, max = 200, value=0),
         numericInput("censQ2", "Censoring Quantile:", min=0, max = 1, value=0.0, step=0.01)
       ),
       
       wellPanel(
         h5("Sample subgroup 3"),
         numericInput("obs3", "Sample Size:", min=0, max = 200, value=0),
         numericInput("censQ3", "Censoring Quantile:", min=0, max = 1, value=0.0, step=0.01)
       ),
       
       actionButton("myValue1", "Select Distribution and RUN Sim")
       
        
        
        )
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      h3("Density Plot of Population Distribution"),
      plotOutput("distGraph"),
      plotOutput("Lambdaview"),
     
     
      plotOutput("UCLgraph"),
      plotOutput("judgementPlot"),
      h3("UCL Statistics"),
      dataTableOutput("UCLcoverage"),
      h3("Lambda Method Recommended UCL"),
      dataTableOutput("UCLrecommend")
#       plotOutput("Bayesgraph"),
#       plotOutput("bayesView")
    )
  )
))

