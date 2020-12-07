library(shiny)


ui <- fluidPage("Sensitivity",
                
      checkboxGroupInput(inputId = "Morpho",label ="Morphological parameter",
                        choiceNames=list("Minimal diameter","Maximal diameter",
                                         "Interbranching distance","Gravitromism intensity"),
                        choiceValues=list("dmin","dmax","ibd","igravi")),
      
      checkboxGroupInput(inputId = "Dicot",label ="Phylogenetic group",
                         choiceNames=list("Dicotyledoneous","Monocotyledoneous"),
                         choiceValues=list(TRUE,FALSE)),
      
      checkboxGroupInput(inputId = "Enzyme",label ="Phosphatase production",
                         choiceNames=list("Yes","No"),
                         choiceValues=list(TRUE,FALSE)),
      
      tags$h3("Sol caracteristics:"),
      
      numericInput(inputId="Psurf",
                   label = "Top soil P",
                   value = "100"),
      
      numericInput(inputId="Pdeep",
                   label = "Deep soil P",
                   value = "10"),
      
      numericInput(inputId="Pavail",
                   label = "P availability",
                   value = "0.15"),
      
      tags$h4("Ready to run the model?"),
      actionButton(inputId = "go",label="Run"),
      
      
      tags$h3("Results of the sensitivity analysis:"),
      
      tags$h4("coefficient of variation"),
      verbatimTextOutput("Coef"),
      tags$h4("Simulations"),
      verbatimTextOutput("simul"),
      
      plotOutput("longueur"),
      plotOutput("phosphore"),
      

)


server <- function(input, output) {

  
result <- eventReactive(input$go,{   
  
  Values<-NULL
  sum1<-NULL
  
  
## dmin
  if(input$Morpho=="dmin"){
  Values<-seq(from=0.05,to=0.5,by=0.05)
  
    for (dmin in Values)
    {
      simul<-Archi(
        #### Crop characteristics
        dmin=dmin,
        dmax=1,
        igravi=0.05,
        ibd=5,
        Dicot=input$Dicot,
        Enzyme=input$Enzyme,
        
        ### Soil characteristics
        P.Surf = input$Psurf,
        P.Depth = input$Pdeep,
        P.Availability = input$Pavail, 
        
        ### Type of analysis 
        Graphic=FALSE)

      
      sum1<-rbind(sum1,simul$Summary)
    }
  }
  
## dmax
  if(input$Morpho=="dmax"){
    Values<-seq(from=0.5,to=2.5,by=0.25)
   
    for (dmax in Values)
    {
      simul<-Archi(
        #### Crop characteristics
        dmin=0.25,
        dmax=dmax,
        igravi=0.05,
        ibd=5,
        Dicot=input$Dicot,
        Enzyme=input$Enzyme,
        
        ### Soil characteristics
        P.Surf = input$Psurf,
        P.Depth = input$Pdeep,
        P.Availability = input$Pavail, 
        
        ### Type of analysis 
        Graphic=FALSE)
      
      sum1<-rbind(sum1,simul$Summary)
    }
  }

## ibd
  if(input$Morpho=="ibd"){
    Values<-seq(from=0.0,to=0.01,by=0.001)

    for (ibd in Values)
    {
      simul<-Archi(
        #### Crop characteristics
        dmin=0.25,
        dmax=3,
        igravi=0.05,
        ibd=ibd,
        Dicot=input$Dicot,
        Enzyme=input$Enzyme,
        
        ### Soil characteristics
        P.Surf = input$Psurf,
        P.Depth = input$Pdeep,
        P.Availability = input$Pavail, 
        
        ### Type of analysis 
        Graphic=FALSE)
      
      sum1<-rbind(sum1,simul$Summary)
    }
  }
  
## igravi
  if(input$Morpho=="igravi"){
    Values<-seq(from=1,to=10,by=1)

    for (igravi in Values)
    {
      simul<-Archi(
        #### Crop characteristics
        dmin=0.25,
        dmax=3,
        igravi=igravi,
        ibd=5,
        Dicot=input$Dicot,
        Enzyme=input$Enzyme,
        
        ### Soil characteristics
        P.Surf = input$Psurf,
        P.Depth = input$Pdeep,
        P.Availability = input$Pavail, 
        
        ### Type of analysis 
        Graphic=FALSE)
      
      sum1<-rbind(sum1,simul$Summary)
    }
  }

  for (k in 1:9)
  {
    sum1["cv",k]<-round(sd(sum1[,k],na.rm=TRUE)/mean(sum1[,k],na.rm=TRUE),2)
  }
  sum1[,1]<-sum1[,1]/100
  
  return(Summary.Sensitivity=sum1)  
  
})

lab1 <- eventReactive(input$go,{ 
  if(input$Morpho=="dmin"){x<-"Minimal diameter"}
  if(input$Morpho=="dmax"){x<-"Maximal diameter"}
  if(input$Morpho=="ibd"){x<-"Interbranching distance"}
  if(input$Morpho=="igravi"){x<-"Gravitropism intensity"}
  return(x=x)
})


output$simul<-renderPrint({
  result()[1:10,]
})

output$Coef<-renderPrint({
  result()[11,1:9]
})

output$longueur<-renderPlot({
  plot(result()[,1]~result()[,input$Morpho],ylab="Total root length (m)",xlab=lab1())
})

output$phosphore<-renderPlot({
  plot(result()[,5]~result()[,input$Morpho],ylab="P export (mg)",xlab=lab1())
})


}


shinyApp(ui = ui, server = server)
