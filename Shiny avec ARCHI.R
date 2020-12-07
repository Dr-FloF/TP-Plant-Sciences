

library(shiny)


ui <- fluidPage(
  
    tags$h1("Ideotype assessment for optimise P use"),
    tags$hr(),
    
    fluidRow(
    column(6,
    tags$h3("Plant caracteristics:"),
                
                checkboxGroupInput(inputId = "Dicot",label ="Phylogenetic group",
                                   choiceNames=list("Dicotyledoneous","Monocotyledoneous"),
                                   choiceValues=list(TRUE,FALSE)),
                
                sliderInput(inputId="dmin",label="Minimal diameter",
                            value = 0.2, min = 0.1, max= 0.5),
                
                sliderInput(inputId="dmax",label="Maximal diameter",
                            value = 1, min = 0.5, max= 2.5),
                
                sliderInput(inputId="ibd",label="Interbranching distance",
                            value = 5, min = 3, max= 10),
                
                sliderInput(inputId="igravi",label="Gravitropism intensity",
                            value = 0.005, min = 0, max= 0.01),

                
                checkboxGroupInput(inputId = "Enzyme",label ="Phosphatase production",
                                   choiceNames=list("Yes","No"),
                                   choiceValues=list(TRUE,FALSE)),),
    column(6,
    tags$h3("Sol caracteristics:"),
    
                numericInput(inputId="Psurf",
                          label = "Top soil P",
                          value = "100"),
                
                numericInput(inputId="Pdeep",
                          label = "Deep soil P",
                          value = "10"),
                
                numericInput(inputId="Pavail",
                          label = "P availability",
                          value = "0.15"),),
    
    tags$h4("Ready to run the model?"),
    actionButton(inputId = "go",label="Run"),
    
    tags$h4("Parameters of the modelisation:"),
    verbatimTextOutput("para")
    
    ),
    
    
    fluidRow(
                
      tags$h3("Results of the modelisation:"),            
                
                verbatimTextOutput("stats"),
      
      tags$h3("Root system projection:"),   
                plotOutput("roots"),
    )
                
)


server <- function(input, output) {
  
  
 
    #############################################
    # Sol homogéne et non limitant pour commencer 
    
    SOLUNC <- data.frame(Croiss=rep(1.0,times=100),Ramif=rep(1.0,times=100),ICMeca=rep(0.02,times=100),DirectionConstraint=rep(0,times=100)) #cr???ation d'un tableau avec pour longueur 100, indice de ramification = 1 ie homog???ne et parall???le, direction de contrainte m???canique: le sol joue nn seulement sur la vitesse de croissance et sur sa direction (contrainte angle) 
    write.table(x=SOLUNC,file="./sol.txt",sep = "\t",row.names = FALSE,col.names = TRUE)
    
    
    ###################################
    # Allocation Biomasse non imitante
    
    biom00 <- rep(10,times=300) 
    
    write(x=biom00,file="./biomrac.txt",sep="\n")  
    
    result <- eventReactive(input$go,{    
    ###################################
    if (input$Dicot==TRUE)
    {
      PAR<- data.frame(
        duration=40,
        emissratesem=0.9,
        propDiamSem=1,
        nbmaxsem=1,
        ageadvemis=100,
        distmaxadv=80,
        emissrateadv=0.1,
        propDiamadv=1.0,
        nbmaxadv=5,
        dmin=input$dmin,
        dmax=input$dmax,
        elongvsdiam=22,
        tgravi=1,
        igravi=input$igravi,
        durdp=5.0,
        ibd=input$ibd,
        dldm=0.30,
        vard=0.20,
        tmd=0.10,
        growthdur=600,
        lifeexp=2500,
        radgrowth=0.10 )
    } else {
      PAR<- data.frame(
        duration=40,
        emissratesem=0.9,
        propDiamSem=0.65,
        nbmaxsem=3,
        ageadvemis=3,
        distmaxadv=4,
        emissrateadv=5,
        propDiamadv=1.0,
        nbmaxadv=35,
        dmin=input$dmin,
        dmax=input$dmax,
        elongvsdiam=22,
        tgravi=1,
        igravi=input$igravi,
        durdp=5.0,
        ibd=input$ibd,
        dldm=0.10,
        vard=0.20,
        tmd=0.10,
        growthdur=250,
        lifeexp=2500,
        radgrowth=0 )
    }
    
    
    
    # Ecriture du fichier des paramétres
    write.table(x=PAR,file="paramarch9.txt",sep = "\n",row.names = FALSE,col.names = FALSE)
    
    # Simulation
    system("archisimp9.exe 0.0 0.0 5.0") # position de la semence (x, y, z) enterre un peu ??? 10 mm
    
    # Lecture du fichier de la structure
    SR <- read.table(file="seg.txt",header=T)
    
    # Calcul des paramétres racines 
    SR$depth <- cut(0.5*(SR$Z1+SR$Z2),breaks=seq(from=0,to=1000,by=50)) # on fait des classes
    PROF <- data.frame(depth=seq(from=0,to=950,by=50))
    PROF$length <- tapply(sqrt(((SR$X1-SR$X2)^2)+((SR$Y1-SR$Y2)^2)+((SR$Z1-SR$Z2)^2)),SR$depth,sum) 
    PROF$surface <- tapply(pi*SR$Diam*sqrt(((SR$X1-SR$X2)^2)+((SR$Y1-SR$Y2)^2)+((SR$Z1-SR$Z2)^2)),SR$depth,sum) # mm???
    PROF$Prop.length <- (PROF$length/sum(PROF$length,na.rm=TRUE))*100
    PROF$Prop.surface <- (PROF$surface/sum(PROF$surface,na.rm=TRUE))*100
    
    # Module sol
    
    SOL <- data.frame(depth=seq(from=0,to=950,by=50))
    
    Vol<-NULL
    
    X<-cbind(tapply(SR$X1,SR$depth,min),tapply(SR$X2,SR$depth,min),tapply(SR$X1,SR$depth,max),tapply(SR$X2,SR$depth,max))
    Y<-cbind(tapply(SR$Y1,SR$depth,min),tapply(SR$Y2,SR$depth,min),tapply(SR$Y1,SR$depth,max),tapply(SR$Y2,SR$depth,max))
    for (i in 1:nrow(X))
    {
      Vol[i]<-((max(X[i,3:4])-min(X[i,1:2]))*(max(Y[i,3:4])-min(Y[i,1:2]))*50)/1000
    }
    SOL$Vol.Explo<-Vol # volume exploré en cm3
    SOL$Masse.Explo<-(SOL$Vol.Explo*1.3)/1000
    
    SOL$C.P.Total<-ifelse (SOL$depth<200,input$Psurf,input$Pdeep) #  mg.kg
    
    SOL$P.Total.Explo<-SOL$C.P.Total*SOL$Masse.Explo  # mg accessible 
    
    if (input$Enzyme==TRUE)
    {SOL$P.Available<-ifelse((SOL$P.Total.Explo*input$Pavail*1.3)<SOL$P.Total.Explo,SOL$P.Total.Explo*input$Pavail*1.3,SOL$P.Total.Explo)
    }else{
      SOL$P.Available<-SOL$P.Total.Explo*input$Pavail}
    
    SOL$P.Export<-ifelse(SOL$P.Available*(PROF$surface/100)*0.0005<SOL$P.Available,SOL$P.Available*(PROF$surface/100)*0.0005,SOL$P.Available) # fraction de P accessible autour de la racine
    
    
    ### Summary
    
    Summary<-data.frame(matrix(nrow=1,ncol=5))
    names(Summary)<-c("Length (m)","Surface (cm²)","Shallow Lenght (%) ","Shallow Surface (%) ","P Export (mg)")
    Summary[1,1]<-sum(PROF$length,na.rm=TRUE)/100
    Summary[1,2]<-sum(PROF$surface,na.rm=TRUE)
    Summary[1,3]<-sum(PROF$Prop.length[1:4],na.rm=TRUE)
    Summary[1,4]<-sum(PROF$Prop.surface[1:4],na.rm=TRUE)
    Summary[1,5]<-sum(SOL$P.Export,na.rm=TRUE)

    para<-data.frame(matrix(nrow=1,ncol=6))
    names( para)<-c("dmin (mm)","dmax (mm)","igravi","ibd","Dicot","Enzyme")
    para[1,1]<-input$dmin
    para[1,2]<-input$dmax
    para$igravi<-input$igravi
    para$ibd<-input$ibd
    para$Dicot<-ifelse(input$Dicot==TRUE,"Yes","No")
    para$Enzyme<-ifelse(input$Enzyme==TRUE,"Yes","No")
    
    ### Sorties
    
    result<-list(summary=Summary,SR=SR,para=para)
    return(result)


    
  })
  
  output$para<-renderPrint({
      result()[[3]]
  })
    
    
  output$stats<-renderPrint({
    result()[[1]]
  })
  
  output$roots<-renderPlot({
    
    DessineSRXZ(result()[[2]][result()[[2]]$Diam>0.4,])
    abline(h=-200,lty=3)
    
  })
  
  
}

shinyApp(ui = ui, server = server)
