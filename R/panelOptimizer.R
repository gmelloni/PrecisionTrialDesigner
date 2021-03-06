.busyIndicator <- function(text = "Processing..."
        , image=NULL
        , wait=1000) {
  tagList(
    singleton(tags$head(
      tags$link(rel = "stylesheet"
        , type = "text/css"
      )))
    ,div(class = "mybusyindicator",p(text),img(src=image))
    ,tags$script(sprintf(
    " setInterval(function(){
       if ($('html').hasClass('shiny-busy')) {
        setTimeout(function() {
          if ($('html').hasClass('shiny-busy')) {
            $('div.mybusyindicator').show()
          }
        }, %d)          
      } else {
        $('div.mybusyindicator').hide()
      }
    },100)
    ",wait)
    )
  ) 
}
setGeneric('panelOptimizer', function(object) standardGeneric('panelOptimizer'))
setMethod('panelOptimizer', 'CancerPanel', function(object){
  ktcga_types <- c(
    "In_Frame_Del", "In_Frame_Ins","Missense_Mutation"
    ,"Frame_Shift_Del", "Frame_Shift_Ins"
    , "Nonsense_Mutation",  "Splice_Site"
    , "Translation_Start_Site", "Nonstop_Mutation", "Silent","3'UTR", "3'Flank"
    , "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region"
  )
  knotTransc <- c("3'UTR", "3'Flank", "5'UTR", "5'Flank"
                  , "IGR", "Intron", "RNA", "Targeted_Region")
  knonsynonymous <- setdiff(ktcga_types , c("Silent" , knotTransc))
  kmiss_type <- ktcga_types[c(1,2,3)]
  ktrunc_type <- ktcga_types[c(4,5,6,7,8,9)]
  panel <- cpArguments(object)$panel
  # Is it possible to select all SNV genes?
  genes <- panel[ panel$alteration=="SNV" &
                  panel$exact_alteration %in% c("mutation_type" , "")
                  , "gene_symbol"] %>% unique
  if(length(genes)==0){
    stop("There are no genes in the panel reqeusted for SNV")
  }
  tumor_type <- cpData(object)$mutations$Samples
  tumor_type <-  names(tumor_type)[!vapply(tumor_type , is.null , logical(1))]
  if(length(tumor_type)==0){
    stop("No tumor types available for mutations in this Cancer Panel object")
  }
  mydata <- cpData(object)$mutations$data
  if(is.null(mydata)|nrow(mydata)==0){
    stop("No mutation data in this CancerPanel object")
  }
  repos <- data.frame(Entrez=mydata$entrez_gene_id
                      ,Gene_Symbol=mydata$gene_symbol
                      ,Amino_Acid_Letter=substr(mydata$amino_acid_change,1,1)
                      ,Amino_Acid_Position=mydata$amino_position
                      ,Amino_Acid_Change=mydata$amino_acid_change
                      ,Mutation_Type=mydata$mutation_type
                      ,Sample=mydata$case_id
                      ,Tumor_Type=mydata$tumor_type
                      , stringsAsFactors=FALSE)
  runApp(list(
    ui = fluidPage(
      sidebarLayout(
        sidebarPanel(
          selectInput("genes", "Choose Your Gene:"
                      , choices=sort(genes)
                      , multiple=FALSE
                      , selected=sort(genes)[1]
                      )
          , selectInput("tumor_type" , "Choose Tumor Types:"
                      , choices=tumor_type
                      , multiple=TRUE
                      , selected=tumor_type[1]
                      )
          , numericInput("bandwidth"
                      , label = "Select a bandwidth:"
                      , value = 0
                      , min=0
                      )
          , helpText("Automatic bandwidth overwrites the slider above")
          , radioButtons("autobw" , "Use Silverman's bandwidth"
                        , c("No" =  "No" , "Yes" = "Yes")
                        )
          , radioButtons("metric" , "Use P-value or Q-value" 
                        , c("Q-Value" = "qvalue" , "P-Value"="pvalue") 
                        )
          , helpText("Select a new gene and click Run")
          , actionButton("goButton", "Run")
          , actionButton("stopandsave" , "Close and save")
          , width=3)
        ,mainPanel(
          tabsetPanel(type="tabs"
            , tabPanel("LowMACA Regions" 
                , bsButton(inputId="storelowmaca"
                , label= "Store LowMACA" 
                , style="warning"
                , size="small"
                , disabled=FALSE)
            , helpText(paste("Click to add the regions suggested"
                  , "by LowMACA or select the rows you prefer and then click"))
            , helpText("The regions will be stored in the Optimize Panel tab")
            , plotOutput("lmPlot" , height=400)
            , hr()
            , verbatimTextOutput("NoSignificanceText")
            , hr()
            , DT::dataTableOutput("lowmacaregions"))
            , tabPanel("LowMACA Stats"  
                , helpText("Check significant positions identified by LowMACA")
                , DT::dataTableOutput("lfm"))
            , tabPanel("Manual Selection"
                , bsButton(inputId="storesingle"
                    , label= "Store Single Mutation" 
                    , style="warning"
                    , size="small"
                    , disabled=FALSE)
                , bsButton(inputId="storeregion"
                    , label= "Store Region" 
                    , style="warning"
                    , size="small"
                    , disabled=FALSE)
                , helpText(paste("Click on a dot and Store a single mutation."
                              ,"Drag a region and then click on Store Region"))
                , plotOutput("barplot" 
                    , click = "plot_click"
                    , brush = "plot_brush"
                    , height = 300)
                , plotOutput("barplotdomains" , height=30)
                , hr()
                , verbatimTextOutput("textMessage")
                , hr()
                , tableOutput("infoclick")
                , tableOutput("infobrush"))
            , tabPanel("Optimize Panel" 
                , helpText(paste("All your selected regions are here." 
                    , "To remove them, select the rows and then click Discard"))
                , bsButton(inputId="discard" 
                                , label= "Discard" 
                                , style="danger"
                                , size="small"
                                , disabled=FALSE)
                       , DT::dataTableOutput("panel"))
          )
        , width=9)
      )
      ,.busyIndicator(wait=1000 
                      , image="http://i.giphy.com/l3V0EQrPMh1nnfbFe.gif")
    )
    ,server = {
      function(input, output , session) {
      # --------------------------------------
      # General Parameter for the server
      # --------------------------------------
      opar <- par("mfrow" , "mar")
      on.exit(par(opar))
      options(DT.options = list(pageLength = 20))
      options(xtable.include.rownames=FALSE)
      # --------------------------------------
      # LowMACA analysis (tab 1 and 2)
      # --------------------------------------

        # Create a LowMACA object that depends on tumor types and gene
        # lmObj is reactive to the Run button 
        # (if you change gene you must press it again)
      lmObj <- eventReactive(input$goButton , {
          repos_red <- repos[ repos$Tumor_Type %in% input$tumor_type &
                              repos$Gene_Symbol==input$genes, ]
          if(nrow(repos_red)==0){
            return(NULL)
          }
          lmObj <- tryCatch(.newLowMACAsimple(genes=input$genes) 
                            , error=function(e) NULL)
          if(is.null(lmObj)){
            return(NULL)
          }
          suppressWarnings( lmObj$arguments$params$mutation_type <- "all" )
          lmObj <- .setupSimple(lmObj , repos=repos_red)
          return(lmObj)
        })
        # lmObjEntropy takes lmObj and apply .entropySimple with the selected bw
        # changing bandwidth will automatically 
        # recalculate this step if run button is pressed
      lmObjEntropy <- eventReactive(input$goButton , {
        if(is.null(lmObj()))
          return(NULL)
        if(input$autobw=="Yes"){
          lmObj <- .entropySimple(lmObj() , bw="auto")  
        } else {
          lmObj <- .entropySimple(lmObj() , bw=as.numeric(input$bandwidth))
        }
        })
      # Visualize lmPlot (tab 1)
      output$lmPlot <- renderPlot({
        if(is.null(lmObjEntropy()))
          return(.emptyPlotter("No data for these tumor types"))
        tmp <- lmObjEntropy()
        .lmPlotSimple(tmp)
      })
      myoutputlist <- reactiveValues(regions=NULL
                                    , lfm=NULL
                                    , NoSignificanceText=NULL)
        # Report the significant mutations (lfm function of LowMACA) (tab 2) 
      observe({
        if(is.null(lmObjEntropy())){
          myoutputlist$lfm <- NULL
          # return(NULL)
        } else {
          tmp <- lmObjEntropy()
          out <- .lfmSimple(tmp , metric=input$metric)
          if(nrow(out)==0){
            myoutputlist$lfm <- out
            # return(out)
          } else {
            # Aesthetics improvements
            out <- out[ , c("Gene_Symbol" , "Amino_Acid_Position" 
                  , "Amino_Acid_Change" , "Tumor_Type" , "metric" , "Sample")]
            colnames(out)[colnames(out)=="metric"] <- input$metric
            agg1 <- aggregate(Sample ~ Gene_Symbol + 
                                Amino_Acid_Position + 
                                Tumor_Type, out , FUN=length)
            agg1$Tumor_Type_and_Freq <- paste(agg1$Tumor_Type 
                                              , agg1$Sample , sep=":")
            agg1 <- aggregate(Tumor_Type_and_Freq ~ Gene_Symbol + 
                                Amino_Acid_Position , agg1 
                    , FUN=function(x) paste(unique(x) , collapse="|"))
            finalOut <- merge(agg1 , unique(out[ , c("Gene_Symbol" 
                            , "Amino_Acid_Position" , input$metric)]))
            finalOut$Amino_Acid_Position <- 
              as.integer(finalOut$Amino_Acid_Position)
            finalOut <- finalOut[ order(finalOut$Amino_Acid_Position) , ]
            # Update element lfm of the output
            myoutputlist$lfm <- finalOut
          }
        }
      })
      output$lfm <- DT::renderDataTable({
        myoutputlist$lfm
        } ,rownames=FALSE
      )
       # Report LowMACA results:
          # 1) the significant regions
          # 2) the number of bp they span
          # 3) the percentage of mutations covered
      observe({
        if(is.null(lmObjEntropy())){
          # return(NULL)
          myoutputlist$regions <- NULL
        } else {
          mydf <- lmObjEntropy()$alignment$df
          aminos <- seq_len(nrow(mydf))
          signAminos <- aminos[which(mydf[[input$metric]]<=0.05)]
          if(length(signAminos)==0){
            # return(NULL)
            myoutputlist$regions <- NULL
          } else {
            signAminosGroups <- .splitter(signAminos)
            out <- data.frame(
            Gene_Symbol=rep(lmObjEntropy()$arguments$genes 
                                , length(signAminosGroups))
                , Span=vapply(signAminosGroups , function(x) {
                      if(length(x)==1){
                        return(paste(as.character(x) , as.character(x) 
                                     , sep="-"))
                      } else {
                        return(paste(x[1] , x[length(x)] , sep="-"))
                      }
                  } , character(1))
                , Region_Width_In_bp=vapply(signAminosGroups , function(x) {
                      if(length(x)==1){
                        return(3)
                      } else {
                        span <- x[length(x)] - x[1] + 1
                        return(span*3)
                      }
                  } , numeric(1))
                , stringsAsFactors=FALSE)
            # Calculate the percentage of aminoacids covered by each region
            Span_Percentage=vapply(signAminosGroups , function(x) {
                      if(length(x)==1){
                        (1/length(aminos))*100
                      } else {
                        span <- x[length(x)] - x[1] + 1
                        (span/length(aminos))*100
                      }
                  } , numeric(1))
            out$Span_Percentage <- paste0(format(round(Span_Percentage 
                                                       , 2), nsmall=2) , "%")
            # Calculate the percentage of mutations caught by the region
            mymuts <- lmObjEntropy()$mutations$data
            coveredMutations <- vapply(signAminosGroups , function(x) {
                regionMuts <- nrow(mymuts[ mymuts$Amino_Acid_Position %in% x, ])
                perc <- 100*(regionMuts/nrow(mymuts))
                return(perc)
                        } , numeric(1))
            out$Percentage_Of_Covered_Mutations <- 
              paste0(format(round(coveredMutations , 2), nsmall=2) , "%")
            totalSummary <- paste(
              paste("Gene Symbol:" , lmObjEntropy()$arguments$genes
                    , "-", length(aminos) , "aminoacids")
              , paste("Identified LowMACA regions:" , nrow(out))
              , paste("Total bp:" , sum(out$Region_Width_In_bp))
              , paste("Percentage of the protein:" 
                      , paste0(format(round(sum(Span_Percentage) , 2)
                                      , nsmall=2) , "%"))
              , paste("Percentage of total mutations:" 
                      , paste0(format(round(sum(coveredMutations) , 2)
                                      , nsmall=2) , "%"))
              , sep="\n")
            myoutputlist$NoSignificanceText <- totalSummary
            # Update element regions in myoutputlist
            myoutputlist$regions <- out
          }
        }
      })
        # Message for no significant mutations from LowMACA
      observe({
        if(is.null(lmObj())){
          myoutputlist$NoSignificanceText <- NULL
        } else {
          if(is.null(myoutputlist$regions)){
            myoutputlist$NoSignificanceText <- 
              "No Significant Regions or Mutations identified by LowMACA"
          }
        }
      })
      output$NoSignificanceText <- renderText(myoutputlist$NoSignificanceText)
        # lowmaca regions table visualization (tab 1)
      output$lowmacaregions <- DT::renderDataTable({
        myoutputlist$regions
        } , server=FALSE
          ,rownames=FALSE)
        # Store LowMACA regions in compound regions
        # When clicking Store LowMACA, an event is triggered:
          # If no rows are selected, take all the regions suggested by LowMACA
          # If the first row is selected, take the regions suggested by LowMACA
          # If specific regions are selected, take those
      observeEvent(input$storelowmaca , {
        if(is.null(input$lowmacaregions_rows_selected)){
          lowmacaadd <- myoutputlist$regions[ -1 , , drop=FALSE]
        } else {
          # print(input$lowmacaregions_rows_selected)
          if(1 %in% as.numeric(input$lowmacaregions_rows_selected)){
            lowmacaadd <- myoutputlist$regions[ -1 , , drop=FALSE]
          } else {
            lowmacaadd <- myoutputlist$regions[ 
              as.numeric(input$lowmacaregions_rows_selected) , , drop=FALSE]
          }
        }
        myregion$compoundRegion <- rbind(myregion$compoundRegion 
                                         , lowmacaadd) %>% unique
      })

      # --------------------------------------
      # Manual selection analysis (tab 3)
      # --------------------------------------
        # Collect mutations by position from LowMACA object
        # This reactive will be used by all the output of tab 3
      mut_aligned_df <- reactive({
        if(is.null(lmObj()))
          return(NULL)
        mut_aligned <- lmObj()$mutations$aligned
        # New Version with trunc and missense
        mut_aligned_data <- lmObj()$mutations$data
        mut_aligned_tab <- table(mut_aligned_data$Amino_Acid_Position 
                                 , mut_aligned_data$Mutation_Type) %>%
          melt %>%
          setNames(. , c("Position" , "Mutation_Type" , "Mutations")) %>%
          .[ .$Mutations!=0 , ]
        mut_aligned_tab$Color <- ifelse(mut_aligned_tab$Mutation_Type %in% 
                                          ktrunc_type , "darkred" , "limegreen")
        mut_aligned_tab_agg1 <- aggregate(Mutations ~ Position 
                                          , mut_aligned_tab , FUN=sum)
        mut_aligned_tab_agg2 <- aggregate(Color ~ Position 
                                          , mut_aligned_tab , FUN=function(x) {
          out <- paste(sort(unique(x)) , collapse="|")
          if(out=="darkred|limegreen")
            return("purple")
          else
            return(out)})
        mut_aligned_df <- merge(mut_aligned_tab_agg1 , mut_aligned_tab_agg2)
        list(mut_aligned , mut_aligned_df)
      })
        # Visualize clickable barplot (tab 3)
      output$barplot <- renderPlot({
        if(is.null(lmObj)){
          return(.emptyPlotter("No Mutations or No protein available"))
        }
        mut_aligned <- mut_aligned_df()[[1]]
        mut_aligned_df <- mut_aligned_df()[[2]]
        # Draw a custom axis
        myaxis <- seq(0 , ncol(mut_aligned) , 10)
        myaxis[1] <- 1
        if(ncol(mut_aligned)!=myaxis[length(myaxis)]) {
          myaxis[length(myaxis)] <- ncol(mut_aligned)
        }
        # Impose custom limits
        myylim <- c(0, max(colSums(mut_aligned)) + 1)
        myxlim <- c(0 , ncol(mut_aligned))
        par(mar=c(1.5,4.1,4.1,2.1))
        plot(mut_aligned_df$Position
          , mut_aligned_df$Mutations
          , cex=1.5
          , main=paste0("Lolliplot of " , rownames(mut_aligned) , " - " 
                        , ncol(mut_aligned) , " aminoacids"
                      , "\nclick on the dots or drag to select a region")
          , xaxs = "i"
          , xaxt="n"
          , yaxt="n"
          , frame.plot=FALSE
          , bg=mut_aligned_df$Color
          , pch=21
          , col="gray"
          , xlab=""
          , ylab="Mutations"
          , xlim=myxlim
          , xaxp=c(myxlim[1] , myxlim[2] , 10)
          , ylim=myylim
          )
        axis(side = 1 , at = round(seq(1 , myxlim[2] , length = 20)) 
             , labels = as.character(round(seq(1 , myxlim[2] , length = 20)))
             , mgp = c(0, 0.5, 0)
             , font.axis=2,font.lab=2,cex.lab=0.8,cex.axis=0.8)
        ymarks <- if(myylim[2]<6){
                    seq(0 , myylim[2] , 1) 
                  } else {
                    seq(0 , myylim[2] , length = 5)
                  }
        axis(side = 2 , at = round(ymarks)
             , labels = as.character(round(ymarks))
             #, cex = 0.1 
             #, font=windowsFont("Comic Sans MS")
             , mgp = c(0, 0.5, 0)
             , font.axis=2,font.lab=2,cex.lab=0.8,cex.axis=0.8)
        legend("topleft" 
               , legend = c("Truncation" , "Missense" , "Both") 
               , col = c("darkred" , "limegreen" , "purple")
               , pch=19
               , bty="n")
        # The distance between the bottom of y axis and x axis is 4% of yaxis
        # We have to elong the segements so that they touch the x axis
        for(i in mut_aligned_df$Position){
            segments( i 
                      , 0 - 0.04*diff(myylim)
                      , i 
                      , mut_aligned_df[mut_aligned_df$Position==i 
                          , "Mutations"] - 0.025*diff(myylim)
                      , col="gray" , lwd=2)
        }
        par <- par(opar)
      })
        # Visualize annotated domains (tab 3)
      output$barplotdomains <- renderPlot({
        if(is.null(lmObj())){
          return(NULL)
        }
        mut_aligned <- mut_aligned_df()[[1]]
        mut_aligned_df <- mut_aligned_df()[[2]]
        myxlim <- c(0 , ncol(mut_aligned))
        myPfam <- LowMACAAnnotation::getMyPfam()
        gene <- rownames(mut_aligned)
        domains <- myPfam[myPfam$Gene_Symbol==gene
          , c("Envelope_Start" , "Envelope_End" , "Pfam_Name")]
        ## create empty plot
        # No margins on top and bottom, but keep 
        # the left and right margins as the plot above
        par(mar=c(0,4.1,0,2.1))
        plot(c(0, 1), c(0, 1)
            , ann = FALSE, bty = 'n', type = 'n'
            , xaxt = 'n', yaxt = 'n' , xaxs = "i"
            , xlim=myxlim , ylim=c(0,0.05) 
            , xaxp=c(myxlim[1] , myxlim[2] , 10))
        # Draw a gray rectangle to represent 
        # the protein, the domains will be plotted over
        rect(xleft = myxlim[1] , xright = myxlim[2] 
             , ytop = 0.0332 , ybottom = 0.01666
             , col="gray" , border=NA)
        ## plot domains
        if(nrow(domains)>0) {
          # Assign a different color for each domain name
          myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral"))
                                        , space="Lab")
          domains$Colors <- factor(domains$Pfam_Name) %>% 
            as.numeric %>% myPalette(max(.))[.]
          for(i in seq_len(nrow(domains))) {
            xleft=as.numeric(domains[i , "Envelope_Start"])
            xright=as.numeric(domains[i , "Envelope_End"])
            ytop=0.05
            ybottom=0
            col=domains[i , "Colors"]
            characters <- nchar(domains[i , "Pfam_Name"])
            rect(xleft=xleft , xright=xright 
              , ytop=ytop , ybottom=ybottom 
              , col=col , border=NA)
            text(x=xleft , y=0.025 , as.character(xleft) 
                 , font=2 , col="red" , cex=1 , srt=90)
            text(x=xright , y=0.025 , as.character(xright) 
                 , font=2 , col="red" , cex=1 , srt=90)
            if(characters<=3){
              text(x=(xright+xleft)/2 , y=0.025 
                , domains[i , "Pfam_Name"] 
                , font=2)
            } else {
              if((xright-xleft)<=100){
                text(x=(xright+xleft)/2 , y=0.025 
                  , domains[i , "Pfam_Name"] 
                  , font=2 , cex=0.8)
              } else {
                text(x=(xright+xleft)/2 , y=0.025 
                  , domains[i , "Pfam_Name"] 
                  , font=2)
              }
            }
          }
        } else {
          text(x=sum(myxlim)/2, y=0.025
            , 'no pfam domains', cex=1.5)
        }
        par <- par(opar)
      })
      # Reactive values for manual selection
      # 1) selectRegion: brush points selected
      # 2) selectSingle: click single mutations selction
      # 2) compoundRegion: incremental brush point + LowMACA regions for tab 4
      # 3) textMessage: tell the user is selecting 
           # an overlapping region or NULL when click on store
      myregion <- reactiveValues(selectRegion=NULL , selectSingle=NULL 
                                 , compoundRegion=NULL , textMessage=NULL)

        # Visualize single mutation on mouse click
      output$infoclick <- shiny::renderTable({
        if(is.null(input$plot_click)){
          return(NULL)
        }
        mut_aligned <- mut_aligned_df()[[1]]
        mut_aligned_df <- data.frame(Position=seq_len(ncol(mut_aligned))
                                    , Mutations=mut_aligned[ 1 , , drop=TRUE])
        mut_aligned_df <- mut_aligned_df[ mut_aligned_df$Mutations!=0 , ]
        mynearPoints <- nearPoints(mut_aligned_df , input$plot_click 
                                , xvar="Position" , yvar="Mutations"
                                , maxpoints=1)
        if(nrow(mynearPoints)==0){
          return()
        }
        mynearPoints$Mutations <- as.integer(mynearPoints$Mutations)
        perc_mutation <- 100*(mynearPoints$Mutations/sum(
          mut_aligned_df$Mutations))
        perc_mutation <- paste0(format(round(perc_mutation , 2)
                                       , nsmall=2) , "%")
        mynearPoints$Percentage_Of_Covered_Mutations <- perc_mutation
        mymuts <- lmObj()$mutations$data
        mydataPosition <- mymuts[ mymuts$Amino_Acid_Position==
                                    mynearPoints$Position ,  ]
        tumtypeTab <- table(mydataPosition$Tumor_Type)
        mynearPoints$Tumor_Type <- paste( paste(names(tumtypeTab) 
                                              , unname(tumtypeTab) 
                                              , sep=":") 
                                          , collapse="|")
        aminochangeTab <- table(mydataPosition$Amino_Acid_Change)
        mynearPoints$Amino_Acid_Change <- paste( paste(names(aminochangeTab) 
                                                      , unname(aminochangeTab) 
                                                      , sep=":")
                                                , collapse="|")
        out <- data.frame(Gene_Symbol = rownames(mut_aligned)
            , Span = paste( mynearPoints$Position 
                      , mynearPoints$Position  , sep="-")
            , Region_Width_In_bp = 3
            , Span_Percentage = 
              paste0(format(round(100/ncol(mut_aligned) , 2), nsmall=2) , "%")
            , Percentage_Of_Covered_Mutations = perc_mutation
            , stringsAsFactors=FALSE)
        myregion$selectSingle <- out
        return(mynearPoints)
      })
        # Show the features of the selcted region with brush points      
      output$infobrush <- shiny::renderTable({
        if(is.null(input$plot_brush)) {
          myregion$selectRegion <- NULL
          return(NULL)
        }
        mut_aligned <- mut_aligned_df()[[1]]
        mut_aligned_df <- data.frame(Position=seq_len(ncol(mut_aligned))
                                    , Mutations=mut_aligned[ 1 , , drop=TRUE])
        start <- ifelse(floor(input$plot_brush$xmin)<1 , 1 
                        , floor(input$plot_brush$xmin))
        end <- ifelse(ceiling(input$plot_brush$xmax)>nrow(mut_aligned_df) 
                    , nrow(mut_aligned_df)
                    , ceiling(input$plot_brush$xmax))
        myinterval <- start:end
        myblock <- mut_aligned_df[ mut_aligned_df$Position %in% myinterval , ]
        perc_space <- end - start + 1
        perc_space <- 100*(perc_space/(nrow(mut_aligned_df)))
        perc_space <- paste0(format(round(perc_space , 2), nsmall=2) , "%")
        perc_mutation <- 100*(sum(myblock$Mutations)/sum(
          mut_aligned_df$Mutations))
        perc_mutation <- paste0(format(round(perc_mutation , 2)
                                       , nsmall=2) , "%")
        out <- data.frame(Gene_Symbol=rownames(mut_aligned)
                          , Span = paste( start , end  , sep="-")
                          , Region_Width_In_bp = as.integer(3*(end - start + 1))
                          , Span_Percentage = perc_space
                          , Percentage_Of_Covered_Mutations = perc_mutation
                          , stringsAsFactors=FALSE)
        myregion$selectRegion <- out
        return(out)
      })
        # Store a single mutation in manual selection
      observeEvent(input$storesingle , {
          if(!is.null(myregion$selectSingle)){
          tmp <- myregion$selectSingle
          if(is.null(myregion$compoundRegion)){
            myregion$compoundRegion <- tmp
          } else {
            if(!tmp$Gene_Symbol %in% myregion$compoundRegion$Gene_Symbol){
              myregion$compoundRegion <- unique(rbind(myregion$compoundRegion
                                                      , tmp))
            } else {
              tomerge <- myregion$compoundRegion[ 
                myregion$compoundRegion$Gene_Symbol==tmp$Gene_Symbol 
                , , drop=FALSE ]
              spantmp <- strsplit(as.character(tmp$Span) , "-") %>% 
                unlist %>% as.integer
              for(i in seq_len(nrow(tomerge))){
                myspan <- as.character(tomerge[i , "Span"])
                intervalAmino <- strsplit(myspan , "-") %>% 
                  unlist %>% as.integer
                intervalAmino <- intervalAmino[1]:intervalAmino[2]
                if(spantmp[1] %in% intervalAmino | spantmp[2] %in% 
                   intervalAmino){
                  myregion$textMessage <- 
                    paste("This region overlaps with one of the previous ones."
                        ,"When you close and save, the regions will be merged"
                        , sep="\n")
                  myregion$compoundRegion <- 
                    unique(rbind(myregion$compoundRegion , tmp))
                } else {
                  myregion$compoundRegion <- 
                    unique(rbind(myregion$compoundRegion , tmp))
                  myregion$textMessage <- NULL
                }
              }
            }
          }
        }
      })
        # Store a region in manual selection
      observeEvent(input$storeregion , {
          if(!is.null(myregion$selectRegion)){
          tmp <- myregion$selectRegion
          if(is.null(myregion$compoundRegion)){
            myregion$compoundRegion <- tmp
          } else {
            if(!tmp$Gene_Symbol %in% myregion$compoundRegion$Gene_Symbol){
              myregion$compoundRegion <- unique(rbind(myregion$compoundRegion 
                                                      , tmp))
            } else {
              tomerge <- myregion$compoundRegion[ 
                myregion$compoundRegion$Gene_Symbol==tmp$Gene_Symbol 
                , , drop=FALSE ]
              spantmp <- strsplit(as.character(tmp$Span) , "-") %>% 
                unlist %>% as.integer
              for(i in seq_len(nrow(tomerge))){
                myspan <- as.character(tomerge[i , "Span"])
                intervalAmino <- strsplit(myspan , "-") %>% 
                  unlist %>% as.integer
                intervalAmino <- intervalAmino[1]:intervalAmino[2]
                if(spantmp[1] %in% intervalAmino | spantmp[2] %in% 
                   intervalAmino){
                  myregion$textMessage <- 
                    paste("This region overlaps with one of the previous ones."
                          ,"When you close and save, the regions will be merged"
                          , sep="\n")
                  myregion$compoundRegion <- 
                    unique(rbind(myregion$compoundRegion , tmp))
                } else {
                  myregion$compoundRegion <- 
                    unique(rbind(myregion$compoundRegion , tmp))
                  myregion$textMessage <- NULL
                }
              }
            }
          }
        }
      })
        # If a new brush event occurs, the text message is put back to NULL
      observeEvent(input$plot_brush , {
        myregion$textMessage <- NULL
      })
        # Discard a region from the panel Optimize Panel (tab 4)
      observeEvent(input$discard , {
        if(!is.null(input$panel_rows_selected)){
          temp <- myregion$compoundRegion[ 
            -as.numeric(input$panel_rows_selected) , , drop=FALSE]
          temp <- unique(temp)
          myregion$compoundRegion <- temp
        }
      })
        # Show the error message when attempting to 
        # select two overlapping regions in manual selection
      output$textMessage <- renderText({
        myregion$textMessage
        })
        # Render the Optimize panel table
      output$panel <- DT::renderDataTable({
        myregion$compoundRegion
      } , server=FALSE
        , rownames=FALSE)

      # ------------------------------------------------------
      # Output Result
      # ------------------------------------------------------

      # Return compoundRegion content in case 
      # of close and save or if the app is closed
      observeEvent(input$stopandsave , {
          stopApp(isolate({
            if(is.null(myregion$compoundRegion)){
              colsToKeep <- c("drug" , "gene_symbol" 
                , "alteration" , "exact_alteration" 
                , "mutation_specification" , "group")
              panelNew <- panel[ , colsToKeep] %>% unique
              list(regions=NULL , mergedRegions=NULL , panel=panelNew)
            } else if(nrow(myregion$compoundRegion)==0){
              colsToKeep <- c("drug" , "gene_symbol" 
                , "alteration" , "exact_alteration" 
                , "mutation_specification" , "group")
              panelNew <- panel[ , colsToKeep] %>% unique
              list(regions=NULL , mergedRegions=NULL , panel=panelNew)
            } else {
              finalRegionSelected <- myregion$compoundRegion
              splitregions <- strsplit(finalRegionSelected$Span , "-")
              finalRegionSelected$start <- 
                as.integer(vapply(splitregions , '[' , character(1) , 1))
              finalRegionSelected$end <- 
                as.integer(vapply(splitregions , '[' , character(1) , 2))
              finalRegionSelectedsplit <- 
                split(finalRegionSelected , finalRegionSelected$Gene_Symbol)
              reducer <- lapply(finalRegionSelectedsplit , function(df){
                  myir <- IRanges(df$start , df$end)
                  myir <- reduce(myir)
                  out <- data.frame(region=paste(myir@start 
                                    , myir@start + myir@width - 1 , sep="-")
                                    ,gene_symbol=unique(df$Gene_Symbol)
                                    ,stringsAsFactors=FALSE)
                  return(out)
                }) %>% do.call("rbind" , .)
              panelNew <- merge(panel , reducer , all.x=TRUE)
              idxToSubstituteWithRegion <- 
                panelNew$alteration=="SNV" & 
                panelNew$exact_alteration %in% c("" , "mutation_type") & 
                !is.na(panelNew$region)
              panelNew$exact_alteration <- ifelse(idxToSubstituteWithRegion 
                            ,"amino_acid_position" ,panelNew$exact_alteration)
              panelNew$mutation_specification <- 
                ifelse(idxToSubstituteWithRegion
                            , panelNew$region ,panelNew$mutation_specification)
              colsToKeep <- c("drug" , "gene_symbol" 
                      , "alteration" , "exact_alteration"
                      , "mutation_specification" , "group")
              panelNew <- panelNew[ , colsToKeep] %>% unique
              list(regions=myregion$compoundRegion 
                   , mergedRegions=reducer , panel=panelNew)
            }
           })
         )
      })
      session$onSessionEnded(function() {
        stopApp(isolate({
            if(is.null(myregion$compoundRegion)){
              colsToKeep <- c("drug" , "gene_symbol" 
                              , "alteration" , "exact_alteration" 
                              , "mutation_specification" , "group")
              panelNew <- panel[ , colsToKeep] %>% unique
              list(regions=NULL , mergedRegions=NULL , panel=panelNew)
            } else if(nrow(myregion$compoundRegion)==0){
              colsToKeep <- c("drug" , "gene_symbol" 
                              , "alteration" , "exact_alteration" 
                              , "mutation_specification" , "group")
              panelNew <- panel[ , colsToKeep] %>% unique
              list(regions=NULL , mergedRegions=NULL , panel=panelNew)
            } else {
              finalRegionSelected <- myregion$compoundRegion
              splitregions <- strsplit(finalRegionSelected$Span , "-")
              finalRegionSelected$start <- 
                as.integer(vapply(splitregions , '[' , character(1) , 1))
              finalRegionSelected$end <- 
                as.integer(vapply(splitregions , '[' , character(1) , 2))
              finalRegionSelectedsplit <- 
                split(finalRegionSelected , finalRegionSelected$Gene_Symbol)
              reducer <- lapply(finalRegionSelectedsplit , function(df){
                myir <- IRanges(df$start , df$end)
                myir <- reduce(myir)
                out <- data.frame(region=paste(myir@start 
                                  , myir@start + myir@width - 1 , sep="-")
                                  ,gene_symbol=unique(df$Gene_Symbol)
                                  ,stringsAsFactors=FALSE)
                return(out)
                }) %>% do.call("rbind" , .)
              panelNew <- merge(panel , reducer , all.x=TRUE)
              idxToSubstituteWithRegion <- panelNew$alteration=="SNV" & 
                              panelNew$exact_alteration %in% 
                c("" , "mutation_type") & 
                              !is.na(panelNew$region)
              panelNew$exact_alteration <- ifelse(idxToSubstituteWithRegion 
                            ,"amino_acid_position" ,panelNew$exact_alteration)
              panelNew$mutation_specification <- 
                ifelse(idxToSubstituteWithRegion ,panelNew$region 
                       ,panelNew$mutation_specification)
              colsToKeep <- c("drug" , "gene_symbol" 
                              , "alteration" , "exact_alteration" 
                              , "mutation_specification" , "group")
              panelNew <- panelNew[ , colsToKeep] %>% unique
              list(regions=myregion$compoundRegion
                   , mergedRegions=reducer , panel=panelNew)
            }
           })
         )
      })
    }}
  ))
})