#' @title shiny-based GUI of NetLCP
#' @description launch shiny-based GUI of NetLCP
#'
#'
#' @examples
#' NetLCPGUI()
#' @export
#---
# title: "NetLCP-GUI"
# author: "Ming-Yu Ran"
# date: "2022/4/26"
#---
NetLCPGUI = function(){
  ###############################################.
  ## Header ----
  ###############################################.

  ui = shiny::fluidPage(
    theme = bslib::bs_theme(version = 4, bootswatch = "minty"), #Theme of the app
    shiny::navbarPage(
      title = shiny::div(shiny::tags$a("NetLCP", href= "https://github.com/rmyhandsome"), style = "font-size:25px;"), # Navigation title
      windowTitle = "NetLCP", #title for browser tab
      collapsible = TRUE, #tab panels collapse into menu in small screens
      ###############################################.
      ## Highlighting regulatory elements
      ###############################################.
      shiny::tabPanel(
        title = "Highlighting regulatory elements",

        shiny::sidebarLayout(
          shiny::sidebarPanel(shiny::fileInput("prioritization_input_file", shiny::h3("File input:")),
                              shiny::textAreaInput("prio_input_text", label = shiny::h3("EntrezID/miRBaseID"), value = NULL),
                              shiny::verbatimTextOutput("prio_input_judge"),
                              shiny::radioButtons("prio_type", label = shiny::h3("Biological Elements"),
                                  choices = list("KEGG" = "KEGG", "Reactome" = "Reactome", "Wikipathway" = "Wikipathway", "LncRNA" = "lncRNA", "CircRNA" = "circRNA"),
                                  selected = "KEGG"),
                              # shiny::radioButtons("empirical_pval", label = shiny::h3("Empirical Pvalue"),
                              #     choices = list("Yes (tips: it could take a long time)" = "TRUE", "No" = "FALSE"),
                              #     selected = "FALSE"),
                     shinyjs::useShinyjs(),
                     shiny::actionButton("prioritize_action", label = "Begin", style = "position: relative; left: 70%;")),
          shiny::mainPanel(
            shiny::tabsetPanel(
              shiny::tabPanel(title = "Prioritization Results",
                              shiny::br(),
                              shiny::br(),
                              shiny::verbatimTextOutput("progress"),
                              shiny::br(),
                              shiny::br(),
                              shiny::tableOutput("prioritization_out")
                    ),
              shiny::tabPanel(title = "Prioritization Output", shiny::br(), shiny::br(), shiny::br(), shiny::downloadButton('download_prio_result', 'Download'))
                  )
                )
        )
     ),
     ###############################################.
     ## local Regulation
     ###############################################
     shiny::tabPanel(
       title = "CREs in local area",

       shiny::sidebarLayout(
         shiny::sidebarPanel(shiny::fileInput("extract_input_file", shiny::h3("File input:")),
                             shiny::textAreaInput("extr_input_text", label = shiny::h3("Element ID"), value = NULL),
                             shiny::verbatimTextOutput("extract_input_judge"),
                             shiny::selectInput("reg_select", label = shiny::h3("Regulation Type"),
                                  choices = list("miRNA-mRNA" = "miRNA-mRNA", "circRNA-miRNA" = "circRNA-miRNA", "lncRNA-miRNA" = "lncRNA-miRNA",
                                                 "lncRNA-mRNA" = "lncRNA-mRNA", "lncRNA-miRNA-mRNA" = "lncRNA-miRNA-mRNA", "circRNA-miRNA-mRNA" = "circRNA-miRNA-mRNA",
                                                 "miRNA-pathway" = "miRNA-pathway",
                                                 "mRNA-pathway" = "mRNA-pathway", "miRNA-mRNA-pathway" = "miRNA-mRNA-pathway",
                                                 "lncRNA-miRNA-mRNA-pathway", "circRNA-miRNA-mRNA-pathway" = "circRNA-miRNA-mRNA-pathway")),
                             shiny::radioButtons("extract_range", label = shiny::h3("Range"),
                                   choices = list("CREs between input elements" = "FALSE", "All associated CREs in storage" = "TRUE"),
                                   selected = "FALSE"),
                      shinyjs::useShinyjs(),
                      shiny::actionButton("extract_action", label = "Inspect", style = "position: relative; left: 83%;")),
         shiny::mainPanel(
           shiny::tabsetPanel(
             shiny::tabPanel(title = "Results",
                             shiny::br(),
                             shiny::tableOutput("extract_out")
             ),
             shiny::tabPanel(title = "Output", shiny::br(), shiny::br(), shiny::br(), shiny::downloadButton('download_extr_result', 'Download'))
           )
         )
       )
     ),

     ###############################################.
     ## Regulation visualization
     ###############################################
     shiny::tabPanel(
       title = "CREs visualization",
       shiny::sidebarLayout(
         shiny::sidebarPanel(
           shiny::sliderInput("PR_node_degree", label = shiny::h3("Node Degree"), min = 0, max = 100, value = 0),
           shiny::radioButtons("PR_Net_layout", label = shiny::h3("Network Layout"),
                        choices = list("Circle" = "layout_in_circle", "Nicely" = "layout_nicely"),
                        selected = "layout_in_circle"),
           shiny::textAreaInput("PR_filter_input_text", label = shiny::h3("Filter by elements or CREs"), value = NULL),
           shiny::actionButton("PR_filter_action", label = "Filter",style = "position: relative; left: 83%;")
         ),
         shiny::mainPanel(
           shiny::tabsetPanel(
             shiny::tabPanel(title = "Statistics", plotly::plotlyOutput("PR_stat", width = "100%", height = "800px", inline = FALSE, reportTheme = TRUE)),
             shiny::tabPanel(title = "Network", visNetwork::visNetworkOutput("PR_net", width = "100%", height = "800px"))
           )
         )
       )
     ),

     ###############################################.
     ## Elements eQTLs regulation
     ###############################################
     shiny::tabPanel(
       title = "Prioritizing CREs",
       shiny::sidebarLayout(
         shiny::sidebarPanel(
           shiny::sliderInput("ER_node_degree", shiny::h3("Node Degree"), min = 0, max = 200, value = 0),
           shiny::sliderInput("ER_top_CREs", shiny::h3("Top CREs"), min = 0, max = 200, value = 10),
           shiny::radioButtons("ER_Net_layout", label = shiny::h3("Network Layout"),
                        choices = list("Circle" = "layout_in_circle", "Nicely" = "layout_nicely"),
                        selected = "layout_in_circle"),
           shiny::textAreaInput("ER_filter_input_text", label = shiny::h3("Filter by elements or CREs"), value = NULL),
           shiny::actionButton("ER_filter_action", label = "Filter",style = "position: relative; left: 83%;")
         ),
         shiny::mainPanel(
           shiny::tabsetPanel(

             shiny::tabPanel(title = "Statistics",
                             shiny::verticalLayout(
                               plotly::plotlyOutput("ER_stat_single", width = "100%", height = "200px", inline = FALSE, reportTheme = TRUE),
                               plotly::plotlyOutput("ER_stat_reg", width = "100%", height = "600px", inline = FALSE, reportTheme = TRUE)
                      )),
             shiny::tabPanel(title = "Network", visNetwork::visNetworkOutput("ER_net", width = "100%", height = "800px")),
             shiny::tabPanel(title = "Prioritized CREs", shiny::br(), shiny::tableOutput("prioritizedreg_out"))
           )
         )
       )
     ),
     ###############################################.
     ## Binding variants regulation
     ###############################################
     shiny::tabPanel(
       title = "CREs Variant 'switches'",
       shiny::sidebarLayout(
         shiny::sidebarPanel(
           shiny::textAreaInput("BR_filter_input_text", label = shiny::h3("Filter by elements or CREs"), value = NULL),
           shiny::actionButton("BR_filter_action", label = "Filter",style = "position: relative; left: 83%;")
         ),
         shiny::mainPanel(
           shiny::tabsetPanel(
             shiny::tabPanel(title = "Statistics", plotly::plotlyOutput("BR_stat_reg", width = "100%", height = "800px", inline = FALSE, reportTheme = TRUE)),
             shiny::tabPanel(title = "Network", visNetwork::visNetworkOutput("BR_net", width = "100%", height = "800px"))

           )
         )
       )
     )
    )
  )
  server = function(input, output, session) {
    ###############################################.
    ## Prioritizaion
    ###############################################
    # file upload
    prioInputFile <- shiny::reactiveValues(file_prio_input = NULL)
    shiny::observeEvent(input$prioritization_input_file, {
      prioInputFile$file_prio_input <- as.character(as.matrix(read.table(input$prioritization_input_file$datapath,
                                 sep = '\t', header = FALSE,
                                 stringsAsFactors = FALSE)))
      output$prio_input_judge = shiny::renderText({
        paste0("You should input at least 50 unique elements, now totally input ", length(c(prioInputFile$file_prio_input, unlist(strsplit(input$prio_input_text, "\\s+")))), " elements......")
      })

    })
    # text input
    shiny::observeEvent(input$prio_input_text, {

      output$prio_input_judge = shiny::renderText({

        paste0("You should input at least 50 unique elements, now totally input ", length(c(prioInputFile$file_prio_input, unlist(strsplit(input$prio_input_text, "\\s+")))), " elements......")
      })

    })


    #
    shiny::observeEvent(c(input$prioritization_input_file, input$prio_input_text), {

      if(length(c(prioInputFile$file_prio_input, unlist(strsplit(input$prio_input_text, "\\s+")))) < 50){
        output$progress = shiny::renderText({
          "Wait for valid input......"
        })
      }else{
        output$progress = shiny::renderText({
          'Prioritization is prepared to begin'
        })
      }

      shinyjs::enable('prioritize_action')

    })


    #
    shiny::observeEvent(c(input$prio_type, input$empirical_pval), {

      shinyjs::enable('prioritize_action')

    })
    #

    prio_input = shiny::eventReactive(input$prioritize_action,{
      c(prioInputFile$file_prio_input, unlist(strsplit(input$prio_input_text, "\\s+")))
    })
    prio_type_input = shiny::eventReactive(input$prioritize_action,{
      input$prio_type
    })
    empirical_pval_input = shiny::eventReactive(input$prioritize_action,{
      input$empirical_pval
    })

    shiny::observeEvent(input$prioritize_action, {
      if(length(c(prioInputFile$file_prio_input, unlist(strsplit(input$prio_input_text, "\\s+")))) >= 50){
        shinyjs::disable('prioritize_action')
      }
    })

    prio_output = shiny::reactiveValues(output = NULL)

    output$prioritization_out = shiny::renderTable({
      prio_output$output = BioRegElePrioritization(transcriptomeList = prio_input(), prioType = prio_type_input(), empiricalPvalue = FALSE)
      prio_output$output
    })

    output$download_prio_result = shiny::downloadHandler(
      filename = function(){
        paste(prio_type_input(), "_prioritization_result_", Sys.Date(), ".csv", sep="")
      },
      content = function(file){
        write.csv(prio_output$output, file, row.names = F)
      }
    )

    ###############################################.
    ## Regulation Extract
    ###############################################
    # file upload
    ExtrInputFile <- shiny::reactiveValues(file_extract_input = NULL)
    shiny::observeEvent(input$extract_input_file, {
      ExtrInputFile$file_extract_input <- as.character(as.matrix(read.table(input$extract_input_file$datapath,
                                                                         sep = '\t', header = FALSE,
                                                                         stringsAsFactors = FALSE)))
      output$extract_input_judge = shiny::renderText({
        paste0("now totally input ", length(c(ExtrInputFile$file_extract_input, unlist(strsplit(input$extr_input_text, "\\s+")))), " elements......")
      })
      shinyjs::enable('extract_action')

    })
    # text input
    shiny::observeEvent(input$extr_input_text, {

      output$extract_input_judge = shiny::renderText({

        paste0("now totally input ", length(c(ExtrInputFile$file_extract_input, unlist(strsplit(input$extr_input_text, "\\s+")))), " elements......")
      })
      shinyjs::enable('extract_action')

    })

    #
    shiny::observeEvent(c(input$reg_select, input$extract_range), {

      shinyjs::enable('extract_action')

    })
    #

    Extr_input = shiny::eventReactive(input$extract_action,{
      c(ExtrInputFile$file_extract_input, unlist(strsplit(input$extr_input_text, "\\s+")))
    })
    Extr_type_input = shiny::eventReactive(input$extract_action,{
      input$reg_select
    })
    Extr_range = shiny::eventReactive(input$extract_action,{
      input$extract_range
    })

    Extr_output = shiny::reactiveValues(output = NULL)

    ############# click extract event ################
    shiny::observeEvent(input$extract_action, {
      if(length(c(ExtrInputFile$file_extract_input, unlist(strsplit(input$extr_input_text, "\\s+")))) > 0 ){
        shinyjs::disable('extract_action')
      }

      if(Extr_type_input() == "circRNA-miRNA-mRNA" | Extr_type_input() == "lncRNA-miRNA-mRNA" | Extr_type_input() == "miRNA-mRNA-pathway" | Extr_type_input() == "lncRNA-miRNA-mRNA-pathway" | Extr_type_input() == "circRNA-miRNA-mRNA-pathway"){
        output$extract_out = shiny::renderTable({
          Extr_output$output = multieleRegulation(elementList = Extr_input(), regulationType = Extr_type_input(),  allRegulation = Extr_range())
          Extr_output$output
        })
      }else{
        output$extract_out = shiny::renderTable({
          Extr_output$output = binaryRegulation(elementList = Extr_input(), regulationType = Extr_type_input(),  allRegulation = Extr_range())
          Extr_output$output
        })
      }

      ######## Statistics ######
      shiny::observe({
          if(!is.character(Extr_output$output)){
            ##### Pure Regulation ####
            PR_filter_input = shiny::eventReactive(input$PR_filter_action,{
              input$PR_filter_input_text
            })
            shiny::observeEvent(input$PR_filter_action,{
              output$PR_stat = plotly::renderPlotly({regStat(regData = Extr_output$output, filterDegree = input$PR_node_degree, selectCREs = unlist(strsplit(PR_filter_input(), "\\s+")))})
              output$PR_net = visNetwork::renderVisNetwork({regNetVis(regData = Extr_output$output, filterDegree = input$PR_node_degree, selectCREs = unlist(strsplit(PR_filter_input(), "\\s+")), netLayout = input$PR_Net_layout)})
            })
            # netvis
            output$PR_stat = plotly::renderPlotly({regStat(regData = Extr_output$output, filterDegree = input$PR_node_degree, selectCREs = NULL)})
            output$PR_net = visNetwork::renderVisNetwork({regNetVis(regData = Extr_output$output, filterDegree = input$PR_node_degree, selectCREs = NULL, netLayout = input$PR_Net_layout)})

            ###### Elements eQTLs regulation ######
            # Extr_output$output
            eqtls_out = shiny::reactiveValues(output = NULL)
            eqtls_out$output = eQTLsDetection(regData = Extr_output$output)

            ER_filter_input = shiny::eventReactive(input$ER_filter_action,{
              input$ER_filter_input_text
            })
            shiny::observeEvent(input$ER_filter_action,{
              output$ER_stat_single = plotly::renderPlotly({eQTLsSingleEleStat(regData = Extr_output$output, eQTLsData = eqtls_out$output, filterDegree = input$ER_node_degree, selectCREs = unlist(strsplit(ER_filter_input(), "\\s+")))})
              reg_out = shiny::reactiveValues(output = NULL)
              reg_out$output = eQTLsRegStat_GUI(regData = Extr_output$output, eQTLsData = eqtls_out$output, regulationType = Extr_type_input(), filterDegree = input$ER_node_degree, topCREs = input$ER_top_CREs, selectCREs = unlist(strsplit(ER_filter_input(), "\\s+")))
              output$ER_stat_reg = plotly::renderPlotly({reg_out$output[[2]]})
              output$prioritizedreg_out = shiny::renderTable({reg_out$output[[1]]})
              output$ER_net = visNetwork::renderVisNetwork({eQTLsNetVis(regData = Extr_output$output, eQTLsData = eqtls_out$output, filterDegree = input$ER_node_degree, selectCREs = unlist(strsplit(ER_filter_input(), "\\s+")), netLayout = input$ER_Net_layout)})
            })
            output$ER_stat_single = plotly::renderPlotly({eQTLsSingleEleStat(regData = Extr_output$output, eQTLsData = eqtls_out$output, filterDegree = input$ER_node_degree, selectCREs = NULL)})
            reg_out = shiny::reactiveValues(output = NULL)
            reg_out$output = eQTLsRegStat_GUI(regData = Extr_output$output, eQTLsData = eqtls_out$output, regulationType = Extr_type_input(), filterDegree = input$ER_node_degree, topCREs = input$ER_top_CREs, selectCREs = NULL)
            output$ER_stat_reg = plotly::renderPlotly({reg_out$output[[2]]})
            output$prioritizedreg_out = shiny::renderTable({reg_out$output[[1]]})
            output$ER_net = visNetwork::renderVisNetwork({eQTLsNetVis(regData = Extr_output$output, eQTLsData = eqtls_out$output, filterDegree = input$ER_node_degree, selectCREs = NULL, netLayout = input$ER_Net_layout)})

            ###### regulation variants ######
            regvar_out = shiny::reactiveValues(output = NULL)
            regvar_out$output = regVarDetection(regData = Extr_output$output,regulationType = Extr_type_input())

            BR_filter_input = shiny::eventReactive(input$BR_filter_action,{
              input$BR_filter_input_text
            })
            shiny::observeEvent(input$BR_filter_action,{
              output$BR_stat_reg = plotly::renderPlotly({regVarStat(regVar = regvar_out$output, regulationType = Extr_type_input(), selectCREs = unlist(strsplit(BR_filter_input(), "\\s+")))})
              output$BR_net = visNetwork::renderVisNetwork({regVarNetVis(regVar = regvar_out$output, regulationType = Extr_type_input(), selectCREs = unlist(strsplit(BR_filter_input(), "\\s+")))})

            })
        }
      })

    })

    output$download_extr_result = shiny::downloadHandler(
      filename = function(){
        paste(Extr_type_input(), "_Extraction_result_", Sys.Date(), ".csv", sep="")
      },
      content = function(file){
        write.csv(Extr_output$output, file, row.names = F)
      }
    )

    ###############################################.
    ## Visual Var
    ###############################################



  }

  shiny::shinyApp(ui, server)
}

