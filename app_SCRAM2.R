#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# A. Gilbert
# Biodiversity Research Institute
# 276 Canco Rd
# Portland, ME 04103
# This tool was modified from one largely developed by Chris Field at the University of Rhode Island and depends on code
# from the Stochastic Band model
# created: 16 Feb 2022
# Modified 06 April 2022 - pre-render movement data for plotting and tables for PIPL, REKN, ROST; prerender summary flt ht for same
#   Add back probability prediction for exceeding collision threshold
#   Create probability of occurrence polygons layers for PIPL, REKN, ROST from baked data files for mapping
# 27 April 22 - recreated output report using DPLYR and ggplot to simplify
# 05 May 22 - fixed issues with report and output tables to tabs in results figure 
# 10 May 22 - 0.74.3 - table reorganization in output
# 12 May 22 - Fixed from issues and rewrote baked movemement files with NAs in months according to when we had movememnt, per Pam Loring
# 19 May 22 - Added download folder for automatic saving of results in case of long model runs, but muted as it only works with local file systems which doesn't means
#   only local file service on remote machines.
# 25 May 22 - update baked movement models
# 13 Jun 22 - fully update movement models
# 15 Jun 22 - remove other species input for later release
# 21 Jul 22 - 0.80 - Upate movement and flight height models, bug fixes and improvements based on reviews
# 20 Sep 22	- 0.81 - update species data 
# 21 Sep 22 also change to report total monthly collision and the prediction interval for that estimate.
# 22 Sep 22 - change def of hub height add to correct to air gap
# 03 Oct 22 - 0.90 - fix movement model - some errors fixed by EMA, slight changes to report and manual. 
# 07 Oct 22 - 0.91 - fix bug where prob. of exceedenc not updating on subsequent runs and figures don't update. Also fixed issue with report names
#  not updating on re-runs, or zip files. 
# 11 Oct 22 - 0.91.1 - fix histogram results and slight change to report output
# 22 Dec 22 - 1.0 - Changed species occup. model maps to show mean monthly, not mean cum. daily (over all months), 
#   also added last position models to compare to first position (of day) models
# 24 Feb 23 - 1.0.1 - fixed reporting issue and made change to be much more like the Band 2012 migration model
# 01 Mar 23 - 1.0.2 - add coastal flag as dialog box and report warning.
# 09 Mar 23 - 1.0.3 - removed coastal flag and modified the number of animals in a model cell by  multiplying by the proportion of transient animals
# to reduce the effect of coastal staging animals on risk.
# 17 Aug 23 - 2.0 - modify SCRAM to use StochLab package and it's functions for calculating crm modified to accept annex 6 migrant flux calculations based on movement data


# load scripts
source("scripts/helpers.R")
source("scripts/utils.R")
source("scripts/get_mig_flux_SCRAM.R")
source("scripts/band_SCRAM.R")
source("scripts/stoch_SCRAM.R")

# SCRAM_version = "1.0.3 - Cathartic Adela" 
SCRAM_version = "2.0.0 - Altruistic Anaheim"   #https://www.cayennediane.com/big-list-of-hot-peppers/

options(shiny.trace = F)

##################################################################################################################################################
##                                                                                                                                              ##
##                                                              DASHBOARD HEADER                                                                ##
##                                                                                                                                              ##
##################################################################################################################################################
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(
    titleWidth = 400,
    
    tags$li(
      class = "dropdown",
      actionLink("appvrsn", label = tags$b(paste("Stochastic Collision Risk Assessment for Movement: v", SCRAM_version), style = "font-size: 14px")),
      style = "float: left"
    ),
    tags$li(
      class = "dropdown",
      a(id = "download_manual",
        icon('fa-solid fa-book', "fa-2x"),
        style = "padding-top: 10px; padding-bottom: 10px",
        target = '_blank',
        href = "SCRAM_manual_v103_031723.pdf"),
      style = "float: left"
    ),
    tags$li(
      class = "dropdown",
      a(
        icon('github', "fa-2x"),
        href = 'https://github.com/Biodiversity-Research-Institute/SCRAM2',
        style = "padding-top: 10px; padding-bottom: 10px; padding-left: 5px; padding-right: 5px",
        target = '_blank',
        id = "lbl_codeLink"
      ),
      style = "float: left"
    ),
    
    tags$li(
      class = "dropdown",
      a(
        icon('bug', "fa-2x"),
        href = 'https://github.com/Biodiversity-Research-Institute/SCRAM2/issues',
        style = "padding-top: 10px; padding-bottom: 10px; padding-left: 5px; padding-right: 5px",
        target = '_blank',
        id = "lbl_issuesLink"
      ),
      style = "float: left"
    ),
    
    tags$li(
      class = "dropdown",
      a(
        img(src = "BRI_color_logo_no_words.png", height = "40px"),
        href = 'https://briwildlife.org',
        style = "padding-top: 5px; padding-bottom: 5px; padding-left: 5px;, padding-right: 0px; margin-right: -5px",
        target = '_blank',
        id = "lbl_BRILogoLink"
      ),
      style = "float: right"
    ),
    tags$li(
      class = "dropdown",
      a(
        img(src = "URI.png", height = "40px"),
        href = 'https://URI.edu',
        style = "padding-top: 5px; padding-bottom: 5px; padding-left: 5px;, padding-right: 0px; margin-right: -15px",
        target = '_blank',
        id = "lbl_URILogoLink"
      ),
      style = "float: right"
    ),
    tags$li(
      class = "dropdown",
      a(
        img(src = "USFWS.png", height = "40px"),
        href = 'https://www.fws.gov/',
        style = "padding-top: 5px; padding-bottom: 5px; padding-left: 5px;, padding-right: 0px; margin-right: -5px",
        target = '_blank',
        id = "lbl_FWSLogoLink"
      ),
      style = "float: right"
    ),
    tags$li(
      class = "dropdown",
      a(
        img(src = "BOEM.png", height = "40px"),
        href = 'https://www.BOEM.gov',
        style = "padding-top: 5px; padding-bottom: 5px; padding-left: 5px;, padding-right: 0px; margin-right: -5px",
        target = '_blank',
        id = "lbl_BOEMLogoLink"
      ),
      style = "float: right"
    )
  ), #dashboardHeader
  
  ##################################################################################################################################################
  ##                                                                                                                                              ##
  ##                                                              DASHBOARD SIDEBAR                                                               ##
  ##                                                                                                                                              ##
  ##################################################################################################################################################
  
  dashboardSidebar(
    collapsed = F,
    width=415,
    # title= 
    
    sidebarMenu(
      id = "sidebar",
      style = "overflow-y:scroll; max-height: 800px; position:relative;",
      tags$a(
        img(src = "SCRAM_logo_400px.png", alt="Stochastic Collision Risk Assessment for Movement", width = "400px", class="header_img"),
        href = 'https://briwildlife.org/SCRAM',
        style = "margin-bottom: 10px;",# padding-right:20px; margin-bottom: -10px; margin-top:-45px;",
        target = '_blank',
        id = "lbl_SCRAMLogoLink"
      ),
      
      h4("1) SCRAM run details:", style = "padding-left: 10px; margin-bottom: 0px"),
      textInput(inputId = "project_name", label = "Project name: ", value = "", width = "380px", placeholder = "Project"),
      div(style = "margin-top: -20px;"),  #reduce space between elements
      textInput(inputId = "modeler", label = "Name of person running SCRAM: ", value = "", width = "380px", placeholder = "Name"),
      
      # conditionalPanel( 
      #   #show only when project names entered
      #   condition = ('input.project_name != "" & input.modeler != ""'),
        
        ################### Input: Select the migration calculation type - movement model or annex 6 type
        radioButtons(inputId = "migration_calc_type",
                     # label ="Select included species data or your own:",
                     label ="Select migration calcuation mode:",
                     choices = c("Movement model" = "occup_model", "Band 2012 annex 6" = "band_2012_annex_6"),
                     selected = "occup_model"),
      # ),
      
      #################Enter wind farm parameters
      conditionalPanel( 
        #show only when species data have been inputted
        # condition = 'input.species_input',
        condition = 'input.project_name != "" & input.modeler != ""',
        
        h4("2) Load wind farm parameters:", 
           style = "padding-left: 10px; margin-bottom: 0px"),
        fileInput(
          "file_wf_param",
          "Upload wind farm data",
          accept = c('text/csv',
                     'text/comma-separated-values,text/plain',
                     '.csv'),
          width = '95%'),
        
        div(style = "margin-top: 0px;"),  #reduce space between elements
      
      conditionalPanel( 
        #show only when project names entered
        # condition = ('input.project_name != "" & input.modeler != ""'),
        condition = ('output.fileUploaded'),
        
        h4("3) Select the species or load species data:", style = "padding-left: 10px; margin-bottom: 0px"),
        
        ################### Input: Select the model type and species to model

        # div(style = "margin-top: -20px;"),  #reduce space between elements
        
        radioButtons(inputId = "species_input",
                     # label ="Select included species data or your own:",
                     label ="Select included species:",
                     # choices = c("Piping Plover" = "Piping_Plover", "Red Knot" = "Red_Knot", "Roseate Tern" = "Roseate_Tern", "Common Tern" = "Common_Tern","Use your own species data" = "Other"),
                     choices = c("Piping Plover" = "Piping_Plover", "Red Knot" = "Red_Knot", "Roseate Tern" = "Roseate_Tern"),
                     selected = character(0)), #start with no items selected
        
        #Deactivated ability to add Other species for now. Mute below until can be fixed for next version
        # conditionalPanel(  #only expand if you selected other to hide the load button
        #   condition = "input.species_input == 'Other'",
        #   # downloadButton("downloadSpeciesExample", "Download example species input",
        #   #                style = "margin-left: 20px; margin-top: 0px; background-color: blue; color: white; font-weight: bold;"),
        #   fileInput("file_spp_param", "Upload species data and flight height distributions", accept = ".csv", multiple = TRUE, width = '90%', options(shiny.maxRequestSize = 25 * 1024^2))
        #   ),

          
          
          #################Upload migration passage data
          conditionalPanel(
            #show only when band_2012_annex_6 has been selected
            condition = 'input.migration_calc_type == "band_2012_annex_6"',
            div(style = "margin-top: -20px;"),  #reduce space between elements
            uiOutput("migr_front_type"),
            div(style = "margin-top: -20px;"),  #reduce space between elements
            
            # Add output for Annex 6 calcs showing the frontal width and closest state
            htmlOutput("migr_front_width", 
                       style = "padding: 10px; margin-top: -10px"),
            htmlOutput("nearest_state", 
                       style = "padding: 10px; margin-top: -10px")
            
            # fileInput("file_migr_passages", "Upload migration passage data", accept = c('text/csv',
            #                                                                             'text/comma-separated-values,text/plain',
            #                                                                             '.csv'), width = '95%'),
            
          )
        ),
        
      ),
      
      #################Enter CRM options
      conditionalPanel( 
        #show only when wind farm data have been inputted
        # condition = "output.fileUploaded",
        condition = 'input.species_input',
        h4("4) Select CRM parameter options:", style = "padding-left: 10px; margin-bottom: -10px"),
        radioButtons("crm_option_radio", "Band (2012) equivalent CRM options:",
                     c("Option 1: faster approximation" = "1", "Option 3: slower but more precise" = "3")),
        div(style = "margin-top: -20px;"),  #reduce space between elements
        sliderInput("iter_slider", label = "Model iterations (rec. min. 1,000)", min = 100, 
                    max = 10000, value = 1000, step=100, width = '95%'), 
        htmlOutput("iter_message", style = "margin-top: -10px; margin-left: 10px"),
        numericInput(
          inputId = "inputthreshold",
          label = "Annual collision threshold",
          value = 1, #round(1/35, 3),
          min = 0,
          max = NA,
          step = NA,
          width = '75%'),
        # hr(),
        h4("5) Run CRM:", style = "padding-left: 10px;"),

        fluidRow(column(12, 
                        uiOutput("runUI")), 
                 # column(6, 
                 #        uiOutput("cancelUI"))
                 ),
        br()
      )
    )
  ),
  
  ##################################################################################################################################################
  ##                                                                                                                                              ##
  ##                                                              DASHBOARD BODY                                                                  ##
  ##                                                                                                                                              ##
  ##################################################################################################################################################
  
  dashboardBody(
    # tags$head(
    #   tags$link(rel = "stylesheet",
    #             type = "text/css",
    #             href = "www/style.css")
    # ),
    
    useShinyjs(),
    
    # adds a indicator on slow load
    add_busy_spinner(spin = "atom", 
                     color = "blue",
                     position = "top-left", 
                     margins = c("50%", "50%"), 
                     onstart = T),
    
    tabsetPanel(
      id = "tabsetpan",
      type = "tabs",
      selected = "start",
      tabPanel(
        "Start Here",
        value = "start",
        fluidRow(
          box(
            # title = "Instructions",
            status = "primary",
            solidHeader = F,
            collapsible = T,
            width = 12,
            style = "margin-top: -20px; padding-bottom: 20px",
            fluidRow(
              column(width = 8, htmlOutput("user_instructions", style = "margin-top: -14px")),
              # Mute download buttons - for species data for now, don't allow uploading of custom species data
              # Also move manual to book button
              column(width = 4, 
                     htmlOutput("downloads", style = "margin-top: -14px"),
                     #        downloadButton("downloadSpeciesExample", "Example species input",
                     #                       style = "margin-bottom: 8px; background-color: blue; color: white; font-weight: bold;"),
                     #        br(),
                     downloadButton("downloadTurbineExample", "Example wind farm input",
                                    style = "margin-bottom: 8px; background-color: orange; color: white; font-weight: bold;"),
                     br(),
                     img(src = "CVOW_turbines.jpg", width = 200)
              )
            )
          ), 
          fluidRow(
            p("This tool was developed by Biodiversity Research Institute, The University of Rhode Island, and 
                    U.S. Fish and Wildlife Service with funding from the Bureau of Ocean Energy Management.", 
              style = "margin-left: 20px; margin-right: 20px; margin-top: -10px; color: steelblue; font-size: 14px")
          )
        )
      ), #tabpanel
      
      
      tabPanel(
        "Wind Farm Inputs",
        value = "wind_farm_panel",
        fluidRow(
          box(
            status = "primary",
            solidHeader = F,
            collapsible = T,
            width = 12,
            style = "margin-top: -20px; padding-bottom: 20px",
            htmlOutput("check_windfarm_instructions", style = "margin-top: -14px")
          )
        ),
        fluidRow(
          column(7,
                 fluidRow(
                   box(
                     title = "Wind Farm Inputs",
                     status = "success",
                     solidHeader = TRUE,
                     width = 12,
                     tabsetPanel(
                       tabPanel("Wind Farm Specs", dataTableOutput("wind_farm_data1"),
                                style = "height:500px; overflow-y: scroll;overflow-x: scroll;"),
                       tabPanel("Turbine Ops Data", dataTableOutput("ops_data"),
                                style = "height:500px; overflow-y: scroll;overflow-x: scroll;"))
                   )
                 )
          ),
          column(5,
                 fluidRow(
                   box(
                     title = "Wind Farm Map",
                     status = "primary",
                     solidHeader = TRUE,
                     width = 12,
                     leafletOutput("studymap", height = "535px", width = "100%")
                   )
                 ),
                 fluidRow(
                   htmlOutput("mean_prob_txt", 
                              style = "margin-left: 20px; margin-right: 20px; margin-top: -10px; color: steelblue; font-size: 10px")) 
                 # p("\u00B9Mean occup. prob.= mean monthly occupancy probability"
          )
        )
      ), #tabpanel wind farm data
      
      tabPanel(
        "Species Input",
        value = "species_panel",
        
        fluidRow(
          #show the species data prior to modeling for checks
          column(6,
                 fluidRow(
                   box(
                     title = "Species Inputs",
                     status = "primary",
                     solidHeader = TRUE,
                     width = 12,
                     dataTableOutput("species_data"),
                     style = "height:577px; overflow-y: scroll;overflow-x: scroll;")
                 )
          ),
          #Show the flight height data raw and as figure to help to make sure user check these before running
          column(6,
                 fluidRow(
                   box(title = "Flight Height Data",
                       status = "primary",
                       solidHeader = TRUE,
                       width = 12,
                       tabsetPanel(
                         id = "tabsetpan",
                         type = "tabs",
                         selected = "flt_ht_plot",
                         tabPanel(title = "Plot", 
                                  value = "flt_ht_plot",
                                  plotOutput("flt_ht_plot"), 
                                  style = "height:420px; overflow-y: scroll;overflow-x: scroll;"),
                         tabPanel(title = "Data", 
                                  value = "flt_ht_data",
                                  dataTableOutput("flt_ht_data"))
                         #style = "overflow-y: scroll;")
                       ), #tabsetPanel
                       column(3, h5("Flight heights displayed (m)", 
                                    style="margin-top:20px; margin-left:5px; margin-right:5px; margin-bottom:5px")),
                       column(9, sliderInput("slider_flt_ht", label = "", min = 0, max = 1000, value=c(0,1000),step=25, width = '100%')))
                 )
          ) #flight height column
        ), #species data fluidRow
        fluidRow(
          #show the species count data prior to modeling for checks
          column(12,
                 fluidRow(
                   box(
                     title = "Population Data",
                     solidHeader = TRUE,
                     width = 12,
                     dataTableOutput("count_data"),
                     htmlOutput("popn_notes", style = "margin-top: 0px; margin-left:20px; margin-right:20px; margin-bottom:10px; font-size: 12px")
                   )))
        )
        
      ), #tabpanel Species Data

      #CRM results tab
      tabPanel("CRM Results", value = "crm_results",
               fluidRow(
                 box(
                   title = "Results",
                   status = "primary",
                   solidHeader = F,
                   collapsible = F,
                   width = 12,
                   textOutput("run_start_txt"),
                   textOutput("run_end_txt"),
                   textOutput("run_success_msg"), 
                   htmlOutput("prob_exceed_threshold_msg")
                 )
               ),
               fluidRow(
                 column(7,
                        h4("Output dashboard", style = " margin-top: -10px;"), 
                        uiOutput("plot_tabs"),
                        uiOutput("plot_results_caption")
                        
                 ), 
                 #buttons for sensitivity Analysis, downloading output, and generating report
                 column(5,
                        h4("Next steps:", style = " margin-top: -10px;"), 
                        # uiOutput("runGSA2"), 
                        # br(),
                        uiOutput("download_output"), 
                        br(),
                        uiOutput("genreport")
                 )
               )
      )# tabpanel CRM results
    ) # tabsetpanel
  ) #dashboard body
) #dashboard page

verbose <- F

##################################################################################################################################################
##                                                                                                                                              ##
##                                                              Define server logic                                                             ##
##                                                                                                                                              ##
##################################################################################################################################################

server <- function(input, output, session) {
  
  # input for the model options
  crm_option_radio_react <- eventReactive(input$crm_option_radio, {input$crm_option_radio})
  
  output$fileUploaded  <- reactive({
    val <- !(is.null(input$file_wf_param))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  # number of iterations
  # iter_slider_react <- reactiveVal()
  iter_slider_react <- eventReactive(input$iter_slider, {as.numeric(input$iter_slider)})
  
  # approximate the run time for 1 species and 1 turbine; reactive with user input for number of iterations
  approx_run_times <- c(0.1, 1.2, 0.8, 5.9, 1.9, 3.1)
  
  # update estimated run time when user chooses a new option
  # below needs updating based on current models - remove for now. 
  # observeEvent(req(input$crm_option_radio), {output$iter_message <- renderText({
  #   if((approx_run_times*as.numeric(input$iter_slider))[as.numeric(input$crm_option_radio)] > 90 & (approx_run_times*as.numeric(input$iter_slider))[as.numeric(input$crm_option_radio)] <= 3600){
  #     time_estimate <- round((approx_run_times*as.numeric(input$iter_slider))[as.numeric(input$crm_option_radio)]/60)
  #     unit_time <- "minutes"
  #   }
  #   if((approx_run_times*as.numeric(input$iter_slider))[as.numeric(input$crm_option_radio)] > 3600){
  #     time_estimate <- round((approx_run_times*as.numeric(input$iter_slider))[as.numeric(input$crm_option_radio)]/3600)
  #     unit_time <- "hour(s)"
  #   }
  #   if((approx_run_times*as.numeric(input$iter_slider))[as.numeric(input$crm_option_radio)] <= 90){
  #     time_estimate <- (approx_run_times*as.numeric(input$iter_slider))[as.numeric(input$crm_option_radio)]
  #     unit_time <- "seconds"
  #   }
  #   paste("~", time_estimate, unit_time, "per wind farm/turbine option")
  # })
  # })
  
  # load labels to display versions of species names without underscores
  SpeciesLabels <- read.csv("data/SpeciesLabels.csv")
  
  # Set study map footnote text
  observeEvent(c(input$species_input), {
    req(input$migration_calc_type == "occup_model" & length(input$species_input) > 0)
    
    output$mean_prob_txt <- renderText("\u00B9Mean occup. prob.= mean monthly occupancy probability")
    
  })

  # # Calculate the corridor width when the species is selected
  # observeEvent(c(input$species_input), {
  #   req(input$migration_calc_type == "band_2012_annex_6" & length(input$species_input) > 0)
  #   
  #   #determine the location to use for 
  #   
  #   annex6_migr_spp <- annex6_migr_data %>% 
  #     dplyr::filter(Species == input$species_input)
  #   # output$migr_front_type <- renderUI({
  #   #   selectInput("migr_front_opt", "Migratory front width option", choices = annex6_migr_spp$migr_front_type)
  #   # })
  #   
  # })
  
  # #show the migration corridor width when a migration front type is chosen
  #  observe({
  #    req(!is.null(input$migr_front_opt))
  #    annex6_migr_width <- annex6_migr_data %>% 
  #      dplyr::filter(Species == input$species_input & migr_front_type == input$migr_front_opt)
  #    updateTextInput(session, "migr_front_width", value=annex6_migr_width$migr_front_km)
  #    
  #  })
  
  # Render and remove buttons -----------------------------------------------
  
  # if all conditions are met (conditions dependent on whether user is using baked-in species or uploading their own data) 
  # could simplify using bothdataup()
  observeEvent(c(input$file_wf_param, input$species_input), {
    # render the "run CRM" button to allow user to run main script
    output$runUI <- renderUI({
      div(style="display: flex; align-items: center; justify-content: center;",
          actionButton("run", "Run CRM", style = "width: 100px; background-color: green; color: white; font-weight: bold;"),
          id="run_div")
    })
    # # render cancel button according to conditions for "run CRM"
    # output$cancelUI <- renderUI({
    #   div(style="display:inline-block; float:left; padding-left: 10px",
    #       actionButton("cancel", "Cancel", style = "width: 100px; background-color: red; color: white; font-weight: bold;"),
    #       id="cancel_div")
    # })
  })
  
  # box appears if main function is canceled
  observeEvent(input$cancel,{
    if(running())interruptor$interrupt("Canceled")
  })
  
  # if conditions above are no longer met, remove "run CRM" button
  observeEvent(req((length(input$file_spp_param$datapath) < 3 & length(which(input$species_input == "Other")) == 1 & !is.null(input$file_spp_param))|
                     length(input$species_input) == 0 | nrow(isolate(windfarm_loc$cell_sf))==0), { 
                       removeUI(
                         selector = "#run_div")
                     })
  
  # remove cancel button according to conditions for "run CRM"
  observeEvent(req((length(input$file_spp_param$datapath) < 3 & length(which(input$species_input=="Other")) == 1 &
                      !is.null(input$file_spp_param))|length(input$species_input) == 0 | nrow(isolate(windfarm_loc$cell_sf))==0), { 
                        removeUI(
                          selector = "#cancel_div")
                      })
  
  
  # Start Here tab ----------------------------------------------------------
  
  output$user_instructions <- renderText(
    "<h4>INSTRUCTIONS:</h4>
   <p>1) Enter the project name and person conducting the analysis. This will be saved in output. <br>
    2) Select the species of interest included with SCRAM. <br>
    3) Check the species data for expected values in the tables and figures in the 'Species Data' tab. If data values are not as
    expected, do NOT run SCRAM. These values are currently fixed and can't be changed. Future updates should allow custom data. <br>
    4) Download the example wind farm inuput data using the button to the right and either modify for your specific use
    or use it directly to demonstrate the use of the tool. <br>
    5) Upon proper loading of the wind farm data, the wind farm data tab will be shown.
    Check the wind farm data for correct values in the include maps and tables. 
    Correct any errors and reload as necessary. <br>
    6) Choose which version of the CRM to run.<br>
    7) Select the number of iterations (100-10,0000).<br>
    8) Set a threshold for the maximum acceptable number of collisions. <br>
    9) Run CRM. <br>
    10) Run sensitivity analyses (optional). <br>
    11) Generate summary report and/or download results for each iteration. <br>
    12) Check the CRM results.</p>")
  
  output$check_windfarm_instructions <- renderText(
    "<h4>INSTRUCTIONS:</h4>
     Check the wind farm data carefully below before running this tool to make sure it's correct.<br>
     Fix any data in your original data file and upload again."
  )
  
  output$downloads <- renderText("<h4>DOWNLOAD INPUT FILE(S):</h4>")
  
  
  # Species data (parameters, flight height, popn data) --------------------
  
  selected_species <- reactive(input$species_input)
  
  # #navigate to species tab once selected
  # observeEvent(input$species_input,
  #              updateTabItems(session, inputId = "tabsetpan", selected = "species_panel")
  # )
  
  #load bird data
  bird_data <- read.csv("data/BirdData_020923.csv", header = T)

  species_params_vals <- reactiveValues(migrate_resident = NULL, )
  
  #reprocess data for input to stochCRM function
  observe({
    req(input$species_input)
    
    #filter to species
    species_params <- bird_data %>%
      filter(Species == input$species_input)
    
    #select flight speed data and rename
    species_params_vals$flt_spd_data <- species_params %>% 
      select(Flight_Speed, Flight_SpeedSD) %>% 
      rename(mean = Flight_Speed, sd = Flight_SpeedSD)
    
    #select body length data and rename
    species_params_vals$body_lt_data <- species_params %>% 
      select(Body_Length, Body_LengthSD) %>% 
      rename(mean = Body_Length, sd = Body_LengthSD)
    
    #select Wingspan data and rename
    species_params_vals$wing_span_data <- species_params %>% 
      select(Wingspan, WingspanSD) %>% 
      rename(mean = Wingspan, sd = WingspanSD)
    
    #select Avoidance data and rename
    species_params_vals$avoid_ext_data <- species_params %>% 
      select(Avoidance, AvoidanceSD) %>% 
      rename(mean = Avoidance, sd = AvoidanceSD)
    
    species_params_vals$flap_glide <- species_params[which(species_params$Species == isolate(input$species_input)), "Flight"]
    
    #select Avoidance data and rename
    species_params_vals$prop_upwind_data <- species_params %>% 
      select(Prop_Upwind, Prop_UpwindSD) %>% 
      rename(mean = Prop_Upwind, sd = Prop_UpwindSD)
    
    #add nocturnal activity dataframe since not an input at this point
    species_params_vals$noct_act_data <- data.frame(mean = 1, sd = 0)
    
  })
  
  # species_params_vals$migrate_resident = NULL
  
  # observeEvent(input$migration_calc_type, {
  #   req(!is.null(input$species_input))
  #   species_params_vals$migrate_resident_reac
  #   #create a parameter for migrant vs. resident for crm calc switch
  #   if (input$species_input %in% c("Piping_Plover","Red_Knot", "Roseate_Tern") & input$migration_calc_type == "occup_model"){
  #     species_params_vals$migrate_resident = "resident"
  #   } else {  
  #     species_params_vals$migrate_resident = "migrant"
  #   }
  # })
  
  observe({
    req(!is.null(input$species_input))
    #create a parameter for migrant vs. resident for crm calc switch
    if (input$species_input %in% c("Piping_Plover","Red_Knot", "Roseate_Tern")) {
      species_params_vals$migrate_resident <- "migrant"
      if (input$migration_calc_type == "occup_model") {
        #model distribution type inputs
        species_params_vals$model_input_dist_type <- "trunc" #"last" - seems to produce worse state estimates, use first daily position (trunc)
      } else { #annex 6 calcs
        species_params_vals$model_input_dist_type <- "mean"
      }
    } else {  
      species_params_vals$migrate_resident <- "resident"
      
    }
  })
  
  output$species_data <-
    DT::renderDataTable({
      req(input$species_input)
      # bird_data <- read.csv("data/BirdData_020923.csv", header = T)
      species_params_defs <- c(
        "Parameter definitions",
        "Mean Proportion of birds that avoid turbines (0-1)",
        "Standard deviation of the avoidance rate",
        "Mean body length of the target species (m)",
        "Standard deviation of the species body length (m)",
        "Mean species wingspan length (m)",
        "Standard deviation of the species wingspan length (m)",
        "Mean species flight speed (m/s)",
        "Standard deviation of the species flight speed (m/s)",
        "Mean proportion of upwind flight for the species (0-1)",
        "Standard deviation of the proportion of upwind flight",
        # "Nocturnal flight activity",
        # "Standard deviation of the nocturnal flight activity",
        "Flight type, either flapping or gliding"
      )
      species_data_row <-
        as.data.frame(cbind(species_params_defs , t(bird_data[which(bird_data$Species == input$species_input),])))
      colnames(species_data_row) <-
        sub("_", " ", species_data_row[1, ])
      return(species_data_row[-1, , drop = FALSE])
    },
    # selection = list(mode = "single", selected = 1),
    options = list(dom = 't', 
                   paging = FALSE,
                   bSort=FALSE,
                   scrollY = "500px",
                   searching = FALSE))
  
  # Species flight height data ----------------------------------------------
  
  # flight height distributions

  #load flight height data
  #Required only for model Options 2 and 3, a data frame with bootstrap samples
  # of flight height distributions (FHD) of the species derived from general (country/
  #regional level) data. FHD provides relative frequency distribution of bird
  # flights at 1-+ -metre height bands, starting from sea surface. The first column
  # must be named as height, expressing the lower bound of the height band 
  # (thus it’s first element must be 0). Each remaining column should provide a bootstrap
  # sample of the proportion of bird flights at each height band, with no column
  # naming requirements.
  species_fhd <- eventReactive(input$species_input,{
    read.csv(paste0("data/", input$species_input,"_ht_dflt_v2.csv")) %>% 
      select(-Species) %>% 
      rename(height = Height_m) %>% 
      mutate(height = height - 1)
  })
  
  # Calculate the expected proportion of bird flights at collision risk height (i.e. at rotor height, between
  # bottom and top of the rotor) based on the bird’s flight height distribution.
  
  prop_chr_fhd <- eventReactive(c(input$species_input, wind_farm_vals), {
    prop_in_fhd_list <- list()
    for (t in 1:wind_farm_vals$num_turb_mods) {
      fhd_rotor_boot_prop <-
        apply(species_fhd()[, -1], 2, function(fhd_boot) {
          #proportion of flights at collision risk height derived from site survey, assumed to be Beta distributed.
          d_y <- suppressWarnings(stochLAB::get_fhd_rotor(
              hub_height = (wind_farm_vals$airgap_data$mean + wind_farm_vals$rtr_radius_data$mean),
              fhd = fhd_boot,
              rotor_radius = wind_farm_vals$rtr_radius_data$mean,
              tidal_offset = 0,
              yinc = 0.05
            ))
          stochLAB::get_prop_crh_fhd(d_y)
        })

      prop_in_fhd_list <-
        append(prop_in_fhd_list, list(data.frame(
          mean = mean(fhd_rotor_boot_prop),
          sd = sd(fhd_rotor_boot_prop)
        )))
    }
    return(prop_in_fhd_list)
  })
    
  #Summary table of flight height data for plotting and tables
  #  Created these for default/included data to speed up process, and otherwise summarize loaded flight height data
  # flt_ht_summary_react <- eventReactive(c(input$file_spp_param, input$species_input), {
  #   if(!is.null(input$file_spp_param)){
  #     flt_ht_boot_table <- spp_flt_ht_react()[[1]]
  #     n_cols <- ncol(flt_ht_boot_table)
  #     #summarize across all boot samples
  #     flt_ht_summary <- flt_ht_boot_table %>%
  #       dplyr::group_by(Height_m) %>%
  #       dplyr::rowwise() %>%
  #       dplyr::summarise(mean_prop = mean(c_across(3:n_cols-1)), min_prop = min(c_across(3:n_cols-1)), max_prop = max(c_across(3:n_cols-1)))
  #   } else {
  #     #return  pre-rendered flt height summaries depending on species
  #     # load(file=paste0("data/", input$species_input, "_ht_dflt_summary.RData"))
  #     load(file=paste0("data/", input$species_input, "_ht_dflt_v2_summary.RData"))
  #   }
  #   return(flt_ht_summary)
  # })
  
  flt_ht_summary_react <- reactive({
    load(file=paste0("data/", input$species_input, "_ht_dflt_v2_summary.RData"))
    return(flt_ht_summary)
  })
  
  #output table showing the min, max, mean flight heights from the bootstrap tables
  output$flt_ht_data <-
    DT::renderDataTable(
      datatable(
        flt_ht_summary_react() %>% 
          filter(Height_m >= input$slider_flt_ht[1] & Height_m <= input$slider_flt_ht[2]),
        selection = list(mode = "single", selected = 1),
        options = list(
          paging = FALSE,
          scrollY = "500px",
          searching = FALSE
        )))
  
  #flight height data plot
  output$flt_ht_plot <- renderPlot({
    ggplot(flt_ht_summary_react()) +
      geom_pointrange(aes(x = Height_m, y = mean_prop, ymin = min_prop, ymax = max_prop)) +
      geom_point(aes(x = Height_m, y = mean_prop), col = "red") +
      xlim(input$slider_flt_ht[1], input$slider_flt_ht[2]) +
      xlab("Flight height (m)") +
      ylab("Proportion") +
      coord_flip() +
      theme_bw()
    
  })
  

  # Band 2012 Annex species data inputs -------------------------------------
  
  #set data for annex 6 calcs
  # observeEvent(c(input$migration_calc_type, input$species_input), {
  observe( {
      
    req(input$migration_calc_type == "band_2012_annex_6" & length(input$species_input) > 0)
    
    # Set study map footnote text
    output$mean_prob_txt <- renderText("")
    
    annex6_migr_spp <- annex6_migr_data %>% 
      dplyr::filter(Species == input$species_input)
    
    
  })
  
  # annex6_migr_data <- read.csv("data/band_2012_annex6_migr_data_091123.csv")
  annex6_migr_data <- read.csv("data/Band_2012_REKN_PIPL_ROST_100223.csv")
  
  spp_corridor <- reactive({
    req(input$migration_calc_type == "band_2012_annex_6" & !is.null(input$species_input))
     if (input$species_input == "Red_Knot"){
       sf::read_sf("data/species_corridors/REKN_BandGeneralized.shp")
     } else if (input$species_input == "Piping_Plover") {
       sf::read_sf("data/species_corridors/PIPL_Band2012.shp")
     } else if (input$species_input == "Roseate_Tern") {
       sf::read_sf(sf::read_sf("data/species_corridors/ROST_Band2012_04142023.shp")) %>% 
         sf::st_transform(4326)
     }
  })
  
  annex6_vals <- reactiveValues(nearest_state_due_W = NULL, 
                                coastal_locations = NULL, 
                                mig_corridor_line = NULL, 
                                migr_front_km = 0, 
                                spp_popn_data = NULL)
  
  
  # Species population data -------------------------------------------
  
  #load population data
  popn_data <- read.csv("data/CountData_USFWS_20220919.csv", header = T)
  
  popn_data_vals <- reactiveValues()
  
  observe({
    req(input$migration_calc_type == "occup_model" & !is.null(input$species_input))
        
    #filter to species data and reformat to month, mean, and sd data columns
    spp_popn_SD <- popn_data %>%
      dplyr::filter(Species == input$species_input) %>%
      dplyr::select(ends_with("SD")) %>%
      t %>%
      as.data.frame()

    popn_data_vals$spp_popn_mean <- popn_data %>%
      dplyr::filter(Species == input$species_input) %>%
      dplyr::select(!ends_with("SD")) %>%
      dplyr::select(!"Species") %>%
      t %>%
      data.frame() %>%
      dplyr::rename("mean" = 1) %>%
      dplyr::mutate(month = month.abb, sd = spp_popn_SD[, 1]) %>%
      dplyr::select(month, mean, sd)

    # output$count_data <-
    #   DT::renderDataTable({
    #     count_tbl <- popn_data %>%
    #       dplyr::filter(Species==input$species_input) %>%
    #       dplyr::select(-Species)
    #     return(count_tbl)},
    #     options = list(dom = 't',
    #                    paging = FALSE,
    #                    bSort=FALSE,
    #                    scrollX = TRUE
    #     ),
    #     rownames= FALSE
    #   )
    
    spp_popn_SD_display <- popn_data %>%
      dplyr::filter(Species == input$species_input) %>%
      dplyr::select(ends_with("SD"))%>% 
      dplyr::mutate(Var = "SD") %>% 
      dplyr::select(c("Var", paste0(month.abb, "SD"))) %>%
      dplyr::rename_with(~ gsub("SD", "", .x, fixed = T))
    
    popn_data_vals$spp_popn_mean_display <- popn_data %>%
      dplyr::filter(Species == input$species_input) %>%
      dplyr::select(!ends_with("SD")) %>%
      dplyr::select(!c("Species")) %>% 
      dplyr::mutate(Var = "mean") %>% 
      dplyr::select(c("Var", month.abb)) %>% 
      rbind(spp_popn_SD_display)
    
    output$count_data <- DT::renderDataTable({
      popn_data_vals$spp_popn_mean_display},
      options = list(dom = 't', 
                     paging = FALSE,
                     bSort=FALSE, 
                     scrollX = TRUE
      ),
      rownames= FALSE
    )
    })
  
  
  
  observe({
    
    #annex 6 migrant data 
    req(input$migration_calc_type == "band_2012_annex_6" & !is.null(input$species_input) & !is.null(annex6_vals$filter_lcn))

    #filter the data and sum over the rows to get total accross states or regions
    annex6_vals$spp_popn_data <- annex6_migr_data %>% 
      dplyr::filter(Species == input$species_input & Location %in% annex6_vals$filter_lcn) %>%
      dplyr::select(-c(Location)) %>% 
      dplyr::group_by(Species) %>% 
      dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

    #filter to species data and reformat to month, mean, and sd data columns
    spp_popn_SD <- annex6_vals$spp_popn_data %>%
      # dplyr::filter(Species == input$species_input & migr_front_km == as.numeric(input$migr_front_width)) %>%
      # dplyr::filter(Species == input$species_input & Location %in% annex6_vals$coastal_locations) %>%
      dplyr::select(ends_with("SD")) %>%
      t %>%
      as.data.frame()
    
    popn_data_vals$spp_popn_mean <- annex6_vals$spp_popn_data %>%
      # dplyr::filter(Species == input$species_input) %>%
      dplyr::select(!ends_with("SD")) %>%
      as.data.frame() %>% 
      dplyr::select(!c(Species)) %>%
      t %>%
      data.frame() %>%
      dplyr::rename("mean" = 1) %>%
      dplyr::mutate(month = month.abb, sd = spp_popn_SD[, 1]) %>%
      dplyr::select(month, mean, sd)
    
    spp_popn_SD_display <- annex6_vals$spp_popn_data %>%
      # dplyr::filter(Species == input$species_input) %>%
      dplyr::select(ends_with("SD"))%>% 
      dplyr::mutate(Var = "SD") %>% 
      dplyr::select(c("Var", paste0(month.abb, "SD"))) %>%
      dplyr::rename_with(~ gsub("SD", "", .x, fixed = T))

    popn_data_vals$spp_popn_mean_display <- annex6_vals$spp_popn_data %>%
      # dplyr::filter(Species == input$species_input) %>%
      dplyr::select(!ends_with("SD")) %>%
      dplyr::select(!c("Species")) %>% 
      dplyr::mutate(Var = "mean") %>% 
      dplyr::select(c("Var", month.abb)) %>% 
      rbind(spp_popn_SD_display)

    output$count_data <- DT::renderDataTable({
      popn_data_vals$spp_popn_mean_display},
      options = list(dom = 't', 
                     paging = FALSE,
                     bSort=FALSE, 
                     scrollX = TRUE
      ),
      rownames= FALSE
    )
 
  })
  

  #Notes about assumptions and limitations for the three species for which population info is provided.
  spp_popn_notes <- reactive({
    req(input$species_input)
    Species <- input$species_input
    if (input$migration_calc_type == "occup_model"){
      if (Species=="Piping_Plover"){
        return(c("Entire Atlantic coast population could be present in area during months listed.",
                 "Occurrence through October to include birds stopping over in mid-Atlantic (e.g. North Carolina). Number of birds still present in Atlantic likely lower.",     
                 "Estimate of HY fledges, uses the 20-year (2002 - 2021) average productivity (unweighted)."))
      }
      else if (Species=="Red_Knot"){
        return(c("All pass through in spring - #s consistent w/Lyons et al super-population estimate for 2020 in DE Bay: 40,444 (95 perc. credible interval: 33,627–49,966).",
                 "Winter population estimates represent the total # of adults and sub-adults (in general); they do not include hatch-year (HY) birds in the fall.",
                 "Southern and northern wintering birds could be present during July - Sept.",
                 "Only northern wintering birds could be present during Oct - Nov.",
                 "Only southeast US and Caribbean birds could be present during Dec.",
                 "Birds from western Gulf population are excluded from totals in Atlantic region due to lack of information on extent to which they use the Atlantic region.",
                 "Numbers do not include HY birds in fall.",
                 "Dec number coming from Lyons et al 2017. Just includes SE US Birds, not Caribbean.",
                 "Issues with double counting addressed because birds may be present in different areas of Atlantic region for weeks to months."))
      }
      else if (Species=="Roseate_Tern"){
        return(c("Entire NW Atlantic pop could be present in area during months listed.",
                 "Average of most recent (2018 and 2019) productivity data from three largest colonies (representing >90 perc. of population) representative of entire population.",
                 "Fledging and post-breeding dispersal period occurs from July through Sept.",
                 "Numbers of non-breeding adults are not included.",
                 "Does not include non-breeding 1 and 2 year old birds that return but do not breed.",
                 "From Gochfeld and Burger (2020): Northeastern birds first arrive at Nantucket and Martha's Vineyard, MA, in large flocks, then disperse north as well as west. They arrive 26 Apr-20 May at Bird I., MA (Nisbet 1980, Nisbet 1981b, Nisbet 1989b), slightly later at Falkner I., CT, and Great Gull I., NY.",
                 "From Gochfeld and Burger (2020): Apparently all birds migrate directly from the staging area around Cape Cod across the w. North Atlantic to the West Indies (Nisbet 1984, C. Mostello). Very small numbers occur at sea off N. Carolina from late Aug to late Sep, with a peak in early Sep; the latest date was 28 Oct (D. Lee)."))
      } else {return("")} 
    } else { #annex 6 text
      
      if (Species=="Piping_Plover"){
        return(c("Population data are from 2021 (USFWS 2021a) and exclude an unknown (but likely small) number of nonbreeding birds.",
                 "The Southern recovery unit population is excluded.",
                 "The SB total includes YOY, calculated as the unweighted mean 20-year productivity rates (2002 - 2021) times the 2021 breeding pair estimate for each state within the Eastern Canada, New England, and NY-NJ recovery units.",
                 "The eastern edge of the migration corridor runs southwest parallel to the general orientation of the coast to account for major migration staging areas in North Carolina. The eastern edge of the corridor south of Cape Hatteras is also constrained westward to account for much larger numbers of piping plovers wintering in the western Bahamas."))
      }
      else if (Species=="Red_Knot"){
        return(c(
          "Population data are from Table 2, above.",
          "Birds from the Western recovery unit population are sometimes documented on the Atlantic coast. However, available tracking and resighting data show that the prevailing migration corridor for these birds is overland across the mid-continent (Perkins 2023, USFWS 2021b, USFWS 2014). On this basis, birds from the Western recovery unit are excluded from this analysis.",
          "In many years, a percentage of northbound birds do not depart the mid-Atlantic until early June. But for the purposes of this analysis, we attribute them all to May.",
          "Some juveniles and nonbreeding adults remain south of the migration front, others cross the migration front once in spring and spend the breeding season just south of the breeding grounds, while still others may remain resident in the mid-Atlantic for prolonged periods and may cross the migration front multiple times. We have no estimate of the total number of nonbreeding adults in a typical year, or their distribution across the species nonbreeding range. However, we do estimate the total number of juveniles. Modeling by Schwarzer (2011) found that the Florida population was stable at around 8.75 percent juveniles among wintering birds, and available data suggest the three populations considered in this analysis are currently stable (USFWS 2021b). Thus, we assume 8.75 percent of the total wintering birds are juveniles (i.e., of the 59,269 total birds, we assume 5,186 are juveniles.) We have little information on the distribution of juveniles across the species’ range during any month. In light of data gaps, we assume all breeding adults, nonbreeding adults, and juveniles cross the migration front twice per year.",
          "The SB total includes YOY, calculated as 1 chick per pair. Number of pairs is calculated as [the total wintering population (59,269) minus juveniles (5,186)] divided by 2. We have no way to estimate nonbreeding adults, so we include them with breeding adults, then attempt to compensate by using a reproductive rate of 1 chick per pair, below the range estimated by Wilson and Morrison (2018) as needed for a stable population."))
      }
      else if (Species=="Roseate_Tern"){
        return(
          c(
            "Migration numbers were generated based on 2021 breeding population numbers and productivity rates from the US and Canada.",
            "Spring migration totals were calculated as the number of breeding pairs in each region multiplied by 2 adults per breeding pair.",
            "Fall migration totals included all adults from spring migration plus the approximate number of YOY.",
            "YOY totals were calculated by multiplying the number of breeding pairs in the US and Canada by the average productivity of these pairs (approximately 1 YOY per pair).",
            "Migration months were determined based on peak migration during the spring and fall migration seasons, as reported by Gochfeld and Burger (2020).",
            "Number of spring and fall migrants were then assumed to be divided evenly across migration months."
          )
        )
      } else {return("")} 
    }

    
  })

  
  #create count notes for rendering.
  output$popn_notes <- renderText(
    c("<h5>Population data assumptions/limitations:</h5>",
      paste0(1:length(spp_popn_notes()),") ", spp_popn_notes(), " <br>"))
  )
  
  
  # Wind farm data ----------------------------------------------------------
  
  observeEvent(input$file_wf_param, {
    Sys.sleep(2)  #wait until species data has rendered and then switch also helps with pre-rendering maps, etc.
    updateTabItems(session, inputId = "tabsetpan", selected = "wind_farm_panel")})
  
  load(file = file.path(data_dir, "BOEM_halfdeg_grid_sf.RData")) #BOEM_halfdeg_grid_sf
  
  wind_farm_df <- eventReactive(input$file_wf_param, {
    infile <- input$file_wf_param
    if (is.null(infile)) {
      return(NULL)
    }
    suppressWarnings(read.csv(infile$datapath, header=T))
  })
  
  #reactive value to hold the selected grid cell for calculating
  windfarm_loc <- reactiveValues(
    center = NULL,
    Latitude = NULL, 
    Longitude = NULL,
    cell_sf = NULL,
    )
  

  #perform lat/long checks prior to creating spatial data and create a spatial sf object holding the wind farm location and data
  observeEvent(wind_farm_df(), {
    
    if(length(which(colnames(isolate(wind_farm_df()))=="Latitude")) > 0 & length(which(colnames(isolate(wind_farm_df())) == "Longitude")) > 0){
      windfarm_lats <- isolate(wind_farm_df())[,c("Latitude", "Longitude")]
      #check to see if multiple lat/long value pairs for inputs - not allowed in this version
      if(nrow(windfarm_lats) >= 1){
        if(nrow(unique(windfarm_lats)) >= 2){
          showModal(modalDialog(
            title = "Check values for Latitude and Longitude",
            footer = modalButton("OK"),
            paste("Latitude and/or Longitude are not the same in the different wind farm option parameters.
              Only one location is allowed per model run.
              Please upload wind farm options for the same location.")
          ))
        } else {
          windfarm_loc$Latitude = isolate(wind_farm_df())[1,"Latitude"]
          windfarm_loc$Longitude = isolate(wind_farm_df())[1,"Longitude"]
          
          
          windfarm_loc$center <- sf::st_as_sf(isolate(wind_farm_df())[1,], coords=c("Longitude", "Latitude" ), crs = sf::st_crs(4326))

          # Gets cell values from lat/longs provided to provide info on low confidence areas as well as grid cell area in sq. km 
          #intersect with the model grid polygon
          windfarm_loc$cell_sf <- BOEM_halfdeg_grid_sf[unlist(sf::st_intersects(windfarm_loc$center, BOEM_halfdeg_grid_sf)), ]

          if(nrow(windfarm_loc$cell_sf)==0){
            showModal(modalDialog(
              title = "Check values for Latitude and Longitude",
              footer = modalButton("OK"),
              paste("Latitude and/or Longitude values fall outside the geographic extent of this tool.
              Please upload appropriate values with wind farm data or consult documentation for more information on geographic scope.")
            ))
          } 
        } 
      } 
    } else {  #missing lat or long fields
      showModal(modalDialog(
        title = "Check wind farm data for Latitude and Longitude",
        footer = modalButton("OK"),
        paste("Latitude and/or Longitude are missing.
              Please upload appropriate values with wind farm data.")
      ))
    }
  })
  
  #set the migration corridor
  observe({
    req(!is.null(windfarm_loc$center) & !is.null(input$species_input))
    
    #create the annex6 migration corridor line and split at migration corridor
    #create an extended line for processsing corridor and getting closest state info
    EW_line <- create_EW_line(sf::st_coordinates(windfarm_loc$center), length_km = 5000)
    
    annex6_vals$mig_corridor_line <- EW_line %>% 
      lwgeom::st_split(spp_corridor()) %>% st_intersection(spp_corridor()) #%>% st_cast()
    
    # corr = spp_corridor()
    # plot(corr$geometry)
    # plot(annex6_vals$mig_corridor_line, add= T, col = "blue")

    #calculate migration corridor width at the wind farm in km, convert to equal area proj first
    annex6_vals$migr_front_km <- round(as.numeric(sf::st_length(st_transform(annex6_vals$mig_corridor_line, crs = 9822))) / 1000, 0) #Albers equal area conic proj
    output$migr_front_width <- renderText(paste("Migratory front width:",  annex6_vals$migr_front_km, "km"))
    # sf::write_sf(annex6_vals$mig_corridor_line, "mig_corridor_line.shp", )
    
    # Split the line on the state boundaries to get a series of split lines for which the centroid of these allow us to determine the closest
    # section to the start of the line and thus the section of state closest west of the wind farm. 
    # line_seg_centers <- lwgeom::st_split(annex6_vals$mig_corridor_line, state_boundaries_wgs84) %>% st_cast() %>% st_centroid()
    line_seg_centers <- lwgeom::st_split(EW_line, state_boundaries_wgs84) %>% st_cast() %>% st_centroid()
    
    nearest_seg_center <- line_seg_centers[(nngeo::st_nn(windfarm_loc$center, line_seg_centers)[[1]]), ]
    # Error in wk_handle.wk_wkb(wkb, s2_geography_writer(oriented = oriented,  : 
    # Loop 3 is not valid: Edge 30 h
    # https://github.com/r-spatial/sf/issues/1902
    sf_use_s2(FALSE) #Spherical geometry (s2) switched off
    
    annex6_vals$nearest_state_due_W <- state_boundaries_wgs84[st_intersects(state_boundaries_wgs84, nearest_seg_center, sparse = F), ]$NAME
    output$nearest_state <- renderText(paste("Nearest coastal state due west:",  annex6_vals$nearest_state_due_W))
    
    #now derive list of coastal locations including and north of the state selected
    annex6_vals$coastal_locations <- get(paste0(state.abb[grep(annex6_vals$nearest_state_due_W, state.name)], "_states_north"))
    
    #set the filter location for population data
    if (input$species_input == "Red_Knot") {
      annex6_vals$filter_lcn <- c("All")
    } else if (input$species_input == "Roseate_Tern") {
      annex6_vals$filter_lcn <- c("All")
    } else if(input$species_input == "Piping_Plover") {
      annex6_vals$filter_lcn <- annex6_vals$coastal_locations
    }
    
  })
  
  #render the wind farm data to the tables on the wind farm data tab
  output$wind_farm_data1 <- DT::renderDataTable({
    req(wind_farm_df())
    wf_t <- isolate(wind_farm_df()) %>%
      dplyr::select(-contains("Op"))
    wf_t <- as.data.frame(t(wf_t))
    colnames(wf_t) <- paste("Run", wf_t[1, ])
    #add parameter defs
    wf1_params_defs <- c("Parameter definitions",
                         "The number of turbines in the wind farm array",
                         "The turbine model used in the analysis",
                         "Number of blades for the turbine model",
                         "Rotor radius (hub to blade tip; m)",
                         "Standard deviation of rotor radius (m)",
                         "Air gap (m)",
                         "Standard deviation of air gap (m)",
                         "Chord width of blade (m)",
                         "Standard deviation of blade width (m)",
                         "Mean wind speed (m/s) between cut-in and cut-out speeds or wind speed rating of the turbine model if wind speed is unavailable",
                         "Standard deviation of wind speed or wind speed rating (m/s)",
                         "Pitch angle of blades (degrees relative to rotor plane)",
                         "Standard deviation of pitch angle of blades",
                         "Wind farm width (km)",
                         "Latitude (decimal degrees)",
                         "Longitude (decimal degrees)")
    wf_t <- cbind("Parameter definitions" = wf1_params_defs, wf_t)
    return(wf_t[-1,, drop = FALSE])
  },
  options = list(
    dom = 't',
    paging = FALSE,
    bSort=FALSE
  ))

  # wind farm parameters
  # wind_farm_params_react <- eventReactive(input$file_wf_param, {
  #   suppressWarnings(read.table(input$file_wf_param$datapath, header=TRUE, sep = ","))
  # })

  #show the Wind Farm operational data as as table for QA/QC
  output$ops_data <-
    DT::renderDataTable({
      req(wind_farm_df())
      ops_data <- isolate(wind_farm_df()) %>%
        dplyr::select(Run, matches("Op", ignore.case=F)) %>%
        t()
      colnames(ops_data) <- paste("Run", ops_data[1, ])
      #add parameter defs
      op_defs <- c("defs", (rep(c("Wind availability (maximum amount of time turbines can be operational/month).",
                                  "Mean time that turbines will not be operational (“down time”)",
                                  "Standard deviation of mean operational time"),12)))
      ops_data <- cbind("Parameter definitions"=op_defs, ops_data)
      return(ops_data[-1, , drop = FALSE])},
      options = list(dom = 't',
                     paging = FALSE,
                     bSort=FALSE)
    )
  
  #format and store windfarm data for stochCRM processing
  wind_farm_vals <- reactiveValues()
  
  observeEvent(wind_farm_df(), {

    wind_farm_vals$num_turb_mods <- nrow(wind_farm_df())

    #select airgap data and rename
    wind_farm_vals$airgap_data <- wind_farm_df() %>%
      select(HubHeightAdd_m, HubHeightAddSD_m) %>%
      rename(mean = HubHeightAdd_m, sd = HubHeightAddSD_m)

    #select rotor radius data and rename
    wind_farm_vals$rtr_radius_data <- wind_farm_df() %>%
      select(RotorRadius_m, RotorRadiusSD_m) %>%
      rename(mean = RotorRadius_m, sd = RotorRadiusSD_m)

    # select blade width data and rename
    wind_farm_vals$bld_width_data <- wind_farm_df() %>%
      select(BladeWidth_m, BladeWidthSD_m) %>%
      rename(mean = BladeWidth_m, sd = BladeWidthSD_m)

    # select blade pitch data and rename
    wind_farm_vals$bld_pitch_data <- wind_farm_df() %>%
      select(Pitch, PitchSD) %>%
      rename(mean = Pitch, sd = PitchSD)

    # select rotation speed data and rename
    wind_farm_vals$wind_speed_data <- wind_farm_df() %>%
      select(WindSpeed_mps, WindSpeedSD_mps) %>%
      rename(mean = WindSpeed_mps, sd = WindSpeedSD_mps)

    # select rotation speed data and rename
    # a single row data frame with
    # columns mean and sd, the mean and standard deviation of the operational rotation
    # speed, in revolutions per minute. Assumed to follow a tnorm-lw0 distribution.

    # get rotor speed according to diameter and mean and sd wind speeds
    # May want to check these wind speed/rotational speed data to ask for more accurate data from developers
    # more in line with what stoch_crm requires - should be discussed further
    # 5.5 is the TSR; 60 is for converting radians/s to rpm
    # browser()
    wind_farm_vals$rtn_speed_data <- data.frame()
    for (t in 1:wind_farm_vals$num_turb_mods){
      RotorSpeed_rpm_mean <-
        (5.5 *  wind_farm_vals$wind_speed_data$mean[t] * 60) / (pi * 2 *  wind_farm_vals$rtr_radius_data$mean[t])
      RotorSpeed_rpm_sd <-
        (5.5 *  wind_farm_vals$wind_speed_data$sd[t] * 60) / (pi * 2 *  wind_farm_vals$rtr_radius_data$sd[t])
      RotorSpeed_rpm_sd[is.nan(RotorSpeed_rpm_sd)] = 0
      RotorSpeed_rpm_sd[is.infinite(RotorSpeed_rpm_sd)] = 0
      wind_farm_vals$rtn_speed_data <- rbind(wind_farm_vals$rtn_speed_data, data.frame(mean = RotorSpeed_rpm_mean, sd = RotorSpeed_rpm_sd))
    }

    # select operational data
    # A data frame with the monthly estimates of operational wind availability. It
    # must contain the columns:
    # - month, (unique) month names,
    # - pctg, the percentage of time wind conditions allow for turbine operation per month.
    pctg = wind_farm_df() %>%
      select(ends_with("Op")) %>%
      t

    wind_farm_vals$trb_wind_avbl_data <- list()
    for(i in 1:ncol(pctg)){
      wind_farm_vals$trb_wind_avbl_data <- append(wind_farm_vals$trb_wind_avbl_data,
                                                  list(data.frame(month = month.abb, pctg = pctg[,i])))
    }

    # select downtime data
    # A data frame with monthly estimates of maintenance downtime, assumed to
    # follow a tnorm-lw0 distribution.
    wind_farm_vals$trb_downtime_data <- list()

    down_mean = wind_farm_df() %>%
      select(ends_with("OpMean")) %>%
      t

    down_sd = wind_farm_df() %>%
      select(ends_with("OpSD")) %>%
      t
    for(i in 1:ncol(down_mean)){
      wind_farm_vals$trb_downtime_data <- append(wind_farm_vals$trb_downtime_data,
                                               list(data.frame(month = month.abb, mean = down_mean[, i], sd = down_sd[, i])))
    }
  })
  
  bladeIcon <- makeIcon(
    iconUrl = "www/outline_wind_power_black_36dp.png",
    iconWidth = 36, iconHeight = 36,
    iconAnchorX = 0, iconAnchorY = 36,
  )
  
  # Species occupancy model processing and mapping ----------------------------------------------------------
  
  #Load movement model file for the grid cell in which the WF lies
  #Depends on the species and the model type which varies with migratory status
  #read from zipped file to reduce number of files to upload to R Shiny and handle generally
  
  bird_flux_vals <- reactiveValues()
  bird_flux_vals$num_birds_cell_perday <- list()

  bird_dens_from_model <- eventReactive(input$run, { 
    
    req(nrow(popn_data_vals$spp_popn_mean) > 0 & input$migration_calc_type == "occup_model" & !is.null(input$species_input))

    if (length(SpeciesLabels$name_underscore == input$species_input) > 0) {
      spp_code <- SpeciesLabels[which(SpeciesLabels$name_underscore == input$species_input), "spp_code"]
    }
    
    bird_flux_vals$trans_prop <- windfarm_loc$cell_sf[[paste0(spp_code, "_prop_trans")]]
    
    movement_model_post <- read.csv(unz(paste0("data/movements/", input$species_input, "_movements_", species_params_vals$model_input_dist_type, ".zip"), 
                                                                            paste0("MovementBaked_", input$species_input, "_", species_params_vals$model_input_dist_type, windfarm_loc$cell_sf$id, ".csv")), header=T)
    
    bird_dens <- data.frame(i=1:length(movement_model_post[,1]))
    num_birds_cell_perday <- data.frame(i=1:length(movement_model_post[,1]))

    #Sample counts each month as tnormal and combine with movment draw - resample 1000 times to create
    #a resample matrix of monthly density (birds/km) which can be resampled in stoch_crm
    popn_df <- popn_data_vals$spp_popn_mean
    
    #loop through months for processing bird density data
    for (i in 1:12){
      #calc samples from mean/sd of population
      m <- month.abb[[i]]
      popn_samples <- as.numeric(stochLAB::rtnorm_dmp(1000, mean = popn_df[i, "mean"], sd = popn_df[i, "sd"], lower = 0))

      #Multiply population sample * the posterior of the movement model * the proportion of transient behavior / grid cell mean width in the model cell 
      #to get a sample for drawing for crm calcs
      num_birds_cell_perday[[m]] <- popn_samples * movement_model_post[[m]] * bird_flux_vals$trans_prop / 
        monthDays(as.Date(paste0("01/", str_pad(i, width = 2, side = "left", pad = "0"), "/2023")))
      
      bird_dens[[m]] <- popn_samples * movement_model_post[[m]] * bird_flux_vals$trans_prop / sqrt(windfarm_loc$cell_sf$area)  #birds/km
      
    }
    
    bird_flux_vals$num_birds_cell_perday <- append(bird_flux_vals$num_birds_cell_perday, list(as.data.frame(num_birds_cell_perday)))

    #remove index field for input to model
    bird_dens <- bird_dens %>% 
      select(-i)
    return(bird_dens)
    
  }) # end bird_dens_from_model
  
  bird_dens_from_annex6 <- eventReactive(input$run, {
    
    req(nrow(popn_data_vals$spp_popn_mean) > 0 & input$migration_calc_type == "band_2012_annex_6" & 
          !is.null(input$species_input) & annex6_vals$migr_front_km > 0)
    
    bird_dens <- data.frame(i=1:1000)

    #Sample counts each month as tnormal and combine with movment draw - resample 1000 times to create
    #a resample matrix of monthly density (birds/km) which can be resampled in stoch_crm
    popn_df <- popn_data_vals$spp_popn_mean

    #loop through months for processing bird density data
    for (i in 1:12){
      #calc samples from mean/sd of population
      m <- month.abb[[i]]
      
      #When mean 0 produce NA as equal to lower bound of 0, make lower bound a very small negative
      popn_samples <- as.numeric(stochLAB::rtnorm_dmp(1000, mean = popn_df[i, "mean"], sd = popn_df[i, "sd"], lower = -0.00000000000001))
      
      bird_dens[[m]] <- popn_samples / as.numeric(annex6_vals$migr_front_km)  #birds/km
      
    }
    
    #remove index field for input to model
    bird_dens <- bird_dens %>% 
      select(-i)
    return(bird_dens)
  })
  
  
  #add appropriate species model output to map in leaflet when species chosen
  #create the color palette functions for coloring
  spp_move_data <- reactive({
      
    req(input$migration_calc_type == "occup_model" & length(input$species_input) > 0)
    
    species <- input$species_input
    grid_layer <-  paste0(species, "_monthly_prob_BOEM_half_deg_", species_params_vals$model_input_dist_type)
    load(file.path(data_dir, paste0(species, "_monthly_prob_BOEM_half_deg_", species_params_vals$model_input_dist_type,".RData")))
    assign(grid_layer, spp_monthly_prob_BOEM_half_deg)
  })
  
  meanpal <- reactive({
    #issue with zero inflated bins - need to generate non-zero quants
    #changed from using cumulative daily to mean monthly (non-culative daily) occupancy
    non_zero_mean <- spp_move_data()$mean_daily_allmonths[spp_move_data()$mean_daily_allmonths>0]
    mean_bins <- c(0, quantile(non_zero_mean, probs=seq(0,1,1/7)))
    colorBin("YlOrRd", spp_move_data()$mean_daily_allmonths, bins = mean_bins)
  })
  
  #modify map to add occupancy data if the model is used for calculation
  observe({
      
    #render the map with the lat/longs given in the study area map panel
    req(input$migration_calc_type == "occup_model" & nrow(spp_move_data())>0)

    pal <- meanpal()
    
    # CIpal <- eventReactive(spp_move_data(), {
    #   #issue with zero inflated bins - need to generate non-zero quants
    #   non_zero_CIrange <- spp_move_data()$mean_CI_range_daily_allmonths[spp_move_data()$mean_CI_range_daily_allmonths>0]
    #   CI_bins <- c(0, quantile(non_zero_CIrange, probs=seq(0,1,1/7)))
    #   colorBin("Purples", spp_move_data()$mean_CI_range_daily_allmonths, bins = CI_bins)
    # })
    output$studymap <- renderLeaflet({
      leaflet(options = leafletOptions(preferCanvas = T, tolerance = 1)) %>%
        addEsriBasemapLayer(esriBasemapLayers$Oceans, autoLabels = TRUE) %>%
        addEsriFeatureLayer(
          #add BOEM renewable lease areas as a WFS
          url = "https://services7.arcgis.com/G5Ma95RzqJRPKsWL/ArcGIS/rest/services/Wind_Lease_Boundaries__BOEM_/FeatureServer/0",
          weight = 1, fill=FALSE, color = "yellow",#"#808080", #gray
          layerId = "BOEM_wind_leases",
          group = "BOEM wind leases") %>%
        addEsriFeatureLayer(
          #add BOEM WPAs as a WFS
          url = "https://services7.arcgis.com/G5Ma95RzqJRPKsWL/arcgis/rest/services/Wind_Planning_Area_Boundaries__BOEM_/FeatureServer/0",
          weight = 1, fill=FALSE, color = "green", #C0C0C0", #silver
          layerId = "BOEM_wpa",
          group = "BOEM wind planning areas") %>%
        # MOTUS antenna data
        addMarkers(
          data = wind_farm_df(),
          lat = ~ Latitude,
          lng = ~ Longitude,
          icon = bladeIcon,
          popup =
            paste0(
              "Run: ",
              wind_farm_df()$Run,
              "<br/>Number of turbines: ",
              wind_farm_df()$Num_Turbines,
              "<br/>Turbine model: ",
              wind_farm_df()$TurbineModel,
              "<br/>Rotor radius (m): ",
              wind_farm_df()$RotorRadius_m,
              "<br/>Air gap (m): ",
              wind_farm_df()$HubHeightAdd_m,
              "<br/>BladeWidth (m): ",
              wind_farm_df()$BladeWidth_m,
              "<br/>Pitch: ",
              wind_farm_df()$Pitch,
              "<br/>Wind farm width (km): ",
              wind_farm_df()$WFWidth_km
            ),
          group = "Wind farm"
        ) %>%
        setView(lat = mean(wind_farm_df()$Latitude), lng = mean(wind_farm_df()$Longitude), zoom = 6) %>%
        addPolygons(data=spp_move_data(), weight = 1, fillOpacity = 0.5, opacity = 1,
                    color = ~pal(mean_daily_allmonths),
                    highlightOptions = highlightOptions(color = "white", weight = 2),
                    label = ~formatC(mean_daily_allmonths, digits = 2, format = "g"),
                    group="Occup. prob.") %>% 
        #custom legend formatting labels
        addLegend(data=spp_move_data(), pal = pal, values = ~mean_daily_allmonths,
                  title = "\u00B9Mean occup. prob.",
                  labFormat = function(type, cuts, p) {
                    n = length(cuts)
                    paste0(formatC(cuts[-n], digits = 2, format = "g"), " &ndash; ", formatC(cuts[-1], digits = 2, format = "g"))
                  },
                  position = "bottomright", group="Occup. prob.") %>% 
        #Layers control
        addLayersControl(
          overlayGroups = c("Wind farm", "BOEM wind leases", "BOEM wind planning areas", "Occup. prob."), #, "CI range"),
          position = "topright", #"topleft",
          options = layersControlOptions(collapsed = TRUE)
        )

    })
  })
  
  
  observe({
  #render the map with the lat/longs given in the study area map panel
    req(input$migration_calc_type == "band_2012_annex_6" & !is.null(spp_corridor()))
    
    output$studymap <- renderLeaflet({
      leaflet(options = leafletOptions(preferCanvas = T, tolerance = 1)) %>%
        addEsriBasemapLayer(esriBasemapLayers$Oceans, autoLabels = TRUE) %>%
        addEsriFeatureLayer(
          #add BOEM renewable lease areas as a WFS
          url = "https://services7.arcgis.com/G5Ma95RzqJRPKsWL/ArcGIS/rest/services/Wind_Lease_Boundaries__BOEM_/FeatureServer/0",
          weight = 1, fill=FALSE, color = "yellow",#"#808080", #gray
          layerId = "BOEM_wind_leases",
          group = "BOEM wind leases") %>%
        addEsriFeatureLayer(
          #add BOEM WPAs as a WFS
          url = "https://services7.arcgis.com/G5Ma95RzqJRPKsWL/arcgis/rest/services/Wind_Planning_Area_Boundaries__BOEM_/FeatureServer/0",
          weight = 1, fill=FALSE, color = "green", #C0C0C0", #silver
          layerId = "BOEM_wpa",
          group = "BOEM wind planning areas") %>%
        # MOTUS antenna data
        addMarkers(
          data = wind_farm_df(),
          lat = ~ Latitude,
          lng = ~ Longitude,
          icon = bladeIcon,
          popup =
            paste0(
              "Run: ",
              wind_farm_df()$Run,
              "<br/>Number of turbines: ",
              wind_farm_df()$Num_Turbines,
              "<br/>Turbine model: ",
              wind_farm_df()$TurbineModel,
              "<br/>Rotor radius (m): ",
              wind_farm_df()$RotorRadius_m,
              "<br/>Air gap (m): ",
              wind_farm_df()$HubHeightAdd_m,
              "<br/>BladeWidth (m): ",
              wind_farm_df()$BladeWidth_m,
              "<br/>Pitch: ",
              wind_farm_df()$Pitch,
              "<br/>Wind farm width (km): ",
              wind_farm_df()$WFWidth_km
            ),
          group = "Wind farm"
        ) %>%
        addPolygons(data = spp_corridor(), opacity = 0.3, color = "red") %>% 
        addPolylines(data = annex6_vals$mig_corridor_line, color = "orange", weight = 2, dashArray = "5, 10") %>% 
        setView(lat = mean(wind_farm_df()$Latitude), lng = mean(wind_farm_df()$Longitude), zoom = 6) %>%
        #Layers control
        addLayersControl(
          overlayGroups = c("Wind farm", "BOEM wind leases", "BOEM wind planning areas"), #, "Occup. prob."), #, "CI range"),
          position = "topright", #"topleft",
          options = layersControlOptions(collapsed = TRUE)
        )
    })
  })


  # Additional data checks prior to running SCRAM --------------------------------------
  
  # an object that indexes whether at least one file has been uploaded for both species and turbine data
  bothdataup <- reactiveVal()
  bothdataup  <- eventReactive(c(input$"file_spp_param", input$"file_wf_param"), {length(input$file_spp_param$datapath) > 0 & length(input$file_wf_param$datapath) > 0})
  
  # dialog box that alerts the user that some data was uploaded, but not every data type was provided
  observeEvent(c(input$file_spp_param, input$species_input), {
    if(!is.null(input$file_spp_param)){
      if(length(which(input$species_input=="Other")) == 0&
         (length(which(sapply(1:length(input$file_spp_param$datapath), function(x) length(suppressWarnings(read.table(input$file_spp_param[[x,"datapath"]], header=TRUE, sep = ","))[1,]))
                       == 10)) == 0|
          length(which(sapply(1:length(input$file_spp_param$datapath), function(x) length(suppressWarnings(read.table(input$file_spp_param[[x,"datapath"]], header=TRUE, sep = ","))[1,]))
                       == 25)) == 0|
          length(which(sapply(1:length(input$file_spp_param$datapath), function(x) length(suppressWarnings(read.table(input$file_spp_param[[x,"datapath"]], header=TRUE, sep = ","))[1,]))
                       == 1002)) == 0)){
        showModal(modalDialog(
          title = "Some data are missing",
          footer = modalButton("OK"),
          paste("Species data, flight heights, or survey data are not provided; default values will be used for any datasets that are not provided.")
        ))
      }}
  })
  
  # dialog box for when user chooses "Other" species option
  observeEvent(c(input$file_spp_param, input$species_input), {
    {if(length(input$file_spp_param$datapath) < 3&length(which(input$species_input=="Other")) == 1&!is.null(input$file_spp_param)){
      showModal(modalDialog(
        title = "Cannot run CRM",
        footer = modalButton("OK"),
        paste("Please provide species data, flight heights, and movement data when choosing the 'Use your own species' option. See documentation ('Start here') for more information.")
      ))
    }}
  })
  

  # Main crm execution code -------------------------

  # source("scripts/GSA.R")

  stoch_SCRAM_out <- reactiveVal()
  
  run_times <- reactiveValues()
  filenames <- reactiveValues()
  
  # primary function that called computational script (using a promise)
  observeEvent(input$run, {

    run_times$start <- Sys.time()
    output$run_start_txt <- renderText(paste0("The model run was started at: ", strftime(run_times$start, "%Y-%m-%d %H:%M:%S %Z", tz = "America/New_York")))
    
    CRSpecies <- input$species_input
    
    if (input$migration_calc_type == "occup_model") {
      bird_density_vals <- bird_dens_from_model()
    } else { #annex 6 calcs instead
      bird_density_vals <- bird_dens_from_annex6()
    }

    model_out <- list()

    # Start of the species loop -----------------------------------------------    
    for (s in 1:length(CRSpecies)){
      for (t in 1:wind_farm_vals$num_turb_mods){
        
        cTurbModel <- paste0("turbModel", wind_farm_df()[t, "TurbineModel_MW"])
        model_out[[CRSpecies[s]]][[cTurbModel]] <- 
          # bind the inputs and the outputs together into a list
          c(            
            #add input parameter, output and sample param draws summaries
            # inputs = list(
              species = CRSpecies[s],
              migration_calc_type = input$migration_calc_type,
              # migr_front_type = input$migr_front_opt,
              migrant_states = paste(annex6_vals$filter_lcn, collapse = ", "),
              migr_front_width = annex6_vals$migr_front_km,
              model_option = crm_option_radio_react(),  
              model_type = species_params_vals$model_input_dist_type,
              model_iterations = iter_slider_react(),
              turbine_model = wind_farm_df()[t, "TurbineModel_MW"],
              num_turbines = wind_farm_df()[t, "Num_Turbines"],
              windfarm_width_km = wind_farm_df()[t, "WFWidth_km"],
              lat = wind_farm_df()[t, "Latitude"],
              long = wind_farm_df()[t, "Longitude"],
              # ),
            stoch_SCRAM(
              model_options = c(crm_option_radio_react()),
              n_iter = iter_slider_react(),
              migrant_resident = species_params_vals$migrate_resident,
              flt_speed_pars = species_params_vals$flt_spd_data,
              body_lt_pars = species_params_vals$body_lt_data,
              wing_span_pars = species_params_vals$wing_span_data,
              avoid_bsc_pars = species_params_vals$avoid_ext_data,
              avoid_ext_pars = species_params_vals$avoid_ext_data,
              noct_act_pars = species_params_vals$noct_act_data,
              prop_crh_pars = isolate(prop_chr_fhd())[[t]],
              #Required only for model Option 1, a single row data frame with columns mean and sd.
              bird_dens_opt = "resample",
              bird_dens_dt = bird_density_vals,
              flight_type = species_params_vals$flap_glide,
              prop_upwind = species_params_vals$prop_upwind_data$mean,
              gen_fhd_boots = isolate(species_fhd()),
              site_fhd_boots = NULL,
              n_blades =  wind_farm_df()$Num_Blades[t],
              air_gap_pars = wind_farm_vals$airgap_data[t, ],
              rtr_radius_pars = wind_farm_vals$rtr_radius_data[t, ],
              bld_width_pars = wind_farm_vals$bld_width_data[t, ],
              bld_chord_prf = chord_prof_5MW,
              rtn_pitch_opt = "probDist",
              bld_pitch_pars = wind_farm_vals$bld_pitch_data[t, ],
              #Only required if rtn_pitch_opt = "probDist"
              rtn_speed_pars = wind_farm_vals$rtn_speed_data[t, ],
              #Only required if rtn_pitch_opt = "probDist"
              windspd_pars = NULL,
              #Only required if rtn_pitch_opt = "windSpeedReltn"
              rtn_pitch_windspd_dt = NULL,
              #Only required if rtn_pitch_opt = "windSpeedReltn"
              trb_wind_avbl = wind_farm_vals$trb_wind_avbl_data[[t]],
              trb_downtime_pars = wind_farm_vals$trb_downtime_data[[t]],
              wf_n_trbs = wind_farm_df()$Num_Turbines[t],
              wf_width = wind_farm_df()$WFWidth_km[t],
              wf_latitude = wind_farm_df()$Latitude[1],
              tidal_offset = 0,
              lrg_arr_corr = T,
              xinc = 0.05,
              yinc = 0.05,
              out_format = c("draws", "summaries"),
              out_sampled_pars = T,
              out_period = "months",
              season_specs = NULL,
              verbose = T,
              log_file = NULL,
              #"data/stoch_crm_SCRAM.log",
              seed = 11
            ))
      }
    }
    stoch_SCRAM_out(model_out)

    # running(FALSE) # done with run
    run_times$end <- Sys.time()
    output$run_end_txt <-
      renderText(paste0(
        "The model run was completed at: ",
        strftime(run_times$end, "%Y-%m-%d %H:%M:%S %Z", tz = "America/New_York")
      ))
    # run_elaps_time <- run_end_time - run_start_time
    updateTabItems(session, inputId = "tabsetpan", selected = "crm_results")

  })  #end observeEvent input$run
  
  # Render results output tab -----------------------------------------------

  # rendertext to report the probability of collisions exceeding a user-specified threshold
  # Divide the total number of collision results exceeding the threshold, divided by the total number of runs
  prob_exceed_threshold <- reactive({
    prob_threshold_list <- c()
    CRSpecies <- isolate(input$species_input)
    n <- 1
    for (s in 1:length(CRSpecies)){
      for (t in 1:wind_farm_vals$num_turb_mods){
        
        cTurbModel <- paste0("turbModel", wind_farm_df()[t, "TurbineModel_MW"])
        threshold_text <- ""
        #add species and turbine input variable to carry the probabilities for each run through
        # sum_collisions <- rowSums(CRM_fun$output[[t]], na.rm = T)
        
        sum_collisions <- rowSums(stoch_SCRAM_out()[[s]][[t]]$collisions[[1]], na.rm = T)
        
        threshold_text <-
          length(which(sum_collisions > isolate(input$inputthreshold)
          )) /
          length(sum_collisions)
        
        if (threshold_text == 1) {
          threshold_text <-
            paste(">", isolate(round(
              1 - 1 / input$iter_slider, log10(input$iter_slider)
            )), sep = " ")
        }
        else if (threshold_text == 0) {
          threshold_text <-
            paste("<", isolate(round(((1 / input$iter_slider)
            ), log10(
              input$iter_slider
            ))), sep = " ")
        }
        else {
          #add some formatting to other probabilities
          threshold_text <-
            formatC(threshold_text, digits = 3, format = "fg")
        }
        prob_threshold_list[n] <-
          paste0(
            "Run ",
            n,
            ": the probability of exceeding specified threshold (",
            isolate(input$inputthreshold),
            ") is ",
            threshold_text,
            "."
          )
        n <-  n + 1
      }
    }
    return(prob_threshold_list)
  })

  #now render in correct format to output to Shiny with newlines as needed.
  output$prob_exceed_threshold_msg <- renderUI({
    HTML(paste(prob_exceed_threshold(), collapse = " <br> "))
  })

  # # after main function is run, but before sensitivity analyses are, create a .csv that alerts the user (will be included in download package if user downloads raw results)
  # observeEvent(input$run,
  #              if(!is.null(CRM_fun())){
  #                write.csv("Sensitivity analyses not run", file = paste0(tempdir(), "/", CRM_fun()[['CRSpecies']][1], "_", paste0("turbModel", CRM_fun()[['Turbines']][1]),"_sensitivity.csv"), row.names = FALSE)
  #              })

  # main plot for annual cumulative collision estimation
  results_plots <- eventReactive(run_times$end, {
    
    num_species <- length(stoch_SCRAM_out())
    num_turb_mods <- length(stoch_SCRAM_out()[[1]])
    
    plot_list = vector("list",num_turb_mods)
    nbins = 20
    for(s in 1:num_species) {
      for (t in 1:num_turb_mods) {
        species <- stoch_SCRAM_out()[[s]][[t]]$species
        cTurbModel <- paste0("turbModel", stoch_SCRAM_out()[[s]][[t]]$turbine_model)
        if(!is.null(stoch_SCRAM_out()[[s]][[t]]$collisions)){
          # ATG - an issue with plots not rendering in dynamic tabs; use local to get output correct
          # https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
          plot_list[[t]] <- local({
            
            NA_index <- which(is.na(stoch_SCRAM_out()[[s]][[t]]$collisions[[1]][1,]))
            
            outvector <-
              round(rowSums(stoch_SCRAM_out()[[s]][[t]]$collisions[[1]], na.rm = TRUE), digits = 3)
            xmin <- round(min(outvector, na.rm = TRUE))
            xmax <- round(max(outvector, na.rm = TRUE))
            hist_outvector <-
              hist(outvector, breaks = nbins, plot = F)
            freq_outvector <- hist_outvector$counts
            bar_outvector <-
              as.data.frame(cbind(
                mids = hist_outvector$mids,
                density = freq_outvector / sum(freq_outvector)
              ))
            
            if (sum(outvector, na.rm = T) == 0) {
              #no collisions provide a modified figure
              plot_xmin = 0
              plot_xmax = 10
              plot_ymin = 0
              plot_ymax = 10
              x_ann = 5
              y_ann = 5
              fig_text = "No collisions predicted"
              
            } else {
              plot_xmin = -1
              if (xmax + 2 < 10) {
                plot_xmax = 10
              } else {
                plot_xmax = xmax + 2
              }
              plot_ymin = 0
              plot_ymax = max(freq_outvector) / sum(freq_outvector) * 1.1  #calculate frequency of first element to set as max y
              x_ann = 0
              y_ann = 0
              fig_text = ""
            }
            
            #main plot title
            # if (length(which(SpeciesLabels[, t] == species)) > 0) {
            if (length(SpeciesLabels$name_underscore == species) > 0) {
              main_label <-
                paste0(
                  "Model option ", stoch_SCRAM_out()[[s]][[t]]$model_option,
                  " for species ", SpeciesLabels[SpeciesLabels$name_underscore == species, "name_space"],
                  " (turbine model ", stoch_SCRAM_out()[[s]][[t]]$turbine_model, ")"
                )
            } else{
              main_label <-
                paste0(
                  "Model option ", stoch_SCRAM_out()[[s]][[t]]$model_option,
                  " for species ",  species,
                  " (turbine model ", stoch_SCRAM_out()[[s]][[t]]$turbine_model, ")"
                )
            }
            
            bold <- rep(2, 12)
            month_col <- rep("dark blue", 12)
            month_col[NA_index] <- rgb(0, 0, 0, 0.8)
            bold[NA_index] <- 1 # make bold those months with data
            p1 <-
              ggplot2::ggplot(bar_outvector, aes(x = mids, y = density)) +
              #center on integers use binwidth = 1 and center = 0; but modified to make bar end at number for thresholding
              # stat_bin(aes(y = stat(count / sum(count))), col="darkgreen", fill="darkgreen", binwidth = 0.5, center = 0.5) +
              #change to geom_col to have more control over bars with addition of limits to 0 and above to threshold
              geom_col(fill = "darkgreen") +
              annotate(
                "rect",
                xmin = 0,
                xmax = input$inputthreshold,
                ymin = 0,
                ymax = plot_ymax,
                fill = rgb(1, 1, 1, 0.65),
                col = rgb(0, 0, 0, 0)
              ) +
              ggtitle(main_label) +
              scale_x_continuous(breaks = scales::pretty_breaks()) +
              # xlim(c(plot_xmin, plot_xmax)) +
              # ylim(c(plot_ymin, plot_ymax)) +
              xlab("Total collisions per year over months highlighted below") +
              ylab("Frequency") +
              annotate("text",
                       x = x_ann,
                       y = y_ann,
                       label = fig_text) +
              geom_vline(
                xintercept = input$inputthreshold,
                color = "red",
                linetype = "dashed"
              ) +
              annotate(
                "text",
                x = input$inputthreshold,
                y = plot_ymax * 0.5,
                col = "red",
                label = "Input threshold",
                angle = 90,
                vjust = -0.75
              ) +
              theme_classic()
            
            p2 <-
              cowplot::ggdraw(
                cowplot::add_sub(
                  p1,
                  label = month.abb,
                  x = seq(0.1, 0.9, 0.8 / 11),
                  color = month_col,
                  size = 10,
                  fontface = bold
                )
              )
            # Add plot to list to render in each tab
            plot(p2)
          })
          # names(plot_list[[n]]) <- main_label
          # n = n + 1
        }
      }
    }
    return(plot_list)
  })

  # ATG - causing flashing in plot - figure out why that is.
  # Started happening with rebuild of baked movemement was changed
  # This has something to do with plots telling user they are updating see:
  # https://github.com/rstudio/shiny/issues/1591
  # Adding  tags$style(type="text/css", ".recalculating {opacity: 1.0;}" removes this feature
  output$plot_tabs = renderUI({
    nTabs = length(results_plots())
    if (nTabs>0){
      histTabs = lapply(1:nTabs, function(x) {
        tabPanel(paste('Run', x),
                 renderPlot(results_plots()[[x]]), height = 400,
                 tags$style(type="text/css", ".recalculating {opacity: 1.0;}"))
      })
      do.call(tabsetPanel, histTabs)
    }
  })

  output$plot_results_caption <-
    renderUI(p("Figure 1: A histogram of the frequency of the total number of collisions per year.
                            The heights of the bars show the relative frequency of each value.
                            Months for which movement data were provided or available are shown in bold;
                            only bold months are shown in histogram of annual collisions.",
               style = "margin-top: 10px; margin-left: 5px; margin-right: 10px; font-size: 12px; font-style: italic;"))


  # # create a reactive object for the sensitivity analyses
  # GSA_fun <- eventReactive(
  #   input$runGSA, {
  #     GSA_approx(CRM_fun$output, input$crm_option_radio)
  #   })

  # # write results of sensitivity analysis to .csv (to be downloaded if user chooses)
  # observeEvent(input$runGSA, {
  #   for(t in 1:length(CRM_fun$output)){
  #       write.csv(GSA_fun()[[t]], file = paste0(tempdir(), "/SCRAM_", CRM_fun()[['CRSpecies']][i], "_", paste0("turbModel", CRM_fun()[['Turbines']][e]),"_sensitivity.csv"), row.names = FALSE)
  #     }
  # })

  option_labels <- c("Option 1: faster approximation", NA,"Option 3: slower but more precise")

  # # dialog box for sensitivity analyses
  # observeEvent(input$runGSA, {
  #   {showModal(modalDialog(
  #     title = "Sensitivity analysis ran successfully",
  #     footer = modalButton("OK"),
  #     paste("Global sensivity analysis for ", option_labels[as.numeric(input$crm_option_radio)], " ran successfully. Results ready to download.", sep="")
  #   ))}
  # })


  # print message alerting when main script is run successfully
  observeEvent(input$run, {
    output$run_success_msg <- renderText({
    # if(!is.null(CRM_fun$output[1])){
    if(!is.null(stoch_SCRAM_out())){
      paste(option_labels[as.numeric(input$crm_option_radio)], " ran successfully.", sep="")
    }
  })
  })

  # # when results are available, render button for running sensitivity analyses
  # observeEvent({req(!is.null(!is.null(stoch_SCRAM_out())))}, {
  #   output$runGSA2 <- renderUI({
  #     actionButton("runGSA", HTML("Run sensitivity <br/> analysis"),
  #                  style = "width: 150px; margin-left: 0px; margin-top: 0px; background-color: navy; color: white; font-weight: bold;")
  #   })
  # })

  # when results are available, render button for downloading raw model run data
  observeEvent({req(!is.null(stoch_SCRAM_out()))}, {
    output$download_output <- renderUI({
      downloadButton("downloadDataRaw", HTML("Download model <br/> results"),
                     style = "width: 150px; margin-left: 0px; margin-top: 0px; background-color: darkkhaki; color: white; font-weight: bold;")
    })
  })

  # if results are available, render button for generating report
  observeEvent({req(!is.null(stoch_SCRAM_out()))}, {
    output$genreport <- renderUI({
      downloadButton("report", HTML("Generate output <br/> report"),
                     style = "width: 150px; margin-left: 0px; margin-top: 0px; background-color: darkviolet; color: white; font-weight: bold;")
    })
  })


  # SCRAM model outputs -----------------------------------------------------

  # download handler for report using R Markdown
  output$report <- downloadHandler(
    # for PDF output, change this to "report.pdf"
    # must use function() of else the report name will not change on re-run
    filename = function() {paste0("SCRAM2_report_", strftime(run_times$end, "%Y%m%d_%H%M%S"), ".pdf")},
    content = function(file) {
      # copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report_SCRAM2_v081823.Rmd")
      img1 <- file.path(tempdir(), "SCRAM_logo_2_4inch.jpg")
      img2 <- file.path(tempdir(), "BRI_color_logo_no_words.png")
      img3 <- file.path(tempdir(), "URI.png")
      img4 <- file.path(tempdir(), "USFWS.png")
      img5 <- file.path(tempdir(), "BOEM.png")

      file.copy("scripts/report_SCRAM2_v081823.Rmd", tempReport, overwrite = TRUE)
      # need to copy images to temp dir otherwise can't be found
      # see: (https://stackoverflow.com/questions/35800883/using-image-in-r-markdown-report-downloaded-from-shiny-app?rq=1)
      file.copy("www/SCRAM_logo_2_4inch.jpg", img1, overwrite = TRUE)
      file.copy("www/BRI_color_logo_no_words.png", img2, overwrite = TRUE)
      file.copy("www/URI.png", img3, overwrite = TRUE)
      file.copy("www/USFWS.png", img4, overwrite = TRUE)
      file.copy("www/BOEM.png", img5, overwrite = TRUE)

      species_params <- bird_data %>%
        filter(Species == input$species_input)

      
      # set up parameters to pass to Rmd document
      if(input$migration_calc_type == "occup_model") {
        params <- list(
          SCRAM_version = SCRAM_version,
          project = input$project_name,
          modeler = input$modeler,
          run_start_time = isolate(run_times$start),
          run_end_time = isolate(run_times$end),
          iterations = input$iter_slider,
          threshold = input$inputthreshold,
          migration_calc_type = input$migration_calc_type,
          migrant_states = "",
          migr_front_width = "",
          migr_front_line = "",
          migr_corridor = "",
          model_option = input$crm_option_radio,
          model_type = species_params_vals$model_input_dist_type,
          model_cell = windfarm_loc$cell_sf,
          trans_prop = bird_flux_vals$trans_prop,
          model_output = stoch_SCRAM_out(),
          daily_cell_flux = bird_flux_vals$num_birds_cell_perday,
          spp_move_data_summary = spp_move_data(),
          prob_exceed = prob_exceed_threshold(),
          species_labels = SpeciesLabels,
          species_popn_data = popn_data[which(popn_data$Species == input$species_input), ],
          species_popn_assumptions = spp_popn_notes()
        )
      } else { #annex 6
        params <- list(
          SCRAM_version = SCRAM_version,
          project = input$project_name,
          modeler = input$modeler,
          run_start_time = isolate(run_times$start),
          run_end_time = isolate(run_times$end),
          iterations = input$iter_slider,
          threshold = input$inputthreshold,
          migration_calc_type = input$migration_calc_type,
          migrant_states = annex6_vals$filter_lcn,
          migr_front_width = annex6_vals$migr_front_km,
          migr_front_line = annex6_vals$mig_corridor_line,
          migr_corridor = spp_corridor(),
          model_option = input$crm_option_radio,
          model_type = species_params_vals$model_input_dist_type,
          model_cell = windfarm_loc$cell_sf,
          trans_prop = "",
          model_output = stoch_SCRAM_out(),
          daily_cell_flux = "",
          spp_move_data_summary = "",
          prob_exceed = prob_exceed_threshold(),
          species_labels = SpeciesLabels,
          species_popn_data = annex6_vals$spp_popn_data,
          species_popn_assumptions = spp_popn_notes()
        )
        
      }

      save(params, file = "scripts/test_params.RData")

      # can't render to PDF - error with latexpdf. Found this solution:
      # https://stackoverflow.com/questions/66056764/knitr-cannot-find-pdflatex-when-creating-pdf-from-shiny-app
      # "You should NOT specify the output_file argument in render() Instead you need to rename the file AFTER rendering."

      out <- rmarkdown::render(tempReport,
                               params = params,
                               envir = new.env(parent = globalenv())
      )
      file.rename(out, file)
    }
    # contentType = "application/pdf"
  )

  # download handler for example for turbine data
  output$downloadTurbineExample <- downloadHandler(
    # filename = "TurbineData_example.zip",
    filename = "TurbineData_inputs_example.zip",
    content = function(file) {
      file.copy("data/TurbineData_inputs_example.zip", file)
    }
  )

# download handler for raw results download
  output$downloadDataRaw <- {
    downloadHandler(
      # must use function() of else the zip name will not change on re-run
      filename = function() {
        paste0('SCRAM2_model_output_', strftime(run_times$end, "%Y%m%d_%H%M%S"), '.zip')
      },
      content = function(fname) {
        tmpdir = tempdir()
        fnames4zip1 <- list()
        
        num_species <- length(stoch_SCRAM_out())
        num_turb_mods <- length(stoch_SCRAM_out()[[1]])
        
        #loop through species and turbines to process and write out data
        for(s in 1:num_species){
          for(t in 1:num_turb_mods){
            model_output <- stoch_SCRAM_out()[[s]][[t]]
            species <- model_output$species
            turbine_type <- paste0("turbModel", model_output$turbine_model)
            model_option <- input$crm_option_radio
            model_type <- model_output$model_type
            pred_monthly_coll <- model_output$collisions[[1]]
            colnames(pred_monthly_coll) <- paste0("crm_pred_", month.abb)
            pred_monthly_coll <- cbind(run=1:iter_slider_react(), pred_monthly_coll)
            write.csv(pred_monthly_coll, file = paste0(tmpdir, "/SCRAM2_crm_pred_monthly_coll_opt", model_option,  "_", species, "_", turbine_type,".csv"), row.names = FALSE)
            
            #process sampled paramaters (turbine and bird) as outputted by stochCRM function
            sampled_pars <- model_output$sampled_pars
            
            turbine_pars_1by5_list <- sampled_pars[c("air_gap", "hub_height", "bld_pitch", "bld_width", "rtr_radius", "rtn_speed")]
            
            turbine_pars_1by5_df <- do.call("rbind", turbine_pars_1by5_list) %>% 
              mutate(turbine_parameter = row.names(.)) %>% 
              select(turbine_parameter, mean, sd, median, pctl_2.5, pctl_97.5)
            
            write.csv(turbine_pars_1by5_df, file = paste0(tmpdir, "/SCRAM2_sampled_turbine_pars_1_opt", model_option,  "_", species, "_", turbine_type,".csv"), row.names = FALSE)
            
            turbine_pars_12by6_list <- sampled_pars[c("prop_oper_mth", "downtime")]
            
            turbine_pars_12by6_df <- do.call("rbind", turbine_pars_12by6_list) %>% 
              mutate(turbine_parameter = row.names(.)) %>% 
              select(turbine_parameter, mean, sd, median, pctl_2.5, pctl_97.5)
            
            write.csv(turbine_pars_12by6_df, file = paste0(tmpdir, "/SCRAM2_sampled_turbine_pars_2_opt", model_option,  "_", species, "_", turbine_type,".csv"), row.names = FALSE)

            if(model_option == "1") {
              bird_pars_1by5_list <- sampled_pars[c("body_lt", "flt_speed", "noct_actv", "wing_span", "avoid_bsc", "avoid_bsc")]
            } else if(model_option == "3") {
              bird_pars_1by5_list <- sampled_pars[c("body_lt", "flt_speed", "noct_actv", "wing_span", "avoid_bsc", "avoid_ext")]
            }
            
            bird_pars_1by5_df <- do.call("rbind", bird_pars_1by5_list) %>% 
              mutate(bird_parameter = row.names(.)) %>% 
              select(bird_parameter, mean, sd, median, pctl_2.5, pctl_97.5)
            
            write.csv(bird_pars_1by5_df, file = paste0(tmpdir, "/SCRAM2_sampled_bird_pars_opt", model_option,  "_", species, "_", turbine_type,".csv"), row.names = FALSE)
            
            #compress files for uploading the wind farm specific data
            fnames4zip1 <- c(fnames4zip1, paste0(tmpdir, "/SCRAM2_crm_pred_monthly_coll_opt", model_option,  "_", species, "_", turbine_type,".csv"))
            fnames4zip1 <- c(fnames4zip1, paste0(tmpdir, "/SCRAM2_sampled_turbine_pars_1_opt", model_option,  "_", species, "_", turbine_type,".csv"))
            fnames4zip1 <- c(fnames4zip1, paste0(tmpdir, "/SCRAM2_sampled_turbine_pars_2_opt", model_option,  "_", species, "_", turbine_type,".csv"))
            fnames4zip1 <- c(fnames4zip1, paste0(tmpdir, "/SCRAM2_sampled_bird_pars_opt", model_option,  "_", species, "_", turbine_type,".csv"))
            
          } #turbine model loop
          
          #write species specific data
          write.csv(popn_data[which(popn_data$Species==species),], file = paste0(tmpdir, "/SCRAM_", species, "_popn_data.csv"), row.names = FALSE)
          write.csv(bird_data[which(bird_data$Species==species),], file = paste0(tmpdir, "/SCRAM_", species, "_character_data.csv"), row.names = FALSE)
          write.csv(bird_flux_vals$num_birds_cell_perday, file = paste0(tmpdir, "/SCRAM2_num_birds_cell_perday_opt", model_option,  "_", species,".csv"), row.names = FALSE)

          fnames4zip1 <- c(fnames4zip1, paste0(tmpdir, "/SCRAM_", species, "_popn_data.csv"))
          fnames4zip1 <- c(fnames4zip1, paste0(tmpdir, "/SCRAM_", species, "_character_data.csv"))
          fnames4zip1 <- c(fnames4zip1, paste0(tmpdir, "/SCRAM2_num_birds_cell_perday_opt", model_option,  "_", species,".csv"))
          fnames4zip1 <- c(fnames4zip1, paste0("data/", species, "_ht_dflt_v2.csv"))
          fnames4zip1 <- c(fnames4zip1, paste0("data/movements/", species, "_movements_", model_output$model_type, ".zip"))
          
        } #spp model loop
        
        file.copy(input$file_wf_param$datapath, to = file.path(tmpdir, input$file_wf_param$name))
        
        # set up parameters to pass to Rmd document
        if(input$migration_calc_type == "occup_model") {
          scram_params <- list(
            SCRAM_version = SCRAM_version,
            project = input$project_name,
            modeler = input$modeler,
            run_start_time = isolate(run_times$start),
            run_end_time = isolate(run_times$end),
            iterations = input$iter_slider,
            threshold = input$inputthreshold,
            migration_calc_type = input$migration_calc_type,
            model_option = input$crm_option_radio,
            model_type = species_params_vals$model_input_dist_type,
            model_cell = windfarm_loc$cell_sf,
            trans_prop = bird_flux_vals$trans_prop,
            model_output = stoch_SCRAM_out(),
            daily_cell_flux = bird_flux_vals$num_birds_cell_perday,
            spp_move_data_summary = spp_move_data(),
            prob_exceed = prob_exceed_threshold(),
            species_labels = SpeciesLabels,
            species_popn_data = popn_data[which(popn_data$Species == input$species_input), ],
            species_popn_assumptions = spp_popn_notes()
          )
        } else { #annex 6
          scram_params <- list(
            SCRAM_version = SCRAM_version,
            project = input$project_name,
            modeler = input$modeler,
            run_start_time = isolate(run_times$start),
            run_end_time = isolate(run_times$end),
            iterations = input$iter_slider,
            threshold = input$inputthreshold,
            migration_calc_type = input$migration_calc_type,
            migrant_states = annex6_vals$filter_lcn,
            migr_front_width = annex6_vals$migr_front_km,
            model_option = input$crm_option_radio,
            model_cell = windfarm_loc$cell_sf,
            model_output = stoch_SCRAM_out(),
            prob_exceed = prob_exceed_threshold(),
            species_labels = SpeciesLabels,
            species_popn_data = annex6_vals$spp_popn_data,
            species_popn_assumptions = spp_popn_notes())
        }
        save(scram_params, file = file.path(tmpdir, paste0('SCRAM2_model_pars_output_', strftime(isolate(run_times$end), "%Y%m%d_%H%M%S"),'.RData')))

        #add additional files for compressing
        fnames4zip1 <- c(fnames4zip1, file.path(tmpdir, input$file_wf_param$name))
        fnames4zip1 <- c(fnames4zip1, file.path(tmpdir, paste0('SCRAM2_model_pars_output_', strftime(isolate(run_times$end), "%Y%m%d_%H%M%S"),'.RData')))
        
        utils::zip(zipfile=fname, files=unlist(c(fnames4zip1)), flags = "-r9Xj")
        
      },
      contentType = "application/zip"
    )}

} # end server function

# Run the application
shinyApp(ui = ui, server = server)


