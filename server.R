#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(Matrix)
library(tidyverse)
library(magrittr)
library(plotly)
library(adegenet)
library(poppr)
library(hierfstat)
library(DT)
library(shinythemes)
library(shinybusy)
library(shinyalert)
library(stringr)
library(leaflet)
library(maptools)
library(sf)
date <- parse_date("2020/10/21")
vers <- 1.0


walleye_gen <- read_rds("./Data/Bowles_walleye_2020_10_20.rds")
# walleye_gen <- read_rds("./Data/Bowles_walleye_2020_10_20_no_outlier.rds")

locations <- tribble( ~location,   ~lat,     ~lon,
                      "CHALIFOUR", 50.511214, -73.845557,
                      "PERCH",     50.336787, -73.764081,
                      "ICON",      50.227803, -73.811791,
                      "TAKWA",     51.347222, -72.782778)


data_fail<- function(){shinyalert(title = "No Data!",
                                  text = "Please select fish species and click 'Get Data'!",
                                  type = "error")}

fail_loc<- function(){shinyalert(title = "Find Fish Missing Loci!",
                                 text = "Please click 'Find Fish Missing Loci'!",
                                 type = "error")}

fail_fish<- function(){shinyalert(title = "Find Loci Missing Fish!",
                                  text = "Please click 'Find Loci Missing Fish'!",
                                  type = "error")}



# points to line function taken from: https://rpubs.com/walkerke/points_to_line
points_to_line <- function(data, long, lat, id_field = NULL, sort_field = NULL) {
    
    # Convert to SpatialPointsDataFrame
    coordinates(data) <- c(long, lat)
    
    # If there is a sort field...
    if (!is.null(sort_field)) {
        if (!is.null(id_field)) {
            data <- data[order(data[[id_field]], data[[sort_field]]), ]
        } else {
            data <- data[order(data[[sort_field]]), ]
        }
    }
    
    # If there is only one path...
    if (is.null(id_field)) {
        
        lines <- SpatialLines(list(Lines(list(Line(data)), "id")))
        
        return(lines)
        
        # Now, if we have multiple lines...
    } else if (!is.null(id_field)) {  
        
        # Split into a list by ID field
        paths <- sp::split(data, data[[id_field]])
        
        sp_lines <- SpatialLines(list(Lines(list(Line(paths[[1]])), "line1")))
        
        # I like for loops, what can I say...
        for (p in 2:length(paths)) {
            id <- paste0("line", as.character(p))
            l <- SpatialLines(list(Lines(list(Line(paths[[p]])), id)))
            sp_lines <- spRbind(sp_lines, l)
        }
        
        return(sp_lines)
    }
}



low_col = "midnightblue"
high_col = "magenta1"
na_col = "grey10"



# FSTsum class to store FST information for easy calling and plotting
FSTsum <- setClass(
    Class = "FSTsum",
    slots = c(
        sum_tibble = "data.frame",
        x_argument = "character",
        y_argument = "character"
    ),
    prototype = list(
        sum_tibble = tibble(
            sites_x = character(),
            sites_y = character(),
            Fst = numeric(),
            date_x = numeric(),
            date_y = numeric(),
            pop_x = character(),
            pop_y = character(),
            popdate_x = character(),
            popdate_y = character(),
            lat_x = numeric(),
            lon_x = numeric(),
            lat_y = numeric(),
            lon_y = numeric()
        ),
        x_argument = "pop_x",
        y_argument = "pop_y"
    ),
    
    validity = function(object){}
    
)



d_FST <- function(d, by.sites = "ALL",
                  by.years = "ALL"){
    if(by.years == "ALL"){
        year_sel <- d %>% 
            select(., date_x) %>% 
            unique %>% unlist %>% as.numeric
        
    }else{year_sel <- by.years}
    
    
    if(by.sites == "ALL"){
        sites_sel <- d %>% 
            select(., site_x) %>% 
            unique %>% unlist %>% as.character
        
    }else{sites_sel <- by.sites}
    
    
    dplot<-d %>% 
        dplyr::filter(., date_x %in% year_sel,
                      date_y %in% year_sel,
                      pop_x %in% sites_sel,
                      pop_y %in% sites_sel)
    
    if(length(sites_sel) > 1 & length(year_sel) > 1){
        x_arg <- "popdate_x"
        y_arg <- "popdate_y"
    }else if(length(sites_sel) > 1 & length(year_sel) == 1){
        x_arg <- "pop_x"
        y_arg <- "pop_y"
    }else if(length(sites_sel) == 1 & length(year_sel) >1){
        x_arg <- "date_x"
        y_arg <- "date_y"
    }else{stop("Error: missing site or year")}
    
    d_fst <- FSTsum(sum_tibble = dplot,
                    x_argument = x_arg,
                    y_argument = y_arg)
    
    # d_fst <- list(dplot, x_arg, y_arg)
    return(d_fst)
    
}



plot_FST <- function(d, x_col, y_col){
    ggplot(d, 
           aes(x = eval(parse(text=x_col)), y= eval(parse(text=y_col)), fill=Fst)) +
        geom_tile() +
        scale_fill_gradient(na.value = na_col, low = low_col, high = high_col) +
        scale_y_discrete(position = "right") +
        theme(panel.grid = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        xlab("") + ylab("")
    
    
    
}

FST_exp <- paste0("Determining genetic relatedness (or difference) among populations is important.",
                  "This allows geneticists to assess diversity within  species among different sites.",
                  " Geneticists can use this information to determine the amount of gene flow between populations",
                  " (by migration and interbreeding). Furthermore, if there is more difference among geographically",
                  " close populations, this may allow for increased ability to adapt to changing conditions because",
                  " diversity is high. FST is a common measure of genetic relatedness among populations.",
                  " FST values range from 0 (completely related) to 1 (completely unrelated).",
                  " It is calculated by comparing percentage of heterozygous genotypes within each population between",
                  " populations.","<br>"," Consider the gene", "<b>","<i>", " A", "</b>","</i>"," with alleles",
                  "<b>","<i>"," A", "</b>","</i>", " and","<b>","<i>"," a","</b>","</i>",". You would compare the proportion of individuals that are",
                  "<b>","<i>"," Aa","</b>","</i>", " between the populations. For population genetics, we use",
                  " individual nucleotides instead of genes. If the same spot (or locus) on a chromosome can be coded",
                  " by different nucleotides on homologous chromosomes, then it is said to be",
                  "<i>"," polymorphic","</i>", ". Geneticists compare the heterozygosity of these Single",
                  " Nucleotide Polymorphisms (SNPs) between populations. By doing this for many different SNPs,",
                  " you can  calculate FST. This app will walk you through how to do this with real genetic data!",
                  " The data are from different breeding populations of fish from sites around the Lake Mistassini area.
                   Because all of these populations in this system are known to return to their hatching site to spawn, 
                   the terms 'population' and 'site' are used interchangeably. Have fun!" )

sp_sel_exp <- paste0("Let's get started! Select a species to work with, the click the 'Get Data' button. This will import",
                     " data that are publicly available from a peer-reviewed study. You can access data sets and publications",
                     " by navigating to the 'Documentation' tab.", "<br>",
                     " Once you import your data, you will see a table that has the number of alleles found at each locus",
                     " for each fish. Each column represents a different SNP locus (a long number code) and each row represents a fish.",
                     " For example, if fish <b> 3 </b> has the number '2' in column one, that means it is heterozygous at",
                     " that locus because it has 2 different alleles, one on each chromosome. One could be the nucleotide",
                     "<b>", " A", "</b>", " and the other", "<b>"," C", "</b>",". If a fish is homozygous at a locus,",
                     " it will have a '1' in that column. If the identity of the locus could not be determined, it is said",
                     " to be 'missing' and is represented by a '0' or is left blank. Note that all of the",
                     " fish have been genotyped at thousands of loci (the plural of locus). This table only lets", 
                     " you see the first 25 loci because we will use other tools to summarise the data better!")

filter_loci_exp <- paste0("As you saw above, not all loci have been successfully genotyped in each fish. In this step,",
                          " you will screen each locus to see how many fish have been successfully genotyped. First,",
                          " click the button 'Find Fish Missing from Loci' to perform the screening. This takes a few",
                          " seconds to compute. Only click the button ONCE! If the fish GIF in the top corner is flopping,",
                          " this means the data are being processed. Once the fish stops flopping, you can click",
                          " 'Filter out loci'. The plot will show you which loci you will keep for analysis. Use the slider",
                          " to control the acceptable minimum number of fish a locus must have been successfully genotyped for
                          in order to be included in the analysis. For example, if the tolerance is set to 2%, only loci for which 
                          genotyping has failed in  <i> less than 2% of fish </i> will be included (i.e. > 98% genotyping
                          success for those SNP loci). The blue points on the graph indicate these loci. You can count the loci by scrolling through the data table and observing how many are
                          'kept'. After adjusting the slider, click 'Filter out loci' again.
                           Loci",
                          " <br> <br> Note: The more loci that you keep, the longer the final computation will take. You want",
                          " to keep at least 50 loci, but if you go too much over 100 loci, the computation may take several minutes.")

filter_fish_exp <- paste0("Each fish generally hasn't been genotyped at all loci in the dataset. We want to get rid of the fish",
                          " with the most incomplete data. Where do you think you should draw the line? Use the slider",
                          " to adjust the number of fish kept in the analysis. Notice that the slider adjsuts the threshold",
                          " of the successfully genotyped loci. So adjusting the slider to '10' means that you will only anlyze",
                          " fish that have at least 90% of there loci genotyped.This doesn't affect computation as much, so", 
                          " keep as many as you'd like! Just keep in mind that you want to use good quality fish! <br>",
                          " Click 'Find Loci Missing from Fish' to complete the computation, then adjust the slider and",
                          " click 'Filter out fish' to update your plot")

filter_exp <-paste0("The above steps have now identified the appropriate filters, it's time to finalize them and calculate",
                    " FST values. Check to make sure that no more than 200 loci are included in the analysis. If there are,",
                    " the computation will be very lengthy. <br> Make sure that you only click the filter button <b> ONCE </b>. Additional clicks will delay computation.",
                    " Remember, if the fish in the top right corner is flopping, a computation is occurring. Don't click anything else!",
                    " a summary table will be produced upon completion.")

FST_plot_exp <- paste0("This plots a color coded plot of the genetic distances among populations on the x & y axes.",
                       " Darker purlple squares (low FST) indicate less genetic distance (more similarity) between populations.",
                       " Brighter pink squares (high FST) indicate more genetic distance (less similarity) between populations (sites).",
                       " Toggle which comparisons you'd like to plot by checking the appropriate boxes. You can even compare the",
                       " same populations to each other throughout time!")

FST_map_exp <- paste0("This map will show either genetic distance of populations over years (click 'Plot Genetic Distance Between Years')",
                      " or genetic distance among populations in the same year (click 'Plot Genetic Distance Among Sites').",
                      " Use the toggles above to control which comparisons to make. The color scheme is the same as for the FST",
                      " grid plot above. For between years comparisons, both circle size and color represent genetic distance (larger circle = more difference).")


# Define server logic required for App
shinyServer(function(input, output, session) {
    session$onSessionEnded(stopApp)
    
    
    
    # fill in text blocks
    output$FST_exp <- renderText({FST_exp})
    output$sp_sel_exp <- renderText({sp_sel_exp})
    output$filter_loci_exp <- renderText({filter_loci_exp}) 
    output$filter_fish_exp <- renderText({filter_fish_exp}) 
    output$filter_exp <- renderText({filter_exp})
    output$FST_plot_exp <- renderText({FST_plot_exp})
    output$FST_map_exp <- renderText({FST_map_exp})
    
    
    #Data import tab
    
    # create placeholders so that there aren't just empty spots!  
    locmiss_tib <-tibble(proportion_scored = rep(0,10), locus = rep(NA,10),
                         rownum = c(1:10))
    output$locmiss_tbl <- renderDataTable({locmiss_tib[,-3]})
    
    output$locmiss_plot <- renderPlot({
        ggplot(data = locmiss_tib, 
               aes(x = rownum, 
                   y = proportion_scored)) + 
            geom_point()+
            xlab("Loci")+
            scale_y_continuous(name = "Percentage of fish with \n locus genotyped",
                               labels = scales::percent)
    })
    
    
    
    indmiss_tib <-tibble(proportion_scored = rep(0,10), fish = rep(NA,10),
                         rownum = c(1:10))
    
    output$indmiss_tbl <- renderDataTable({indmiss_tib[,-3]})
    
    output$indmiss_plot <- renderPlot({
        ggplot(data = indmiss_tib, 
               aes(x = rownum, 
                   y = proportion_scored)) + 
            geom_point()+
            xlab("Fish")+
            scale_y_continuous(name = "Percentage of loci with \n fish genotyped",
                               labels = scales::percent)
    })
    
    
    #Initialize reactive values to be updated by user
    vals <- reactiveValues(
        gen_obj = NULL,
        
        locmiss_tib = NULL,
        indmiss_tib = NULL,
        
        locmiss_update = NULL,
        indmiss_update = NULL,
        
        gen_obj_filtered = NULL,
        pwise_fst = NULL,
        pwise_fst_matrix = NULL,
        
        d_fst = FSTsum(),
        d_fst_yr = NULL,
        d_fst_site = NULL,
        lds = NULL,
        map = NULL
        
    )
    
    #Import the data
    observeEvent(input$GetData, {
        
        vals$gen_obj <- switch(input$species,
                               "Walleye" = walleye_gen
        )
        
        gen_tbl <- as_tibble(vals$gen_obj$tab)
        
        output$tbl_head <- renderDataTable(gen_tbl[, c(1:25)], 
                                           options = list(scrollX = TRUE))
        
        
        
    })
    
    #calculate the missing # fish missing from  each locus
    observeEvent(input$locmiss_eval, {
        if(is.null(vals$gen_obj)==TRUE){data_fail()}else{
            # starting locmiss_tib
            locmiss_obj <- vals$gen_obj %>% 
                propTyped(., by ="loc")
            
            vals$locmiss_tib <- locmiss_obj%>% 
                as_tibble %>% 
                mutate(locus = names(locmiss_obj)) %>%
                rename(proportion_scored = value) %>% 
                arrange(., desc(proportion_scored)) %>% 
                bind_cols(., rownum = c(1:nrow(.)))
            
            output$locmiss_tbl <- renderDataTable({
                vals$locmiss_tib %>% 
                    select(!rownum) %>% 
                    datatable %>% 
                    formatRound(columns = "proportion_scored", 4)})
        }
        
    })  
    
    
    #Create plot for the loci being filtered in/out 
    
    observeEvent(input$make_loc_filter_plots,{
        if(is.null(vals$gen_obj)){data_fail()}else if(is.null(vals$locmiss_tib)==TRUE){fail_loc()} else{
            
            vals$locmiss_update <- vals$locmiss_tib %>% 
                mutate(selected_loci = cut(.$proportion_scored,
                                           breaks = c(0,(1 - (input$locmiss_threshold/100)),1e7),
                                           labels = c("Rejected", "Kept")))
            
            output$locmiss_tbl <- renderDataTable({
                vals$locmiss_update %>%
                    select(!rownum) %>% 
                    datatable %>% 
                    formatRound(columns = "proportion_scored", 4)}) 
            
            
            output$locmiss_plot <-renderPlot({
                ggplot(data = vals$locmiss_update, 
                       aes(x = rownum, y = proportion_scored,
                           color = selected_loci)) + 
                    geom_point() +
                    xlab("Loci") +
                    scale_y_continuous(name = "Percentage of fish with \n locus genotyped",
                                       labels = scales::percent) +
                    theme(legend.position = "top", legend.title = element_blank())
 
            })
        }
        
        
    })
    
    #calculate the missing # loci missing from  each fish
    
    
    observeEvent(input$indmiss_eval, {
        if(is.null(vals$gen_obj)==TRUE){data_fail()} else{
            # starting indmiss_tib
            indmiss_obj<-vals$gen_obj %>% 
                propTyped(., by ="ind")
            
            vals$indmiss_tib <- indmiss_obj %>%    
                as_tibble %>% 
                mutate(fish = names(indmiss_obj)) %>%
                rename(proportion_scored = value) %>% 
                arrange(., desc(proportion_scored)) %>% 
                bind_cols(., rownum = c(1:nrow(.)))
            
            output$indmiss_tbl <- renderDataTable({
                vals$indmiss_tib %>% 
                    select(!rownum) %>% 
                    datatable %>% 
                    formatRound(columns = "proportion_scored", 4)})
        }
    }) 
    
    #Create plot for the fish being filtered in/out   
    observeEvent(input$make_ind_filter_plots,{
        if(is.null(vals$gen_obj)==TRUE){data_fail()}else if(is.null(vals$indmiss_tib)==TRUE){fail_fish()} else{
            
            vals$indmiss_update <- vals$indmiss_tib %>% 
                mutate(selected_fish = cut(.$proportion_scored,
                                           breaks = c(0,(1 - (input$indmiss_threshold/100)),1e7),
                                           labels = c("Rejected", "Kept")))
            
            output$indmiss_tbl <- renderDataTable({
                vals$indmiss_update %>% 
                    select(!rownum) %>% 
                    datatable %>% 
                    formatRound(columns = "proportion_scored", 4)})
            
            output$indmiss_plot <- renderPlot({
                ggplot(data = vals$indmiss_update, 
                       aes(x = rownum, y = proportion_scored,
                           color = selected_fish)) + 
                    geom_point()+
                    xlab("Fish")+
                    scale_y_continuous(name = "Percentage of loci with \n fish genotyped",
                                       labels = scales::percent) +
                    theme(legend.position = "top", legend.title = element_blank())
                
                
            })
            
            
        }
        
        
        
    })
    
    
    # finalize filters. Also, calculate and average pairwise FST values for all
    # comparisons in remaining data
    
    observeEvent(input$filter_gen_obj,{
        if(is.null(vals$gen_obj)==TRUE){data_fail()}else{
            vals$gen_obj_filtered <- vals$gen_obj %>% 
                missingno(., type = "loci", cutoff = input$locmiss_threshold/100) %>% 
                missingno(., type = "geno", cutoff = input$indmiss_threshold/100) %>% 
                `[`(loc = which(isPoly(.)==TRUE))
            
            
            vals$pwise_fst <- genet.dist(vals$gen_obj_filtered, method = "WC84")
            
            ##organize fst information
            
            vals$pwise_fst_matrix<- vals$pwise_fst %>% 
                as.matrix %>% 
                replace(0, .00001) %>%  #this correction will allow for some plotting shennanigans later
                
                tril %>% ##keep lower triangular portion of matrix
                as.matrix %>% 
                
                #### make ordering vector based on FST similarity
                #replace negative vals with 0s
                
                as_tibble %>% 
                replace(., . < 0 , 0.00001) %>% 
                bind_cols(sites_x = names(.),.) %>% 
                pivot_longer(., cols = everything()[-1],
                             names_to = "sites_y",
                             values_to = "Fst")%>% 
                dplyr::mutate(.,date_x = 
                                  str_extract_all(sites_x, "(\\d){4}"),
                              date_y = str_extract_all(sites_y, "(\\d){4}"),
                              pop_x = str_extract(sites_x, "[^_]+"),
                              pop_y = str_extract(sites_y, "[^_]+")) %>% 
                dplyr::mutate(., across(date_x:pop_y, unlist)) %>% 
                dplyr::mutate(., across(date_x:date_y, as.numeric)) %>% 
                dplyr::mutate(., popdate_x = str_c(pop_x, " ", date_x),
                              popdate_y = str_c(pop_y, " ", date_y)) %>% 
                mutate(Fst = na_if(Fst, 0)) %>%
                left_join(., locations, by = c("pop_x" = "location")) %>% 
                rename(lat_x = lat, lon_x = lon) %>% 
                left_join(., locations, by = c("pop_y" = "location")) %>% 
                rename(lat_y = lat, lon_y = lon)
            
            
            # print(head(vals$pwise_fst_matrix))   #this was a quality control check.
            
            
            
            
            summary_tib <- tibble(attribute = c("Number of Loci Included",
                                                "Tolerance of % Fish Missing from Loci",
                                                "Number of Fish Included",
                                                "Tolerance of % Loci Missing from Fish",
                                                "Site Names",
                                                "Date Time"),
                                  values = c(length(vals$gen_obj_filtered$loc.n.all),
                                             input$locmiss_threshold,
                                             length(vals$gen_obj_filtered$pop),
                                             input$indmiss_threshold,
                                             str_flatten(as.character(levels(vals$gen_obj_filtered$pop)), collapse =", "),
                                             format(Sys.time(), "%a %b %d %X %Y")))
            
            output$summary_tbl <- renderTable({summary_tib})
            
        }
        
    })
    
    
    
    
    # The checkbox options for the plots have to be made reactively from the data site. This means, that in the future,
    # when different data sets are used, the fields will match what is in the data set.
    fst_plot_data <- observe({
        req(vals$pwise_fst_matrix)
        
        updateCheckboxGroupInput(session = session, inputId = "location",
                                 label = "Select Locations",
                                 choices = unique(vals$pwise_fst_matrix$pop_x), selected = NULL)
        
        updateCheckboxGroupInput(session = session, inputId = "year",
                                 label = "Select Years",
                                 choices = unique(c(vals$pwise_fst_matrix$date_x, 
                                                    vals$pwise_fst_matrix$date_y)), selected = NULL)
    })
    
    
    
    
    observeEvent(input$make_plot, {
        output$fst_heat <- renderPlot({
            if(is.null(vals$pwise_fst_matrix)== T){
                
                shinyalert(title = "No Data!",
                           text = "Please finalize filters!",
                           type = "error")
                
            }else if(((length(input$year) <= 1 & length(input$location) <=1))|
                     is.null(input$year)==T | is.null(input$location)==T){
                # print(c(input$year, length(input$year), input$location, length(input$location)))
                
                shinyalert(title = "Not Enough Comparisons",
                           text = "Please select at least one year or site and two of the other variable",
                           type = "error")
                
            }else { 
                
                vals$d_fst <- d_FST(d = vals$pwise_fst_matrix, by.sites = input$location, by.years = input$year)
                plot_FST(d= vals$d_fst@sum_tibble, x_col = vals$d_fst@x_argument, y_col = vals$d_fst@y_argument)
                
            }
            
            
            
            
        })
    })
    
    # instantiate leaflet plot
    map <- leaflet() %>% fitBounds(-74, 50, -72.5, 51.7) %>% addTiles()
    output$map <- renderLeaflet({
        map
    })
    
    # make reactive color palette
    pal_cont <- reactive({
        colorNumeric(c(low_col, high_col), vals$d_fst@sum_tibble$Fst, na.color = na_col) 
    })
    
    # make reactive data objects for plotting leaflet
    
    #among years comparisons
    d_fst_yr <- reactive({
        vals$d_fst@sum_tibble %>% 
            filter(Fst !=0 ) %>% 
            mutate(Fst = .$Fst %>% replace_na(0)) %>% 
            filter(pop_x == pop_y)
    })
    
    #among sites comparisons
    d_fst_site <- reactive({
        vals$d_fst@sum_tibble %>% 
            filter(Fst != 0) %>% 
            mutate(Fst = .$Fst %>% replace_na(0)) %>% 
            filter(pop_x != pop_y)
    })
    
    observeEvent(input$plot_yrs,{
        
        # print(d_fst_yr())
        dy <- d_fst_yr() %>% as_tibble
        
        # print(nrow(dy)) 
        pal <- pal_cont()
        
        if(length(input$location) < 1){
            shinyalert(title = "Invalid number of sites",
                       text = "Please select at least one site when mapping genetic distance between years!",
                       type = "error")
        }else if(length(input$year) !=2){
            # print(length(input$year))
            shinyalert(title = "Invalid number of years",
                       text = "Please select 2 years (ONLY) to compare between",
                       type = "error")
        }else if(nrow(dy)>0){
            
            # update leaflet plot with points for years
            leafletProxy("map", data = d_fst_yr()) %>% 
                clearShapes() %>% 
                addCircles(lat = ~dy$lat_x, lng = ~dy$lon_x,
                           radius = ~(dy$Fst^2)*1e9, 
                           fillColor = ~pal(dy$Fst),
                           fillOpacity = .8, stroke = FALSE) %>% 
                addPopups(lat = ~dy$lat_x, lng = ~dy$lon_x, popup = ~dy$pop_x)
        }
        
    })
    
    observeEvent(input$plot_sites,{
        
        
        # further wrangling required to get data in format appropriate for points_to_line() function
        
        d_fst_site_x <- d_fst_site()%>% as_tibble %>%  
            filter(pop_x != pop_y) %>% 
            mutate(pair = as.factor(c(1:nrow(.)))) %>% 
            select(ends_with("_x"), Fst, pair) 
        
        d_fst_site_y <- d_fst_site() %>% as_tibble %>% 
            filter(pop_x != pop_y) %>% 
            mutate(pair = as.factor(c(1:nrow(.)))) %>% 
            select(ends_with("_y"), Fst, pair)
        
        names(d_fst_site_y) <- names(d_fst_site_x)
        
        ds <- d_fst_site_x %>% 
            bind_rows(., d_fst_site_y)
        
        pal <- pal_cont()
        
        # print(ds)
        
        if(length(input$year)!=1){
            shinyalert(title = "Invalid number of years",
                       text = "Please select ONLY one year when mapping genetic distance among sites!",
                       type = "error")
        }else if(length(input$location) < 1){
            shinyalert(title = "Invalid number of sites",
                       text = "Please select multiple sites to compare!",
                       type = "error")
        }else if(nrow(ds)>3){                                  # If more than 3 rows in ds, use plot_to_lines helper function
            ds_col <- pal(ds$Fst[1:(nrow(ds)/2)])                # otherwise, use regular addPolylines funcionality
            dss <- points_to_line(ds, "lon_x", "lat_x", "pair")
            # print(nrow(ds))
            leafletProxy("map", data = d_fst_site()) %>%
                clearShapes() %>%
                addPolylines(data= dss, color = ds_col )  
        }else if(nrow(ds)<3)
            
            leafletProxy("map", data = d_fst_site()) %>%
            clearShapes() %>%
            addPolylines(lat = ~ds$lat_x, lng = ~ds$lon_x, color= ~pal(ds$Fst[1]), 
                         opacity = .8) %>% 
            addPopups(lat = ~ds$lat_x, lng = ~ds$lon_x, popup = ~ds$pop_x)
        
        
    })

})
