#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(tidyverse)

# craft the table
gene.go_ko <- read_delim("data/all.proteins.emapper.out.emapper.annotations", delim = "\t", col_names = FALSE, comment = "#") %>%
    # select(c(1,6,7,8,9,16,22))
    select(c(1,6,7,8,9,10,11,12,13,14,15,16,22))

gene.go_ko <- data.frame(base::lapply(gene.go_ko, function(x) {
    gsub(",", ", ", x)
}))

colnames(gene.go_ko) <- c("gene",
                          "preferred_name",
                          "GO",
                          "EC",
                          "KEGG_ko",
                          "KEGG_Pathway",
                          "KEGG_Module",
                          "KEGG_Reaction",
                          "KEGG_rclass",
                          "BRITE",
                          "KEGG_TC",
                          "CAZy",
                          "Free_text_description")

bins_from_annotations <- gene.go_ko %>%
    dplyr::select(gene) %>%
    tidyr::separate(gene, into = c("AS", "TE", "serial"), sep = "_") %>%
    tidyr::unite(bin, AS, TE, sep = "_") %>% select(-serial)

gene.go_ko <- dplyr::bind_cols(bins_from_annotations, gene.go_ko)

# Define UI for application that shoots this table
ui <- fluidPage(

    # Application title
    titlePanel("TE Biogas: bin annotations"),
    tags$img(src = "SLUBI_logo_smaller.png", align="right"),

    # Sidebar with checkboxes for displaying columns 
    sidebarLayout(
        sidebarPanel(
            checkboxGroupInput("show_vars", "Columns/annotations to show:",
                               names(gene.go_ko), selected = names(gene.go_ko))
            # sliderInput("bins",
            #             "Number of bins:",
            #             min = 1,
            #             max = 50,
            #             value = 30)
        ),

        DT::dataTableOutput("gene.go_ko.table")
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$gene.go_ko.table = DT::renderDataTable({
        DT::datatable(gene.go_ko[, input$show_vars, drop=FALSE],
                      filter = list(position="top",
                                    clear=TRUE,
                                    plain=FALSE),
                      extensions = c("Buttons",
                                     "ColReorder",
                                     "FixedHeader"),
                      options = list(
                          dom = "Blfrtip",
                          colReorder = TRUE,
                          fixedHeader = TRUE,
                          pageLength = 50,
                          lengthMenu = c(5, 10, 25, 50, 100, 200, 500, 1000),
                          buttons = list(
                              # list(extend = "colvis", columns = 1:ncol(.)),
                              c('copy', 'csv', 'excel'))
                      )
        )
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
