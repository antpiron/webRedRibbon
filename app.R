library(shiny)
library(RedRibbon)
library(data.table)
library(stringr)

maxFileUploadMiB <- 50 

ui <- fluidPage(
    tags$head(tags$link(rel = "stylesheet", type="text/css", href="RedRibbon.css")),
    titlePanel("RedRibbon online generator"),
    sidebarLayout(

        ## Sidebar panel for inputs ----
        sidebarPanel(
            h1("Citation"),
            p(tags$cite("Anthony Piron, Florian Szymczak, Theodora Papadopoulou, Maria Inês Alvelos, Matthieu Defrance, Tom Lenaerts, Décio L Eizirik, Miriam Cnop. RedRibbon: A new rank-rank hypergeometric overlap for gene and transcript expression signatures. Life Science Alliance. 2023 Dec 8;7(2):e202302203."),
              a("doi: 10.26508/lsa.202302203.", href = "https://doi.org/10.26508/lsa.202302203")),
            
            hr(),
            h1("Data"),
            selectInput(
                inputId = "DataSource",
                label = NULL,
                choices = list( "File" = "file", "Demo" = "demo", "DESeq2" = "DESeq2"),
                selected = "file",
                multiple = FALSE,
                selectize = TRUE,
                width = NULL,
                size = NULL
            ),
            
            conditionalPanel(
                condition = "input.DataSource == 'file'",
                div(title = paste("The file should be a three columns file. `id` column is the identifier, `a` the first statistic (e.g., DESeq2 stat or direction signed P-Value), `b` the second statistic. The column are separated either by spaces, tabulations, commas, pipe, semicolon, or double points. Tabulations are preferred. Example:\n\nid     \ta      \tb\ngene1 \t-0.2 \t-0.3 \ngene42\t10.6 \t-11.6\n...\n\nMaximal file size is", maxFileUploadMiB, "MiB."),
                    fileInput("File", NULL))
            ),
            
            conditionalPanel(
                condition = "input.DataSource == 'demo'",
               
                    sliderInput(inputId = "N",
                                label = "List size:",
                                min = 100,
                                max = 2000,
                                value = 500)
            ),

            conditionalPanel(
                condition = "input.DataSource == 'DESeq2'",
                fileInput("DESeq2A", "a:"),
                fileInput("DESeq2B", "b:"),
                selectInput(
                    inputId = "DESeq2Column",
                    label = "Statistic:",
                    ## TODO: if pvalue: do not split in quadrants, do whole
                    choices = list( "log Fold Change (column named `log2FoldChange`)" = "log2FoldChange",
                                   "Statistic (column named `stat`)" = "stat",
                                   "P-value (column named `pvalue`)" = "pvalue"),
                    selected = "log2FoldChange",
                    multiple = FALSE,
                    selectize = TRUE,
                    width = NULL,
                    size = NULL
                )
            ),

            
            hr(),
            h1("Algorithm"),
            div(title = "Evolutionnary algorithm is fast and highly accurate. Classic is as implemented in Plaisier et al., 2010.",
                selectInput(
                    inputId = "Algorithm",
                    label = NULL,
                    choices = list("Evolutionnary (Recommended)" = "ea", "Classic" = "classic"),
                    selected = "ea",
                    multiple = FALSE,
                    selectize = TRUE,
                    width = NULL,
                    size = NULL
                )),
            
            div(title = "Compute the adjusted P-value. This can take some time to generate depending on the number of permutations. We advise to experiments first without permutation and once, you have a satifactory result to activate the option.",
                checkboxInput(inputId = "Permutation", "Permute", value = FALSE, width = NULL)),
            checkboxInput(inputId = "Quadrants", "Split into quadrants", value = TRUE, width = NULL),
            checkboxInput(inputId = "ShowData", "Show data", value = TRUE, width = NULL),

            hr(),
            h1("Plot"),
            textInput(inputId = "ListA", label = "Name of the first list:", value = "a"),
            textInput(inputId = "ListB", label = "Name of the second list:", value = "b"),
            sliderInput(inputId = "FontSize",
                        label = "Font size:",
                        min = 5,
                        max = 100,
                        value = 20),
            sliderInput("Width", "Plot width (Points):", min = 512, max = 8192, value = 1024)
            ## sliderInput("dpi", "DPI:", min = 48, max = 1200, value = 96)
            
        ),
        
        ## Main panel for displaying outputs ----
        mainPanel(
            
            ## Output: Histogram ----
            plotOutput(outputId = "ggRedRibbon", width = "100%", height = "100%"),
            tableOutput(outputId = "Results"),
            dataTableOutput(outputId = "table")
            
        )
    )
)

gen.dt <- function (n)
{
    set.seed(42)
    n2 <- n %/% 2
    n4 <- n %/% 4
    a <- (1:n) - n2
    b <- a
    a[n4:(n4+n2-1)] <- rnorm(n2, sd=100)
    b[n4:(n4+n2-1)] <- rnorm(n2, sd=100)
    df <- data.table(
        id = paste0("gene", 1:n),
        a = a,
        b = b)

    return(df)
}

server <- function(input, output, session)
{
    dpi <- 96 ## reactive({ as.integer(input$dpi) } )

    width <- debounce(reactive({ as.integer(input$Width) } ), 1000)
    height <- reactive({ as.integer(768 * width() / 1024)} )
    size_factor <- reactive({
        width() / 1024
    })


    base_size <- debounce(reactive({ max(5, as.integer(input$FontSize * size_factor() )) }), 1000)

    demo_n <- debounce( reactive({ as.integer(input$N)  }), 1000)

    ## repel_force  <- reactive({ as.integer(0.5 * sqrt(width()^2 + height()^2) * size_factor())  })

    observe({
      # periodically collect
      invalidateLater(1000 * 30,session)
      gc()
    })
 
    genDT <- reactive({
        if ("demo" == input$DataSource)
        {
            gen.dt(demo_n())
        } else if ("file" == input$DataSource)
        {
            req(input$File)
            dt <- fread(file=input$File$datapath)
            ## TODO: check the columns
            dt
        } else if ("DESeq2" == input$DataSource)
        {
            req(input$DESeq2A)
            req(input$DESeq2B)

            ## print(input$DESeq2A)
            ## print(input$DESeq2B)

            DESeq2A <- fread(file=input$DESeq2A$datapath)
            DESeq2A <- DESeq2A[,c(1, grep(input$DESeq2Column, colnames(DESeq2A))), with=FALSE]
            colnames(DESeq2A) <- c("id", "a")
            
            DESeq2B <- fread(file=input$DESeq2B$datapath)
            DESeq2B <- DESeq2B[,c(1, grep(input$DESeq2Column, colnames(DESeq2B))), with=FALSE]
            colnames(DESeq2B) <- c("id", "b")
            
            dt <- merge(DESeq2A, DESeq2B, by="id")
            
            dt
        } else
            req(FALSE)
    })
    
    rr <- reactive({
        dt <- genDT()
        rr <- RedRibbon(dt, enrichment_mode="hyper-two-tailed")
        
        rr
    })

    quad <- reactive({
        ## print(input$Algorithm)
        quadrants(rr(), algorithm=input$Algorithm, permutation=input$Permutation, whole=! input$Quadrants)
    })

    output$ggRedRibbon <- renderPlot({

        ggRedRibbon(rr(), quadrants=quad(), show.quadrants = input$Quadrants, labels = c(input$ListA, input$ListB), base_size = base_size(), repel.force = 200) +
            coord_fixed(ratio = 1, clip = "off")  +
            theme(legend.key.height = unit(size_factor() * 0.5, 'cm'),
                  legend.key.width = unit(size_factor() * 1, 'cm'))
        
    }, res = dpi, width = width, height = height)


    output$Results <- renderTable({
        dt <- genDT()
        qd <- quad()

        ## print(names(qd$downdown))
        df <- t(sapply(names(qd),
                       function(dir)
                       {
                           c(
                               dir,
                               qd[[dir]]$pvalue,
                               str_wrap(paste(dt[qd[[dir]]$positions,]$id, sep = ",", collapse = ", "))
                           )
                       }))
        colnames(df) <- c("Direction", "P-value", "IDs")
        df
        ## data.table(
        ##     Direction=c("Down-Down", "Up-Up", "Down-Up", "Up-Down"),
        ##     "- log P-value" = c( qd$downdown$pvalue, qd$upup$pvalue, qd$downup$pvalue, qd$updown$pvalue ),
        ##     IDs=c(str_wrap(paste(dt[qd$downdown$positions,]$id, sep = ",", collapse = ", ")),
        ##           str_wrap(paste(dt[qd$upup$positions,]$id, sep = ",", collapse = ", ")),
        ##           str_wrap(paste(dt[qd$downup$positions,]$id, sep = ",", collapse = ", ")),
        ##           str_wrap(paste(dt[qd$updown$positions,]$id, sep = ",", collapse = ", ")))
        ## )
        
    })
    
    output$table <- renderDataTable({
        if (input$ShowData)
        {
            dt  <- genDT()
            dt
        }
    }, options = list(pageLength = 10))
}

options(shiny.mxRequestSize = maxFileUploadMiB * 1024^2)
shinyApp(ui, server)
