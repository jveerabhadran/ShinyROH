# This is the main file for shinyROH application.
# This app visualizes and analyzes Run of Homozygosity (ROH) data for the user. 
# It generates violin plots of ROH lengths, categorized by chromosome and ROH class, and allows for the comparison of data between the two users. 
# The app also combines the datasets for a side-by-side analysis and displays the results in interactive data tables. 
# Author: Jyothishree
# Version: 1.00

# Load necessary libraries
library(shiny)
library(data.table)
library(gtools)
library(tidyr)
library(dplyr)
library(mclust)
library(MASS)
library(ggplot2)
library(plotly)
library(DT)
library(future.apply)

# Set the file size limits
# Set max file upload size to 200MB
options(shiny.maxRequestSize = 200 * 1024^2) 
# Set the maximum size of global variables for parallel processing
options(future.globals.maxSize =1024^3)

# Set up parallel processing with multiple workers
plan(multisession, workers = parallel::detectCores() - 2)

# Define UI layout
ui <- fluidPage(
  titlePanel("ShinyROH"), # Name of the app
  
  sidebarLayout(
    sidebarPanel(
      # File upload inputs and numeric inputs for user parameters
      h3("Upload SNP Data"),
      fileInput("file1", "Choose SNP Data File", accept = c(".txt")),
      fileInput("file2", "Choose SNP Data File (Optional)", accept = c(".txt")),
      
      h3("Parameters"),
      numericInput("n_boot", "Number of Bootstrap Iterations", value = 10, min = 1),
      numericInput("num_snps", "Number of SNPs per Window", value = 50, min = 1),
      numericInput("step_snps", "Step SNPs for Sliding Window", value = 10, min = 1),
      numericInput("error_rate", "Genotyping Error Rate", value = 0.01, min = 0, max = 1),
      actionButton("process_btn", "Process Data"),
      
      hr(),
      
      h3("Download Processed Data"),
      downloadButton("download_lod_data_t1", "LOD Data (User 1)"),
      downloadButton("download_snp_roh_freq_t1", "SNP ROH Freq (User 1)"),
      downloadButton("download_chr_roh_freq_t1", "Chromosome ROH Freq (User 1)"),
      downloadButton("download_lod_data_t2", "LOD Data (User 2)"),
      downloadButton("download_snp_roh_freq_t2", "SNP ROH Freq (User 2)"),
      downloadButton("download_chr_roh_freq_t2", "Chromosome ROH Freq (User 2)"),
      
      hr(),
      
      h3("Worldwide ROH Plotting"),
      # Worldwide Data File Input and Button
      fileInput("file1_worldwide", "SNP-ROH Freq 1", accept = c(".txt")),
      fileInput("file2_worldwide", "SNP-ROH Freq 2 (Optional)", accept = c(".txt")),
      actionButton("process_worldwide_btn", "Process Worldwide Data")
    ),
    
    
    mainPanel(
      tabsetPanel(
        tabPanel("Individual Plots", 
                 plotlyOutput("plot_t1"), 
                 plotlyOutput("plot_t2")),
        tabPanel("Combined Plot", plotlyOutput("plot_combined")),
        tabPanel("Worldwide ROH Plot", plotlyOutput("plot_worldwide")),
        tabPanel("Data Table", DTOutput("data_table"))
      )
    )
  )
)

# Define server logic for handling file input, data processing, and rendering outputs
server <- function(input, output, session) {
  
  # Worldwide ROH Distribution Plot
  # Define function to reshape the data for plotting
  reshape_data <- function(data, class_name) {
    data_long <- data %>%
      select(-chromosome, -position, -marker) %>%
      pivot_longer(cols = AFRICA:All, names_to = "Region", values_to = "ROH_Frequency") %>%
      mutate(Class = class_name)
    return(data_long)
  }
  
  # Load inbuilt data files for Classes A, B, and C
  classA_path <- "C:/Users/jyoth/OneDrive/Desktop/PG/BINP29/Population_Genetics/App/data/classA.txt"
  classB_path <- "C:/Users/jyoth/OneDrive/Desktop/PG/BINP29/Population_Genetics/App/data/classB.txt"
  classC_path <- "C:/Users/jyoth/OneDrive/Desktop/PG/BINP29/Population_Genetics/App/data/classC.txt"
  
  # Check if files exist at the specified location
  if (!file.exists(classA_path) | !file.exists(classB_path) | !file.exists(classC_path)) {
    stop("One or more class files (classA.txt, classB.txt, classC.txt) are missing in the data directory!")
  }
  
  # Load the inbuilt data into data frames
  classA <- read.table(classA_path, header = TRUE, sep = "\t")
  classB <- read.table(classB_path, header = TRUE, sep = "\t")
  classC <- read.table(classC_path, header = TRUE, sep = "\t")
  
  # Reshape the inbuilt data for plotting
  classA_long <- reshape_data(classA, "Class A")
  classB_long <- reshape_data(classB, "Class B")
  classC_long <- reshape_data(classC, "Class C")
  
  # Process User input data when the button is clicked
  observeEvent(input$process_worldwide_btn, {
    
    req(input$file1_worldwide)  # Ensure file1 is uploaded
    
    tryCatch({
      
      # Process User Data 1
      user_worldwide_data1 <- future_lapply(input$file1_worldwide$datapath, function(file) {
        read.table(file, header = TRUE, sep = "\t") %>%
          mutate(ROH_Frequency = roh_freq) %>%
          select(Class = roh_class, ROH_Frequency) %>%
          mutate(Class = gsub("Class A \\(Short ROH\\)", "Class A", Class),
                 Class = gsub("Class B \\(Intermediate ROH\\)", "Class B", Class),
                 Class = gsub("Class C \\(Long ROH\\)", "Class C", Class),
                 Region = "User Data 1")
      })
      
      # Process User Data 2 (if provided)
      user_worldwide_data2 <- NULL
      if (!is.null(input$file2_worldwide)) {
        user_worldwide_data2 <- future_lapply(input$file2_worldwide$datapath, function(file) {
          read.table(file, header = TRUE, sep = "\t") %>%
            mutate(ROH_Frequency = roh_freq) %>%
            select(Class = roh_class, ROH_Frequency) %>%
            mutate(Class = gsub("Class A \\(Short ROH\\)", "Class A", Class),
                   Class = gsub("Class B \\(Intermediate ROH\\)", "Class B", Class),
                   Class = gsub("Class C \\(Long ROH\\)", "Class C", Class),
                   Region = "User Data 2")
        })
      }
      
      # Check if User Data 1 is valid
      if (is.null(user_worldwide_data1) || length(user_worldwide_data1) == 0) {
        stop("User Data 1 is missing or invalid.")
      }
      
      # Combine both user data sets (only if file2 exists and is not NULL)
      if (is.null(user_worldwide_data2) || length(user_worldwide_data2) == 0) {
        user_worldwide_data_combined <- user_worldwide_data1  # Only User Data 1
      } else {
        user_worldwide_data_combined <- bind_rows(user_worldwide_data1, user_worldwide_data2)  # Both User Data 1 and User Data 2
      }
      
      # Combine reshaped inbuilt data with user data after user input
      all_worldwide_data <- bind_rows(classA_long, classB_long, classC_long, user_worldwide_data_combined)
      all_worldwide_data <- as.data.frame(all_worldwide_data)
      all_worldwide_data_clean <- na.omit(all_worldwide_data)
      
      # Compute max ROH frequency for scaling y-axis in the plot
      max_value <- max(all_worldwide_data$ROH_Frequency, na.rm = TRUE)
      
      # Render the worldwide ROH frequency distribution plot
      output$plot_worldwide <- renderPlotly({
        plot <- ggplot(all_worldwide_data_clean, aes(x = Region, y = ROH_Frequency, fill = Region)) +
          geom_violin(trim = TRUE, scale = "width", adjust = 2, alpha = 0.8) + 
          geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", alpha = 0.6) +
          facet_grid(rows = vars(Class), scales = "free_y") +
          scale_y_continuous(limits = c(0, max_value * 1.5)) +  
          labs(title = "ROH Frequency Distribution by Geographical Location and Class", 
               x = "Geographical Location", 
               y = "ROH Frequency") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1), 
                text = element_text(size = 9), 
                legend.position = "bottom", 
                legend.text = element_text(size = 9), 
                legend.key.size = unit(0.5, "cm")) +
          guides(fill = guide_legend(nrow = 1))
       
        ggplotly(plot) # Convert to interactive plot
        
      })
      
    }, error = function(e) {
      showNotification(paste("Error: ", e$message), type = "error")
      })
    })
  
  # Function to process SNP data
  process_snp_data <- function(file) {
    data <- fread(file$datapath, header = TRUE, sep = "\t")
    data <- na.omit(data)
    data[, genotype := fifelse(allele1 == allele2, 0, 1)]
    data <- data[chromosome %in% 1:22]
    
    # Bootstrap function for allele frequencies
    bootstrap_allele_freq <- function(data, n_boot = 10) {
      boot_res <- future_lapply(1:n_boot, function(i) {
        sample_data <- data[sample(.N, .N, replace = TRUE)]  # Resample data with replacement
        sample_data[, .(
          pA = (sum(allele1 == "A") + sum(allele2 == "A")) / (2 * .N),
          pG = (sum(allele1 == "G") + sum(allele2 == "G")) / (2 * .N),
          pC = (sum(allele1 == "C") + sum(allele2 == "C")) / (2 * .N),
          pT = (sum(allele1 == "T") + sum(allele2 == "T")) / (2 * .N)
        ), by = "rsid"]
      }, future.seed = TRUE)
      
      rbindlist(boot_res)[, lapply(.SD, mean), by = "rsid"]  # Aggregate the mean allele frequencies
    }
    
    allele_freqs <- bootstrap_allele_freq(data)
    
    # Return both data and bootstrapped allele frequencies
    return(list(data = data, allele_freqs = allele_freqs))
  }
  
  # Function to compute LOD scores
  compute_lod_scores_parallel <- function(data, allele_freqs, num_snps = 50, step_snps = 10, error_rate = 0.01) {
    results <- future_lapply(unique(data$chromosome), function(chr) {
      
      chr_data <- data[chromosome == chr][order(position)]
      chr_data <- merge(chr_data, allele_freqs, by = "rsid", all.x = TRUE)
      
      # Precompute allele frequencies for this chromosome
      chr_data[, `:=`(
        pA = fifelse(is.na(pA) | pA == 0, 0.01, pA),
        pG = fifelse(is.na(pG) | pG == 0, 0.01, pG),
        pC = fifelse(is.na(pC) | pC == 0, 0.01, pC),
        pT = fifelse(is.na(pT) | pT == 0, 0.01, pT)
      )]
      
      start_indices <- seq(1, nrow(chr_data) - num_snps, by = step_snps)
      lod_list <- vector("list", length(start_indices))
      
      for (i in seq_along(start_indices)) {
        segment_data <- chr_data[start_indices[i]:(start_indices[i] + num_snps - 1),]
        segment_data[, `:=`(
          p1 = fifelse(allele1 == "A", pA, 
                       fifelse(allele1 == "G", pG, 
                               fifelse(allele1 == "C", pC, pT))),
          p2 = fifelse(allele2 == "A", pA, 
                       fifelse(allele2 == "G", pG, 
                               fifelse(allele2 == "C", pC, pT)))
        )]
        
        segment_data[, `:=`(
          P_ROH = fifelse(allele1 == allele2,  
                          (1 - error_rate) * p1 + error_rate * p1^2,  
                          2 * error_rate * p1 * p2),
          P_noROH = fifelse(allele1 == allele2, p1^2, 2 * p1 * p2)
        )]
        
        segment_data[, lod := log10(P_ROH / P_noROH)]
        segment_data[, lod := fifelse(is.nan(lod), 0, lod)]
        
        lod_score <- sum(segment_data$lod, na.rm = TRUE)
        lod_list[[i]] <- data.table(
          chromosome = chr,
          start_snp = segment_data$rsid[1], 
          end_snp = segment_data$rsid[nrow(segment_data)], 
          lod_score = lod_score,
          roh_length = as.numeric(abs(segment_data$position[nrow(segment_data)] - segment_data$position[1]))
        )
      }
      
      return(rbindlist(lod_list, use.names = TRUE, fill = TRUE))
    }, future.seed = TRUE)
    
    return(rbindlist(results, use.names = TRUE, fill = TRUE))
  }
  
  # Process the SNP data when the button is clicked
  observeEvent(input$process_btn, {
    req(input$file1)
    
    # Validate file format
    file_extension <- tools::file_ext(input$file1$name)
    if (file_extension != "txt") {
      showModal(modalDialog(
        title = "Invalid file format",
        "Please upload a .txt file.",
        easyClose = TRUE,
        footer = NULL
      ))
      return(NULL)
    }
    
    # Process User 1 (t1) SNP data
    withProgress(message = 'Processing SNP data for User 1...', {
      data1 <- process_snp_data(input$file1)  # Process SNP data for User 1
      lod_data_t1 <- compute_lod_scores_parallel(data1$data, data1$allele_freqs)  # Compute LOD scores
    })
    
    
    # Process SNP data for User 2 (if exists)
    lod_data_t2 <- NULL
    output_t2 <- NULL  # Initialize output for User 2
    if (!is.null(input$file2)) {
      req(input$file2)
      
      # Validate file format for User 2
      file_extension2 <- tools::file_ext(input$file2$name)
      if (file_extension2 != "txt") {
        showModal(modalDialog(
          title = "Invalid file format",
          "Please upload a .txt file for User 2.",
          easyClose = TRUE,
          footer = NULL
        ))
        return(NULL)
      }
      
      withProgress(message = 'Processing SNP data for User 2...', {
        data2 <- process_snp_data(input$file2)
        lod_data_t2 <- compute_lod_scores_parallel(data2$data, data2$allele_freqs)
      })
      
      # Store processed results for User 2 in output_t2
      output_t2 <- list(
        snp_roh_freq_t2 = snp_roh_freq_t2,
        chr_roh_freq_t2 = chr_roh_freq_t2
      )
    }
    
    # Add Download Handlers for LOD Data, SNP, and Chromosome ROH Frequency
    output$download_lod_data_t1 <- downloadHandler(
      filename = function() {
        paste("lod_data_t1_", Sys.Date(), ".txt", sep = "")
      },
      content = function(file) {
        write.table(lod_data_t1, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    output$download_snp_roh_freq_t1 <- downloadHandler(
      filename = function() {
        paste("snp_roh_freq_t1_", Sys.Date(), ".txt", sep = "")
      },
      content = function(file) {
        write.table(snp_roh_freq_t1, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    output$download_chr_roh_freq_t1 <- downloadHandler(
      filename = function() {
        paste("chr_roh_freq_t1_", Sys.Date(), ".txt", sep = "")
      },
      content = function(file) {
        write.table(chr_roh_freq_t1, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    if (!is.null(lod_data_t2)) {
      output$download_lod_data_t2 <- downloadHandler(
        filename = function() {
          paste("lod_data_t2_", Sys.Date(), ".txt", sep = "")
        },
        content = function(file) {
          write.table(lod_data_t2, file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
      )
      
      output$download_snp_roh_freq_t2 <- downloadHandler(
        filename = function() {
          paste("snp_roh_freq_t2_", Sys.Date(), ".txt", sep = "")
        },
        content = function(file) {
          write.table(snp_roh_freq_t2, file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
      )
      
      output$download_chr_roh_freq_t2 <- downloadHandler(
        filename = function() {
          paste("chr_roh_freq_t2_", Sys.Date(), ".txt", sep = "")
        },
        content = function(file) {
          write.table(chr_roh_freq_t2, file, sep = "\t", row.names = FALSE, quote = FALSE)
        }
      )
    }
    
    # Compute LOD scores for both users
    lod_threshold_t1 <- estimate_lod_threshold(lod_data_t1)
    lod_threshold_t2 <- if (!is.null(lod_data_t2)) estimate_lod_threshold(lod_data_t2) else NA
    
    # Apply KDE threshold to classify ROH
    lod_data_t1[, roh_class := ifelse(lod_score >= lod_threshold_t1, "ROH", "Non-ROH")]
    if (!is.null(lod_data_t2)) {
      lod_data_t2[, roh_class := ifelse(lod_score >= lod_threshold_t2, "ROH", "Non-ROH")]
    }
    
    # Compute ROH frequencies per chromosome for User 1 (T1)
    autosomal_genome_size <- 2.8e9  # Human autosome size in base pairs
    roh_total_length_t1 <- sum(lod_data_t1$roh_length, na.rm = TRUE)
    roh_frequency_t1 <- roh_total_length_t1 / autosomal_genome_size
    
    # If User 2 data exists, compute ROH frequencies for User 2 (T2)
    if (!is.null(lod_data_t2)) {
      roh_total_length_t2 <- sum(lod_data_t2$roh_length, na.rm = TRUE)
      roh_frequency_t2 <- roh_total_length_t2 / autosomal_genome_size
    }
    
    # Gaussian Mixture Model (GMM) for ROH classification for User 1 (T1)
    gmm_model_t1 <- Mclust(lod_data_t1$roh_length, G = 3)
    
    # For User 2 (T2), run GMM only if data exists
    if (!is.null(lod_data_t2)) {
      gmm_model_t2 <- Mclust(lod_data_t2$roh_length, G = 3)
    }
    
    # Extract and sort component means for t1 and t2
    component_means_t1 <- sort(gmm_model_t1$parameters$mean)
    
    # For User 2 (T2), extract and sort means only if model exists
    if (!is.null(gmm_model_t2)) {
      component_means_t2 <- sort(gmm_model_t2$parameters$mean)
    }
    
    # Calculate cutoffs for User 1 (T1)
    cutoff_A_B_t1 <- (component_means_t1[1] + component_means_t1[2]) / 2
    cutoff_B_C_t1 <- (component_means_t1[2] + component_means_t1[3]) / 2
    
    # Calculate cutoffs for User 2 (T2) only if component_means_t2 is available
    if (!is.null(component_means_t2)) {
      cutoff_A_B_t2 <- (component_means_t2[1] + component_means_t2[2]) / 2
      cutoff_B_C_t2 <- (component_means_t2[2] + component_means_t2[3]) / 2
    }
    
    # Classify ROH length based on GMM for User 1 (T1)
    lod_data_t1[, roh_class := fifelse(roh_length < cutoff_A_B_t1, "Class A (Short ROH)",
                                       fifelse(roh_length < cutoff_B_C_t1, "Class B (Intermediate ROH)", 
                                               "Class C (Long ROH)"))]
    
    # For User 2 (T2), classify ROH length only if data exists
    if (!is.null(lod_data_t2)) {
      lod_data_t2[, roh_class := fifelse(roh_length < cutoff_A_B_t2, "Class A (Short ROH)",
                                         fifelse(roh_length < cutoff_B_C_t2, "Class B (Intermediate ROH)", 
                                                 "Class C (Long ROH)"))]
    }
    
    # Prepare output for displaying
    output_t1 <- list(
      lod_data_t1 = lod_data_t1,
    )
    
    # Show data to the user as data table
    output$lod_data_t1_table <- renderDT({
      datatable(output_t1$lod_data_t1)
    })
    
    # Show LOD data for User 2 (if exists)
    if (!is.null(lod_data_t2)) {
      output$lod_data_t2_table <- renderDT({
        datatable(lod_data_t2)
      })
    }
    
    # Plots
    # Render plot for User 1 (lod_data_t1)
    output$plot_t1 <- renderPlotly({
      # Create a violin plot with a boxplot overlay for User 1's data
      ggplotly(
        ggplot(lod_data_t1, aes(x = as.factor(chromosome), y = abs(roh_length) / 1e6, fill = roh_class)) +
          geom_violin(trim = TRUE, scale = "width", adjust = 2, alpha = 0.8) +  
          geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 0.5) +
          scale_fill_manual(values = c("Class A (Short ROH)" = "lightblue", "Class B (Intermediate ROH)" = "pink", "Class C (Long ROH)" = "lightgreen")) +  
          labs(title = "Violin Plot of ROH Lengths for User 1", x = "Chromosome", y = "Total Length of ROH (Mb)", fill = "ROH Class") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
    })
    
    # Render plot for User 2 (lod_data_t2), same as User 1 plot with different data
    output$plot_t2 <- renderPlotly({
      ggplotly(
        ggplot(lod_data_t2, aes(x = as.factor(chromosome), y = abs(roh_length) / 1e6, fill = roh_class)) +
          geom_violin(trim = TRUE, scale = "width", adjust = 2, alpha = 0.8) +  
          geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 0.5) +
          scale_fill_manual(values = c("Class A (Short ROH)" = "lightblue", "Class B (Intermediate ROH)" = "pink", "Class C (Long ROH)" = "lightgreen")) +  
          labs(title = "Violin Plot of ROH Lengths for User 2", x = "Chromosome", y = "Total Length of ROH (Mb)", fill = "ROH Class") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
    })
    
    # Combined plot: combining data from both User 1 and User 2
    output$plot_combined <- renderPlotly({
      # Combine User 1 and User 2 data into one data frame
      combined_lod_data <- rbindlist(list(lod_data_t1, lod_data_t2), fill = TRUE)
      # Add a 'user' column to identify the source of each dataset (User 1 or User 2)
      combined_lod_data[, user := rep(c("User 1 (T1)", "User 2 (T2)"), c(nrow(lod_data_t1), nrow(lod_data_t2)))]
      
      # Convert chromosome to factor for proper ordering
      combined_lod_data[, chromosome := factor(chromosome, levels = mixedsort(unique(chromosome)))]
      
      # Define color palette for the users
      user_colors <- c("User 1 (T1)" = "green", "User 2 (T2)" = "red")
      
      # Create a violin plot with boxplot overlays for both users, facet by ROH class
      p_combined <- ggplot(combined_lod_data, aes(x = as.factor(chromosome), y = abs(roh_length) / 1e6, fill = user)) +
        geom_violin(trim = TRUE, scale = "width", adjust = 2, alpha = 0.8) +
        geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 0.5) +
        facet_grid(rows = vars(roh_class), scales = "free_y") +
        scale_fill_manual(values = user_colors) +
        labs(title = "Violin Plot of ROH Lengths by Class and User", 
             x = "Chromosome", 
             y = "Total Length of ROH (Mb)",
             fill = "User") +
        
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid.major.x = element_blank(),
              panel.background = element_rect(fill = "white", color = NA)) +
        
        
        # Add vertical dashed lines to separate chromosomes visually
        geom_vline(xintercept = seq(1.5, length(unique(combined_lod_data$chromosome)) - 0.5, by = 1),
                   linetype = "dashed", color = "grey50", alpha = 0.5)
      
      # Convert ggplot to plotly object for interactivity
      ggplotly(p_combined)
    })
    
    # Render data table with added 'user' column added
    output$data_table <- renderDT({
      # Copy data for User 1 and add a 'user' column
      user_data_t1 <- copy(lod_data_t1)
      user_data_t1$user <- "User 1"
      
      # If User 2 data exists, process and combine it with User 1 data
      if (!is.null("lod_data_t2") && nrow(lod_data_t2) > 0) {
        user_data_t2 <- copy(lod_data_t2)
        user_data_t2$user <- "User 2"
        combined_data <- rbind(user_data_t1, user_data_t2) # Combine both datasets
      } else {
        combined_data <- user_data_t1 # Only User 1 data if no User 2 data
      }
      
      # Render the combined data table with a page length of 5 rows per page
      datatable(combined_data, options = list(pageLength = 5))
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)