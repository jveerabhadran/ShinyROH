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
options(future.globals.maxSize = 1024^3)

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
        tabPanel("Worldwide ROH Plot", plotOutput("plot_worldwide")),
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
  
  # URL of the file on GitHub
  classA_url <- "https://raw.githubusercontent.com/jveerabhadran/BINP29_Population_Genetics/main/data/class_data_files/classA.txt"
  classB_url <- "https://raw.githubusercontent.com/jveerabhadran/BINP29_Population_Genetics/main/data/class_data_files/classB.txt"
  classC_url <- "https://raw.githubusercontent.com/jveerabhadran/BINP29_Population_Genetics/main/data/class_data_files/classC.txt"
  
  # Define the local path to save the downloaded file
  classA_path <- "classA.txt"
  classB_path <- "classB.txt"
  classC_path <- "classC.txt"
  
  # Download the files
  download.file(classA_url, classA_path)
  download.file(classB_url, classB_path)
  download.file(classC_url, classC_path)
  
  # Read the downloaded files into R
  classA <- read.table(classA_path, header = TRUE, sep = "\t")
  classB <- read.table(classB_path, header = TRUE, sep = "\t")
  classC <- read.table(classC_path, header = TRUE, sep = "\t")
  
  # Reshape the data for plotting
  classA_long <- reshape_data(classA, "Class A")
  classB_long <- reshape_data(classB, "Class B")
  classC_long <- reshape_data(classC, "Class C")
  
  # Process User input data when the button is clicked
  observeEvent(input$process_worldwide_btn, {
    
    req(input$file1_worldwide)  # Ensure file1 is uploaded
    
    tryCatch({
      
      # Process User Data 1
      user_worldwide_data1 <- future_lapply(input$file1_worldwide$datapath, function(file) {
        data <- read.table(file, header = TRUE, sep = "\t")
        
        # Check if the file is empty
        if (nrow(data) == 0) {
          showNotification("Error: User Data 1 file is empty.", type = "error")
          return(NULL)
        }
        
        data %>%
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
          data <- read.table(file, header = TRUE, sep = "\t")
          
          # Check if the file is empty
          if (nrow(data) == 0) {
            showNotification("Error: User Data 2 file is empty.", type = "error")
            return(NULL)
          }
          
          data %>%
            mutate(ROH_Frequency = roh_freq) %>%
            select(Class = roh_class, ROH_Frequency) %>%
            mutate(Class = gsub("Class A \\(Short ROH\\)", "Class A", Class),
                   Class = gsub("Class B \\(Intermediate ROH\\)", "Class B", Class),
                   Class = gsub("Class C \\(Long ROH\\)", "Class C", Class),
                   Region = "User Data 2")
        })
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
      output$plot_worldwide <- renderPlot({
        ggplot(all_worldwide_data_clean, aes(x = Region, y = ROH_Frequency, fill = Region)) +
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
      })
      
    }, error = function(e) {
      showNotification(paste("Error: ", e$message), type = "error")
    })
  })
  
  # Function to process SNP data
  process_snp_data <- function(file) {
    if (is.null(file) || !file.exists(file$datapath)) {
      showModal(modalDialog(
        title = "Error",
        "File not found! Please upload a valid file.",
        easyClose = TRUE
      ))
      return(NULL)  # Stop processing but keep the app running
    }
    
    # Attempt to read the file
    data <- tryCatch({
      fread(file$datapath, header = TRUE, sep = "\t")
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        "Incorrect file format! See the GitHub README.md file for the required format.",
        easyClose = TRUE
      ))
      return(NULL)  # Stop processing but keep the app running
    })
    
    # Ensure the file was read successfully
    if (is.null(data)) return(NULL)
    
    # Check for required columns
    required_cols <- c("rsid", "allele1", "allele2", "chromosome", "position")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      showModal(modalDialog(
        title = "Error",
        paste("Incorrect file format! Missing required columns:", paste(missing_cols, collapse = ", "), 
              "<br> See the GitHub README.md file for details."),
        easyClose = TRUE,
        footer = modalButton("OK")
      ))
      return(NULL)
    }
    
    # Convert columns to correct data types
    data[, `:=`(
      chromosome = as.integer(chromosome),
      position = as.numeric(position)
    )]
    
    # Remove rows with missing values
    data <- na.omit(data)
    
    # Validate autosomal chromosome range (1-22)
    if (any(!data$chromosome %in% 1:22)) {
      showModal(modalDialog(
        title = "Error",
        "Invalid chromosome values detected! Chromosome values must be between 1 and 22.",
        easyClose = TRUE
      ))
      return(NULL)
    }
    
    # Calculate genotype
    data[, genotype := fifelse(allele1 == allele2, 0, 1)]
    
    # Bootstrap function for allele frequencies
    bootstrap_allele_freq <- function(data, n_boot = 10) {
      if (nrow(data) < 10) {
        showModal(modalDialog(
          title = "Error",
          "Not enough data for bootstrapping. Please upload a larger dataset.",
          easyClose = TRUE
        ))
        return(NULL)
      }
      
      boot_res <- future_lapply(1:n_boot, function(i) {
        sample_data <- data[sample(.N, .N, replace = TRUE)]
        sample_data[, .(
          pA = (sum(allele1 == "A") + sum(allele2 == "A")) / (2 * .N),
          pG = (sum(allele1 == "G") + sum(allele2 == "G")) / (2 * .N),
          pC = (sum(allele1 == "C") + sum(allele2 == "C")) / (2 * .N),
          pT = (sum(allele1 == "T") + sum(allele2 == "T")) / (2 * .N)
        ), by = "rsid"]
      }, future.seed = TRUE)
      
      result <- rbindlist(boot_res, use.names = TRUE, fill = TRUE)
      
      if (nrow(result) == 0) {
        showModal(modalDialog(
          title = "Error",
          "Bootstrap failed to compute allele frequencies.",
          easyClose = TRUE
        ))
        return(NULL)
      }
      
      result[, lapply(.SD, mean), by = "rsid"]
    }
    
    # Compute allele frequencies
    allele_freqs <- tryCatch({
      bootstrap_allele_freq(data)
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("Bootstrap computation failed:", e$message),
        easyClose = TRUE
      ))
      return(NULL)
    })
    
    # Ensure allele_freqs is not NULL
    if (is.null(allele_freqs)) return(NULL)
    
    return(list(data = data, allele_freqs = allele_freqs))
  }
  
  
  # Function to compute LOD scores with error handling
  compute_lod_scores_parallel <- function(data, allele_freqs, num_snps = 50, step_snps = 10, error_rate = 0.01) {
    if (nrow(data) == 0) {
      stop("Error: Input data is empty.")
    }
    
    if (nrow(allele_freqs) == 0) {
      stop("Error: Allele frequency data is empty.")
    }
    
    results <- future_lapply(unique(data$chromosome), function(chr) {
      chr_data <- data[chromosome == chr][order(position)]
      
      if (nrow(chr_data) < num_snps) {
        warning(paste("Warning: Chromosome", chr, "has fewer SNPs than required. Skipping."))
        return(NULL)
      }
      
      chr_data <- merge(chr_data, allele_freqs, by = "rsid", all.x = TRUE)
      
      # Handle missing allele frequencies
      chr_data[, `:=`(
        pA = fifelse(is.na(pA) | pA == 0, 0.01, pA),
        pG = fifelse(is.na(pG) | pG == 0, 0.01, pG),
        pC = fifelse(is.na(pC) | pC == 0, 0.01, pC),
        pT = fifelse(is.na(pT) | pT == 0, 0.01, pT)
      )]
      
      lod_list <- list()
      start_indices <- seq(1, nrow(chr_data) - num_snps, by = step_snps)
      
      # Loop over the start indices to calculate LOD scores for each segment
      for (i in seq_along(start_indices)) {
        segment_data <- chr_data[start_indices[i]:(start_indices[i] + num_snps - 1), ]
        
        segment_data[, `:=`(
          p1 = fifelse(allele1 == "A", pA, 
                       fifelse(allele1 == "G", pG, 
                               fifelse(allele1 == "C", pC, pT))),
          p2 = fifelse(allele2 == "A", pA, 
                       fifelse(allele2 == "G", pG, 
                               fifelse(allele2 == "C", pC, pT)))
        )]
        
        # Calculate probabilities for ROH and non-ROH
        segment_data[, `:=`(
          P_ROH = fifelse(allele1 == allele2,  
                          (1 - error_rate) * p1 + error_rate * p1^2,  
                          2 * error_rate * p1 * p2),
          P_noROH = fifelse(allele1 == allele2, p1^2, 2 * p1 * p2)
        )]
        
        # Handle zero division and NaN issues
        segment_data[, lod := fifelse(P_noROH == 0, 0, log10(P_ROH / P_noROH))]
        segment_data[, lod := fifelse(is.nan(lod) | is.infinite(lod), 0, lod)]
        
        lod_score <- sum(segment_data$lod, na.rm = TRUE)
        
        lod_list[[i]] <- data.table(
          chromosome = chr,
          start_snp = segment_data$rsid[1], 
          end_snp = segment_data$rsid[nrow(segment_data)], 
          lod_score = lod_score,
          roh_length = abs(segment_data$position[nrow(segment_data)] - segment_data$position[1])
        )
      }
      
      return(rbindlist(lod_list, use.names = TRUE, fill = TRUE))
    }, future.seed = TRUE)
    
    final_results <- rbindlist(results, use.names = TRUE, fill = TRUE)
    
    if (nrow(final_results) == 0) {
      stop("Error: No valid LOD scores were computed.")
    }
    
    return(final_results)
  }
  
  # Function to calculate LOD threshold using KDE
  estimate_lod_threshold <- function(lod_scores) {
    tryCatch({
      if (is.null(lod_scores) || !("lod_score" %in% colnames(lod_scores))) {
        stop("Error: Invalid LOD scores data.")
      }
      kde_fit <- density(lod_scores$lod_score, kernel = "gaussian")
      lod_threshold <- kde_fit$x[which.max(kde_fit$y)]  # Mode of KDE as threshold
      return(lod_threshold)
    }, error = function(e) {
      showModal(modalDialog(title = "Error", paste("LOD threshold estimation failed:", e$message), easyClose = TRUE))
      return(NA)
    })
  }
  
  # Function to compute LOD scores
  compute_lod_scores_for_user <- function(data) {
    tryCatch({
      if (is.null(data) || !"data" %in% names(data) || !"allele_freqs" %in% names(data)) {
        stop("Invalid input data for LOD score calculation.")
      }
      lod_data <- compute_lod_scores_parallel(data$data, data$allele_freqs)
      return(lod_data)
    }, error = function(e) {
      showModal(modalDialog(title = "Error", paste("LOD score computation failed:", e$message), easyClose = TRUE))
      return(NULL)
    })
  }
  
  # Function to classify ROH based on LOD threshold
  classify_roh <- function(lod_data, lod_threshold) {
    tryCatch({
      if (is.null(lod_data) || !("lod_score" %in% colnames(lod_data))) {
        stop("Invalid LOD data for classification.")
      }
      lod_data[, roh_class := ifelse(lod_score >= lod_threshold, "ROH", "Non-ROH")]
      return(lod_data)
    }, error = function(e) {
      showModal(modalDialog(title = "Error", paste("ROH classification failed:", e$message), easyClose = TRUE))
      return(NULL)
    })
  }
  
  # Function to fit GMM model and classify ROH length
  fit_gmm_and_classify <- function(lod_data) {
    tryCatch({
      if (is.null(lod_data) || !"roh_length" %in% colnames(lod_data)) {
        stop("Invalid data for GMM classification.")
      }
      gmm_model <- Mclust(lod_data$roh_length, G = 3)
      
      component_means <- sort(gmm_model$parameters$mean)
      cutoff_A_B <- (component_means[1] + component_means[2]) / 2
      cutoff_B_C <- (component_means[2] + component_means[3]) / 2
      
      lod_data[, roh_class := fifelse(roh_length < cutoff_A_B, "Class A (Short ROH)",
                                      fifelse(roh_length < cutoff_B_C, "Class B (Intermediate ROH)", 
                                              "Class C (Long ROH)"))]
      return(lod_data)
    }, error = function(e) {
      showModal(modalDialog(title = "Error", paste("GMM classification failed:", e$message), easyClose = TRUE))
      return(NULL)
    })
  }
  
  # Function to compute chromosome-level ROH frequency
  compute_chr_roh_freq <- function(lod_data, autosomal_genome_size) {
    tryCatch({
      if (is.null(lod_data) || !("chromosome" %in% colnames(lod_data))) {
        stop("Invalid LOD data for chromosome ROH frequency.")
      }
      chr_roh_freq <- lod_data[, .(roh_count = .N, roh_class = unique(roh_class)), by = chromosome]
      chr_roh_freq[, chr_freq := roh_count / sum(roh_count), by = chromosome]
      return(chr_roh_freq)
    }, error = function(e) {
      showModal(modalDialog(title = "Error", paste("Chromosome ROH frequency computation failed:", e$message), easyClose = TRUE))
      return(NULL)
    })
  }
  
  # Function to compute SNP-based ROH frequency per chromosome and ROH class
  compute_snp_roh_freq <- function(lod_data) {
    if (!inherits(lod_data, "data.table")) {
      stop("Error: Input 'lod_data' must be a data.table.")
    }
    
    required_cols <- c("chromosome", "roh_class")
    missing_cols <- setdiff(required_cols, names(lod_data))
    
    if (length(missing_cols) > 0) {
      stop(paste("Error: Missing required column(s):", paste(missing_cols, collapse = ", ")))
    }
    
    if (nrow(lod_data) == 0) {
      warning("Warning: Input data.table is empty. Returning empty result.")
      return(data.table(chromosome = character(), roh_class = character(), snp_count = integer(), roh_freq = numeric()))
    }
    
    tryCatch({
      snp_roh_freq <- lod_data[, .(snp_count = .N, roh_class = unique(roh_class)), by = .(chromosome, roh_class)]
      snp_roh_freq[, roh_freq := snp_count / sum(snp_count), by = chromosome]
      return(snp_roh_freq)
    }, error = function(e) {
      stop(paste("Error while computing SNP ROH frequency:", e$message))
    })
  }
  
  # Process the SNP data when the button is clicked
  observeEvent(input$process_btn, {
    req(input$file1)
    
    # Validate file format for User 1
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
    
    # Error-handling wrapper function
    safe_process <- function(expr, error_msg) {
      tryCatch(
        expr,
        error = function(e) {
          showModal(modalDialog(
            title = "Error",
            paste(error_msg, "Details:", e$message),
            easyClose = TRUE
          ))
          return(NULL)
        }
      )
    }
    
    # Process SNP Data for User 1
    lod_data_t1 <- safe_process({
      withProgress(message = 'Processing SNP data for User 1...', {
        data1 <- process_snp_data(input$file1)
        compute_lod_scores_for_user(data1)
      })
    }, "Failed to process SNP data for User 1.")
    
    if (is.null(lod_data_t1)) return(NULL)
    
    # Process User 2 (if file exists)
    lod_data_t2 <- NULL
    if (!is.null(input$file2)) {
      req(input$file2)
      
      file_extension2 <- tools::file_ext(input$file2$name)
      if (file_extension2 != "txt") {
        showModal(modalDialog(
          title = "Invalid File Format",
          "Please upload a .txt file for User 2.",
          easyClose = TRUE
        ))
        return(NULL)
      }
      
      lod_data_t2 <- safe_process({
        withProgress(message = 'Processing SNP data for User 2...', {
          data2 <- process_snp_data(input$file2)
          compute_lod_scores_for_user(data2)
        })
      }, "Failed to process SNP data for User 2.")
      
      if (is.null(lod_data_t2)) return(NULL)
    }
    
    # Compute LOD Thresholds
    lod_threshold_t1 <- safe_process(estimate_lod_threshold(lod_data_t1),
                                     "Failed to compute LOD threshold for User 1.")
    
    lod_threshold_t2 <- if (!is.null(lod_data_t2)) {
      safe_process(estimate_lod_threshold(lod_data_t2),
                   "Failed to compute LOD threshold for User 2.")
    } else {
      NA
    }
    
    # Classify ROH
    lod_data_t1 <- safe_process(classify_roh(lod_data_t1, lod_threshold_t1),
                                "Failed to classify ROH for User 1.")
    
    if (!is.null(lod_data_t2)) {
      lod_data_t2 <- safe_process(classify_roh(lod_data_t2, lod_threshold_t2),
                                  "Failed to classify ROH for User 2.")
    }
    
    # Compute ROH frequencies
    autosomal_genome_size <- 2.8e9  # Human autosomal genome size
    chr_roh_freq_t1 <- safe_process(compute_chr_roh_freq(lod_data_t1, autosomal_genome_size),
                                    "Failed to compute chromosome ROH frequency for User 1.")
    
    chr_roh_freq_t2 <- if (!is.null(lod_data_t2)) {
      safe_process(compute_chr_roh_freq(lod_data_t2, autosomal_genome_size),
                   "Failed to compute chromosome ROH frequency for User 2.")
    } else {
      NULL
    }
    
    # Fit GMM models and classify ROH lengths
    lod_data_t1 <- safe_process(fit_gmm_and_classify(lod_data_t1),
                                "Failed to fit GMM model for User 1.")
    
    if (!is.null(lod_data_t2)) {
      lod_data_t2 <- safe_process(fit_gmm_and_classify(lod_data_t2),
                                  "Failed to fit GMM model for User 2.")
    }
    
    # Compute SNP frequencies inside ROHs
    snp_roh_freq_t1 <- safe_process(compute_snp_roh_freq(lod_data_t1),
                                    "Failed to compute SNP ROH frequency for User 1.")
    
    snp_roh_freq_t2 <- if (!is.null(lod_data_t2)) {
      safe_process(compute_snp_roh_freq(lod_data_t2),
                   "Failed to compute SNP ROH frequency for User 2.")
    } else {
      NULL
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
      filename = function() { paste("SNP_ROH_Frequency_User1.txt") },
      content = function(file) {
        write.table(snp_roh_freq_t1, file, sep = "\t", row.names = FALSE)
      }
    )
    
    output$download_chr_roh_freq_t1 <- downloadHandler(
      filename = function() { paste("Chromosome_ROH_Frequency_User1.txt") },
      content = function(file) {
        write.table(chr_roh_freq_t1, file, sep = "\t", row.names = FALSE)
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
        filename = function() { paste("SNP_ROH_Frequency_User2.txt") },
        content = function(file) {
          write.table(snp_roh_freq_t2, file, sep = "\t", row.names = FALSE)
        }
      )
      
      output$download_chr_roh_freq_t2 <- downloadHandler(
        filename = function() { paste("Chromosome_ROH_Frequency_User2.txt") },
        content = function(file) {
          write.table(chr_roh_freq_t2, file, sep = "\t", row.names = FALSE)
        }
      )
    }
    
    # Prepare output for displaying
    # Ensure lod_data_t1 is not NULL or empty
    if (!is.null(lod_data_t1) && nrow(lod_data_t1) > 0) {
      output_t1 <- list(
        lod_data_t1 = lod_data_t1
      )
      
      # Show data to the user as a data table
      output$lod_data_t1_table <- renderDT({
        datatable(output_t1$lod_data_t1)
      })
    } else {
      output$lod_data_t1_table <- renderDT({
        datatable(data.frame(Message = "No data available"))
      })
    }
    
    # Show LOD data for User 2 (if exists and not empty)
    if (!is.null(lod_data_t2) && nrow(lod_data_t2) > 0) {
      output$lod_data_t2_table <- renderDT({
        datatable(lod_data_t2)
      })
    } else {
      output$lod_data_t2_table <- renderDT({
        datatable(data.frame(Message = "No data available"))
      })
    }
    
    # Plots
    # Render plot for User 1 (lod_data_t1)
    output$plot_t1 <- renderPlotly({
      if (exists("lod_data_t1") && nrow(lod_data_t1) > 0) {
        # Create a violin plot with a boxplot overlay for User 1's data
        ggplotly(
          ggplot(lod_data_t1, aes(x = as.factor(chromosome), y = abs(roh_length) / 1e6, fill = roh_class)) +
            geom_violin(trim = TRUE, scale = "width", adjust = 2, alpha = 0.8) +  
            geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 0.5) +
            scale_fill_manual(values = c("Class A (Short ROH)" = "lightblue", 
                                         "Class B (Intermediate ROH)" = "pink", 
                                         "Class C (Long ROH)" = "lightgreen")) +  
            labs(title = "Violin Plot of ROH Lengths for User 1", 
                 x = "Chromosome", 
                 y = "Total Length of ROH (Mb)", 
                 fill = "ROH Class") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        )
      } else {
        # If no data for User 1, don't render any plot or return a message
        return(NULL)
      }
    })
    
    # Only render plot for User 2 (lod_data_t2) if data exists
    output$plot_t2 <- renderPlotly({
      if (!is.null(lod_data_t2) && nrow(lod_data_t2) > 0) {
        # Create a violin plot with a boxplot overlay for User 2's data
        ggplotly(
          ggplot(lod_data_t2, aes(x = as.factor(chromosome), y = abs(roh_length) / 1e6, fill = roh_class)) +
            geom_violin(trim = TRUE, scale = "width", adjust = 2, alpha = 0.8) +  
            geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 0.5) +
            scale_fill_manual(values = c("Class A (Short ROH)" = "lightblue", 
                                         "Class B (Intermediate ROH)" = "pink", 
                                         "Class C (Long ROH)" = "lightgreen")) +  
            labs(title = "Violin Plot of ROH Lengths for User 2", 
                 x = "Chromosome", 
                 y = "Total Length of ROH (Mb)", 
                 fill = "ROH Class") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        )
      } else {
        # If no data for User 2, don't render any plot
        return(NULL)
      }
    })
    
    # Combined plot for User 1 and User 2 (only if both datasets exist)
    output$plot_combined <- renderPlotly({
      # Check if both lod_data_t1 and lod_data_t2 exist
      if (exists("lod_data_t1") && nrow(lod_data_t1) > 0 && 
          !is.null(lod_data_t2) && nrow(lod_data_t2) > 0) {
        
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
      } else {
        # If no data for User 2, don't render combined plot
        return(NULL)
      }
    })
    
    output$data_table <- renderDT({
      # Copy data for User 1 and add a 'user' column
      user_data_t1 <- copy(lod_data_t1)
      user_data_t1$user <- "User 1"
      
      # If User 2 data exists, process and combine it with User 1 data
      if (!is.null(lod_data_t2) && nrow(lod_data_t2) > 0) {
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