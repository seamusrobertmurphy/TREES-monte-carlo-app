# Diagnostic script to identify the exact problem
# Run this in R console where your app.R is located

diagnose_csv_issues <- function() {
  cat("=====================================\n")
  cat("CSV FILE DIAGNOSTIC\n")
  cat("=====================================\n\n")
  
  # 1. Check current working directory
  cat("1. Current working directory:\n")
  cat("   ", getwd(), "\n\n")
  
  # 2. List all files in current directory
  cat("2. Files in current directory:\n")
  files <- list.files()
  for(f in files) {
    cat("   -", f, "\n")
  }
  cat("\n")
  
  # 3. Check if CSV files exist
  cat("3. CSV File Existence Check:\n")
  csv_files <- c("AllBiomassData_Prepared.csv", "LIF_Prepared.csv", "LDF_Prepared.csv")
  
  for(csv in csv_files) {
    if(file.exists(csv)) {
      cat("   ✓", csv, "EXISTS\n")
      
      # Try to read it
      tryCatch({
        data <- read.csv(csv, stringsAsFactors = FALSE)
        cat("     Successfully read:", nrow(data), "rows,", ncol(data), "columns\n")
        cat("     Column names:", paste(names(data), collapse=", "), "\n")
        
        # Check for specific columns we need
        if(csv == "AllBiomassData_Prepared.csv") {
          expected <- c("AGB_tC_ha", "BGB_tC_ha", "Saplings_tC_ha", "Litter_tC_ha", 
                       "Standing_Dead_tC_ha", "Lying_Dead_tC_ha", "Soil_tC_ha")
          found <- expected[expected %in% names(data)]
          missing <- expected[!expected %in% names(data)]
          
          if(length(found) > 0) {
            cat("     Found columns:", paste(found, collapse=", "), "\n")
          }
          if(length(missing) > 0) {
            cat("     MISSING columns:", paste(missing, collapse=", "), "\n")
          }
        }
        
      }, error = function(e) {
        cat("     ERROR reading file:", e$message, "\n")
      })
      
    } else {
      cat("   ✗", csv, "NOT FOUND\n")
    }
    cat("\n")
  }
  
  # 4. Check for case sensitivity issues
  cat("4. Case sensitivity check:\n")
  all_files <- list.files(pattern = "\\.csv$", ignore.case = TRUE)
  cat("   All CSV files (case-insensitive):\n")
  for(f in all_files) {
    cat("     -", f, "\n")
  }
  cat("\n")
  
  # 5. Try alternate reading methods
  cat("5. Testing alternate reading methods:\n")
  if(file.exists("AllBiomassData_Prepared.csv")) {
    # Method 1: Default read.csv
    tryCatch({
      data1 <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
      cat("   Method 1 (read.csv): Success,", ncol(data1), "columns\n")
    }, error = function(e) {
      cat("   Method 1 (read.csv): Failed -", e$message, "\n")
    })
    
    # Method 2: read.csv with header check
    tryCatch({
      data2 <- read.csv("AllBiomassData_Prepared.csv", header = TRUE, stringsAsFactors = FALSE)
      cat("   Method 2 (explicit header): Success,", ncol(data2), "columns\n")
    }, error = function(e) {
      cat("   Method 2 (explicit header): Failed -", e$message, "\n")
    })
    
    # Method 3: read.table
    tryCatch({
      data3 <- read.table("AllBiomassData_Prepared.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
      cat("   Method 3 (read.table): Success,", ncol(data3), "columns\n")
    }, error = function(e) {
      cat("   Method 3 (read.table): Failed -", e$message, "\n")
    })
  }
  
  cat("\n=====================================\n")
  cat("DIAGNOSTIC COMPLETE\n")
  cat("=====================================\n")
}

# Run the diagnostic
diagnose_csv_issues()

# Additional check: Show first few lines of CSV if it exists
if(file.exists("AllBiomassData_Prepared.csv")) {
  cat("\n\nFirst 3 lines of AllBiomassData_Prepared.csv:\n")
  cat("----------------------------------------\n")
  lines <- readLines("AllBiomassData_Prepared.csv", n = 3)
  for(i in 1:length(lines)) {
    cat("Line", i, ":", lines[i], "\n")
  }
}
