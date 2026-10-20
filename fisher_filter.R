library(tidyverse)
df <- read_csv("ABE-ALL-nuc.csv")

# Automatically extract sample names
samples <- unique(
  sub("^[ATCG]_Count_", "", 
      grep("^[ATCG]_Count_", names(df), value = TRUE))
)

samples

for (s in samples) {
  
  base_cols <- paste0(c("A", "T", "C", "G"), "_Count_", s)
  
  counts <- as.matrix(df[, base_cols])
  colnames(counts) <- c("A", "T", "C", "G")
  
  # Exclude reference base
  ref_idx <- match(toupper(df$Ref), colnames(counts))
  valid <- !is.na(ref_idx)
  counts[cbind(which(valid), ref_idx[valid])] <- NA
  
  # Highest non-reference count
  alt_count <- apply(counts, 1, max, na.rm = TRUE)
  
  # Base corresponding to highest non-reference count
  alt <- apply(
    counts, 1,
    function(x) names(x)[which.max(x)]
  )
  
  # Strand-specific rules
  force_G <- (df$Strand == "+" & toupper(df$Ref) == "A") %in% TRUE
  force_C <- (df$Strand == "-" & toupper(df$Ref) == "T") %in% TRUE
  
  alt[force_G] <- "G"
  alt[force_C] <- "C"
  
  # If ALT is forced, ALT count should also be the count
  # corresponding to that forced nucleotide
  alt_count[force_G] <- df[[paste0("G_Count_", s)]][force_G]
  alt_count[force_C] <- df[[paste0("C_Count_", s)]][force_C]
  
  # Save results
  df[[paste0("ALT_", s)]] <- alt
  df[[paste0("ALT_Count_", s)]] <- alt_count
}

control_samples <- c("C3", "C4", "C5")
trt_samples     <- c("S1", "S2", "S3")

# Function to get REF count for one sample
get_ref_count <- function(df, sample) {
  
  sapply(seq_len(nrow(df)), function(i) {
    
    ref <- toupper(df$Ref[i])
    
    if (is.na(ref) || !ref %in% c("A", "T", "C", "G")) {
      return(NA_real_)
    }
    
    df[[paste0(ref, "_Count_", sample)]][i]
  })
}


# Create REF_Count columns
samples <- c(control_samples, trt_samples)

for (s in samples) {
  df[[paste0("REF_Count_", s)]] <- get_ref_count(df, s)
}

for (s in samples) {
  
  df[[paste0("Pct_", s)]] <-
    df[[paste0("ALT_Count_", s)]] /
    (df[[paste0("ALT_Count_", s)]] +
       df[[paste0("REF_Count_", s)]]) * 100
}

df$Control_Mean_Pct <- rowMeans(
  df[paste0("Pct_", control_samples)],
  na.rm = TRUE
)

df$Treatment_Mean_Pct <- rowMeans(
  df[paste0("Pct_", trt_samples)],
  na.rm = TRUE
)



df$ALT_Control_Pooled <- rowSums(
  df[paste0("ALT_Count_", control_samples)],
  na.rm = TRUE
)

df$REF_Control_Pooled <- rowSums(
  df[paste0("REF_Count_", control_samples)],
  na.rm = TRUE)  


results <- list()

for (s in trt_samples) {
  
  pval <- rep(NA_real_, nrow(df))
  
  for (i in seq_len(nrow(df))) {
    
    tbl <- matrix(
      c(
        df$ALT_Control_Pooled[i],
        df$REF_Control_Pooled[i],
        df[[paste0("ALT_Count_", s)]][i],
        df[[paste0("REF_Count_", s)]][i]
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    if (!anyNA(tbl) && sum(tbl) > 0) {
      pval[i] <- fisher.test(tbl, alternative = "less")$p.value
    }
  }
  
  # BH correction within each treatment sample
  pval_adj <- p.adjust(pval, method = "BH")
  
  results[[s]] <- data.frame(
    index = seq_len(nrow(df)),
    Region = df$Region,
    POSITION = df$Pos,
    Sample = s,
    Ref = df$Ref,
    ALT = df[[paste0("ALT_", s)]],
    ALT_Control = df$ALT_Control_Pooled,
    REF_Control = df$REF_Control_Pooled,
    ALT_Treatment = df[[paste0("ALT_Count_", s)]],
    REF_Treatment = df[[paste0("REF_Count_", s)]],
    pval = pval,
    pval_adj = pval_adj,
    Control_Mean_Pct = df$Control_Mean_Pct,
    Treatment_Pct = df[[paste0("Pct_", s)]]
  )
}

# Combine all treatment samples
result_tbl <- do.call(rbind, results)

rownames(result_tbl) <- NULL

write_csv(result_tbl, "results_fisher_09292026.csv")
