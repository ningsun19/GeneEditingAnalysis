
# Automatically extract sample names
samples <- unique(
  sub("^[ATCG]_Count_", "", 
      grep("^[ATCG]_Count_", names(df), value = TRUE))
)

samples

control_samples <- c("C3", "C4", "C5")
trt_samples     <- c("S1", "S2", "S3")

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
      pval[i] <- fisher.test(tbl)$p.value
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

write_csv(result_tbl, "results_fisher.csv")
