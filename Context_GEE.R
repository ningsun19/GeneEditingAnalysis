
library(stringr)
library(dplyr)
library(tidyr)

my_function_eff2 <- function(x, y) {
  # x: original data frame
  # y: pattern to locate in first cell
  
  # locate pattern in first cell
  location <- unlist(str_locate_all(x[1, 10], y))
  
  # get target sequence column
  VEC <- x[["TargetSequence"]]
  
  # cut out the 20-ish bp region
  cut_string <- substr(VEC, location[1] - 1, location[2])
  cut_string <- toupper(cut_string)
  
  # split to columns (one column per base)
  temp <- as.data.frame(do.call("rbind", strsplit(cut_string, "")),
                        stringsAsFactors = FALSE)
  
  # add Reads column at the end
  temp$Reads <- x$Reads
  
  ref <- unlist(strsplit(toupper(substr(x[1, 10], location[1]-1, location[2])), ""))
  
  # find C positions in the first row
  location_c <- which(ref == "C") 
  # keep only within (1, 12) 
  location_c <- location_c[location_c > 1 & location_c < 12]
  
  # 1) get the C bases + Reads
  result_c <- temp %>%
    select(all_of(location_c), ncol(temp))
  names(result_c) <- c(paste0(location_c, "_c"), "Reads_c")
  
  n_rows <- nrow(temp)
  
  # 2) m1: take from row 1 only, then repeat
  m1_df <- as.data.frame(
    lapply(location_c, function(loc_c) {
      rep(ref[loc_c - 1], n_rows)  # repeat reference base
    })
  )
  names(m1_df) <- paste0(location_c, "_m1")
  
  # 3) p1: take from row 1 only, then repeat
  p1_df <- as.data.frame(
    lapply(location_c, function(loc_c) {
      rep(ref[loc_c + 1], n_rows)  # repeat reference base
    })
  )
  names(p1_df) <- paste0(location_c, "_p1")
  
  # 4) bind them and pivot longer like you did
  final <- bind_cols(result_c, m1_df, p1_df) %>%
    pivot_longer(
      -Reads_c,
      names_to = c("location", ".value"),
      names_sep = "_"
    ) %>%
    mutate(
      # keep your y indicator
      y = ifelse(c == "T", 1, 0)
    )
  
  final
}

# A8
data_names <- c("A8_1", "A8_2", "A8_3", "A8_4")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object

  temp <- my_function_eff2(data_obj, "CCCAACTAAATACTACCGTA")
    
  # Create new variable name based on dataset name and suffix
  new_name <- paste0(nm, "_result")
    
  assign(new_name, temp)
  
}

data_names <- c("N1_1", "N1_2", "N1_3", "N1_4", "N1_5", "ND1_3_1", "ND1_3_2", "ND1_3_3", "ND1_3_4",
                "ND1_3_5", "ND1_3_6", "PN1_1", "PN1_2", "PN1_3", "PN1_4", "PN1_5", "PN1_6", "PN1_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "CTATCAACATTACTAATAAG")
    
  new_name <- paste0(nm, "_result")
    
  assign(new_name, temp)
}


data_names <- c("N2_1", "N2_2", "N2_3", "N2_4", "N2_5", "ND2_3_1", "ND2_3_2", "ND2_3_3", "ND2_3_4",
                "ND2_3_5", "ND2_3_6", "PN2_1", "PN2_2", "PN2_3", "PN2_4", "PN2_5", "PN2_6", "PN2_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "TCCATCATAGCAGGCAGTTG")
  
  new_name <- paste0(nm, "_result")
  
  assign(new_name, temp)
}


data_names <- c("N3_1", "N3_2", "N3_3", "N3_4", "ND3_3_1", "ND3_3_2", "ND3_3_3", 
                "ND3_3_5", "ND3_3_6", "ND3_3_7", "PN3_1", "PN3_2", "PN3_3", "PN3_4", "PN3_5", "PN3_6", "PN3_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "AAATCCACCCCTTACGAGTG")
  
  new_name <- paste0(nm, "_result")
  
  assign(new_name, temp)
}


data_names <- c("N4_1", "N4_2", "N4_3", "N4_4", "N4_5", "N4_6", "ND4_3_1", "ND4_3_2", "ND4_3_3", "ND4_3_4",
                "ND4_3_5", "ND4_3_6",  "PN4_1", "PN4_2", "PN4_3", "PN4_4", "PN4_5", "PN4_6", "PN4_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "CGCATCATAATCCTCTCTCA")
  
  new_name <- paste0(nm, "_result")
  
  assign(new_name, temp)
}

data_names <- c("N5_1", "N5_2", "N5_3", "N5_4")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "GCAGCCGGAAGCCTATTCGC")
  
  new_name <- paste0(nm, "_result")
  
  assign(new_name, temp)
}


data_names <- c("C2_1", "C2_2", "C2_3", "C2_4")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "ACCTACGAGTACACCGACTA")
  
  new_name <- paste0(nm, "_result")
  
  assign(new_name, temp)
}

data_names <- c("C3_1", "C3_2", "C3_3", "C3_4", "C3_5")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "CAGCCCATGACCCCTAACAG")
  
  new_name <- paste0(nm, "_result")
  
  assign(new_name, temp)
}


data_names <- c("OX1_3_1", "OX1_3_2", "OX1_3_3", "OX1_3_4", "OX1_3_5", "OX1_3_6")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_eff2(data_obj, "ATAATCATCGCTATCCCCAC")
  
  new_name <- paste0(nm, "_result")
  
  assign(new_name, temp)
}

results_files <- ls(pattern = "_result$")

results_list <- mget(results_files)

df_combined <- dplyr::bind_rows(results_list, .id = "source")

m <- df_combined %>% 
  mutate(cell = case_when(
    grepl("A8", source) ~ "A8",
    grepl("OX1", source) ~ "C1",
    grepl("C2", source) ~ "C2",
    grepl("C3", source) ~ "C3",
    grepl("N1|ND1", source) ~ "N1",
    grepl("N2|ND2", source) ~ "N2",
    grepl("N3|ND3", source) ~ "N3",
    grepl("N4|ND4", source) ~ "N4",
    grepl("N5", source) ~ "N5"
  ))

m <- m %>% 
  mutate(group = if_else(location %in% c("7", "6"), "1", if_else(location %in% c("5", "4", "9", "8"), "2", "3")))

m2 <- m %>% 
  select( -c) %>% 
  group_by(y, m1, p1, group, cell) %>% 
  summarise(n = sum(Reads_c)) 

m_wide <- m2 %>% 
  pivot_wider(names_from = y, values_from = n, names_prefix = "y")

m_wide$m1   <- factor(m_wide$m1)
m_wide$p1   <- factor(m_wide$p1)
m_wide$group <- factor(m_wide$group)
m_wide$cell  <- factor(m_wide$cell)

model <- geepack::geeglm(cbind(y1, y0) ~ m1 + p1 + group, id = as.factor(cell) , corstr ="exchangeable" , family = "binomial", data = m_wide)
summary(model)

m_wide$fitted <- predict(model, type = "response")

plot(m_wide$y1/(m_wide$y0 + m_wide$y1), m_wide$fitted)
abline(0, 1)

ggplot(m_wide, aes(x = y1/(y0+y1), y = fitted, color = cell)) +
  geom_point()

library(lsmeans)

LRWeights_2 <- as.data.frame(matrix(NA, nrow = 4, ncol = 2))
LRWeights_2[, 1] <- summary(lsmeans(model, "m1", type = "response"))[2]
LRWeights_2[, 2] <- summary(lsmeans(model, "p1", type = "response"))[2]


logodds <- function(p){
  log(p/(1-p))
}

LRWeights_2_res <-  LRWeights_2 %>%
  mutate_all(logodds)

LRWeights_2_res <- LRWeights_2_res - mean(data.matrix(LRWeights_2_res))

LRWeights_2_res$C <- rep(NA, 4)
LRWeights_2_res <- LRWeights_2_res[, c(1, 3, 2)]
  
row.names(LRWeights_2_res) <- c("A", "C", "G", "T")
LRWeights_2_matrix <- data.matrix(LRWeights_2_res) 

library(ggseqlogo)


tiff("logistic weights C.tif", width = 1000, height = 900, res = 300)

ggseqlogo(LRWeights_2_matrix, method = "custom", seq_type = "dna" ) +
  scale_x_continuous(breaks = 1:3, labels = c(-1, 0, 1)) +
  geom_abline(slope = 0) + 
  annotate('text', x = 2, y = 0, size = 30, label = "C", color = "grey") +
  theme(axis.title.y = element_text(size = 12, face = "bold", color = "black", margin = margin(r = 5)),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 5))
        ) + 
  ylab ("Relative Context Effect") +
  xlab ("Position Relative to C") 
  
dev.off()


# purity ------------------------------------------------------------------


my_function_purity <- function(x, y) {
  # x: original data frame
  # y: pattern to locate in first cell
  
  # locate pattern in first cell
  location <- unlist(str_locate_all(x[1, 10], y))
  
  # get target sequence column
  VEC <- x[["TargetSequence"]]
  
  # cut out the 20-ish bp region
  cut_string <- substr(VEC, location[1] - 1, location[2])
  cut_string <- toupper(cut_string)
  
  # split to columns (one column per base)
  temp <- as.data.frame(do.call("rbind", strsplit(cut_string, "")),
                        stringsAsFactors = FALSE)
  
  # add Reads column at the end
  temp$Reads <- x$Reads
  
  ref <- unlist(strsplit(toupper(substr(x[1, 10], location[1]-1, location[2])), ""))
  
  # find C positions in the first row
  location_c <- which(ref == "C") 
  # keep only within (1, 12) 
  location_c <- location_c[location_c > 1 & location_c < 12]
  
  # 1) get the C bases + Reads
  result_c <- temp %>%
    select(all_of(location_c), ncol(temp))
  names(result_c) <- c(paste0(location_c, "_c"), "Reads_c")
  
  n_rows <- nrow(temp)
  
  # 2) m1: take from row 1 only, then repeat
  m1_df <- as.data.frame(
    lapply(location_c, function(loc_c) {
      rep(ref[loc_c - 1], n_rows)  # repeat reference base
    })
  )
  names(m1_df) <- paste0(location_c, "_m1")
  
  # 3) p1: take from row 1 only, then repeat
  p1_df <- as.data.frame(
    lapply(location_c, function(loc_c) {
      rep(ref[loc_c + 1], n_rows)  # repeat reference base
    })
  )
  names(p1_df) <- paste0(location_c, "_p1")
  
  # 4) bind them and pivot longer like you did
  final <- bind_cols(result_c, m1_df, p1_df) %>%
    pivot_longer(
      -Reads_c,
      names_to = c("location", ".value"),
      names_sep = "_"
    ) %>%
    mutate(y = ifelse(c == "T", 1, 0)) %>% 
    filter(!c == "C") 
  
  final
}

# A8
data_names <- c("A8_1", "A8_2", "A8_3", "A8_4")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "CCCAACTAAATACTACCGTA")
  
  # Create new variable name based on dataset name and suffix
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
  
}

data_names <- c("N1_1", "N1_2", "N1_3", "N1_4", "N1_5", "ND1_3_1", "ND1_3_2", "ND1_3_3", "ND1_3_4",
                "ND1_3_5", "ND1_3_6", "PN1_1", "PN1_2", "PN1_3", "PN1_4", "PN1_5", "PN1_6", "PN1_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "CTATCAACATTACTAATAAG")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}


data_names <- c("N2_1", "N2_2", "N2_3", "N2_4", "N2_5", "ND2_3_1", "ND2_3_2", "ND2_3_3", "ND2_3_4",
                "ND2_3_5", "ND2_3_6", "PN2_1", "PN2_2", "PN2_3", "PN2_4", "PN2_5", "PN2_6", "PN2_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "TCCATCATAGCAGGCAGTTG")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}


data_names <- c("N3_1", "N3_2", "N3_3", "N3_4", "ND3_3_1", "ND3_3_2", "ND3_3_3", 
                "ND3_3_5", "ND3_3_6", "ND3_3_7", "PN3_1", "PN3_2", "PN3_3", "PN3_4", "PN3_5", "PN3_6", "PN3_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "AAATCCACCCCTTACGAGTG")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}


data_names <- c("N4_1", "N4_2", "N4_3", "N4_4", "N4_5", "N4_6", "ND4_3_1", "ND4_3_2", "ND4_3_3", "ND4_3_4",
                "ND4_3_5", "ND4_3_6",  "PN4_1", "PN4_2", "PN4_3", "PN4_4", "PN4_5", "PN4_6", "PN4_7")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "CGCATCATAATCCTCTCTCA")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}

data_names <- c("N5_1", "N5_2", "N5_3", "N5_4")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "GCAGCCGGAAGCCTATTCGC")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}


data_names <- c("C2_1", "C2_2", "C2_3", "C2_4")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "ACCTACGAGTACACCGACTA")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}

data_names <- c("C3_1", "C3_2", "C3_3", "C3_4", "C3_5")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "CAGCCCATGACCCCTAACAG")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}


data_names <- c("OX1_3_1", "OX1_3_2", "OX1_3_3", "OX1_3_4", "OX1_3_5", "OX1_3_6")  

for (nm in data_names) {
  data_obj <- get(nm)  # get the actual data object
  
  temp <- my_function_purity(data_obj, "ATAATCATCGCTATCCCCAC")
  
  new_name <- paste0(nm, "_result_p")
  
  assign(new_name, temp)
}

results_files2 <- ls(pattern = "_result_p$")

results_list2 <- mget(results_files2)

df_combined2 <- dplyr::bind_rows(results_list2, .id = "source")

m2 <- df_combined2 %>% 
  mutate(cell = case_when(
    grepl("A8", source) ~ "A8",
    grepl("OX1", source) ~ "C1",
    grepl("C2", source) ~ "C2",
    grepl("C3", source) ~ "C3",
    grepl("N1|ND1", source) ~ "N1",
    grepl("N2|ND2", source) ~ "N2",
    grepl("N3|ND3", source) ~ "N3",
    grepl("N4|ND4", source) ~ "N4",
    grepl("N5", source) ~ "N5"
  ))

m2 <- m2 %>% 
  mutate(group = if_else(location %in% c("7", "6"), "1", if_else(location %in% c("5", "4", "9", "8"), "2", "3")))

m2 <- m2 %>% 
  select(-c) %>% 
  group_by(y, m1, p1, group, cell, source) %>% 
  summarise(n = sum(Reads_c)) 

m_wide2 <- m2 %>% 
  pivot_wider(names_from = y, values_from = n, names_prefix = "y")

model2 <- geepack::geeglm(cbind(y1, y0) ~ m1 + p1 + group, id = as.factor(source) , corstr ="exchangeable" , family = "binomial", data = m_wide2)
summary(model2)

LRWeights_p <- as.data.frame(matrix(NA, nrow = 4, ncol = 2))
LRWeights_p[, 1] <- summary(lsmeans(model2, "m1", type = "response"))[2]
LRWeights_p[, 2] <- summary(lsmeans(model2, "p1", type = "response"))[2]


logodds <- function(p){
  log(p/(1-p))
}

LRWeights_p_res <-  LRWeights_p %>%
  mutate_all(logodds)

LRWeights_p_res <- LRWeights_p_res - mean(data.matrix(LRWeights_p_res))

LRWeights_p_res$C <- rep(NA, 4)
LRWeights_p_res <- LRWeights_p_res[, c(1, 3, 2)]

row.names(LRWeights_p_res) <- c("A", "C", "G", "T")
LRWeights_p_matrix <- data.matrix(LRWeights_p_res) 

tiff("purity_logistic weights.tif", width = 1000, height = 900, res = 300)

ggseqlogo(LRWeights_p_matrix, method = "custom", seq_type = "dna" ) +
  scale_x_continuous(breaks = 1:3, labels = c(-1, 0, 1)) +
  geom_abline(slope = 0) + 
  annotate('text', x = 2, y = 0, size = 30, label = "C", color = "grey") +
  theme(axis.title.y = element_text(size = 12, face = "bold", color = "black", margin = margin(r = 5)),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 5))
  ) + 
  ylab ("Relative Context Effect") +
  xlab ("Position Relative to C") 

dev.off()


