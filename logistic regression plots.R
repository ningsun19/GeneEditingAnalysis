
# efficiency --------------------------------------------------------------


my_function_eff <- function(x, y, number){ # x is the name of original file
  location <- unlist(str_locate_all(x[1, 1], y))
  VEC <- x[['TargetSequence']]
  cut_string <- substr(VEC, location[1]-1, location[2]) # the length is always 20
  cut_string <- toupper(cut_string)
  temp <- as.data.frame(do.call("rbind", strsplit(cut_string, "")))
  temp$Reads <- x$Reads
    location_c <- which(temp[1, ] == "C") 
  location_c <- subset(location_c, location_c > 1 & location_c < 12)
  location_x <- location_c + number 
  result_c <- temp %>% 
    select(c(location_c, 22))
  result_x_value <- temp[1,] %>% 
    select(location_x)
  result_x <- as.data.frame(matrix(rep(result_x_value, each = nrow(temp)), ncol = length(location_x)))
  names(result_x) <- paste0(names(result_c)[- ncol(result_c)], "_x")
  names(result_c) <- paste0(names(result_c), "_c")
  cbind(result_c, result_x) %>% 
    pivot_longer(
      -Reads_c,
      names_to = c("location", ".value"),
      names_sep = "_")%>% 
    mutate(y = ifelse(c == "T", 1, 0)) %>% 
    mutate(x = unlist(x)) 
}


for (i in c(-1, 1)){
  temp <- my_function_eff(atp8_combined, "CCCAACTAAATACTACCGTA", i)
  assign(paste0("atp8_m", i, sep = ""), temp)
} 

for (i in c(-1, 1)){
  temp <- my_function_eff(nd1_combined, "CTATCAACATTACTAATAAG", i)
  assign(paste0("nd1_m", i, sep = ""), temp)
} 


for (i in c( -1, 1)){
  temp <- my_function_eff(nd2_combined, "TCCATCATAGCAGGCAGTTG", i)
  assign(paste0("nd2_m", i, sep = ""), temp)
} 

for (i in c( -1, 1)){
  temp <- my_function_eff(nd3_combined, "AAATCCACCCCTTACGAGTG", i)
  assign(paste0("nd3_m", i, sep = ""), temp)
} 

for (i in c(-1, 1)){
  temp <- my_function_eff(nd4_combined, "CGCATCATAATCCTCTCTCA", i)
  assign(paste0("nd4_m", i, sep = ""), temp)
} 

for (i in c( -1, 1)){
  temp <- my_function_eff(nd5_combined, "GCAGCCGGAAGCCTATTCGC", i)
  assign(paste0("nd5_m", i, sep = ""), temp)
} 


for (i in c(-1, 1)){
  temp <- my_function_eff(cox2_combined, "ACCTACGAGTACACCGACTA", i)
  assign(paste0("cox2_m", i, sep = ""), temp)
} 

for (i in c(-1, 1)){
  temp <- my_function_eff(cox3_combined, "CAGCCCATGACCCCTAACAG", i)
  assign(paste0("cox3_m", i, sep = ""), temp)
} 


m1 <- bind_rows(`atp8_m-1`, `nd1_m-1`, `nd2_m-1`, `nd3_m-1`, `nd4_m-1`,`nd5_m-1`, `cox2_m-1`, `cox3_m-1`)

p1 <- bind_rows(atp8_m1, nd1_m1, nd2_m1, nd3_m1, nd4_m1,nd5_m1, cox2_m1, cox3_m1)

m <- cbind(m1, p1)

names(m)[c(4, 9)] <- c("m1", "p1")

m <- m %>% 
  select(1, 2, 4, 5, 9)

m <- m %>% 
  mutate(group = if_else(location %in% c("V7", "V6"), "1", if_else(location %in% c("V5", "V4", "V9", "V8"), "2", "3")))

m <- m %>% 
  group_by(y, m1, p1, group) %>% 
  summarise(n = sum(Reads_c)) 

m_wide <- m %>% 
  pivot_wider(names_from = y, values_from = n, names_prefix = "y")


model_position <- glm(cbind(y1, y0) ~ m1 + p1 + group, family = "binomial", data = m_wide)

summary(model_position)

library(lsmeans)

LRWeights_2 <- as.data.frame(matrix(NA, nrow = 4, ncol = 2))
LRWeights_2[, 1] <- summary(lsmeans(model_position, "m1", type = "response"))[2]
LRWeights_2[, 2] <- summary(lsmeans(model_position, "p1", type = "response"))[2]


logodds <- function(p){
  log(p/(1-p))
}

LRWeights_2_res <-  LRWeights_2 %>%
  mutate_all(logodds)

LRWeights_2_res <- LRWeights_2_res - median(data.matrix(LRWeights_2_res))

LRWeights_2_res$C <- rep(NA, 4)
LRWeights_2_res <- LRWeights_2_res[, c(1, 3, 2)]
  
row.names(LRWeights_2_res) <- c("A", "C", "G", "T")
LRWeights_2_matrix <- data.matrix(LRWeights_2_res) 

library(ggseqlogo)


tiff("logistic weights.tif", width = 1000, height = 900, res = 300)

ggseqlogo(LRWeights_2_matrix, method = "custom", seq_type = "dna" ) +
  scale_x_continuous(breaks = 1:3, labels = c(-1, 0, 1)) +
  geom_abline(slope = 0) + 
  annotate('text', x = 2, y = 0, size = 30, label = "C", color = "grey") +
  theme(axis.title.y = element_text(size = 12, face = "bold", color = "black", margin = margin(r = 5)),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 5))
        ) + 
  ylab ("Logistic Regression Weights") +
  xlab ("Position Relative to C") 
  
dev.off()


# purity ------------------------------------------------------------------

library(tidyverse)

my_function_purity <- function(x, y, number){ # x is the name of original file
  location <- unlist(str_locate_all(x[1, 1], y))
  VEC <- x[['TargetSequence']]
  cut_string <- substr(VEC, location[1]-1, location[2]) # the length is always 20
  cut_string <- toupper(cut_string)
  temp <- as.data.frame(do.call("rbind", strsplit(cut_string, "")))
  temp$Reads <- x$Reads
    location_c <- which(temp[1, ] == "C") 
  location_c <- subset(location_c, location_c > 1 & location_c < 12)
  location_x <- location_c + number 
  result_c <- temp %>% 
    select(c(location_c, 22))
  result_x_value <- temp[1,] %>% 
    select(location_x)
  result_x <- as.data.frame(matrix(rep(result_x_value, each = nrow(temp)), ncol = length(location_x)))
  names(result_x) <- paste0(names(result_c)[- ncol(result_c)], "_x")
  names(result_c) <- paste0(names(result_c), "_c")
  cbind(result_c, result_x) %>% 
    pivot_longer(
      -Reads_c,
      names_to = c("location", ".value"),
      names_sep = "_")%>% 
    mutate(y = ifelse(c == "T", 1, 0)) %>% 
    filter(!c == "C") %>% 
    mutate(x = unlist(x)) 
}


for (i in c(-1, 1)){
  temp <- my_function_purity(atp8_combined, "CCCAACTAAATACTACCGTA", i)
  assign(paste0("purity_atp8_m", i, sep = ""), temp)
} 


for (i in c(-1, 1)){
  temp <- my_function_purity(nd1_combined, "CTATCAACATTACTAATAAG", i)
  assign(paste0("purity_nd1_m", i, sep = ""), temp)
} 


for (i in c( -1, 1)){
  temp <- my_function_purity(nd2_combined, "TCCATCATAGCAGGCAGTTG", i)
  assign(paste0("purity_nd2_m", i, sep = ""), temp)
} 

for (i in c( -1, 1)){
  temp <- my_function_purity(nd3_combined, "AAATCCACCCCTTACGAGTG", i)
  assign(paste0("purity_nd3_m", i, sep = ""), temp)
} 

for (i in c(-1, 1)){
  temp <- my_function_purity(nd4_combined, "CGCATCATAATCCTCTCTCA", i)
  assign(paste0("purity_nd4_m", i, sep = ""), temp)
} 

for (i in c( -1, 1)){
  temp <- my_function_purity(nd5_combined, "GCAGCCGGAAGCCTATTCGC", i)
  assign(paste0("purity_nd5_m", i, sep = ""), temp)
} 


for (i in c(-1, 1)){
  temp <- my_function_purity(cox2_combined, "ACCTACGAGTACACCGACTA", i)
  assign(paste0("purity_cox2_m", i, sep = ""), temp)
} 

for (i in c(-1, 1)){
  temp <- my_function_purity(cox3_combined, "CAGCCCATGACCCCTAACAG", i)
  assign(paste0("purity_cox3_m", i, sep = ""), temp)
} 


purity_m1 <- bind_rows(`purity_atp8_m-1`, `purity_nd1_m-1`, `purity_nd2_m-1`, `purity_nd3_m-1`, `purity_nd4_m-1`,`purity_nd5_m-1`, `purity_cox2_m-1`, `purity_cox3_m-1`)

purity_p1 <- bind_rows(purity_atp8_m1, purity_nd1_m1, purity_nd2_m1, purity_nd3_m1, purity_nd4_m1,purity_nd5_m1, purity_cox2_m1, purity_cox3_m1)

purity_m <- cbind(purity_m1, purity_p1)

names(purity_m)[c(4, 9)] <- c("m1", "p1")

purity_m <- purity_m %>% 
  select(1, 2, 4, 5, 9)

purity_m <- purity_m %>% 
  mutate(group = if_else(location %in% c("V7", "V6"), "1", if_else(location %in% c("V5", "V4", "V9", "V8"), "2", "3")))

purity_m <- purity_m %>% 
  group_by(y, m1, p1, group) %>% 
  summarise(n = sum(Reads_c)) 

purity_m_wide <- purity_m %>% 
  pivot_wider(names_from = y, values_from = n, names_prefix = "y")

purity_m_wide <- purity_m_wide %>% 
  mutate(y0 = if_else(is.na(y0), 0, y0))

purity_model_position <- glm(cbind(y1, y0) ~ m1 + p1 + group, family = "binomial", data = purity_m_wide)

summary(purity_model_position)

library(lsmeans)

LRWeights_p <- as.data.frame(matrix(NA, nrow = 4, ncol = 2))
LRWeights_p[, 1] <- summary(lsmeans(purity_model_position, "m1", type = "response"))[2]
LRWeights_p[, 2] <- summary(lsmeans(purity_model_position, "p1", type = "response"))[2]

logodds <- function(p){
  log(p/(1-p))
}


LRWeights_p_res <-  LRWeights_p %>%
  mutate_all(logodds)


LRWeights_p_res <- LRWeights_p_res - median(data.matrix(LRWeights_p_res))

LRWeights_p_res$C <- rep(NA, 4)
LRWeights_p_res <- LRWeights_p_res[, c(1, 3, 2)]

row.names(LRWeights_p_res) <- c("A", "C", "G", "T")
LRWeights_p_matrix <- data.matrix(LRWeights_p_res) 

library(ggseqlogo)


tiff("logistic weights_purity.tif", width = 1000, height = 900, res = 300)

ggseqlogo(LRWeights_p_matrix, method = "custom", seq_type = "dna" ) +
  scale_x_continuous(breaks = 1:3, labels = c(-1, 0, 1)) +
  # scale_y_continuous(limits = c(-0.4, 0.4)) + 
  geom_abline(slope = 0) + 
  annotate('text', x = 2, y = 0, size = 30, label = "C", color = "grey") +
  theme(axis.title.y = element_text(size = 12, face = "bold", color = "black", margin = margin(r = 5)),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 5))
  ) + 
  ylab ("Logistic Regression Weights") +
  xlab ("Position Relative to C") 

dev.off()

