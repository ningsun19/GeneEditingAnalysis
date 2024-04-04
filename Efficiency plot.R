library(gridExtra)
library(grid)
library(stringr)
library(readxl)
library(tidyverse)

my_function <- function(x, y){ # x = file name, y = wild sequence
  location <- unlist(str_locate_all(x[1, 1], y))
  VEC <- x[['TargetSequence']]
  cut_string <- substr(VEC, location[1], location[2]) 
  cut_string <- toupper(cut_string)
  temp <- as.data.frame(do.call("rbind", strsplit(cut_string, "")))
  temp$Reads <- x$Reads
  temp <- temp %>% 
    filter(Reads >= 100)
  Sum <- sum(temp$Reads)
  
  result <- data.frame(type = c("A", "C", "G", "T", "-", "N"))
  
  for (i in 1:20){
    temp_prop <- temp %>% 
      group_by_at(i) %>% 
      summarise(sum_read = sum(Reads)) %>% 
      mutate(prop = sum_read/Sum*100) %>% 
      select(-sum_read)
    
    names(temp_prop)[1] <- "type"
    
    result <- left_join(result, temp_prop, by = "type")
  }
  
  names(result)[2:21] <- unlist(strsplit(y, ''))
  
  result
}



# CALCULATION -------------------------------------------------------------

ATP8_1_result <- my_function(ATP8_1, "CCCAACTAAATACTACCGTA")

ATP8_2_result <- my_function(ATP8_2, "CCCAACTAAATACTACCGTA")

ATP8_3_result <- my_function(ATP8_3, "CCCAACTAAATACTACCGTA") 

ATP8_4_result <- my_function(ATP8_4, "CCCAACTAAATACTACCGTA") 

ND1_1_result <- my_function(ND1_1, "CTATCAACATTACTAATAAG")

ND1_2_result <- my_function(ND1_2, "CTATCAACATTACTAATAAG")

ND1_3_result <- my_function(ND1_3, "CTATCAACATTACTAATAAG")

ND1_4_result <- my_function(ND1_4, "CTATCAACATTACTAATAAG")

ND1_5_result <- my_function(ND1_5, "CTATCAACATTACTAATAAG")


ND2_1_result <- my_function(ND2_1, "TCCATCATAGCAGGCAGTTG")

ND2_2_result <- my_function(ND2_2, "TCCATCATAGCAGGCAGTTG")

ND2_3_result <- my_function(ND2_3, "TCCATCATAGCAGGCAGTTG")

ND2_4_result <- my_function(ND2_4, "TCCATCATAGCAGGCAGTTG")

ND2_5_result <- my_function(ND2_5, "TCCATCATAGCAGGCAGTTG")


ND3_1_result <- my_function(ND3_1, "AAATCCACCCCTTACGAGTG")

ND3_2_result <- my_function(ND3_2, "AAATCCACCCCTTACGAGTG")

ND3_3_result <- my_function(ND3_3, "AAATCCACCCCTTACGAGTG")

ND3_4_result <- my_function(ND3_4, "AAATCCACCCCTTACGAGTG")



ND4_1_result <- my_function(ND4_1, "CGCATCATAATCCTCTCTCA")

ND4_2_result <- my_function(ND4_2, "CGCATCATAATCCTCTCTCA")

ND4_3_result <- my_function(ND4_3, "CGCATCATAATCCTCTCTCA")

ND4_4_result <- my_function(ND4_4, "CGCATCATAATCCTCTCTCA")

ND4_5_result <- my_function(ND4_5, "CGCATCATAATCCTCTCTCA")

ND4_6_result <- my_function(ND4_6, "CGCATCATAATCCTCTCTCA")


ND5_1_result <- my_function(ND5_1, "GCAGCCGGAAGCCTATTCGC")

ND5_2_result <- my_function(ND5_2, "GCAGCCGGAAGCCTATTCGC")

ND5_3_result <- my_function(ND5_3, "GCAGCCGGAAGCCTATTCGC")

ND5_4_result <- my_function(ND5_4, "GCAGCCGGAAGCCTATTCGC")

COX2_1_result <- my_function(COX2_1, "ACCTACGAGTACACCGACTA")

COX2_2_result <- my_function(COX2_2, "ACCTACGAGTACACCGACTA")

COX2_3_result <- my_function(COX2_3, "ACCTACGAGTACACCGACTA")

COX2_4_result <- my_function(COX2_4, "ACCTACGAGTACACCGACTA")


COX3_1_result <- my_function(COX3_1, "CAGCCCATGACCCCTAACAG")

COX3_2_result <- my_function(COX3_2, "CAGCCCATGACCCCTAACAG")

COX3_3_result <- my_function(COX3_3, "CAGCCCATGACCCCTAACAG")

COX3_4_result <- my_function(COX3_4, "CAGCCCATGACCCCTAACAG")

COX3_5_result <- my_function(COX3_5, "CAGCCCATGACCCCTAACAG")

# EXPORT ------------------------------------------------------------------
my_list <- Filter(function(x) is(x, "data.frame"), mget(ls()))
list_names <- names(my_list)
location_list <- str_which(list_names, "result")
my_list_use <- my_list[location_list]

lapply(1:length(my_list_use), 
       function(i) 
         write.csv(my_list_use[[i]], 
                   file = paste0(names(my_list_use[i]), ".csv"),
                   row.names = FALSE
         )
)





# Combine by gene type  -------------------------------------------------



ATP8 <- rbind(ATP8_1_result, ATP8_2_result, ATP8_3_result, ATP8_4_result)
ATP8$group <- rep(c("1", "2", "3", "4"), each = 6)
names(ATP8)[2:21] <- paste(names(ATP8)[2:21], 1:20, sep = "_")

COX2 <- rbind(COX2_1_result, COX2_2_result, COX2_3_result, COX2_4_result)
COX2$group <- rep(c("1", "2", "3", "4"), each = 6)
names(COX2)[2:21] <- paste(names(COX2)[2:21], 1:20, sep = "_")

ND3 <- rbind(ND3_1_result, ND3_2_result, ND3_3_result, ND3_4_result)
ND3$group <- rep(c("1", "2", "3", "4"), each = 6)
names(ND3)[2:21] <- paste(names(ND3)[2:21], 1:20, sep = "_")


ND5 <- rbind(ND5_1_result, ND5_2_result, ND5_3_result, ND5_4_result)
ND5$group <- rep(c("1", "2", "3", "4"), each = 6)
names(ND5)[2:21] <- paste(names(ND5)[2:21], 1:20, sep = "_")

COX3 <- rbind(COX3_1_result, COX3_2_result, COX3_3_result, COX3_4_result, COX3_5_result)
COX3$group <- rep(c("1", "2", "3", "4", "5"), each = 6)
names(COX3)[2:21] <- paste(names(COX3)[2:21], 1:20, sep = "_")

ND1 <- rbind(ND1_1_result, ND1_2_result, ND1_3_result, ND1_4_result, ND1_5_result)
ND1$group <- rep(c("1", "2", "3", "4", "5"), each = 6)
names(ND1)[2:21] <- paste(names(ND1)[2:21], 1:20, sep = "_")

ND2 <- rbind(ND2_1_result, ND2_2_result, ND2_3_result, ND2_4_result, ND2_5_result)
ND2$group <- rep(c("1", "2", "3", "4", "5"), each = 6)
names(ND2)[2:21] <- paste(names(ND2)[2:21], 1:20, sep = "_")

ND4 <- rbind(ND4_1_result, ND4_2_result, ND4_3_result, ND4_4_result, ND4_5_result, ND4_6_result)
ND4$group <- rep(c("1", "2", "3", "4", "5", "6"), each = 6)
names(ND4)[2:21] <- paste(names(ND4)[2:21], 1:20, sep = "_")


# plot function -----------------------------------------------------------


efficiency <- function(x){
  
  WT_label <- substr(names(x)[2:21], 1, 1)
  title <- deparse(substitute(x))
  
  df_long <- pivot_longer(x, !c(type, group), names_to = "WT", values_to = "pct")
  
  
  df <- df_long %>%
    mutate(WT_type = substr(WT, 1, 1),
           position = substr(WT, 3, nchar(WT))) %>% 
    filter(type == "T")  %>% 
    mutate(pct = if_else(WT_type == "C" & is.na(pct), 0L, pct)) %>% 
    mutate(pct_new = if_else(WT_type == "C", pct, NA_real_)) 
  
  n <- length(df$group)/20
    
  error_bar <- df %>% 
    group_by(position) %>% 
    summarise(mean_pct = mean(pct_new),
              sd_pct = sd(pct_new)) |> 
    mutate(se_pct = sd_pct / sqrt(n) ) 

  
  ggplot(error_bar, aes(x = fct_reorder(position, as.numeric(position)), y = mean_pct)) +
    geom_bar(stat = "identity", alpha = 1, fill = "darkgrey") +
    geom_errorbar(aes (x = position, ymax = mean_pct + se_pct, ymin = mean_pct - se_pct)) +
    geom_point(data = df, aes(x = position, y = pct_new), size = 1.3) +
    scale_x_discrete(name = NULL,
                     labels = WT_label) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), limits = c(0, 50)) + 
    labs(title = title, y = NULL) +
    theme_classic()+
    theme(plot.title = element_text(size = 10))
  
}



plot_atp <- efficiency(ATP8)
plot_nd1 <- efficiency(ND1)
plot_nd2 <- efficiency(ND2)
plot_nd3 <- efficiency(ND3)
plot_nd4 <- efficiency(ND4)
plot_nd5 <- efficiency(ND5)
plot_cox2 <- efficiency(COX2)
plot_cox3 <- efficiency(COX3)



purity <- function(x){
  
  WT_label <- substr(names(x)[2:21], 1, 1)
  
  title = deparse(substitute(x))
  
  df_long <- pivot_longer(x, cols = 2:21, names_to = "WT", values_to = "pct")
  
  
  df <- df_long %>%
    mutate(WT_type = substr(WT, 1, 1),
           position = substr(WT, 3, nchar(WT))) %>% 
    filter(WT_type != type) |> 
    group_by(group, position)  %>% 
    mutate(pct_new = pct/sum(pct, na.rm = T)*100) |>
    mutate(pct_new = if_else(WT_type == "C" & type == "T" & is.na(pct_new) & sum(pct, na.rm = T) > 1, 0L, pct_new)) %>%  
    mutate(pct_only = if_else(WT_type == "C" & type == "T", pct_new, NA_real_)) 
  

  n <- length(df$group)/100

  error_bar <- df %>%
    group_by(position) %>%
    summarise(mean_pct = mean(pct_only, na.rm = T),
              sd_pct = sd(pct_only, na.rm = T)) |>
    mutate(se_pct = sd_pct / sqrt(n) )


  ggplot(error_bar, aes(x = fct_reorder(position, as.numeric(position)), y = mean_pct)) +
    geom_bar(stat = "identity", alpha = 1, fill = "darkgrey") +
    geom_errorbar(aes (x = position, ymax = mean_pct + se_pct, ymin = mean_pct - se_pct)) +
    geom_point(data = df, aes(x = position, y = pct_only), size = 1.3) +
    scale_x_discrete(name = NULL,
                   labels = WT_label) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.1))) + 
    labs(title = title, y = NULL) +
    theme_classic()+
    theme(plot.title = element_text(size = 10))
  
}

p_atp8 <- purity(ATP8)
p_nd1 <- purity(ND1)
p_nd2 <- purity(ND2)
p_nd3 <- purity(ND3)
p_nd4 <- purity(ND4)
p_nd5 <- purity(ND5)
p_cox2 <- purity(COX2)
p_cox3 <- purity(COX3)



# control group -----------------------------------------------------------
A8_C_result <- my_function(A8_C1, "CCCAACTAAATACTACCGTA")
A8_C2_result <- my_function(A8_C2, "CCCAACTAAATACTACCGTA")
A8_C3_result <- my_function(A8_C3, "CCCAACTAAATACTACCGTA")
A8_C4_result <- my_function(A8_C4, "CCCAACTAAATACTACCGTA")

ATP8_C <- rbind(A8_C_result, A8_C2_result, A8_C3_result, A8_C4_result)
ATP8_C$group <- rep(c("1", "2", "3", "4"), each = 6)
names(ATP8_C)[2:21] <- paste(names(ATP8_C)[2:21], 1:20, sep = "_")


ND1_C_result <- my_function(N1_C1, "CTATCAACATTACTAATAAG")
ND1_C2_result <- my_function(N1_C2, "CTATCAACATTACTAATAAG")
ND1_C3_result <- my_function(N1_C3, "CTATCAACATTACTAATAAG")
ND1_C4_result <- my_function(N1_C4, "CTATCAACATTACTAATAAG")
ND1_C5_result <- my_function(N1_C5, "CTATCAACATTACTAATAAG")

ND1_C <- rbind(ND1_C_result, ND1_C2_result, ND1_C3_result, ND1_C4_result, ND1_C5_result)
ND1_C$group <- rep(c("1", "2", "3", "4", "5"), each = 6)
names(ND1_C)[2:21] <- paste(names(ND1_C)[2:21], 1:20, sep = "_")


C2_C1_result <- my_function(C2_C1, "ACCTACGAGTACACCGACTA")
C2_C2_result <- my_function(C2_C2, "ACCTACGAGTACACCGACTA")
C2_C3_result <- my_function(C2_C3, "ACCTACGAGTACACCGACTA")


COX2_C <- rbind(C2_C1_result, C2_C2_result, C2_C3_result)
COX2_C$group <- rep(c("1", "2", "3"), each = 6)
names(COX2_C)[2:21] <- paste(names(COX2_C)[2:21], 1:20, sep = "_")


C3_C1_result <- my_function(C3_C1, "CAGCCCATGACCCCTAACAG")
C3_C2_result <- my_function(C3_C2, "CAGCCCATGACCCCTAACAG")
C3_C3_result <- my_function(C3_C3, "CAGCCCATGACCCCTAACAG")
C3_C4_result <- my_function(C3_C4, "CAGCCCATGACCCCTAACAG")

COX3_C <- rbind(C3_C1_result, C3_C2_result, C3_C3_result, C3_C4_result)
COX3_C$group <- rep(c("1", "2", "3", "4"), each = 6)
names(COX3_C)[2:21] <- paste(names(COX3_C)[2:21], 1:20, sep = "_")


N2_C1_result <- my_function(N2_C1, "TCCATCATAGCAGGCAGTTG")
N2_C2_result <- my_function(N2_C2, "TCCATCATAGCAGGCAGTTG")
N2_C3_result <- my_function(N2_C3, "TCCATCATAGCAGGCAGTTG")
N2_C4_result <- my_function(N2_C4, "TCCATCATAGCAGGCAGTTG")
N2_C5_result <- my_function(N2_C5, "TCCATCATAGCAGGCAGTTG")


ND2_C <- rbind(N2_C1_result, N2_C2_result, N2_C3_result, N2_C4_result, N2_C5_result)
ND2_C$group <- rep(c("1", "2", "3", "4", "5"), each = 6)
names(ND2_C)[2:21] <- paste(names(ND2_C)[2:21], 1:20, sep = "_")

N3_C1_result <- my_function(N3_C1, "AAATCCACCCCTTACGAGTG")
N3_C2_result <- my_function(N3_C2, "AAATCCACCCCTTACGAGTG")
N3_C3_result <- my_function(N3_C3, "AAATCCACCCCTTACGAGTG")
N3_C4_result <- my_function(N3_C4, "AAATCCACCCCTTACGAGTG")

ND3_C <- rbind(N3_C1_result, N3_C2_result, N3_C3_result, N3_C4_result)
ND3_C$group <- rep(c("1", "2", "3", "4"), each = 6)
names(ND3_C)[2:21] <- paste(names(ND3_C)[2:21], 1:20, sep = "_")


N4_C1_result <- my_function(N4_C1, "CGCATCATAATCCTCTCTCA")
N4_C2_result <- my_function(N4_C2, "CGCATCATAATCCTCTCTCA")
ND4_C <- rbind(N4_C1_result, N4_C2_result)
ND4_C$group <- rep(c("1", "2"), each = 6)
names(ND4_C)[2:21] <- paste(names(ND4_C)[2:21], 1:20, sep = "_")


N5_C1_result <- my_function(N5_C1, "GCAGCCGGAAGCCTATTCGC")
N5_C2_result <- my_function(N5_C2, "GCAGCCGGAAGCCTATTCGC")
N5_C3_result <- my_function(N5_C3, "GCAGCCGGAAGCCTATTCGC")
N5_C4_result <- my_function(N5_C4, "GCAGCCGGAAGCCTATTCGC")

ND5_C <- rbind(N5_C1_result, N5_C2_result, N5_C3_result, N5_C4_result)
ND5_C$group <- rep(c("1", "2", "3", "4"), each = 6)
names(ND5_C)[2:21] <- paste(names(ND5_C)[2:21], 1:20, sep = "_")



plot_atp_C <- efficiency(ATP8_C)
plot_nd1_C <- efficiency(ND1_C)
plot_nd2_C <- efficiency(ND2_C)
plot_nd3_C <- efficiency(ND3_C)
plot_nd4_C <- efficiency(ND4_C)
plot_nd5_C <- efficiency(ND5_C)
plot_cox2_C <- efficiency(COX2_C)
plot_cox3_C <- efficiency(COX3_C)



# output plot -------------------------------------------------------------



tiff("purity.tif", width = 2400, height = 1200, res = 300)
grid.arrange(p_atp8, p_cox2, p_cox3, p_nd1, p_nd2, p_nd3, p_nd4, p_nd5, 
             ncol = 4, nrow = 2,
             left = textGrob("Purity (%)", rot = 90))
dev.off()


tiff("efficacy.tif", width = 2400, height = 1200, res = 300)
grid.arrange(plot_atp, plot_cox2, plot_cox3, plot_nd1, plot_nd2, plot_nd3, plot_nd4, plot_nd5, 
             ncol = 4, nrow = 2,
             left = textGrob("Efficiency (%)", rot = 90))

dev.off()


tiff("efficacy_c.tif", width = 2400, height = 1200, res = 300)
grid.arrange(plot_atp_C, plot_cox2_C, plot_cox3_C, plot_nd1_C, plot_nd2_C, plot_nd3_C, plot_nd4_C, plot_nd5_C, 
             ncol = 4, nrow = 2,
             left = textGrob("Efficiency (%)", rot = 90))
dev.off()
