library(tidyverse)
library(readxl)

DNA <- read_excel("C:/Users/y.liu/Downloads/DNA.xlsx")

result <- vector()

for (i in 1:nrow(DNA)){
  tbl <- matrix(as.numeric(DNA[i, c("alt_control", "ref_control", "alt_trt", "ref_trt")]), nrow = 2, byrow = T)
  temp <- fisher.test(tbl)
  result[i] <- temp$p.value
}

which(result < 0.05)
which(result < 0.05/49)

pvalues <- p.adjust(result, "BH")  # FDR correction
names(pvalues) <- paste0(DNA$Region, "_", DNA$POSITION) 

df <- data.frame(
  index = seq_along(pvalues),
  pval  = pvalues, 
  id = names(pvalues),
  Targeted = DNA$Targeted
)

df <- df %>% 
  mutate(Targeted = case_when(Targeted == "Y" ~ "Targeted",
                               Targeted == "N" ~ "Not targeted"))

plots <- ggplot(df, aes(y = pval, x = "")) +
  geom_jitter(width = 0.2, height = 0, size = 1, aes(color = Targeted)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  theme_classic() +
  labs(
    x = NULL,
    y = expression(FDR~corrected~p~value),
    color = NULL
  ) +
  geom_text(
    data = subset(df, pval < 0.05),
    aes(label = id),
    nudge_x = 0,
    nudge_y = -0.05,
    size = 3.5
  ) 

ggsave(
  filename = "pvalues_plot.tiff",
  plot = plots,
  device = "tiff",
  dpi = 600,
  width = 6,
  height = 4.5,
  units = "in",
  compression = "lzw"
)
