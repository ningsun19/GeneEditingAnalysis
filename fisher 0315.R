
control <- CT |> 
  select(c(2:4, 9, 10, 15, 16))

control <- control |> 
  rowwise() |> 
  mutate(ref_c = sum(c_across(c(2, 4, 6))),
         alt_c = sum(c_across(c(3, 5, 7)))
  ) |> 
  select(ref_c, alt_c)

df_new <- cbind(Editing, control)

df_new <- df_new |> 
  select(c(2:4, 9, 10, 15, 16, 24, 25))

names(df_new)[1] <- "position"

prop_function <- function(df){
  
  fisher.test(matrix(unlist(df[c(2, 3, 8, 9)]), 2, 2), alternative = "less")$p.value
  
}

result_prop_1 <- apply(df_new, 1, prop_function)


prop_function2 <- function(df){
  
  fisher.test(matrix(unlist(df[c(4, 5, 8, 9)]), 2, 2), alternative = "less")$p.value
  
}

result_prop_2 <- apply(df_new, 1, prop_function2)


prop_function3 <- function(df){
  
  fisher.test(matrix(unlist(df[6:9]), 2, 2), alternative = "less")$p.value
  
}

result_prop_3 <- apply(df_new, 1, prop_function3)


result_reduced_prop <- data.frame(cbind(df_new[,1], result_prop_1, result_prop_2, result_prop_3))

result_reduced_prop <- result_reduced_prop |> 
  filter_at(vars(starts_with("result_")), any_vars(. <0.025))

names(result_reduced_prop)[1] <- "position"


result_reduced_prop <- result_reduced_prop |> 
  left_join(df_new, by = "position")


result_reduced_prop <- result_reduced_prop |> 
  mutate(pct1 = ALT_COUNT...4/(ALT_COUNT...4 + REF_COUNT...3),
         pct2 = ALT_COUNT...10/(ALT_COUNT...10 + REF_COUNT...9),
         pct3 = ALT_COUNT...16/(ALT_COUNT...16 + REF_COUNT...15),
         pct_avg = (pct1 + pct2 + pct3)/3)

CT <- CT |> 
  mutate(pct1_c = ALT_COUNT...4/(ALT_COUNT...4 + REF_COUNT...3),
         pct2_c = ALT_COUNT...10/(ALT_COUNT...10 + REF_COUNT...9),
         pct3_c = ALT_COUNT...16/(ALT_COUNT...16 + REF_COUNT...15),
         pct_avg_c = (pct1_c + pct2_c + pct3_c)/3) |> 
  rename(position = POSITION...2) |> 
  select(position, pct1_c, pct2_c, pct3_c, pct_avg_c)

final <- result_reduced_prop |> 
  left_join(CT, by = "position")

write.csv(final, "fisher final2.csv", row.names = F)

