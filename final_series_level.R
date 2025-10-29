# Filter series that start on or before December 1996
data_long_all <- rbind(data_long, data_long_stock_and_watson, data_long_new_fast_variables)


final_series_data <- data_long_all %>%
  group_by(variable) %>%
  filter(min(date) <= as.Date("1996-12-31")) %>%
  ungroup()


