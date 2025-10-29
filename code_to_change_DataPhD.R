library(stringr)

transform_fred_blocks <- function(code) {
  # ----- Pattern 1: Aggregate + Index -----
  pattern_agg_index <- paste0(
    "(?s)",
    "(\\w+)\\s*<-\\s*fredr\\('([^']+)'\\).*?",
    "\\1\\s*<-\\s*\\2.*?",
    # Match aggregate with mean or sum (and handle spaces/newlines)
    "\\1\\s*<-\\s*aggregate\\(.*?(mean|sum).*?\\).*?",
    "index\\(\\1\\)\\s*<-\\s*as.yearmon\\(index\\(\\1\\)\\).*?",
    "data\\s*<-\\s*merge\\(\\s*data\\s*,\\s*\\1\\s*\\).*?",
    "colnames\\(data\\)\\[ncol\\(data\\)\\]\\s*<-\\s*\"([^\"]+)\""
  )
  
  # ----- Pattern 2: Aggregate only -----
  pattern_agg <- paste0(
    "(?s)",
    "(\\w+)\\s*<-\\s*fredr\\('([^']+)'\\).*?",
    "\\1\\s*<-\\s*\\2.*?",
    "\\1\\s*<-\\s*aggregate\\(.*?(mean|sum).*?\\).*?",
    "data\\s*<-\\s*merge\\(\\s*data\\s*,\\s*\\1\\s*\\).*?",
    "colnames\\(data\\)\\[ncol\\(data\\)\\]\\s*<-\\s*\"([^\"]+)\""
  )
  
  # ----- Pattern 3: Index only -----
  pattern_index <- paste0(
    "(?s)",
    "(\\w+)\\s*<-\\s*fredr\\('([^']+)'\\).*?",
    "\\1\\s*<-\\s*\\2.*?",
    "index\\(\\1\\)\\s*<-\\s*as.yearmon\\(index\\(\\1\\)\\).*?",
    "data\\s*<-\\s*merge\\(\\s*data\\s*,\\s*\\1\\s*\\).*?",
    "colnames\\(data\\)\\[ncol\\(data\\)\\]\\s*<-\\s*\"([^\"]+)\""
  )
  
  # ----- Pattern 4: Index + na.approx -----
  pattern_index_approx <- paste0(
    "(?s)",
    "(\\w+)\\s*<-\\s*fredr\\('([^']+)'\\).*?",
    "\\1\\s*<-\\s*\\2.*?",
    "index\\(\\1\\)\\s*<-\\s*as.yearmon\\(index\\(\\1\\)\\).*?",
    "\\1\\s*<-\\s*na\\.approx\\(\\1\\).*?",
    "data\\s*<-\\s*merge\\(\\s*data\\s*,\\s*\\1\\s*\\).*?",
    "colnames\\(data\\)\\[ncol\\(data\\)\\]\\s*<-\\s*\"([^\"]+)\""
  )
  
  # ----- Pattern 5: Simple merge -----
  pattern_simple <- paste0(
    "(?s)",
    "(\\w+)\\s*<-\\s*fredr\\('([^']+)'\\).*?",
    "\\1\\s*<-\\s*\\2.*?",
    "data\\s*<-\\s*merge\\(\\s*data\\s*,\\s*\\1\\s*\\).*?",
    "colnames\\(data\\)\\[ncol\\(data\\)\\]\\s*<-\\s*\"([^\"]+)\""
  )
  
  # ----- Replacement for aggregation / index-aggregation variants -----
  # We use a conditional logic: if it was "sum" in the match, replace with sum; otherwise use mean.
  replacement_agg <- paste0(
    "\\1 <- fredr('\\2')\n\n",
    "\\1 <- \\1 %>% select(date, value)\n",
    "# Aggregate to monthly mean or sum depending on original function\n",
    "\\1 <- \\1 %>%\n",
    "  mutate(year_month = floor_date(date, 'month')) %>%\n",
    "  group_by(year_month) %>%\n",
    "  summarise(value = \\3(value, na.rm = FALSE)) %>%\n",  # <-- this dynamically picks mean/sum
    "  mutate(date = as.Date(year_month)) %>%\n",
    "  select(date, value)\n\n",
    "data <- full_join(data, \\1, by = 'date') %>%\n",
    "  arrange(date)\n",
    "colnames(data)[ncol(data)] <- \"\\4\""
  )
  
  # ----- Replacement for index / index + na.approx variants -----
  replacement_index <- paste0(
    "\\1 <- fredr('\\2')\n\n",
    "\\1 <- \\1 %>% select(date, value)\n",
    "# Convert to monthly frequency (interpolated if necessary)\n",
    "\\1 <- \\1 %>%\n",
    "  mutate(year_month = floor_date(date, 'month')) %>%\n",
    "  group_by(year_month) %>%\n",
    "  summarise(value = mean(value, na.rm = FALSE)) %>%\n",
    "  mutate(date = as.Date(year_month)) %>%\n",
    "  select(date, value)\n\n",
    "data <- full_join(data, \\1, by = 'date') %>%\n",
    "  arrange(date)\n",
    "colnames(data)[ncol(data)] <- \"\\3\""
  )
  
  # ----- Replacement for simple merge -----
  replacement_simple <- paste0(
    "\\1 <- fredr('\\2')\n\n",
    "\\1 <- \\1 %>% select(date, value)\n\n",
    "data <- full_join(data, \\1, by = \"date\") %>%\n",
    "  arrange(date)\n\n",
    "colnames(data)[ncol(data)] <- \"\\3\""
  )
  
  # ----- Apply patterns in order of complexity -----
  code <- str_replace_all(code, regex(pattern_agg_index), replacement_agg)
  code <- str_replace_all(code, regex(pattern_agg), replacement_agg)
  code <- str_replace_all(code, regex(pattern_index_approx), replacement_index)
  code <- str_replace_all(code, regex(pattern_index), replacement_index)
  code <- str_replace_all(code, regex(pattern_simple), replacement_simple)
  
  return(code)
}




script <- readLines("DataPhD.R") |> paste(collapse = "\n")
new_script <- transform_fred_blocks(script)
writeLines(new_script, "DataPhD_modified_for_API.R")