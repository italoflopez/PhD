plotly_simple <- function(series,name){
plot <- plot_ly(x = index(series), y = series, type = 'scatter', mode = 'lines') %>%
  layout(title = 'Time Series Plot',
         xaxis = list(title = 'Date'),
         yaxis = list(title = name))
return (plot)}

stationarity_test <- function(series){
return(list(series%>%log()%>%diff()%>%adf.test(),series%>%log()%>%diff()%>%kpss.test(),series%>%log()%>%diff()%>%pp.test()))}