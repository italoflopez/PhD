zoo_bloomberg_monthly <- function(dataframe){
  dataframe$FECHA <- ifelse(as.numeric(substr(dataframe$FECHA, 8, 9)) > 24, gsub("-(\\d{2})$", "-19\\1", dataframe$FECHA), gsub("-(\\d{2})$", "-20\\1", dataframe$FECHA))
  
  
  # Create a zoo object
  dataframe <- zoo(dataframe$PX_LAST, order.by = as.yearmon(dataframe$FECHA, format = "%d-%b-%Y"))
  
  # Convert the index to Date format with the desired day (01)
  index(dataframe) <- as.Date(format(index(dataframe), "%Y-%m-01"))
  
return(dataframe)}