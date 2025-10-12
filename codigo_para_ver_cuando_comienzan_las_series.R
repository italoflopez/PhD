#Código para ver en qué fecha comienzza cada variable

index(data_stock_and_watson_stationary) <- as.yearmon(index(data_stock_and_watson_stationary))

# Function to find the start date of each series
start_dates <- apply(data_stock_and_watson_stationary, 2, function(col) index(data_stock_and_watson_stationary)[which.min(is.na(col))])

# Display the start dates
print(start_dates)