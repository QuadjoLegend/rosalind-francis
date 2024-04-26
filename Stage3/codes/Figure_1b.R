# First, create a data frame with the data
data <- data.frame(
  Name = c("CCR5", "CCR5", "BRCA1 #1", "BRCA1 #1", "BRCA1 #2", "BRCA1 #2"),
  Replicate = c("Rep_1", "Rep_2", "Rep_1", "Rep_2", "Rep_1", "Rep_2"),
  Day_0_total_reads = c(412679, 173697, 128994, 96675, 74681, 43433),
  Day_0_indels = c(7401, 7246, 7309, 11815, 4681, 5479),
  Day_0_indels_rate = c("33.40%", "32.60%", "19.80%", "25.50%", "13.40%", "14.30%"),
  Day_7_total_reads = c(138860, 127238, 201727, 97051, 218164, 97023),
  Day_7_indels = c(6045, 5746, 8722, 6992, 4893, 3805),
  Day_7_indels_rate = c("34.40%", "34.90%", "10.70%", "10.60%", "5.70%", "5.70%"),
  Day_14_total_reads = c(152553, 127525, 93831, 77133, 88892, 79076),
  Day_14_indels = c(8561, 7953, 1172, 1156, 503, 366),
  Day_14_indels_rate = c("32.80%", "33.50%", "5.70%", "6.40%", "2.50%", "2.00%")
)

# Convert percentage columns to numeric
data$Day_0_indels_rate <- as.numeric(sub("%", "", data$Day_0_indels_rate)) / 100
data$Day_7_indels_rate <- as.numeric(sub("%", "", data$Day_7_indels_rate)) / 100
data$Day_14_indels_rate <- as.numeric(sub("%", "", data$Day_14_indels_rate)) / 100

# Calculate means and standard errors for each day separately
means_day_0 <- aggregate(Day_0_indels_rate ~ Name, data = data, FUN = mean)
means_day_7 <- aggregate(Day_7_indels_rate ~ Name, data = data, FUN = mean)
means_day_14 <- aggregate(Day_14_indels_rate ~ Name, data = data, FUN = mean)
se_day_0 <- aggregate(Day_0_indels_rate ~ Name, data = data, FUN = function(x) sd(x) / sqrt(length(x)))
se_day_7 <- aggregate(Day_7_indels_rate ~ Name, data = data, FUN = function(x) sd(x) / sqrt(length(x)))
se_day_14 <- aggregate(Day_14_indels_rate ~ Name, data = data, FUN = function(x) sd(x) / sqrt(length(x)))

# Merge the means and standard errors for each day
means <- merge(merge(means_day_0, means_day_7, by = "Name"), means_day_14, by = "Name")
se <- merge(merge(se_day_0, se_day_7, by = "Name"), se_day_14, by = "Name")

# Rename columns
colnames(means) <- c("Name", "Mean_Day_0", "Mean_Day_7", "Mean_Day_14")
colnames(se) <- c("Name", "SE_Day_0", "SE_Day_7", "SE_Day_14")

# Merge means and standard errors
means_se <- merge(means, se, by = "Name")

# Subset the data to include only rows for replicate 1 and 2 of the specified genes
subset_data <- subset(data, Name %in% c("CCR5", "BRCA1 #1", "BRCA1 #2") & Replicate %in% c("Rep_1", "Rep_2"))

# Save the means and standard errors data frame as a CSV file
write.csv(means_se, "means_se_data.csv", row.names = FALSE)


# Load the necessary library for plotting
library(plotly)

# Subset the data to include only rows for replicate 1 and 2 of the specified genes
subset_data <- subset(data, Name %in% c("CCR5", "BRCA1 #1", "BRCA1 #2") & Replicate %in% c("Rep_1", "Rep_2"))

# Select relevant columns for plotting
melted_data <- subset_data[, c("Name", "Day_0_indels_rate", "Day_7_indels_rate", "Day_14_indels_rate")]
errors <- means_se[, c("Name", "SE_Day_0", "SE_Day_7", "SE_Day_14")]

# Melt the data for easier plotting
melted_data <- reshape2::melt(melted_data, id.vars = "Name")
errors <- reshape2::melt(errors, id.vars = "Name")

# Define a gradient color scale
colors <- c("#1b9e77", "#d95f02", "#7570b3")

# Plot using Plotly
plot_ly(melted_data, x = ~Name, y = ~value, color = ~variable, colors = colors, type = "bar",
        error_y = list(type = "data", array = errors$value, visible = TRUE)) %>%
  layout(xaxis = list(title = "Gene"),
         yaxis = list(title = "Relative Indel Frequency"),
         barmode = "group",
         legend = list(title = "Day"),
         showlegend = TRUE)
