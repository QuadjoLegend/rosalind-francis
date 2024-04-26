# First, create a data frame with the data
data <- data.frame(
  Name = c("c.5252G>A (R1751Q)", "c.5252G>A (R1751Q)", "c. 191G>A (C64Y)", "c. 191G>A (C64Y)", "81-1G>A", "81-1G>A"),
  Replicate = c("Rep_1", "Rep_2", "Rep_1", "Rep_2", "Rep_1", "Rep_2"),
  Day_0_total_reads = c(85810, 37776, 127404, 218164, 56560, 67735),
  Day_0_sub = c(31630, 34995, 20631, 85059, 12359, 13054),
  Day_0_sub_rate = c("36.90%", "92.63%", "16.20%", "38.98%", "21.90%", "19.27%"),
  Day_7_total_reads = c(158490, 119701, 97051, 236437, 47113, 37316),
  Day_7_sub = c(11509, 12916, 65419, 27495, 19443, 17975),
  Day_7_sub_rate = c("7.26%", "10.79%", "67.40%", "11.62%", "41.27%", "48.17%"),
  Day_21_total_reads = c(123778, 128616, 155508, 158490, 71854, 45202),
  Day_21_sub = c(13211, 13198, 16830, 0, 27286, 17775),
  Day_21_sub_rate = c("10.67%", "10.26%", "10.82%", "0", "37.97%", "39.32%")
)

# Convert percentage columns to numeric
data$Day_0_sub_rate <- as.numeric(sub("%", "", data$Day_0_sub_rate)) / 100
data$Day_7_sub_rate <- as.numeric(sub("%", "", data$Day_7_sub_rate)) / 100
data$Day_21_sub_rate <- as.numeric(sub("%", "", data$Day_21_sub_rate)) / 100

# Calculate means and standard errors for each day separately
means_day_0 <- aggregate(Day_0_sub_rate ~ Name, data = data, FUN = mean)
means_day_7 <- aggregate(Day_7_sub_rate ~ Name, data = data, FUN = mean)
means_day_21 <- aggregate(Day_21_sub_rate ~ Name, data = data, FUN = mean)
se_day_0 <- aggregate(Day_0_sub_rate ~ Name, data = data, FUN = function(x) sd(x) / sqrt(length(x)))
se_day_7 <- aggregate(Day_7_sub_rate ~ Name, data = data, FUN = function(x) sd(x) / sqrt(length(x)))
se_day_21 <- aggregate(Day_21_sub_rate ~ Name, data = data, FUN = function(x) sd(x) / sqrt(length(x)))

# Merge the means and standard errors for each day
means <- merge(merge(means_day_0, means_day_7, by = "Name"), means_day_21, by = "Name")
se <- merge(merge(se_day_0, se_day_7, by = "Name"), se_day_21, by = "Name")

# Rename columns
colnames(means) <- c("Name", "Mean_Day_0", "Mean_Day_7", "Mean_Day_21")
colnames(se) <- c("Name", "SE_Day_0", "SE_Day_7", "SE_Day_21")

# Merge means and standard errors
means_se <- merge(means, se, by = "Name")

# Subset the data to include only rows for specified genes and replicates
subset_data <- subset(data, Name %in% c("c.5252G>A (R1751Q)", "c. 191G>A (C64Y)", "81-1G>A") & Replicate %in% c("Rep_1", "Rep_2"))

# Save the means and standard errors data frame as a CSV file
write.csv(means_se, "means_se_substitution_data.csv", row.names = FALSE)

library(plotly)

# Select relevant columns for plotting
melted_data <- subset_data[, c("Name", "Day_0_sub_rate", "Day_7_sub_rate", "Day_21_sub_rate")]
errors <- means_se[, c("Name", "SE_Day_0", "SE_Day_7", "SE_Day_21")]

# Melt the data for easier plotting
melted_data <- reshape2::melt(melted_data, id.vars = "Name")
errors <- reshape2::melt(errors, id.vars = "Name")

# Define a gradient color scale
colors <- c("#3498DB", "#5DADE2", "#85C1E9")

# Create the bar plot using Plotly
plot_ly(melted_data, x = ~Name, y = ~value, color = ~variable, colors = colors, type = "bar",
        error_y = list(type = "data", array = errors$value, visible = TRUE)) %>%
  layout(xaxis = list(title = "Gene"),
         yaxis = list(title = "Relative Substitution Frequency"),
         barmode = "group",
         legend = list(title = "Day"),
         showlegend = TRUE)

                               