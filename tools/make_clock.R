# Function to draw a simple clock with customized labels and additional text
draw_clock <- function() {
  # Set the size of the output device
  png("school_clock.png", width = 128, height = 200)
  par(mar = c(0, 0, 0, 0))
  
  # Setup plotting area
  plot.new()
  plot.window(xlim = c(-1, 1), ylim = c(-1.5, 1))
  
  # Draw the clock circle
  symbols(0, 0, circles = 1, inches = FALSE, add = TRUE, fg = "black", lwd = 2)
  
  # Add custom labels
  text(cos(seq(0, 2 * pi, length.out = 13))[1:12] * 0.85, 
       sin(seq(0, 2 * pi, length.out = 13))[1:12] * 0.85, 
       labels = rep("splinetime", 12), cex = 0.5)
  
  # Draw the hands
  arrows(0, 0, 0.5, 0, lwd = 2)  # minute hand
  arrows(0, 0, 0, 0.5, lwd = 2)  # hour hand
  
  # Add text below the clock
  text(0, -1.4, "Good Heavens, just look at the time", cex = 0.7)
  
  # Close the plotting device
  dev.off()
}

# Call the function to create the image
draw_clock()
