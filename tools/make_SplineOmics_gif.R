library(magick)

# Function to draw a simple, more stylized clock
draw_clock <- function(time) {
  par(mar = c(2, 2, 2, 2))  # Set margins to be minimal
  plot.new()
  plot.window(xlim = c(-1, 1), ylim = c(-1.3, 1))
  
  # Draw an open circle for the clock
  theta <- seq(0, 2 * pi, length.out = 100)
  circle_x <- cos(theta)
  circle_y <- sin(theta)
  lines(circle_x, circle_y, lwd = 2)

  # Add custom labels for 'spline time'
  text(cos(seq(0, 2 * pi, length.out = 13))[1:12] * 0.85, 
       sin(seq(0, 2 * pi, length.out = 13))[1:12] * 0.85, 
       labels = rep("SplineOmics", 16), cex = 0.7, col = "navy")
  
  # Draw static arrows
  # arrows(0, 0, 0.75, 0, angle = 20, length = 0.15, lwd = 6, col = "gray")
  # arrows(0, 0, 0.5, 0, angle = 20, length = 0.15, lwd = 8, col = "darkgray")
  
  static_arrow1_angle <- pi / 4  # 45 degrees
  static_arrow2_angle <- 3 * pi / 4  # 135 degrees
  
  arrows(0, 0, 0.75 * cos(static_arrow1_angle), 0.75 * sin(static_arrow1_angle), 
         angle = 20, length = 0.15, lwd = 6, col = "black")
  arrows(0, 0, 0.5 * cos(static_arrow2_angle), 0.5 * sin(static_arrow2_angle), 
         angle = 20, length = 0.15, lwd = 8, col = "black")
  
  # Draw the hand as a thicker arrow
  arrow_x <- cos(time * 2 * pi / 60)
  arrow_y <- sin(time * 2 * pi / 60)

  arrows(0, 0, 0.9 * arrow_x, 0.9 * arrow_y, angle = 20, length = 0.15, lwd = 4)
  
  points(0, 0, pch = 16, cex = 1.5, col = "black")
  
  # Add static text below the clock
  text(0, -1.2, paste0("Good Heavens, just look at the time!\nIt's time to ",
                       "use SplineOmics for time-series analysis"), 
       cex = 1.2, 
       col = "black", xpd = TRUE)
}

# Create frames with more granularity for smooth animation
frames <- lapply(59:0, function(t) {
  img <- image_graph(width = 500, height = 500, res = 96)
  draw_clock(t)
  dev.off()  # Close the graphics device to finalize the image
  img
})

# Combine frames into a single GIF
gif <- image_join(frames)
gif <- image_animate(gif, fps = 1)  # Increase FPS for smoother animation
image_write(gif, "SplineOmics_clock.gif")
