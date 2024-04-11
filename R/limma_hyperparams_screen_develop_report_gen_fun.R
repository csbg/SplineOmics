# libraries --------------


load("combo_pair_plots.RData")


# Define subfuntioncs --------


convertPlotToBase64ImgTag <- function(plot, fn5_plot_nrows, width = 7, 
                                      base_height_per_row = 2.5, units = "in", 
                                      html_img_width="100%") {
  
  additional_height_per_row <- 2.1
  height <- base_height_per_row + (fn5_plot_nrows - 1) * 
    additional_height_per_row
  
  # Temporarily save the plot as an image
  img_file <- tempfile(fileext = ".png")
  ggsave(img_file, plot=plot, width = width, height = height, units = units, 
         limitsize = FALSE)
  
  # Convert the image to a Base64 string
  img_base64 <- base64enc::dataURI(file = img_file, mime = "image/png")
  
  # Delete the temporary image file
  unlink(img_file)
  
  # Return the HTML img tag with the Base64 string and a fixed width
  return(sprintf('<img src="%s" alt="Plot" style="width:%s;">', 
                 img_base64, html_img_width))
}


buildPlotReportHtml <- function(header_section, plots, rowcounts, path) {
  
  
  
  index <- 1
  for (plot in plots) {
    plot_nrows <- rowcounts[[index]]
    index <- index + 1
    img_tag <- convertPlotToBase64ImgTag(plot, plot_nrows)
    header_section <- paste(header_section, img_tag, sep="\n")
  }
  
  
  # Close the HTML document
  html_content <- paste(header_section, "</body></html>", sep="\n")
  
  
  dir_path <- dirname(path)
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  writeLines(html_content, path)
}



# Working section -----



