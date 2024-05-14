#' @importFrom ggplot2 ggsave
#' @importFrom base64enc dataURI
#' 
plot2base64 <- function(plot, 
                        plot_nrows, 
                        width = 7, 
                        base_height_per_row = 2.5, 
                        units = "in", 
                        html_img_width="100%") {
  additional_height_per_row <- 2.1
  height <- base_height_per_row + (plot_nrows - 1) * 
    additional_height_per_row
  
  # Temporarily save the plot as an image
  img_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(img_file, plot=plot, width = width, height = height, units = units, 
                  limitsize = FALSE)
  
  # Convert the image to a Base64 string
  img_base64 <- base64enc::dataURI(file = img_file, mime = "image/png")
  
  # Delete the temporary image file
  unlink(img_file)
  
  # Return the HTML img tag with the Base64 string and a fixed width
  return(sprintf('<img src="%s" alt="Plot" style="width:%s;">', 
                 img_base64, html_img_width))
}


build_plot_report_html <- function(header_section, 
                                   plots, 
                                   rowcounts, 
                                   path) {
  index <- 1
  for (plot in plots) {
    plot_nrows <- rowcounts[[index]]
    index <- index + 1
    img_tag <- plot2base64(plot, plot_nrows)
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


#' Make dynamic HTML text descriptions, such as header- and design_text.
#' Future task
#' 
#' @importFrom here here
#' 
generate_report_html_new <- function(plots, 
                                 plot_len, 
                                 report_dir,
                                 filename, 
                                 timestamp) {
  header_text <- paste("limma hyperparams screen", timestamp, sep=" | ")
  
  design_text <- ""
  
  html_content <- paste(
    "<html><head><title>My Plots</title></head><body>",
    "<h1 style='color:red;'>", header_text, "</h1>",
    "<h2></h2>",
    "<p>", design_text, "</p>",
    "</body></html>",
    sep=""
  )
  
  file_name <- sprintf("%s_%s.html", filename, timestamp)
  output_file_path <- here::here(report_dir, file_name)
  
  build_plot_report_html(html_content, 
                         plots, 
                         plot_len, 
                         output_file_path)
}