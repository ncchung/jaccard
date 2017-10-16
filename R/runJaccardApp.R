#' Launch an interactive Shiny app on a local network 
#'
#' @importFrom shiny runApp
#' @export runJaccardApp
runJaccardApp <- function() {
  appDir <- system.file("shinyapps", "jaccard-shiny", package = "jaccard")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `jaccard`.", call. = FALSE)
  }

  runApp(appDir, display.mode = "normal")
}