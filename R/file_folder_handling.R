#' @title Create output directory using the currently-running file name
#' @description \code{outdir_create} Create output directory using the currently-running file name
#' @param save_dir_name Logical. If TRUE, output directory name is saved.
#' @param str_end Integer. Trim `str_end` characters from the file name.
#' @return output_dir Output directory name
#' @export
#' @examples
#' # outdir_create()
outdir_create <- function(save_dir_name = TRUE, str_end = 3) {
  # Extract output directory name
  output_dir <- rstudioapi::getSourceEditorContext()$path %>%
    basename %>% stringr::str_sub(end = - str_end) %>%
    paste0("Out")
  # Create output directory
  dir.create(output_dir)
  # Retrieve output directory name
  if (save_dir_name) return(output_dir)
}


#' @title Save workspace
#' @description \code{save_workspace} xxx
#' @param output_dir xxx.
#' @param list xxx.
#' @export
#' @examples
#' # output_dir <- outdir_create()
#' # save_workspace(output_dir)
save_workspace <- function(output_dir, list = ls(all.names = TRUE)) {
  save(list = ls(all.names = TRUE),
     file = paste0(output_dir, "/", output_dir, ".RData"))
}


#' @title Save session information
#' @description \code{save_session_info} xxx
#' @param output_dir xxx.
#' @param session_info_dir xxx.
#' @export
#' @examples
#' # output_dir <- outdir_create()
#' # save_session_info(output_dir)
save_session_info <- function(output_dir, session_info_dir = "00_SessionInfo") {
  writeLines(utils::capture.output(utils::sessionInfo()),
             paste0(session_info_dir, output_dir, "_", substr(Sys.time(), 1, 10), ".txt"))
}
