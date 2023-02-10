#' @title Create output directory using the currently-running file name
#' @description \code{outdir_create} automatically detect the file name of the current srcipt and create output directory using the file name.
#' @param save_dir_name Logical. If `TRUE`, output directory name is saved to an object.
#' @param str_end Integer. Trim `str_end` characters from the file name.
#' @param suffix Character. Suffix that will be appended to the file name.
#' @return Output directory name.
#' @export
#' @examples
#' # outdir_create()
outdir_create <- function(save_dir_name = TRUE, str_end = 3, suffix = "Out") {
  # Extract output directory name
  output_dir <- rstudioapi::getSourceEditorContext()$path %>%
    basename %>% stringr::str_sub(end = - str_end) %>%
    paste0(suffix)
  # Create output directory
  dir.create(output_dir)
  # Retrieve output directory name
  if (save_dir_name) return(output_dir)
}


#' @title Save workspace using an output directory name
#' @description \code{save_workspace} saves current workspace in a user-specified output directory.
#' @param workspace_dir Character. Workspace will be saved in this directory.
#' @param workspace_file_name Character (optional). Workspace file name.
#' @export
#' @examples
#' # output_dir <- outdir_create()
#' # save_workspace(output_dir)
save_workspace <- function(workspace_dir, workspace_file_name = NULL) {
  # Create output_dir if it does not exists
  if (!dir.exists(workspace_dir)) dir.create(workspace_dir)
  # Specify workspace file name
  if (is.null(workspace_file_name)) workspace_file_name <- workspace_dir
  # Save workspace
  save(list = ls(all.names = TRUE),
     file = paste0(workspace_dir, "/", workspace_file_name, ".RData"))
}


#' @title Save session information using an output directory name
#' @description \code{save_session_info} saves session information and append date information automatically.
#' @param session_info_dir Character. Name of the session information directory.
#' @param create_session_info_dir Logical. If `TRUE`, the function automatically creates an output directory for the session information
#' @param str_end Integer. Trim `str_end` characters from the file name.
#' @param suffix Character. Suffix that will be appended to the file name.
#' @export
#' @examples
#' # output_dir <- outdir_create()
#' # save_session_info(output_dir)
save_session_info <- function(session_info_dir = "00_SessionInfo", create_session_info_dir = FALSE, str_end = 3, suffix = NULL) {
  # Extract the running srcript file name
  current_file_name <- rstudioapi::getSourceEditorContext()$path %>%
    basename %>% stringr::str_sub(end = - str_end)
  if(!is.null(suffix)) current_file_name <- paste0(current_file_name, suffix)
  # Create session_info_dir if it does not exists
  if (!dir.exists(session_info_dir) & create_session_info_dir) {
    dir.create(session_info_dir)
  } else if (dir.exists(session_info_dir) & create_session_info_dir) {
    message("session_info_dir already exist. The directory is NOT overwrited.")
  } else if (!dir.exists(session_info_dir) & !create_session_info_dir){
    stop("session_info_dir does not exist. Please create it.")
  }
  # Save session information
  writeLines(utils::capture.output(utils::sessionInfo()),
             paste0(session_info_dir, "/", current_file_name, "_", substr(Sys.time(), 1, 10), ".txt"))
}
