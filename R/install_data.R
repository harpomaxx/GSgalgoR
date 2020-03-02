#' Title
#'
#' @return
#' @export
#'
#' @examples
use_rna_luad <- function(){
  dest_dir <- paste(find.package("galgoR"),"/data",sep="")
  dest_file <- paste(dest_dir,"/luad_data.rds",sep="")
  if (! file.exists(dest_file))
      utils::download.file("https://github.com/harpomaxx/galgo/raw/0aaa1f070715b2dea541eb1a80691c982e49984f/data/luad_data.rds",dest_file)
  data<-readRDS(dest_file)
}
