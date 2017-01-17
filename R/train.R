#' Trenovani na referencni data
#'
#' Provede shlukovani kmeans, vraci informace o shlucich.
#' @param fold slozka kde jsou umisteny soubory referencnich dat
#' @keywords clustering, kmeans
#' @export
#' @examples
#' cl <- train("slozka/")
#' cl <- train("slozka")

train<-function(fold){
	# nacte potrebny soubor
	mat<-read.table(paste(fold,"data",sep = "/"), header = TRUE, row.names = NULL)
	k<-round(sqrt(nrow(mat)/2))
	# kmeans clustering
	cl<-kmeans(mat, k, algorithm="Forgy")
	return(cl)
}