#' Vzorkovani typu 2
#'
#' Vybere urcity procento fragmentu z kazdeho genomu, vraci vhodne hodnoty prahu. Vytvari referencni soubory ve slozce classDNAlib_sampl2.
#' @param fold slozka kde jsou umisteny soubory referencnich dat
#' @param samples procento fragmetu z genomu
#' @keywords random, sampling
#' @export
#' @examples
#' param <- sampl_2("slozka/", 0.3)
#' param <- sampl_2("slozka", 0.09)

sampl_2<-function(fold, samples){
	# kontorola jestli je samples v rozsahu 0 - 1
	if(samples > 1 || samples < 0) {
		print("Error: Samples musi byt od 0 do 1.")
		return(-1)
	}
	# nacteni potrebnych souboru
	mat2<-read.table(paste(fold,"data",sep = "/"), header = TRUE, row.names = NULL)
	rl<-read.table(paste(fold,"data_list",sep = "/"), header = FALSE, row.names = NULL)
	rn<-read.table(paste(fold,"data_name",sep = "/"), header = FALSE, row.names = NULL)

	#vytvoreni pomocnych promennych
	mat<-matrix(, nrow = 0, ncol = 256)
	rns<-list()
	eu<-0;ar<-0;pr<-0

	# vzorkovani kazdeho genomu
	s<-1;e<-0
	for(i in rl){
		e<-e+i
		# kontrola jestli uzivatel nechce vic vzorku nez je mozne
		if(samples > e-s){
			print("Error: Samples musi byt mensi.")
			return(-1)
		}
		sam<-sample(s:e, round(samples*i), replace = FALSE)
		# pocitani pomocnych hodnot pro prahy
		for(ii in sam){
			mat<-rbind(mat, mat2[ii,])
			rns<-c(rns,rn[ii])
			if(length(grep("eu-",rn[1,ii])) > 0) {eu<-1+eu}
			if(length(grep("ba-",rn[1,ii])) > 0) {pr<-1+pr}
			if(length(grep("ar-",rn[1,ii])) > 0) {ar<-1+ar}
		}
		s<-s+i
	}
	# vystupni soubory
	dir.create("classDNAlib_sampl2")
	write.table(mat, file="classDNAlib_sampl2/data", row.names = FALSE, col.names = TRUE)
	write.table(rns, file="classDNAlib_sampl2/data_name", row.names = FALSE, col.names = FALSE)

	# hodnoty prahu
	z<-ar+eu+pr
	return(c(eu/z, ar/z, pr/z))
}