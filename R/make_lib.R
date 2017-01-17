#' Vyvoreni referencni sady priznakovych vektoru
#'
#' Funkce vezme FASTA soubory ze slozky fold, rozdeli je na fragmenty a spocita priznakove vektory. Vysledna data jsou ve slozce classDNAlib. Potrebuje k cinnosti knihovnu Biostrings.
#' @param fold slozka kde jsou umisteny FASTA soubory
#' @keywords fragments, FASTA, oligonucleotide frequency
#' @export
#' @examples
#' make_lib("slozka/")
#' make_lib("slozka")

make_lib<-function(fold) {
	# seznam FASTA souboru
	list_name<-list.files(fold)

	mat<-matrix(, nrow = 0, ncol = 256);t<-numeric()
	rn<-list() # jmena sekvenci
	rl<-list() # pocet sekvenci vzniklych z jednoho .fa

	# rozdeluenije sekvence na fragmenty, pocitani frekvenci oligonukleotidu
	for (file_name in list_name){
		act<-readDNAStringSet(file=paste(fold,file_name,sep = "/"))

		# upravi sekvence na vhodny tvar
		w<-sum(width(act))
		seq<-DNAStringSet(paste(act,collapse=""))

		s<-1;len<-10000;sublist<-list() # pomocne promenne	
		x<-0;num<-1
		# v cyklu delime sekvence
		while(s < w) {
			x<-x+1
			# deleni na fragmenty
			if((s+len) >= w){temp<-subseq(seq,start=s,width=(w-s+1))} # celych 10 000 bp
			else {temp<-subseq(seq,start=s,width=len)} # zbytek sekvence na konci

			# spocitani frekvenci tetraoligonukleotidu
			t<-oligonucleotideFrequency(temp, 4)/length(temp)
			mat<-rbind(mat,t)
			rn<-c(rn,paste(formatC(num, width=5, flag="0000"),file_name,sep="_"))
			num<-num+1

			s<-s+len
			sublist<-c(sublist,temp)
		}
		rl<-c(rl,x)
	}
	# vystupni soubory
	dir.create("classDNAlib")
	write.table(mat, file="classDNAlib/data", col.names = TRUE, row.names = FALSE)
	write.table(rl, file="classDNAlib/data_list", col.names = FALSE, row.names = FALSE)
	write.table(rn, file="classDNAlib/data_name", col.names = FALSE, row.names = FALSE)
}
