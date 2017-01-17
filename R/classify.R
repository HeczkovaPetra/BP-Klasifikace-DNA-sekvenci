#' Klasifikace DNA sekvence
#'
#' Zpracuje vstupni soubor, najde k nemu nejblizsi shluky, zjistuje pocet zastupcu nadtridy v nejblizsich shlucich, klasifiku podle prahovych hodnot. Potrebuje k cinnosti knihovnu Biostrings. 
#' @param name jmeno a cesta k FASTA souboru s hledanou sekvenci
#' @param fold slozka kde jsou umisteny soubory referencnich dat
#' @param param seznam prahovych hodnot, vystup funkce sampl_1 nebo sampl_2
#' @param csc pocet nejblizsich fragmentu, zadava uzivatel
#' @param cl informace o shlucich, vystup funkce train
#' @keywords classification, archea, prokaryote, eukaryote
#' @export
#' @examples
#' classify("ar-ea-015.fa","slozka/",c(0.5,0.4,0.3),100, cl)
#' classify("cesta/eu-ne-208.fa","slozka",p,30, cl)

classify<-function(name, fold, param, csc, cl){
	# nacteni potrebnych soubory
	mat<-read.table(paste(fold,"data",sep = "/"), header = TRUE, row.names = NULL)
	rn<-read.table(paste(fold,"data_name",sep = "/"), header = FALSE, row.names = NULL)
	k<-nrow(cl$centers)
	# zpracovani hledane sekvence na priznakovy vektor
	t<-numeric()
	act<-readDNAStringSet(file=name)
	tt<-DNAStringSet(paste(act,collapse=""))
	t<-oligonucleotideFrequency(tt, 4)/sum(width(tt))

	# pocitani pomocne hodnoty pro prahy
	trl<-matrix(0, nrow = k, ncol = 3) # V1 eu + V2 pr + V3 arch, radky shluky
	for(j in 1:nrow(mat)){
		if(length(grep("eu-",rn[1,j])) > 0) {trl[cl$cluster[j],1]<-1+trl[cl$cluster[j],1]}
		if(length(grep("ba-",rn[1,j])) > 0) {trl[cl$cluster[j],2]<-1+trl[cl$cluster[j],2]}
		if(length(grep("ar-",rn[1,j])) > 0) {trl[cl$cluster[j],3]<-1+trl[cl$cluster[j],3]}
	}

	# pocitani vzdalenosti hledane sekvence od center clusteru
	vec<-list()
	for(ii in 1:k) {
		vec<-c(vec,sum(mapply(function(x,y) {abs(x-y)}, t, cl$centers[ii,])))
	}
	
	# pocet nejblizsich fragmentu by nemel byt vetsi nez pocet fragmentu v referencni sade
	if(csc > nrow(mat)){
		print("Warning: Referencni data by mela byt vetsi, klasifikator nemusi fungovat spravne.")
		csc<-nrow(mat)
	}
	suma<-0;tr1_e<-0;tr1_p<-0;tr1_a<-0;tr2<-0;vec2<-list();vec3<-list()
	while(suma < csc){
		min<-which.min(vec) # hledani nejblizsich shluku
		suma<-suma+cl$size[min]
		vec[min]<-NaN

		tr1_e<-tr1_e+trl[min,1] # pocet eu fragmentu ve shluku
		tr1_p<-tr1_p+trl[min,2] # pocet proc fragmentu ve shluku
		tr1_a<-tr1_a+trl[min,3] # pocet arch fragmentu ve shluku
		tr2<-tr2+cl$size[min] # celkovy pocet fragmentu v danem shluku

		# orezani na presny pocet nejblizsich fragmentu
		if(suma > csc){
			for(xx in 1:length(rn)) {
				if(cl$cluster[xx] == min){
					vec2<-c(vec2,sqrt(sum(mapply(function(x,y) {(x-y)^2}, t, mat[xx,]))))
					vec3<-c(vec3, rn[xx])
				}
			} # end for
			for(yy in 1:(suma-csc)){
				max<-which.max(vec2)
				if(length(grep("eu-",vec3[[max]])) > 0){tr1_e<-tr1_e-1}
				if(length(grep("ba-",vec3[[max]])) > 0){tr1_p<-tr1_p-1}
				if(length(grep("ar-",vec3[[max]])) > 0){tr1_a<-tr1_a-1}
				vec2[max]<-NaN
			}
		} # end if
	} # end while

	# klasifikace podle prahovych hodnot, kdyz je nad prahem vice hodnot vybere se nejvyzsi
	if(tr1_e/tr2>=param[1] && tr1_a/tr2>=param[2] && tr1_p/tr2>=param[3]){
		mm<-which.max(c(tr1_e, tr1_a, tr1_p))
	} else if(tr1_e/tr2>=param[1] && tr1_a/tr2>=param[2]){
		mm<-which.max(c(tr1_e, tr1_a, NULL))
	} else if(tr1_e/tr2>=param[1] && tr1_p/tr2>=param[3]){
		mm<-which.max(c(tr1_e, NULL, tr1_p))
	} else if(tr1_a/tr2>=param[2] && tr1_p/tr2>=param[3]){
		mm<-which.max(c(NULL, tr1_a, tr1_p))
	} else if(tr1_e/tr2 >= param[1]){mm<-1
	} else if(tr1_a/tr2 >= param[2]){mm<-2
	} else if(tr1_p/tr2 >= param[3]){mm<-3
	} else {mm<-0}		

	# vystup
	if(mm == 0){print("unknown")}
	if(mm == 1){print("eukaryota")}
	if(mm == 2){print("archea")}
	if(mm == 3){print("prokaryota")}
}


