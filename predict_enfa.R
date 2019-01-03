
Predict.enfa <- function (object.enfa, baseline.climate, new.climate, nf, ...){
	## object.enfa <<< should be a object of class "enfa"
	## baseline.climate <<<< object of class "spatialpixelsdataframe" with the climate data used to fit enfa model
	## new.climate <<<< object of class "spatialpixelsdataframe" with the NEW climate data to predict environmental suitability
	
#    if (!inherits(object.enfa, "enfa")){
#        stop("should be an object of class \"enfa\"")
#    warning("the enfa is not mathematically optimal for prediction:\n
#please consider the madifa instead")
#	}
    
    m.baseline <- apply(slot(baseline.climate, "data"), 2, mean)
    sd.baseline <- apply(slot(baseline.climate, "data"), 2, sd) 
    
#    if ((missing(nf)) || (nf > object.enfa$nf)){
        nf <- object.enfa$nf
 #       }
    Zli <- object.enfa$li[, 1:nf]
    f1 <- function(x) {rep(x, object.enfa$pr)}
    Sli <- apply(Zli, 2, f1) #extrai os valores de Zli para os pontos de ocorrencia da sp
    m <- apply(Sli, 2, mean) #media nos pts de presenca
    cov <- t(as.matrix(Sli)) %*% as.matrix(Sli)/nrow(Sli) #cov nos pts de presenca
    if (!missing("new.climate")){
    	new.climate.scale<- sweep(slot(new.climate, "data"),2, m.baseline) #diminui cada coluna de 'new.climate' da sua respectiva media em 'm.baseline'
  		new.climate.scale<- as.matrix(new.climate.scale) %*% diag(1/sd.baseline)

        Zli <- new.climate.scale %*% as.matrix(object.enfa$co)
        }
    maha <- mahalanobis(Zli, center = m, cov = cov)    
    map <- rasterFromXYZ(data.frame(new.climate@coords, maha))*-1        

    return(invisible(map))

}# ends function "Predict.enfa"


