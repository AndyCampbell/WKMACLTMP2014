


################# FOURIER SURROGATES ############################
# This function was originally written by Tristan Rouyer at Ifremer-Sète (tristan.rouyer@ifremer.fr) and is inspired from the following works:
# Theiler J, Eubank S, Longtin A, Galdrikian B, Farmer JD (1992) Testing for nonlinearity in time-series : the method of surrogate data. Physica D 58: 77-94
# Schreiber T, Schmitz, A (2000) Surrogate time series. Physica D 142: 346-382
# From a time-series vector y, the function returns a matrix of n time series with mean, variance and power spectrum identical to the original time series
FourierSurrogates<-function(y,n) {

   L<-length(y)
   z<-fft(y)
   S=matrix(NA,nrow=L,ncol=n)
   for (i in 1:n)
   {
	   if ((length(z)-trunc(length(z)/2)*2)==0) {
	      ph<-2*pi*runif((L-1)/2,0,1)
	      ph<-c(0,ph,0,-ph[length(ph):1])
	   } else {
	      ph<-2*pi*round(runif(floor(length(z)/2),0,1),4);
	      ph<-c(0,ph,-ph[length(ph):1])
	   }
   	   z<-z*exp(1i*ph)
   	   s<-Re(fft(z,inverse=T)/L)
   	   S[,i]=s
   }
return(S)
}



