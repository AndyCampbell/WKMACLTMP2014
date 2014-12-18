

armaSearch = function(
   xx,
   minOrder=c(0,0),
   maxOrder=c(5,5),
   trace=FALSE )
{
   bestAic = 1e9
   len = NROW( xx )
   for( p in minOrder[1]:maxOrder[1] ) for( q in minOrder[2]:maxOrder[2] )
   {
      if( p == 0 && q == 0 )
      {
         next
      }

      formula = as.formula( paste( sep="", "xx ~ arma(", p, ",", q, ")" ) )

      fit = tryCatch( armaFit( formula, data=xx ),
                      error=function( err ) FALSE,
                      warning=function( warn ) FALSE )
      if( !is.logical( fit ) )
      {
         fitAic = fit@fit$aic
         if( fitAic < bestAic )
         {
            bestAic = fitAic
            bestFit = fit
            bestModel = c( p, q )
         }

         if( trace )
         {
            ss = paste( sep="", "(", p, ",", q, "): AIC = ", fitAic )
            print( ss )
         }
      }
      else
      {
         if( trace )
         {
            ss = paste( sep="", "(", p, ",", q, "): None" )
            print( ss )
         }
      }
   }

   if( bestAic < 1e9 )
   {
      return( list( aic=bestAic, fit=bestFit, model=bestModel ) )
   }

   return( FALSE )
}




