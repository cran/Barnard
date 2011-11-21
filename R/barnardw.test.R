`barnardw.test` <-
function(n1,n2,n3,n4,dp=0.001) {
  ret=.C("BarnardW",
    as.integer(n1),
    as.integer(n2),
    as.integer(n3),
    as.integer(n4),
    dp = as.numeric(dp),
    nuisance.vector.x = as.double(vector("double",1.0+1.0/dp)),
    nuisance.vector.y0 = as.double(vector("double",1.0+1.0/dp)),
    nuisance.vector.y1 = as.double(vector("double",1.0+1.0/dp)),
    wald.statistic = as.double(0.0),
    nuisance.parameter = as.double(vector("double",2.0)),
    p.value = as.double(vector("double",2.0)))

  nuisance.matrix<-matrix(cbind(ret$nuisance.vector.x,ret$nuisance.vector.y0,ret$nuisance.vector.y1),ncol=3)
  
  return(list(nuisance.matrix = nuisance.matrix,
              dp = ret$dp,
              contingency.matrix = matrix(c(ret[[1]],ret[[2]],ret[[3]],ret[[4]]),ncol=2,byrow=TRUE),
              alternative = c("One Sided","Two Sided"),
              wald.statistic = ret$wald.statistic,
              nuisance.parameter = ret$nuisance.parameter,
              p.value = ret$p.value))
}

