`barnardw.test` <-
function(a,b,c,d,one.sided=FALSE,dp=0.001) {
  ret=.C("BarnardW",
    as.integer(a),
    as.integer(b),
    as.integer(c),
    as.integer(d),
    one.sided = as.integer(one.sided),
    dp = as.numeric(dp),
    nuisance.vector.x = as.double(vector("double",1.0+1.0/dp)),
    nuisance.vector.y = as.double(vector("double",1.0+1.0/dp)),
    wald.statistic = as.double(0.0))

  nuisance.matrix<-matrix(cbind(ret$nuisance.vector.x,ret$nuisance.vector.y),ncol=2)
  rett<-nuisance.matrix[nuisance.matrix[,2]==max(nuisance.matrix[,2]),]
  
  return(list(nuisance.matrix = nuisance.matrix,
              dp = ret$dp,
              contingency.matrix = matrix(c(ret[[1]],ret[[2]],ret[[3]],ret[[4]]),ncol=2,byrow=TRUE),
              alternative = c("Two Sided","One Sided")[as.integer(ret$one.sided==1)+1],
              wald.statistic = ret$wald.statistic,
              nuisance.parameter = rett[[1]],
              p.value = rett[[2]]))
}

