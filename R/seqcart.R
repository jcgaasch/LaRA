seqcart <- function(
  dataimp,
  misord,
  INDOBS,
  INDMIS,
  control1,
  control2
){

  for(k in misord){
    yobs <- dataimp[INDOBS[, k], k]
    Xobs <- dataimp[INDOBS[, k], !(names(dataimp)%in%k), drop = FALSE]
    Xmis <- dataimp[INDMIS[, k], !(names(dataimp)%in%k), drop = FALSE]
    rpartmethod <- ifelse(is.factor(yobs), "class", "anova")
    tree <- rpart(yobs ~., data = cbind(yobs, Xobs), method = rpartmethod,
      control = rpart.control(minbucket = control1, cp = control2))
    leafdonor <- floor(as.numeric(row.names(tree$frame[tree$where, ])))
    tree$frame$yval <- as.numeric(row.names(tree$frame))
    leafmis <- predict(tree, Xmis, "vector")
    dataimp[INDMIS[, k], k] <- sapply(leafmis, function(x){
      donorpool <- yobs[leafdonor == x]
      if(length(donorpool) == 1){
        obs <- donorpool
      }else{
        di <- sort(runif(length(donorpool) - 1))
        obs <- sample(x = donorpool, size = 1, prob = (c(di, 1) - c(0, di)))
      }
      return(obs)
    })
  }
  return(dataimp)

}
