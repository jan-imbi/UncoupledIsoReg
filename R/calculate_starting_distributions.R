calculateNearestNeighbourDistribution <- function(outputSpace.vals, comparativeSpace.vals){
  ret <- double(length(outputSpace.vals))
  for (val in comparativeSpace.vals){
    ind <- which.min(abs(val - outputSpace.vals))
    ret[ind] <- ret[ind] + 1
  }
  ret <- ret /sum(ret)
  ret
}

calculateMeanValueDistribution <- function(outputSpace.vals, comparativeSpace.vals){
  ret <- double(length(outputSpace.vals))
  val <- sum(comparativeSpace.vals) / length(comparativeSpace.vals)
  ind <- which.min(abs(val - outputSpace.vals))
  ret[ind] <- 1
  return(ret)
}

calculateMedianValueDistribution <- function(outputSpace.vals, comparativeSpace.vals){
  ret <- double(length(outputSpace.vals))
  if (length(comparativeSpace.vals) %% 2 == 0){
    val <- (comparativeSpace.vals[length(comparativeSpace.vals)/2] + comparativeSpace.vals[length(comparativeSpace.vals)/2 + 1])/2
  }
  else{
    val <- comparativeSpace.vals[ceiling(length(comparativeSpace.vals)/2)]
  }
  ind <- which.min(abs(val - outputSpace.vals))
  ret[ind] <- 1
  return(ret)
}

calculateLargestQuantileDistribution <- function(outputSpace.vals, comparativeSpace.vals){
  ret <- double(length(outputSpace.vals))
  indices <- double(length(outputSpace.vals))
  for (val in comparativeSpace.vals){
    ind <- which.min(abs(val - outputSpace.vals))
    indices[ind] <- indices[ind] + 1
  }
  ret[which.max(indices)] <- 1
  ret
}
