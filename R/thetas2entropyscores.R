#' Convert theta scores to entropy bit scores
#'
#' This function takes a vector of theta scores from a fitted model and converts them to entropy bit scores.
#'
#' @param model An fitted TestGardener model object returned by the Analyze function.
#' @param thetas A numeric vector of theta scores.
#' @param items A numeric vector indicating which items to use for computation. By default, it uses all items.
#' @param grid_size An integer specifying the size of the theta grid used for entropy score computation. A higher value leads to improved accuracy. Default is 10000.
#' @param return_grid Whether or not to return the entropy score for each value in the grid used for computation or only the entropy scores for the input thetas. Default is FALSE.
#'
#' @return A matrix with columns: 'theta' (theta scores) and 'entropy_score' (theta corresponding entropy scores).
#'
#' @examples
#' # Read data and fit a model
#' U <- as.matrix(read.table("data/Quant_13B_problem_U.txt"))
#' key <- scan("data/Quant_13B_problem_key.txt","o")
#' key <- as.numeric(unlist(stringr::str_split(key,"")))
#' noption <- rep(5, ncol(U)) # 5 options per item
#' 
#' ScoreList <- list() # option scores
#' for (item in 1:ncol(U)){
#'   scorei <- rep(0,noption[item])
#'   scorei[key[item]] <- 1
#'   ScoreList[[item]] <- scorei
#' }
#' 
#' optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)
#' Math_dataList <- TestGardener::make.dataList(U, key, optList)
#' 
#' AnalyzeResult <- TestGardener::Analyze(
#'   Math_dataList$percntrnk, 
#'   Math_dataList$thetaQnt, 
#'   Math_dataList, 
#'   ncycle = 10, 
#'   itdisp=FALSE, 
#'   verbose=TRUE
#' ) 
#' 
#' # Get the theta scores in the data and compute their corresponding entropy bit scores
#' thetas <- AnalyzeResult$parList[[10]]$theta
#' entropy_scores <- thetas2entropyscores(AnalyzeResult, thetas)
#' entropy_scores
#' 
thetas2entropyscores <- function(model, thetas, items = 1:length(model$parList[[iter]]$WfdList), grid_size = 10000, return_grid = FALSE) {
  theta_grid <- setdiff(seq(0, 100, length.out = grid_size), thetas) |> # remove potential duplicates from grid
    c(thetas) |> # merge with sample thetas
    sort()
  iter <- length(model$parList)
  entropy_scores <- vector("numeric", length = length(theta_grid))
  for (item in items) {
    surprisals <- eval.surp(theta_grid, model$parList[[iter]]$WfdList[[item]]$Wfd)
    categories <- ncol(surprisals)
    surprisals <- surprisals / log(2, categories) # uses 2bit surprisal for all items
    entropies <- rowSums((1 / 2^surprisals) * surprisals)
    
    total_item_dist <- 0
    for (i in 2:length(entropies)) {
      # add up total item entropy change to theta corresponding to surp[i, ]
      total_item_dist <- total_item_dist + sum(abs(entropies[i] - entropies[i - 1]))
      # add this to the total for the same theta
      entropy_scores[i] <- entropy_scores[i] + total_item_dist
    }
  }
  if (return_grid) {
    return(cbind(theta = theta_grid, entropy_score = entropy_scores))
  } else {
    # Extract the entropy scores from the grid
    matching_indices <- which(theta_grid %in% thetas)
    ordered_entropy_scores <- entropy_scores[matching_indices]
    # Order the results based on the original order of the 'thetas' argument
    order_vec <- order(order(thetas))
    return(data.frame(theta = thetas, entropy_score = ordered_entropy_scores[order_vec]))
  }
}
