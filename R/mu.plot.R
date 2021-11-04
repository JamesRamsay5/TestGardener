mu.plot <- function(mufine, scrrng, titlestr=NULL) {
  
  indfine <- seq(0,100,len=101)
  df <- data.frame(indfine=indfine, mufine=mufine)
  p <- ggplot(df, aes(indfine, mufine)) +
    geom_line(size=2) +
    geom_abline(intercept = scrrng[1], slope = (scrrng[2]-scrrng[1])/100, linetype=2) +
    xlim(c(0,100))           + ylim(scrrng) +
    xlab("Percentile Index") + ylab("Expected Test Score") +
    labs(title=titlestr) +
    theme(axis.title=element_text(size=16,face="bold"))
  print(p)
  return(p)
}
  
  