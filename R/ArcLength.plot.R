ArcLength.plot <- function(arclength, arclengthvec, titlestr=NULL) {
  # Last modified 18 May 2022 by Jim Ramsay
  indfine <- seq(0,100,len=101)
  df <- data.frame(indfine=indfine, arclengthvec=arclengthvec)
  pd <- ggplot2::ggplot(df, aes(indfine, arclengthvec)) +
  geom_line(size=2) +
  geom_abline(intercept = 0, slope = arclength/100, linetype=2) +
  xlim(c(0,100)) +
  ylim(c(0,arclength)) +
  xlab("Percentile Index") +
  ylab("arc length (M-bits)") +
  labs(title=titlestr) +
  theme(axis.title=element_text(size=16,face="bold"))
  return(pd)
}
