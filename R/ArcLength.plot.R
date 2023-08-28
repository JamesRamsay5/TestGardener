ArcLength.plot <- function(arclength, arclengthvec, titlestr=NULL) {
  # Last modified 3 June 2022 by Jim Ramsay
  indfine <- seq(0,100,len=101)
  df <- data.frame(arclengthvec=arclengthvec, indfine=indfine)
  pd <- ggplot2::ggplot(df, aes(arclengthvec, indfine)) +
  geom_line(size=2) +
  geom_abline(intercept = 0, slope = 100/arclength, linetype=2) +
  xlim(c(0,arclength)) +
  ylim(c(0,100)) +
  xlab("Arc Length or Information (M-bits)") +
  ylab("Score Index") +
  labs(title=titlestr) +
  theme(axis.title=element_text(size=16,face="bold"))
  return(pd)
}
