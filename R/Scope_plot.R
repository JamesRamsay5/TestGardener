Scope_plot <- function(infoSurp, infoSurpvec, titlestr=NULL) {
  # Last modified 1 November 23 by Jim Ramsay
  indfine <- seq(0,100,len=101)
  df <- data.frame(indfine=indfine, infoSurpvec=infoSurpvec)
  pd <- ggplot2::ggplot(df, aes(indfine, infoSurpvec)) +
  geom_line(linewidth=2) +
  geom_abline(intercept = 0, slope = infoSurp/100, linetype=2) +
  xlim(c(0,100)) +
  ylim(c(0,infoSurp)) +
  xlab("Score Index") +
  ylab("Scope or Information (M-bits)") +
  labs(title=titlestr) +
  theme(axis.title=element_text(size=16,face="bold"))
  return(pd)
}
