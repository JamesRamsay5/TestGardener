\name{NEWS}
\title{Recent changes to the TestGardner package}
\section{Changes in TestGardener version 3.2.7 (2023-22-11)}{
\itemize{
  \item{The demo codes for the Quant_13B_problem data and the
  SDS data have been removed.  They are extensive analyses rather than
  the single topic examples that demos need to be. But the examples
  in the .Rd files have been extended, revised and supplied with
  more comments so that they can serve the demo process.}
  \item{Likewise the vignette for the Quant_13B_problem data has
  been renewed in order to speed up the checking and installing
  processes.  These vignette RMarkdown codes will be provided instead
  in gitHub and possibly CRAN along with other sample analyses in
  RMarkdown format.
  }
  \item{We made the code more readable and interpretable by changing
  object names derived from mathematics to names that describe what
  the purpose of the object is.  For example, "theta" was changed to
  "index" since these values are score index values, and "arclength"
  was changed to "infoSurp".  Likewise, the many objects names beginning
  with capital "W" all pertain to surprisal values, matrices, lists and
  other structures.  For example, a functional data basis object is
  names "Sbasis", a functional data curve is "Sfd" and a functional
  fdPar object to "SfdPar".}
}
}
\section{Changes in TestGardener version 3.2.5 (2023-04-10)}{
\itemize{
  \item{Extensive work on the test vignette Quant_13B_problem.Rmd and on the 
  clinical rating scale vignette hads_dep_3.Rmd (not included in vignette 
  folder, but in inst/HADS) led to detecting multiple problems with the 
  plotting codes using ggplot2 and plotly packages.  
  The argument lists in ICC.plot.R and Wpca.plot.R have been changed.  
  Hopefully these are now safe to use.}
  }
}
\section{Changes in TestGardener version 3.2.4 (2023-23-08)}{
\itemize{
  \item{TestGardener has evolved substantially over 2023 and since 
  the last news item.  There are now 38 functions, only 9 of which have
  not been altered over this period.  The major changes have been
  in the graphical displays because we have switched to ggplot2 and
  plotly to construct them.  Juan Li has become an expert in these
  two packages, and demonstrated a real talent for statistical 
  graphics.  We transported all functions having to do with surprisal
  curves in the fda package to versions in TestGardener.  Unfortunately,
  we have beenb reduced to havingonly only one vignette,
  Quant_13B_problem.Rmd, to CRAN restrictions on time taken to 
  install TestGardener.  We have two large papers under review, one
  for testing data, and the other for rating scales.}
  }
}
\section{Changes in TestGardener version 3.1.3 (2022-20-12)}{
\itemize{
  \item{We observe that the function used to define the best estimate of 
    the score index for an examinee about 10 to 15 percent of the time
    has large nearly flat intervals of low values.  Moreover, it is
    often the case that there are multiple minimizing locations, and
    function thetafun does not always identify the lowest of these.
    New function thetasearch uses a fast grid search to identify the
    lowest point, and replaces the location identified by thetafun where
    required by the integer closest to the global minimum.  This is 
    usually done after the cycles are completed as a final step.}
  }
}
\section{Changes in TestGardener version 3.1.2 (2022-20-12)}{
\itemize{
  \item{The test analyst no longer needs to construct a "grbg"" object in order
    handle missing of illegitimate responses by treating them as an
    additional option, the garbage option as we call it.  This is now
    done automatically within function make.dataList.  Function make.dataList
    no longer needs an argument "grbg"."
    In order to makethis work, the user does need to adhere to a discipline.  
    For each item, the labels or indices of its options must be the first 
    positive integers.  
    For the quantitative SweSAT test, for example, 68 of the items have option 
    indices 1,...,4, and 12 have indices 1,...,5. For this large test given
    to over 55,000 examinees, the garbage option is essential for all items.
    It is assigned within make.dataList to an option with the index one larger 
    than required for the legitimate options, so that after make_dataList
    has done its work, the 68 items have indices 1,...,5 and the others
    indices 1,...,6.  
    Any integer within the data matrix U that exceeds the highest legitimate 
    index is re-assigned to the corresponding garbage index.  
    Do be aware that this process can change the U matrix that you set up 
    initially.
    Users often want to have non-integer labels for options.  This is still
    available by setting these up in the "optList" list vector, and these
    will appear in the most of the plots.
    }
  \item{This version of TestGardener also has various code improvements
    that speed up computation and improve the information in plots.
    }
  }
}
\section{Changes in TestGardener version 3.0.1 (2022-22-07)}{
\itemize{
  \item{Plots using the score index as the abscissa now use arc length along
    along the test information curve as the abscissa, since this interval is
    invariant under warpings of the index and therefore on a ratio scale. The 
    codes in the demo and vignette folders have been modified accordingly.
    }
  }
}
\section{Changes in TestGardener version 3.1.1 (2024-09-01)}{
\itemize{
  \item{A major clean-up and re-organization of the code was required to
  deal with the fact that a plot of the last item was drastically and 
  unbelievably distorted for all data analyses.  
  
  As a part of the cleanup I experimented with adding a new item designed to 
  track  item choices constructed to deal with running out of time or having 
  no clue to what answer was correct.  This experiment meant there were n + 1 
  items technically, even if only the first n were regular choices.  In the 
  end, it didn't seem that this added anything, but this codestill remains  
  because it may prove useful later.}
  }
}
