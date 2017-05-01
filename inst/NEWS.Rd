\name{NEWS}
\title{News for Package \pkg{WRS2}}
\encoding{UTF-8}

\section{Changes in Version 0.9-2}{
  \itemize{
    \item Warnings in vignette removed.
    \item Bug in \code{med2way} fixed.
    } 
}
\section{Changes in Version 0.9-1}{
  \itemize{
    \item Vignette updated.
    } 
}
\section{Changes in Version 0.9-0}{
  \itemize{
    \item \code{Qanova} now allows for a vector of quantiles and applies multiple testing adjustment. 
    \item \code{t3way} stops when design is incomplete. 
    \item \code{ZYmediate} added for robust mediation test.
    \item \code{essays} data added. 
    } 
}
\section{Changes in Version 0.8-0}{
  \itemize{
    \item Bugs in \code{mcp2a} and \code{mcp2atm} fixed. 
    \item Argument \code{pro.dis} added in \code{pbad2way} for projected distance computation.
    \item One-sided testing in \code{yuenbt} added. 
    \item Functions for robust correlation coefficients added: \code{pbcor}, \code{wincor}, \code{pball}, \code{winall}. 
    \item Functions for two-sample correlation tests added: \code{twocor}, \code{twopcor}.
    \item Bug in comparison point determination in \code{ancova} fixed.
    \item Outputs in \code{sppba}, \code{sppbb}, and \code{sppbi} reorganized. 
    \item Several functions for comparing discrete distributions added: \code{binband}, \code{discANOVA}, \code{discmcp}, \code{discstep}
    \item Convenience function \code{medpb2} added for comparing two independent samples medians.
    \item Functions for quantile comparisons added: \code{qcomhd}, \code{Dqcomhd}, \code{Qanova}. 
    } 
}
\section{Changes in Version 0.4-0}{
  \itemize{
    \item Vignette added.
    \item Sanity check for incomplete designs in \code{t2way}. 
    \item \code{mest} for M-estimator and \code{winmean} for Winsorized mean computation added.
    \item Bugs in \code{pbad2way} for unbalanced designs and onestep M-estimator fixed.
    \item Some datasets added.
    \item Yuen's dependent samples t-test added (\code{yuend}).
    \item Running interval smoothers added.
    \item Bug in \code{bwtrim} removed.
    } 
}
\section{Changes in Version 0.3-1}{
  \itemize{
    \item bugfixes in \code{tsplit}, \code{sppbb}, \code{sppba}, \code{sppbi}.
    } 
}
\section{Changes in Version 0.3-0}{
  \itemize{
    \item Post hoc tests  with corresponding S3 methods added: \code{rmmcp}, \code{pairdepb}, \code{mcppb20}, \code{pbad2way}, \code{mcp2a}, \code{mcp2atm}.
    \item Yuen's tests added: \code{yuen}, \code{yuenbt}.
    } 
}
\section{Changes in Version 0.2-0}{
  \itemize{
    \item \code{med2way} added.
    \item \code{ancova}, \code{ancboot} including invisibility dataset added.
    \item \code{tsplit}, \code{sppba}, \code{sppbb}, \code{sppbi} for mixed anova designs added.
   } 
}
\section{Changes in Version 0.1-0}{
  \itemize{
    \item The first version includes the following functions: \code{med1way}, \code{t1way}, \code{t2way}, \code{t3way}, \code{t1waybt}.
  }
}

