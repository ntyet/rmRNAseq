#' Reparameterized Design Matrix
#'
#' This is a reparameterized  matrix from the original design matrix that
#' includes \code{time, line, time x line} such that each column of the
#' reparameterized design matrix corresponds to a test of interest such as
#' line main effect, \code{time2 - time0}, \code{time6 - time0}, \code{time24 -
#' time0}, and interaction effect of \code{timei - time0} with \code{line}. When
#' using this reparameterized design matrix in the analysis, for example,
#' regression coefficients corresponding to its second column \code{(line2)} represents
#' the \code{Line} main effect.
#' @format A design matrix with 8 columns \code{Intercept, line2, time2, time6,
#'   time24, linetime2, linetime6, linetime24}.
#' @usage design
#' @examples
#' data(design)
#' colnames(design)
#' head(design)
"design"

#'Covariate Set Associated with RFI RNA-seq
#'
#'This is a covariate set containing variables/measurements of the pigs that
#'RNA-seq data were obtained.
#'
#'@format \code{covset} is a dataframe with 32 rows and 11 columns,
#'  correspinding to 32 RNA-seq sample from 8 pigs at 4 times. 11 columns
#'  include \describe{ \item{time}{a continuous vector of mapping time points}
#'  \item{timef}{a factor of original time points} \item{line}{a factor of
#'  original line of each pig} \item{ear}{ear id of each pig} \item{line2,
#'  time2, ...}{the columns corresponding to design matrix transformation such that
#'  each column is a corresponding test statistics.} }
#'@examples
#'data(covset)
#'dim(covset)
#'colnames(covset)
#'head(covset)
"covset"

#' RFI RNA-seq Data
#'
#' This is the RFI RNA-seq data motivating our paper.
#'
#' @format The dataset has 11911 rows and 32 columns corresponding to 32 RNA-seq
#'   samples from 8 pigs, each pig has 4 RNA-seq samples measured repeatedly at
#'   0, 2, 6, 24 hours after LPS treatments.
#'
#' @examples
#' data(dat)
#' dim(dat)
#' head(dat)
"dat"

#'Data Containing Results of Our Proposed Method Applying to RFI RNA-seq data
#'
#'This data set contains outputs of \code{\link[limma]{voom}},
#'\code{\link{glsCAR1}} and pqvalue of the tests for interested contrasts when analyzing the RFI RNA-seq
#'data using our proposed method.
#'
#'@format A list with 3 components \describe{
#'\item{ori.res}{a list consisting
#'  of 2 components: \code{v}, which is a voom output including several
#'  components; \code{newlm}, which is the output from \code{glsCAR1}}
#'  \item{pqvalue}{a list of 2 components \code{pv} and
#'  \code{qv}. \code{pv} is a data frame of pvalues for each test such as
#'  \code{line2
#'    time}; \code{qv} is a data frame of corresponding q-values of the above
#'  p-values, calculated using the method by Nettleton 2006.} }
#' @examples
#' data(res)
#' names(res)
#'dim(res$ori.res$v)
#'colnames(res$ori.res$v)
#' colnames(res$ori.res$newlm)
#' colnames(res$pqvalue$pv)
#' colnames(res$pqvalue$qv)
"res"

#' Data Containing results when analyzing RFI using TC_Symm
#'
#' This dataset contains the output of TC_Symm for 1000 genes
#' when applying the TC_Symm to the RFI dataset.
#' @examples
#' data(resSymm)
#' names(resSymm)
#' dim(resSymm[[1]]$v)
#' colnames(resSymm[[1]]$newlm)
#' colnames(resSymm[[2]]$pv)
#' colnames(resSymm[[2]]$qv)

"resSymm"
