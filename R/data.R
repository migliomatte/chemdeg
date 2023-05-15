#' First order kinetic data
#'
#' Synthetic data from a first-order kinetic model with k=0.7
#' @usage ord1
#' @format
#' A data frame with 6 rows and 3 columns:
#' \describe{
#'   \item{t}{time}
#'   \item{concentration}{simulated concentration data at each time point}
#'   \item{std.error}{simulated experimental erro}
#' }
"ord1"

#' Total CQA with ascorbic acid (FOMT model)
#'
#' Degradation data of 1.2 mM 5-caffeoylquinic acid (5-CQA) in the presence of
#' 1.2 mM of ascorbic acid at 37°C. The data refer to total CQA concentration.
#'
#' @usage fomtdata
#' @format
#' A data frame with 8 rows and 2 columns:
#' \describe{
#'   \item{time_h}{Time in hours}
#'   \item{tCQA_AA}{Normalized concentration of total CQA measured at each
#'   time point}
#' }
#' @source Yusaku N. and Kuniyo I. (2013)
#' *Degradation Kinetics of Chlorogenic
#' Acid at Various pH Values and Effects of Ascorbic Acid and Epigallocatechin
#' Gallate on Its Stability under Alkaline Conditions*,
#' Journal of Agricultural and Food Chemistry, \doi{10.1021/jf304105w},
#'  fig 2D, solid diamonds
"fomtdata"

#' Urfa pepper ascorbic acid degradation (2-nd order)
#'
#' Data describing the degradation kinetics of ascorbic acid during dehydration
#' of Urfa peppers. The peppers were treated with hot air at 55, 65 and 75 °C.
#'
#' @usage urfa
#' @format
#' A data frame with 8 rows and 4 columns:
#' \describe{
#'   \item{time_min}{time in minutes}
#'   \item{AA_55}{normalized concentration of ascorbic acid of Urfa peppers
#'   dehydrated at 55°C}
#'   \item{AA_65}{normalized concentration of ascorbic acid of Urfa peppers
#'   dehydrated at 65°C}
#'   \item{AA_75}{normalized concentration of ascorbic acid of Urfa peppers
#'   dehydrated at 75°C}}
#' @source Ş. Dağhan, A. Yildirim, F. Mehmet Yilmaz, H. Vardin and
#' M. Karaaslan (2018)
#' *The effect of temperature and method of drying on isot (Urfa pepper)
#' and its Vitamin C degradation kinetics*,
#' Italian Journal of Food Science,
#' \doi{10.14674/IJFS-1070},
#' fig 5, hot-air
"urfa"
