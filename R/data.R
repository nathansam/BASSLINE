#' VA Lung Cancer Trial Dataset
#'
#' Data from a trial in which a therapy (standard or test chemotherapy)
#'  was randomly applied to 137 patients who were diagnosed with inoperable
#'  lung cancer. The survival times of the patients were measured in days
#'  since treatment.
#'
#'
#' @format A matrix with 137 rows and 8 variables:
#' \describe{
#'   \item{Time}{Survival time (in days)}
#'   \item{Cens}{0 or 1. If 0 the observation is right censored}
#'   \item{Intercept}{The intercept}
#'   \item{Treat}{The treatment applied to the patient (0: standard, 1: test)}
#'   \item{Type.1}{The histological type of the tumor (1: type 1, 0: otherwise)}
#'   \item{Type.2}{The histological type of the tumor (1: type 2, 0: otherwise)}
#'   \item{Type.3}{The histological type of the tumor (1: type 3, 0: otherwise)}
#'   \item{Status}{A continuous index representing the status of the patient:
#'       10—30 completely hospitalized, 40—60 partial confinement, 70—90
#'       able to care for self.}
#'   \item{MFD}{The time between the diagnosis and the treatment (in months)}
#'   \item{Age}{Age (in years)}
#'   \item{Prior}{Prior therapy, 0 or 10}
#' }
#' @source Appendix I of Kalbfleisch and Prentice (1980).
"cancer"
