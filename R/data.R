#' VA Lung Cancer Trial Datatset
#'
#' Data from a trial in which a therapy (standard or test chemotherapy)
#'  was randomly applied to 137 patients who were diagnosed with inoperable
#'  lung cancer. The survival times of the patients were measured in days
#'  since treatment.
#'
#'
#' @format A data frame with 137 rows and 8 variables:
#' \describe{
#'   \item{treat}{The treatment applied to the patient (0: standard, 1: test)}
#'   \item{type}{The histological type of the tumor}
#'   \item{time}{Survival time (in days)}
#'   \item{censor}{0 or 1. If 0 the observation is right censored}
#'   \item{status}{A continuous index representing the status of the patient:
#'       10—30 completely hospitalized, 40—60 partial confinement, 70—90
#'       able to care for self.}
#'   \item{mfd}{The time between the diagnosis and the treatment (in months)}
#'   \item{age}{Age (in years)}
#'   \item{prior}{Prior therapy, 0 or 10}
#' }
#' @source \url{http://doi.org/10.7717/peerj.2555}
"cancer"
