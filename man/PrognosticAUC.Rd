% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrognosticModel.R
\name{PrognosticAUC}
\alias{PrognosticAUC}
\title{Calculate Prognostic Area Under the Curve (AUC)}
\usage{
PrognosticAUC(model, newx, s, acture.y)
}
\arguments{
\item{model}{A survival model object from which predictions will be made.
This model should be capable of generating risk scores, such as a Cox proportional hazards model.}

\item{newx}{A matrix or data frame containing new input data for which predictions are to be made.
This data should have the same features as used to train the model.}

\item{s}{The value of the penalty parameter 'lambda' at which predictions are requested.
This can be a specific value or a character string specifying a criterion,
such as 'lambda.min' or 'lambda.1se', commonly used in models like those from `glmnet`.}

\item{acture.y}{A data frame containing actual survival data, expected to have at least two columns:
'time', which contains the survival time, and 'status', which is the event indicator (1 if the event occurred, 0 otherwise).}
}
\value{
A data frame with two columns containing the AUC values for the 30th and 90th percentile survival times.
        These are named 'probs.3' and 'probs.9' respectively.
}
\description{
This function evaluates the prognostic ability of a survival model by calculating the Area Under the
Curve (AUC) of time-dependent Receiver Operating Characteristic (ROC) curves at specified time points.
The function uses predictions made by the model to compute ROC statistics and AUC at the 30th and 90th
percentiles of survival time, which are commonly used thresholds to assess short-term and long-term risk.
}
\examples{
# Assuming 'fit' is a fitted survival model such as from `coxph` or `glmnet`
new_data <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
actual_outcome <- data.frame(time = rexp(100, rate = 0.1), status = rbinom(100, size = 1, prob = 0.5))
auc_results <- PrognosticAUC(fit, newx = new_data, s = "lambda.min", acture.y = actual_outcome)
}
