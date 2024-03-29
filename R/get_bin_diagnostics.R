#' Helper: get top effects
#' @param prob A probability representing a percentile
#' @param data A data frame
#' @param outcome The name of the outcome column
#' @param exposure The name of the exposure column
#' @param covs A vector of covariates
#' @param pl Whether to use profile likelihood (default: FALSE)
#' @importFrom logistf logistf logistpl.control
#' @export
getTopEffects <- function(prob, data, outcome = "case", exposure, covs = NULL, pl = FALSE) {
    riskBin              <- paste0("Top_", prob)
    data[[riskBin]] <- ifelse(
        data[[exposure]] >= quantile(data[data[[outcome]] == 0, ][[exposure]], probs = 1 - prob),
        1, 0)
    if (!is.null(covs)) {
        tmp_f <- paste0(outcome, " ~ ", riskBin, " + ", paste(covs, collapse = " + "))
    } else {
        tmp_f <- paste0(outcome, " ~ ", riskBin)
    }
    vars <- c(outcome, riskBin)
    tmp_data <- na.omit(subset(data, select = vars))
    if (!is.null(covs)) {
        vars <- c(outcome, riskBin, covs)
    } else {
        vars <- c(outcome, riskBin)
    }
    fitAsso2 <- logistf::logistf(tmp_f, data = tmp_data, plcontrol = logistf::logistpl.control(maxit = 1E5), pl = pl)
    aimTwo::getValues(fitAsso2, riskBin)
}

#' Helper: get top% to middle% effects
#' @param prob A probability representing a percentile
#' @param middle A vector of two probabilities representing a percentile range
#' @param data A data frame
#' @param outcome The name of the outcome column
#' @param exposure The name of the exposure column
#' @param covs A vector of covariates
#' @param pl Whether to use profile likelihood (default: FALSE)
#' @importFrom logistf logistf logistpl.control
#' @export
getTopMidEffects <- function(prob = 0.1, middle = c(0.4, 0.6), data, outcome = "case", exposure, covs = NULL, pl = FALSE) {

    riskBin              <- paste0("Top", prob, "_Mid", middle[2] - middle[1])

    data[[riskBin]] <- fcase(
        data[[exposure]] >= quantile(data[data[[outcome]] == 0, ][[exposure]], probs = 1 - prob), 1,
        data[[exposure]] >= quantile(data[data[[outcome]] == 0, ][[exposure]], probs = middle[1])  &
            data[[exposure]] <= quantile(data[data[[outcome]] == 0, ][[exposure]], probs = middle[2]), 0,
        default = NA
    )

    if (!is.null(covs)) {
        tmp_f <- paste0(outcome, " ~ ", riskBin, " + ", paste(covs, collapse = " + "))
    } else {
        tmp_f <- paste0(outcome, " ~ ", riskBin)
    }

    vars <- c(outcome, riskBin)
    tmp_data <- na.omit(subset(data, select = vars))

    if (!is.null(covs)) {
        vars <- c(outcome, riskBin, covs)
    } else {
        vars <- c(outcome, riskBin)
    }

    fitAsso2 <- logistf::logistf(tmp_f, data = tmp_data, plcontrol = logistf::logistpl.control(maxit = 1E5), pl = pl)
    aimTwo::getValues(fitAsso2, riskBin)

}

#' Helper: get top power
#' @param prob A probability representing a percentile
#' @param data A data frame
#' @param outcome The name of the outcome column
#' @param exposure The name of the exposure column
#' @importFrom pwr EH.h pwr.2p2n.test
#' @export
getTopPower <- function(prob, data, outcome, exposure) {
    riskBin <- glue("Top_{prob}")
    data[[riskBin]] <- ifelse(
        data[[exposure]] >= quantile(data[data[[outcome]] == 0, ][[exposure]], probs = 1 - prob),
        1, 0)
    ptable <- table("Top" = data[[riskBin]], "trait" = data[[outcome]])
    if (nrow(ptable) != 2) {
        return(c("power" = NA, "h" = NA))
    }
    ntop <- sum(ptable[, 2])
    nrest <- sum(ptable[, 1])
    ptop <- ptable[2, 2] / ntop
    prest <- ptable[2, 1] / nrest
    h <- pwr::ES.h(ptop, prest)
    unlist(pwr::pwr.2p2n.test(h = h, n1 = ntop, n2 = nrest, sig.level = 0.05, power = NULL)[c("power", "h")])
}

#' Helper: convert log10P to P
#' @param log10P A log10P value
#' @export
log10toP <- function(log10P) {
    log10P <- abs(as.numeric(log10P))
    if (is.na(log10P)) {
        return(NA)
    }
    if (log10P > 300) {
        part1 <- log10P %/% 100 * 100
        part2 <- log10P - part1
        if (part2 != 0) {
            P <- format(signif(10^-part2, 3), scientific = T)
            P <- paste(as.numeric(gsub("e-.+", "", P)), "e-", as.numeric(gsub(".+-", "", P), sep = "") + part1, sep = "")
        } else {
            P <- paste("1e-", part1, sep = "")
        }
    } else {
        P <- signif(10^-log10P, 3)
    }
    return(as.character(P))
}

#' Helper: get value
#' @param fitAsso A logistf object; used within `getTopEffects`
#' @param exposure The name of the exposure column
#' @export
getValues <- function(fitAsso, exposure) {
    sebeta <- sqrt(diag(vcov(fitAsso)))
    names(sebeta) <- fitAsso$terms
    BETA <- round(fitAsso$coefficient[exposure], 4)
    SEBETA <- round(sebeta[exposure], 4)

    vars <- diag(fitAsso$var)
    names(vars) <- fitAsso$terms
    LOG10P <- -pchisq((fitAsso$coefficient[exposure]^2 / vars[exposure]), 1, lower.tail = F, log = T) / log(10)
    P <- log10toP(LOG10P)

    OR <- round(exp(fitAsso$coefficient[exposure]), 4)
    CI1 <- round(exp(fitAsso$ci.lower[exposure]), 4)
    CI2 <- round(exp(fitAsso$ci.upper[exposure]), 4)
    out <- c(BETA, SEBETA, P, OR, CI1, CI2, signif(LOG10P, 4),
        paste0(
            format(round(OR, 2), nsmall = 2, big.mark = ","),
            " (",
            format(round(CI1, 2), nsmall = 2, big.mark = ","),
            ", ",
            format(round(CI2, 2), nsmall = 2, big.mark = ","),
            ")")
        )
    names(out) <- paste0(exposure, "_", c("beta", "se_beta", "p_value", "or", "lower", "upper", "log10p", "print"))
    return(out)
}

#' Get suite of diagnostics for a risk score and a binary outcome
#' @param data A data frame
#' @param outcome The name of the outcome column
#' @param exposure The name of the exposure column
#' @param covs A vector of covariates
#' @param beta_mod The type of model to use for beta (default: "logistf")
#' @param pctile_or Whether to calculate percentile-based ORs (default: TRUE)
#' @param pl Whether to use profile likelihood (default: FALSE)
#' @param middle A vector of two probabilities representing a percentile range (default: c(0.4, 0.6)) for `getTopMidEffects` OR
#' @importFrom pROC roc
#' @importFrom ResourceSelection hoslem.test
#' @importFrom DescTools BrierScore
#' @importFrom logistf logistf
#' @return A data table with the following columns:
#' \itemize{
#'  \item \code{stat}: The name of the statistic
#'  \item \code{value}: The value of the statistic
#' }
#' @export
get_bin_diagnostics <- function(data, outcome, exposure, covs = NULL, beta_mod = "logistf", pctile_or = TRUE, pl = FALSE, middle = c(0.4, 0.6)) {
    out <- list()
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    data <- data[!is.na(data[[outcome]]) & !is.na(data[[exposure]]), ]
    # formulas
        # without covs
        f_nocovs <- as.formula(paste0(outcome, " ~ ", exposure))
        # with covs
        if (!is.null(covs)) {
            f_covs <- paste0(outcome, " ~ ", exposure, " + ", paste(covs, collapse = " + "))
            f_covs <- as.formula(f_covs)
        } else {
            f_covs <- f_nocovs
        }

    # OR (glm)
    if (beta_mod == "glm") {
        glm_mod <- glm(f_covs, data = data, family = binomial())
        out$beta    <- glm_mod$coefficient[exposure]
        out$se_beta <- sqrt(diag(vcov(glm_mod))[exposure])
        out$or      <- exp(out$beta)
        glm_mod_confint <- suppressMessages(confint(glm_mod))
        out$or_lower <- exp(glm_mod_confint[exposure, 1])
        out$or_upper <- exp(glm_mod_confint[exposure, 2])
        out$or_print <- paste0(
            format(round(out$or, 2), nsmall = 2),
            " (",
            format(round(out$or_lower, 2), nsmall = 2),
            ", ",
            format(round(out$or_upper, 2), nsmall = 2),
            ")"
        )
        out$log10p <- tryCatch(
            {
                -log10(summary(glm_mod)$coefficients[exposure, 4])
            },
            error = function(e) {
                NA
            }
        )
    }

    # OR (logistf)
    if (beta_mod == "logistf") {
        logistf_mod <- logistf::logistf(f_covs, data = data)
        out$beta <- logistf_mod$coefficients[[exposure]]
        out$se_beta <- sqrt(diag(vcov(logistf_mod)))[[exposure]]
        out$or <- exp(out$beta)
        logistf_logistf_mod_confint <- suppressMessages({
            confint(logistf_mod)
        })
        out$or_lower <- exp(logistf_logistf_mod_confint[exposure, 1])
        out$or_upper <- exp(logistf_logistf_mod_confint[exposure, 2])
        out$or_print <- paste0(
            format(round(out$or, 2), nsmall = 2),
            " (",
            format(round(out$or_lower, 2), nsmall = 2),
            ", ",
            format(round(out$or_upper, 2), nsmall = 2),
            ")"
        )
        out$log10p <- tryCatch(
            {
                -log10(logistf_mod$prob[[exposure]])
            },
            error = function(e) {
                NA
            }
        )
    }


    # AUC (pROC::roc)
    tmp_auc <- suppressMessages(pROC::roc(
            response  = data[[outcome]],
            predictor = data[[exposure]],
            family    = binomial(),
            ci        = TRUE
    ))
    out$auc       <- tmp_auc$auc[1]
    out$auc_lower <- tmp_auc$ci[1]
    out$auc_upper <- tmp_auc$ci[3]
    out$auc_print <- paste0(
        format(round(out$auc, 3), nsmall = 3),
        " (",
        format(round(out$auc_lower, 3), nsmall = 3),
        ", ",
        format(round(out$auc_upper, 3), nsmall = 3),
        ")"
    )

    # R2
    nagel_out <- aimTwo::nagelkerke_r2(f_nocovs, data = data, restrictNobs = FALSE)
    out[["R2 (McFadden)"]]                  <- nagel_out$mcfadden_r2
    out$"R2 (Cox and Snell [ML])"           <- nagel_out$cox_snell_r2
    out$"R2 (Nagelkerke [Cragg and Uhler])" <- nagel_out$nagelkerke_r2

    # HL test (ResourceSelection::hoslem.test)
    glm_no_cov_mod <- glm(f_nocovs, data = data, family = binomial())
    hl_out <- ResourceSelection::hoslem.test(glm_no_cov_mod$y, fitted(glm_no_cov_mod), g = 10)
    out$HosmerLemeshow_ChiSq <- hl_out$statistic[1]
    out$HosmerLemeshow_P     <- hl_out$p.value
    out$hl_print <- paste0(
        format(round(out$HosmerLemeshow_ChiSq, 2), nsmall = 2),
        " (",
        formatC(out$HosmerLemeshow_P, format = "e", digits = 2),
        ")"
    )
    # Brier score (DescTools)
    out$Brier_score <- DescTools::BrierScore(glm_no_cov_mod)
    out$bs_print <- format(round(out$Brier_score, 5), nsmall = 5)

    # %ile-based OR
    if (pctile_or) {
        out <- c(out, aimTwo::getTopMidEffects(prob = 0.1, middle = middle, data = data, outcome = outcome, exposure = exposure, covs = covs, pl = pl))
        probs <- c(0.01, 0.02, 0.05, 0.1, 0.25)
        for (pr in probs) {
            out <- c(out, aimTwo::getTopEffects(prob = pr, data = data, outcome = outcome, exposure = exposure, covs = covs, pl = pl))
        }
        
        # determine percentile with at least 80% power (or the most powered percentiles)
        ptest <- sapply(seq(0.005, 0.5, by = 0.005), \(x) aimTwo::getTopPower(x, data = data, outcome = outcome, exposure = exposure))
        presults <- data.table("Top" = seq(0.005, 0.5, by = 0.005), t(ptest))
        #  presults <- data.table('Top' = seq(0.005, 0.5, by=0.005), t(ptest))
        maxp <- presults$Top[which(presults$h == max(presults$h, na.rm = TRUE))[1]]

        p80 <- presults$Top[which(presults$power >= 0.80)]
        underpowered <- F
        if (length(p80) == 0) {
            p80 <- maxp
            underpowered <- T
        } else {
            p80 <- p80[1]
        }

        if (is.na(p80)) {
            prs.p80 <- NA
        } else {
            prs.p80 <- aimTwo::getTopEffects(prob = p80, data = data, outcome = outcome, exposure = exposure, covs = covs, pl = pl)
            names(prs.p80) <- gsub(p80, "MinPower80", names(prs.p80))
        }

        out <- c(out, "MinPower80" = p80, prs.p80, "Top_Underpowered" = underpowered)
    } 

    data.table(
        stat  = names(out),
        value = unlist(out)
    )
}
