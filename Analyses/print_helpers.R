print_result <- function(label, value) {
    cat("\n", label, "\n", sep = "")
    print(value)
    invisible(value)
}

print_model_estimates <- function(label, term, estimate, se) {
    result <- data.frame(
        Term = term,
        HR = exp(estimate),
        CI_low = exp(estimate + qnorm(0.025) * se),
        CI_high = exp(estimate + qnorm(0.975) * se),
        p_value = pnorm(-abs(estimate / se)) * 2,
        check.names = FALSE
    )
    cat("\n", label, "\n", sep = "")
    print(result, digits = 4, row.names = FALSE)
    invisible(result)
}
