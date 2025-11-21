# ------------------------------------------------------------
# Set current script location as working directory + packages
# ------------------------------------------------------------
library(rstudioapi)
a <- getActiveDocumentContext()
setwd(dirname(a$path))
remove(a)

library(dplyr)
library(knitr)
library(kableExtra)

storage <- list()
datasetsDirectory <- "datasets scenario 1/"

# ------------------------------------------------------------
# Helper function
# ------------------------------------------------------------

analyse_dataset <- function(file_path) {
    data <- readRDS(file_path)

    death_rate <- data$subjectData %>%
        summarise(
            death_rate = mean(died),
            death_rate_percent = mean(died) * 100
        )

    mean_recurrence <- data$longData %>%
        group_by(id) %>%
        summarise(nRecurrence = sum(recurrence)) %>%
        summarise(
            meanRecurrence = mean(nRecurrence),
            minRecurrence = min(nRecurrence),
            maxRecurrence = max(nRecurrence)
        )

    percent_no_recEvents <- data$longData %>%
        group_by(id) %>%
        summarise(nRecurrence = sum(recurrence)) %>%
        summarise(
            subjects_no_events = sum(nRecurrence == 0),
            total_subjects = n(),
            percent_no_recEvents = (sum(nRecurrence == 0) / n()) * 100
        )

    percent_no_events <- data$longData %>%
        group_by(id) %>%
        filter(n() == 1, status == 0) %>%
        ungroup() %>%
        summarise(
            subjects_no_events = n(),
            total_subjects = n_distinct(data$longData$id),
            percent_no_events = (n() / n_distinct(data$longData$id)) * 100
        )

    list(
        death_rate = death_rate,
        mean_recurrence = mean_recurrence,
        percent_no_recEvents = percent_no_recEvents,
        percent_no_events = percent_no_events
    )
}


files <- list.files(datasetsDirectory, pattern = "\\.rds$", full.names = TRUE)
for (file in files) {
    result <- analyse_dataset(file)
    storage[[file]] <- result
}

finalRes <- lapply(storage, function(x) {
    list(
        death_rate = x$death_rate,
        mean_recurrence = x$mean_recurrence,
        percent_no_recEvents = x$percent_no_recEvents,
        percent_no_events = x$percent_no_events
    )
})

final_df <- do.call(
    rbind,
    lapply(finalRes, function(x) {
        data.frame(
            death_rate = x$death_rate$death_rate,
            death_rate_percent = x$death_rate$death_rate_percent,
            mean_recurrence = x$mean_recurrence$meanRecurrence,
            min_recurrence = x$mean_recurrence$minRecurrence,
            max_recurrence = x$mean_recurrence$maxRecurrence,
            percent_no_recEvents = x$percent_no_recEvents$percent_no_recEvents,
            percent_no_events = x$percent_no_events$percent_no_events
        )
    })
)
final_df <- as.data.frame(final_df)

summary_stats <- final_df %>%
    summarise(
        across(
            .cols = c(
                death_rate_percent,
                mean_recurrence,
                min_recurrence,
                max_recurrence,
                percent_no_recEvents,
                percent_no_events
            ),
            .fns = list(Mean = mean, SD = sd, Min = min, Max = max),
            .names = "{.col}_{.fn}"
        )
    ) %>%
    tidyr::pivot_longer(
        cols = everything(),
        names_to = c("Variable", ".value"),
        names_pattern = "(.+)_(Mean|SD|Min|Max)"
    ) %>%
    mutate(
        Variable = recode(
            Variable,
            "death_rate_percent" = "Death rate (%)",
            "mean_recurrence" = "Mean recurrence",
            "min_recurrence" = "Min n° recurrence",
            "max_recurrence" = "Max n° recurrence",
            "percent_no_recEvents" = "No rec. event rate (%)",
            "percent_no_events" = "Zero event rate (%)"
        )
    )

kable(
    summary_stats,
    format = "markdown",
    digits = 2,
    caption = "Summary statistics across all datasets"
)
