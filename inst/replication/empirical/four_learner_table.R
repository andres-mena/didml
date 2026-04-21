#!/usr/bin/env Rscript
j_rows <- lapply(c("J01","J02","J03","J04"), function(sid)
  read.csv(sprintf("explorations/empirical_ext/summary/%s_row.csv", sid), stringsAsFactors = FALSE))
k_rows <- lapply(c("K01","K02","K03","K04"), function(sid)
  read.csv(sprintf("explorations/empirical_ext/summary/%s_row.csv", sid), stringsAsFactors = FALSE))
jdf <- do.call(rbind, j_rows); kdf <- do.call(rbind, k_rows)
lbl <- c(lasso = "Lasso", ridge = "Ridge",
         elastic = "Elastic-net ($\\alpha = 0.5$)", rf = "Random forest")
fmt <- function(x) sprintf("%.3f", x); fmt_se <- function(x) sprintf("(%.3f)", x)
lines <- c("\\begin{tabular}{lccccc}",
  "\\toprule",
  "& \\multicolumn{2}{c}{Primary completion ($D \\ge 6$)} & \\multicolumn{2}{c}{High-school completion ($D \\ge 12$)} & \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}",
  " & DML-Wald & DML-TC & DML-Wald & DML-TC & \\\\",
  "Learner & ($\\hat\\alpha$-rule) & ($\\hat\\alpha$-rule) & ($\\hat\\alpha$-rule) & ($\\hat\\alpha$-rule) & \\\\",
  "\\midrule")
for (i in seq_len(nrow(jdf))) {
  line <- sprintf("%s & %s %s & %s %s & %s %s & %s %s & \\\\",
    lbl[[jdf$ml[i]]],
    fmt(jdf$drwald_hat[i]), fmt_se(jdf$drwald_hat_se[i]),
    fmt(jdf$drtc_hat[i]),   fmt_se(jdf$drtc_hat_se[i]),
    fmt(kdf$drwald_hat[i]), fmt_se(kdf$drwald_hat_se[i]),
    fmt(kdf$drtc_hat[i]),   fmt_se(kdf$drtc_hat_se[i]))
  lines <- c(lines, line)
}
lines <- c(lines, "\\bottomrule", "\\end{tabular}")
writeLines(lines, "output/tables/empirical_4learner.tex")
message("Wrote output/tables/empirical_4learner.tex")
