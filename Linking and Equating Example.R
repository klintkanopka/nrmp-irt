# generates figure for linking and equating example
library(mirt)
library(tidyverse)
library(irw)
library(equateIRT)

setwd(here::here())
out_dir   <- "output_2pl"
if (!dir.exists(out_dir)) dir.create(out_dir)

## 1. Load item-response matrix ####
## Note that this requires Redivis authentication. As a backup
## data can be downloaded manually from https://itemresponsewarehouse.org/
dat <- irw_fetch('previc_bohn2023')

set.seed(20260803)

dat_filt <- dat |> 
  group_by(item) |> 
  mutate(p = mean(resp)) |> 
  ungroup() |> 
  filter(p < 0.8 & p > 0.2)

dat_filt_form <- dat_filt |> 
  filter(item %in% sample(unique(dat_filt$item), size = 100))

# 100 total items
length(unique(dat_filt$item))

# high vs low ability groups, interchangeable forms
dat_filt_grp <- dat_filt_form %>% 
  group_by(id) |> 
  mutate(
    ss = sum(resp),
    grp = ifelse(
      ss < 60,
      sample(c('a', 'b'), 1, prob = c(.8, .2)),
      sample(c('a', 'b'), 1, prob = c(.2, .8))
    )
  )|> 
  ungroup()

ggplot(dat_filt_grp |> select(id, grp, ss) |> distinct(), aes(x = ss, y= after_stat(density), fill = grp)) + 
  geom_histogram(position = 'dodge') + geom_density(alpha = 0.2)

all_items <- unique(dat_filt_grp$item)          # 100 items

## split into anchors + two form-specific sets
anchor_items <- sample(all_items, 40)                       # common to both forms
remaining    <- setdiff(all_items, anchor_items)
a_only_items <- sample(remaining, 30)                       # Form A only
b_only_items <- setdiff(remaining, a_only_items)            # Form B only (30)

form_a_items <- c(anchor_items, a_only_items)               # 70-item Form A
form_b_items <- c(anchor_items, b_only_items)               # 70-item Form B

group_a_dat <- dat_filt_grp |>
  filter(grp == 'a', item %in% form_a_items) |>
  pivot_wider(names_from = item, values_from = resp, id_cols = c(id, grp))

group_b_dat <- dat_filt_grp |>
  filter(grp == 'b', item %in% form_b_items) |>
  pivot_wider(names_from = item, values_from = resp, id_cols = c(id, grp))

fit_a <- mirt(group_a_dat |> dplyr::select(-id, -grp), 1, itemtype = "2PL", SE = TRUE)
fit_b <- mirt(group_b_dat |> dplyr::select(-id, -grp), 1, itemtype = "2PL", SE = TRUE)

theta <- seq(-4, 4, length.out = 500)
p2pl  <- function(a, b, th) 1 / (1 + exp(-a * (th - b)))
tcc   <- function(a, b, th) rowSums(mapply(function(aj, bj) p2pl(aj, bj, th), a, b))

ipar_a <- coef(fit_a, simplify = TRUE, IRTpars = TRUE)$items[, c("a", "b")]
ipar_b <- coef(fit_b, simplify = TRUE, IRTpars = TRUE)$items[, c("a", "b")]

tcc_a_pre <- tcc(ipar_a[, "a"], ipar_a[, "b"], theta)
tcc_b_pre <- tcc(ipar_b[, "a"], ipar_b[, "b"], theta)

est_a <- import.mirt(fit_a, display = FALSE)
est_b <- import.mirt(fit_b, display = FALSE)

mods <- modIRT(
  coef  = list(A = est_a$coef, B = est_b$coef),
  var   = list(A = est_a$var,  B = est_b$var),
  names = c("A", "B"), display = FALSE
)

sl <- direc(mods = mods, which = c("B", "A"), method = "Stocking-Lord")
summary(sl)                       # A (slope), B (intercept)

A_const <- sl$A                   # theta_A = A * theta_B + B
B_const <- sl$B

a_b_star  <- ipar_b[, "a"] / A_const
b_b_star  <- A_const * ipar_b[, "b"] + B_const
tcc_b_post <- tcc(a_b_star, b_b_star, theta)

plot_df <- rbind(
  data.frame(Theta = theta, ExpectedScore = tcc_a_pre,  Form = "Form A", Stage = "Before linking"),
  data.frame(Theta = theta, ExpectedScore = tcc_b_pre,  Form = "Form B", Stage = "Before linking"),
  data.frame(Theta = theta, ExpectedScore = tcc_a_pre,  Form = "Form A", Stage = "After linking"),
  data.frame(Theta = theta, ExpectedScore = tcc_b_post, Form = "Form B", Stage = "After linking")
)
plot_df$Stage <- factor(plot_df$Stage, levels = c("Before linking", "After linking"))

ggplot(plot_df, aes(Theta, ExpectedScore, color = Form, linetype = Form)) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~ Stage) +
  labs(
    title    = "Expected Total Score by Ability Level Before and After Linking",
    x = expression(theta), y = "Expected Total Score"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue4", "orange4")) +
  scale_linetype_manual(values = c(1, 8))
