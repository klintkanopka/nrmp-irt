## This script produces the output for the GRM example in the NRMP primer,
## using a pre-selected subset of the open PROMIS Anger dataset
## selected explicitly to serve as a teaching example.

library(tidyverse)
library(mirt)
library(irw)
library(see)

set.seed(20260901)

setwd(here::here())
out_dir <- "output_grm"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

## 1. Load item-response matrix ####
## Note that this requires Redivis authentication. As a backup
## data can be downloaded manually from https://itemresponsewarehouse.org/
d <- irw_fetch('promis1wave1_anger')

# select the focal items

d_resp <- d |>
  filter(
    item %in%
      c(
        'EDANG02',
        'EDANG06',
        'EDANG07',
        'EDANG14',
        'EDANG15',
        'EDANG24',
        'EDANG26',
        'EDANG39',
        'EDANG46'
      )
  ) |>
  pivot_wider(
    id_cols = 'id',
    names_from = item,
    values_from = resp
  ) |>
  na.omit() |>
  select(-id)

m_grm <- mirt(d_resp, 1, itemtype = 'graded', SE = TRUE)


d |>
  na.omit() |>
  count(item, resp) |>
  group_by(item) |>
  mutate(tot_resp = sum(n)) |>
  mutate(p = n / tot_resp) |>
  filter(p <= 0.05)

summary(m_grm)

m2_grm <- M2(m_grm)
itemfit_grm <- itemfit(m_grm)

print(m2_grm)
print(itemfit_grm)

m2_grm |>
  as_tibble() |>
  write_csv(file.path(out_dir, "irt_grm_model_fit.csv"))

itemfit_grm |>
  write_csv(file.path(out_dir, "irt_grm_item_fit.csv"))

pars <- coef(m_grm, simplify = TRUE, IRTpars = TRUE)$items |>
  as_tibble() |>
  mutate(
    item = colnames(d_resp),
    text = c(
      'I was angry when people were unfair',
      'I made myself angry about something just by thinking about it',
      'I tried to get even when I was angry with someone',
      'It was difficult to let people know I was angry',
      'I felt like I was ready to explode',
      'I got into physical fights',
      'I held grudges towards others',
      'When I lost my temper, I hit someone',
      'I spread gossip about people I did not like'
    )
  ) |>
  select(item, text, a, b1, b2, b3, b4)

print(pars)

pars |>
  write_csv(file.path(out_dir, "irt_grm_pars_table.csv"))


pars |>
  pivot_longer(starts_with('b'), names_to = 'par', values_to = 'est') |>
  ggplot(aes(y = est, x = reorder(item, a), color = par)) +
  geom_point(size = 3) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = 'Higher response categories correspond to higher satisfaction; items ordered by discrimination (low to high)',
    x = NULL
  )

pars |>
  pivot_longer(starts_with('b'), names_to = 'par', values_to = 'est') |>
  ggplot(aes(y = a, x = reorder(item, a))) +
  geom_point(size = 3, color = okabeito_colors(7)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = 'Items ordered by discrimination (low to high)',
    x = NULL
  )

plot(m_grm, type = 'trace')
