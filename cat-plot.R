library(meow)
library(mirt)
library(tidyverse)
library(irw)
library(see)


setwd(here::here())
out_dir <- "output_cat"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

d_raw <- irw_fetch('previc_bohn2023')

set.seed(20260903)

good_people <- d_raw |>
  group_by(id) |>
  summarize(p = mean(resp)) |>
  filter(p < 1) |>
  filter(p > 0) |>
  pull(id)

good_items <- d_raw |>
  group_by(item) |>
  summarize(p = mean(resp)) |>
  filter(p < 0.95) |>
  pull(item)

length(good_people)
length(good_items)

id_map <- data.frame(id = good_people, new_id = 1:length(good_people))
item_map <- data.frame(item = good_items, new_item = 1:length(good_items))

d <- d_raw |>
  filter(id %in% good_people) |>
  filter(item %in% good_items) |>
  left_join(id_map, by = 'id') |>
  left_join(item_map, by = 'item') |>
  select(id = new_id, item = new_item, resp)

length(unique(d$id))
# 1186 total, select 100?

oos_id <- 1:100

resp <- d |>
  filter(!id %in% oos_id) |>
  select(id, item, resp) |>
  pivot_wider(id_cols = id, names_from = item, values_from = resp) |>
  arrange(id) |>
  select(-id)

m <- mirt(resp, model = 1, itemtype = "2PL")


oos_resp <- d |>
  filter(id %in% oos_id) |>
  select(id, item, resp) |>
  arrange(id, item)

write_csv(oos_resp, file.path(out_dir, 'oos_resp.csv'))

oos_resp_mat <- oos_resp |>
  pivot_wider(id_cols = id, names_from = item, values_from = resp) |>
  arrange(id) |>
  select(-id)


oos_theta <- fscores(
  m,
  method = 'ML',
  response.pattern = oos_resp_mat
)

oos_pers <- data.frame(id = oos_id, theta = oos_theta[, 1])


write_csv(oos_pers, file.path(out_dir, 'oos_pers.csv'))

oos_item <- coef(m, simplify = TRUE, IRTpars = TRUE)$items |>
  as.data.frame() |>
  rownames_to_column('item') |>
  select(item, a, b)

write_csv(oos_item, file.path(out_dir, 'oos_item.csv'))


cat_sim <- meow(
  select_fun = select_max_info,
  update_fun = update_theta_mle,
  data_loader = data_existing,
  data_args = list(
    resp_path = file.path(out_dir, 'oos_resp.csv'),
    pers_path = file.path(out_dir, 'oos_pers.csv'),
    item_path = file.path(out_dir, 'oos_item.csv')
  ),
  fix = 'item'
)

results_cat <- cat_sim$results |>
  mutate(control = 'CAT')


random_sim <- meow(
  select_fun = select_random,
  update_fun = update_theta_mle,
  data_loader = data_existing,
  data_args = list(
    resp_path = file.path(out_dir, 'oos_resp.csv'),
    pers_path = file.path(out_dir, 'oos_pers.csv'),
    item_path = file.path(out_dir, 'oos_item.csv')
  ),
  fix = 'item'
)

results_random <- random_sim$results |>
  mutate(control = 'Random')


seq_sim <- meow(
  select_fun = select_sequential,
  update_fun = update_theta_mle,
  data_loader = data_existing,
  data_args = list(
    resp_path = file.path(out_dir, 'oos_resp.csv'),
    pers_path = file.path(out_dir, 'oos_pers.csv'),
    item_path = file.path(out_dir, 'oos_item.csv')
  ),
  fix = 'item'
)

results_seq <- seq_sim$results |>
  mutate(control = 'Alphabetical')


results <- bind_rows(results_cat, results_seq)

results |>
  select(iter, control, starts_with("pers_")) |>
  select(iter, control, ends_with("_bias")) |>
  pivot_longer(ends_with("_bias"), names_to = "person", values_to = "bias") |>
  group_by(iter, control) |>
  summarize(rmse = sqrt(mean(bias^2)), .groups = "drop") |>
  filter(rmse < 0.3) |>
  group_by(control) |>
  slice_min(iter)

results |>
  select(iter, control, starts_with("pers_")) |>
  select(iter, control, ends_with("_bias")) |>
  pivot_longer(ends_with("_bias"), names_to = "person", values_to = "bias") |>
  group_by(iter, control) |>
  summarize(rmse = sqrt(mean(bias^2)), .groups = "drop") |>
  ggplot(aes(x = iter, y = rmse, color = control)) +
  geom_line(linewidth = 1) +
  geom_hline(aes(yintercept = 0.3), linewidth = 1, lty = 2) +
  geom_vline(
    aes(xintercept = 14),
    linewidth = 1,
    lty = 2,
    color = okabeito_colors(2)
  ) +
  geom_vline(
    aes(xintercept = 46),
    linewidth = 1,
    lty = 2,
    color = okabeito_colors(1)
  ) +
  labs(
    x = "Number of Items Delivered",
    y = "Vocabulary Measurement Precision",
    color = "Item Selection Method"
  ) +
  scale_color_okabeito() +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.8)) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, 'cat-5B.png'), height = 6, width = 8)


results |>
  select(iter, control, starts_with("pers_")) |>
  select(iter, control, ends_with("_bias")) |>
  pivot_longer(ends_with("_bias"), names_to = "person", values_to = "bias") |>
  group_by(iter, control) |>
  summarize(rmse = sqrt(mean(bias^2)), .groups = "drop") |>
  filter(rmse < 0.3) |>
  group_by(control) |>
  slice_min(iter)


# now the distributions

n_items <- 14
focal_pers <- 1

info_2pl <- function(theta, a, b) {
  p <- plogis(a * (theta - b))
  info <- a^2 * p * (1 - p)
  return(info)
}

focal_resp <- oos_resp_mat[rep(focal_pers, (n_items - 3)), ]

item_mask <- matrix(
  FALSE,
  ncol = ncol(focal_resp),
  nrow = (n_items - 3)
)

item_mask[11, ] <- TRUE
item_mask[, 1:5] <- TRUE


tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA
fscores(m, method = 'EAP', response.pattern = tmp_resp)


oos_item |>
  mutate(info = info_2pl(0.387, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[2:10, 28] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(0.628, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[3:10, 252] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(0.830, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[4:10, 250] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(0.641, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[5:10, 133] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(0.746, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[6:10, 156] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(0.843, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[7:10, 122] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(0.945, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[8:10, 163] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(1.065, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[9:10, 266] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

fscores(m, method = 'EAP', response.pattern = tmp_resp)

oos_item |>
  mutate(info = info_2pl(0.984, a, b)) |>
  arrange(-info) |>
  head(10)

item_mask[10, 150] <- TRUE

tmp_resp <- focal_resp
tmp_resp[!item_mask] <- NA

theta_cat <- fscores(m, method = 'EAP', response.pattern = tmp_resp)

theta_dist <- data.frame(
  item = c(5:14, ncol(focal_resp)),
  theta = theta_cat[, 1],
  se = theta_cat[, 2]
)


ggplot() +
  xlim(-1.25, 2.25) +
  geom_function(
    fun = dnorm,
    aes(col = '5'),
    args = list(mean = theta_dist$theta[1], sd = theta_dist$se[1]),
    lty = 3,
    linewidth = 1
  ) +
  geom_function(
    fun = dnorm,
    aes(col = '8'),
    args = list(mean = theta_dist$theta[4], sd = theta_dist$se[4]),
    lty = 4,
    linewidth = 1
  ) +
  geom_function(
    fun = dnorm,
    aes(col = '11'),
    args = list(mean = theta_dist$theta[7], sd = theta_dist$se[7]),
    lty = 2,
    linewidth = 1
  ) +
  geom_function(
    fun = dnorm,
    aes(col = '14'),
    linewidth = 1,
    args = list(mean = theta_dist$theta[10], sd = theta_dist$se[10]),
  ) +
  geom_vline(
    aes(xintercept = theta_dist$theta[11]),
    lty = 2,
    linewidth = 1,
    color = 'black'
  ) +
  scale_color_manual(
    values = c(
      '5' = "#56B4E9",
      '8' = "#56B4E9",
      '11' = "#56B4E9",
      '14' = 'black'
    ),
    breaks = c('5', '8', '11', '14')
  ) +
  labs(
    x = expression(theta),
    y = 'Estimate Density',
    color = 'Items Adminstered'
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, 'cat-5A.png'), height = 6, width = 8)
