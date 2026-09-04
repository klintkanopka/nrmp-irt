## This script produces the output for the 2PLM example in the NRMP primer,
## using a pre-selected subset of the open Bohn et al. PREVIC dataset
## selected explicitly to serve as a teaching example.
library(mirt)
library(tidyverse)
library(irw)

set.seed(20260803)

setwd(here::here())
out_dir   <- "output_2pl"
if (!dir.exists(out_dir)) dir.create(out_dir)

## 1. Load item-response matrix ####
## Note that this requires Redivis authentication. As a backup
## data can be downloaded manually from https://itemresponsewarehouse.org/
dat <- irw_fetch('previc_bohn2023')

# pull out the 12 selected words
numword <- 12
wide <- dat %>%
  filter(
    item %in% c(
      "Filmstar",
      "Produkt",
      "Versager",
      "Tourist",
      "Bereich",
      "Drama",
      "Erholung",
      "Pflicht",
      "Kunde",
      "Ahorn",
      "Scherz",
      "Ampulle"
    )
  ) %>%
  select(id, item, resp) %>%
  pivot_wider(names_from = item, values_from = resp)

resp_mat <- wide %>%
  select(-id) %>%
  as.data.frame()

rownames(resp_mat) <- wide$id

mod2pl <- mirt(resp_mat, model = 1, itemtype = "2PL", SE = T)

## Model fit: M2 (limited-information overall fit) & item fit (S-X2) ####

m2_fit   <- M2(mod2pl)
item_fit <- itemfit(mod2pl)

print(m2_fit)
print(item_fit)

m2_fit %>% as_tibble() %>% write_csv(file.path(out_dir, "irt_2pl_model_fit.csv"))
item_fit %>% write_csv(file.path(out_dir, "irt_2pl_item_fit.csv"))

## Item parameter estimates ####

irt_params <- coef(mod2pl, IRTpars = TRUE, simplify = TRUE)$items %>%
  as_tibble(rownames = "word")
print(irt_params)

## save item parameters *with standard errors* ####
format_mirt_coefs <- function(mod, digits = 2, IRTpars = FALSE, ...) {
  cfs <- coef(mod, printSE = TRUE, IRTpars = IRTpars, ...)
  item_names <- extract.mirt(mod, 'itemnames')
  cfs_items  <- cfs[names(cfs) %in% item_names]
  
  if (length(cfs_items) == 0) {
    cfs_items <- cfs[!names(cfs) %in% c("GroupPars", "group")]
  }
  
  rows <- lapply(names(cfs_items), function(item_name) {
    mat <- cfs_items[[item_name]]
    
    estimates <- mat["par", ]
    ses       <- mat["SE", ]
    
    formatted <- sapply(seq_along(estimates), function(i) {
      est <- estimates[i]
      se  <- ses[i]
      
      if (is.na(se) || is.nan(se)) {
        sprintf(paste0("%.", digits, "f"), est)
      } else {
        sprintf(paste0("%.", digits, "f (%.", digits, "f)"), est, se)
      }
    })
    
    df_row <- as.data.frame(as.list(formatted), stringsAsFactors = FALSE)
    colnames(df_row) <- colnames(mat)
    df_row$Item <- item_name
    return(df_row)
  })
  
  if (requireNamespace("dplyr", quietly = TRUE)) {
    res <- dplyr::bind_rows(rows)
    res <- dplyr::select(res, Item, dplyr::everything())
  } else {
    all_cols <- unique(unlist(lapply(rows, colnames)))
    all_cols <- c("Item", setdiff(all_cols, "Item"))
    
    aligned <- lapply(rows, function(df) {
      missing_cols <- setdiff(all_cols, colnames(df))
      for (col in missing_cols) df[[col]] <- NA
      df[, all_cols, drop = FALSE]
    })
    res <- do.call(rbind, aligned)
  }
  
  return(res)
}

tbl <- format_mirt_coefs(mod2pl, IRTpars = TRUE, digits = 3) %>% arrange(Item)
tbl %>% write_csv(file.path(out_dir,'irt_pars_table.csv'))

## ICC plot with overlaid bins ####
get_pval <- function(itemnum, theta){
  itemrow <- irt_params[itemnum,]
  return(
    exp(itemrow$a*(theta-itemrow$b))/(1+exp(itemrow$a*(theta-itemrow$b)))
  )
}

iccs <- data.frame(
  name = irt_params$word,
  theta = rep(seq(-3, 3, length.out = 1000), each = numword),
  item = rep(1:numword, times = 1000)
) %>% 
  mutate(p = get_pval(item, theta))

resp_long <- bind_cols(
  wide[,1],
  resp_mat
) %>% 
  pivot_longer(cols = -id, names_to = 'name', values_to = 'resp')

# note that EAP scores are being computed here; also used below
scores <- fscores(mod2pl, full.scores.SE = T) 

pdat_raw <- data.frame(scores) %>% 
  mutate(id = wide$id) %>% 
  mutate(bin = cut(F1, 12)) %>% 
  group_by(bin) %>% 
  mutate(
    mean = mean(F1),
    n = n()
  ) %>% 
  slice(rep(1:n(), each = numword)) %>% 
  ungroup() %>% 
  mutate(
    itemnum = rep(1:numword, times = nrow(resp_mat)),
    name = rep(irt_params$word, times = nrow(resp_mat))
  ) %>% 
  left_join(
    resp_long
  )

pdat_binned <- pdat_raw %>% 
  group_by(bin, mean, name, n) %>% 
  summarise(p = mean(resp)) %>% 
  rename(theta = mean)

ggplot(iccs, aes(x = theta, y = p, color = name)) +
  facet_wrap(.~name) +
  geom_point(size = 0.2) +
  geom_point(data = pdat_binned, color = 'black', shape = 0) +
  scale_color_discrete(name = "Word") +
  theme(legend.position = 'none') +
  labs(x = "θ", y = "Probability of Affirmative Response", title = "Item Characteric Curves + Empirical Probabilities") +
  theme_minimal()

## information-SE plot ####
tseq <- seq(-4, 4, length.out=1000)
info <- testinfo(mod2pl, Theta = tseq)
infodf <- data.frame(
  info = info,
  theta = tseq,
  Item = "Test Information"
)

sedf <- infodf %>% mutate(
  info = 1/sqrt(info), # just to make the fig work; this is actually the SE
  Item = "Conditional SEM"
)

for(i in 1:numword){
  extr <- extract.item(mod2pl, i)
  inf <- iteminfo(extr, tseq)
  infodf <- rbind(
    infodf, 
    data.frame(
      info = inf,
      theta = tseq,
      Item = irt_params[i, 1]$word
    )
  )
}

infodf <- infodf %>% bind_rows(sedf)
infodf$Item <- factor(infodf$Item, levels = c(sort(irt_params$word)[c(1,4:13)], "Test Information", "Ampulle", "Bereich", "Conditional SEM"))
infodf$color = factor(
    case_match(
    infodf$Item,
    "Test Information" ~ 'Test Information',
    "Ampulle" ~ 'Ampulle Item Info',
    "Bereich" ~ "Bereich Item Info",
    "Conditional SEM" ~ "Cond. SEM",
    .default = "Other Items' Item Info"
  ),
  levels = c("Other Items' Item Info", 'Test Information', 'Ampulle Item Info', "Bereich Item Info", "Cond. SEM")
)

ggplot(infodf, aes(x = theta, y = info, color = color, group = Item)) +
  geom_line(aes(linetype = color), linewidth = 1) +
  labs(x = "θ", y = "Fisher Information") +
  # guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(name = "Key", values = c("gray", "#E69F00", "red", "black", "#56B4E9")) +
  scale_linetype_manual(name = "Key", values = c(6, 1, 2, 4, 5)) +
  scale_y_continuous(sec.axis = dup_axis(name = "Conditional Standard Error of Measurement")) +
  theme_minimal()

## histogram of theta estimate distribution
thetadist <- data.frame(scores)
DescTools::Desc(thetadist$F1)  

ggplot(thetadist, aes(x = F1, y = after_stat(density))) +
  geom_histogram(fill = '#56B4E9', color = 'black', bins = 12) +
  geom_density() +
  labs(x = bquote(hat(theta)), y = "Density") +
  theme_minimal()

## person fit ####
pfit <- personfit(x = mod2pl, method = "EAP")
zh <- pfit$Zh
DescTools::Desc(zh)
min(zh)
worst_fit <- which(pfit$Zh==min(pfit$Zh))


resps <- resp_mat[worst_fit,] %>% pivot_longer(names_to = 'word', values_to = 'resp', cols = names(resp_mat))
compare_resps_pars <- left_join(resps, irt_params)

## DIF: Drama and Erholung ####
tseq <- seq(-4, 4, length.out=1000)

wide_gend <- dat %>%
  filter(
    item %in% c(
      "Filmstar",
      "Produkt",
      "Versager",
      "Tourist",
      "Bereich",
      "Drama",
      "Erholung",
      "Pflicht",
      "Kunde",
      "Ahorn",
      "Scherz",
      "Ampulle"
    )
  ) %>%
  select(id, cov_sex, item, resp) %>%
  pivot_wider(names_from = item, values_from = resp)

mg2pl <- multipleGroup(
  data = wide_gend[,3:14], 
  group = wide_gend$cov_sex, 
  model = 1, 
  itemtype = "2PL", 
  SE = TRUE, 
  invariance = c('free_mean', 'free_var', 'intercepts', 'slopes')
)

coef(mg2pl)

difres <- mirt::DIF(mg2pl, scheme = 'drop', which.par = c('d'))
extract.mirt(difres, "DIF_coefficients")

## plot for Drama as sole example ####
difplotdat <- data.frame(
  theta = tseq,
  a = 1.580766,
  Female = -0.1538278,
  Male = -1.0318555
) %>%
  pivot_longer(
    names_to = 'Sex', values_to = 'd', cols = c('Male', 'Female')
  ) %>%
  mutate(
    p = plogis(a*theta + d)
  )

ggplot(difplotdat, aes(x = theta, y = p, group = Sex, color = Sex)) +
  theme_minimal() + 
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c('green', 'purple4')) +
  labs(x = "θ", y = "Probability of Affirmative Response", title = "Male-Female DIF for Word \"Drama\"")
