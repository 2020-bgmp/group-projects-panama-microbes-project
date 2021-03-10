``` r
# read in normalized ASV table and metadata
metadata <- read_tsv("../nodrysamples/panama_metadata_full_v1.txt")
ASV_tab_rar <- read_tsv("../nodrysamples/PanamaPrecip_ITS_ASVs_r10574.txt")

meta_sub <- subset(metadata, select = c("Sample_ID", "Plot", "PlotCreated", "pH.water", 
                                        "ResinP.mg_kg", "MAP", "DON_mg.N_kg"))

# create binary presence/absence matrix from ASV table
ASV_tab_pa <- subset(ASV_tab_rar, select = -c(Kingdom:Species))
ASV_tab_pa[,2:346][ASV_tab_pa[,2:346] >= 1] <- 1

ASV_tab_pa <- rowid_to_column(ASV_tab_pa, var = "ASV_index")
ASV_tab_pa <- rename(ASV_tab_pa, ASV_seq = X)

# transpose -> rows = samples, cols = ASVs
ASV_tab_pa_t <- as_tibble(cbind(Sample_ID = names(ASV_tab_pa)[c(-1,-2)], t(ASV_tab_pa[c(-1,-2)])))

# format sample IDs
colnames(ASV_tab_pa_t) <- c("Sample_ID", paste0("ASV", ASV_tab_pa$ASV_index))
ASV_tab_pa_t$Sample_ID <- gsub("\\.", "_", ASV_tab_pa_t$Sample_ID)

# merge transposed ASV table with metadata
PA <- right_join(meta_sub, ASV_tab_pa_t, by = "Sample_ID")

PA[,8:ncol(PA)] <- sapply(PA[,8:ncol(PA)], as.numeric)

# remove rows with NA 
PA <- drop_na(PA)
```

``` r
# going to omit ASVs w/ < 10 occurrences
PA_mod <- PA[8:ncol(PA)]
PA_mod <- PA_mod[,colSums(PA_mod) > 10]
PA_mod <- cbind(PA[1:7], PA_mod)

# initialize dataframes to hold ASVs w/significant effect sizes
lowP <- data.frame(PA_mod$Sample_ID)
highP <- data.frame(PA_mod$Sample_ID)
unclassified_P <- data.frame(PA_mod$Sample_ID)

lowMAP <- data.frame(PA_mod$Sample_ID)
highMAP <- data.frame(PA_mod$Sample_ID)
unclassified_MAP <- data.frame(PA_mod$Sample_ID)

lowN <- data.frame(PA_mod$Sample_ID)
highN <- data.frame(PA_mod$Sample_ID)
unclassified_N <- data.frame(PA_mod$Sample_ID)

# Resin P
for (asv in colnames(PA_mod)[8:ncol(PA_mod)]){
  form <- as.formula(paste0(asv, "~ log(ResinP.mg_kg + 1) + (1|Plot)"))
  model <- glmer(form, family = binomial, data = PA_mod)
  if (summary(model)[10]$coefficients[8] < .05){ # p-value < 0.05
    if (fixef(model)[[2]] < -.8){ # coefficient/effect size > abs(0.8)
      lowP <- cbind(lowP, PA_mod[asv])
    }
    else if (fixef(model)[[2]] > .8){
      highP <- cbind(highP, PA_mod[asv])
    }
    else{
      unclassified_P <- cbind(unclassified_P, PA_mod[asv])
    }
  }
  else{
    unclassified_P <- cbind(unclassified_P, PA_mod[asv])
  }
}

# MAP
for (asv in colnames(PA_mod)[8:ncol(PA_mod)]){
  form <- as.formula(paste0(asv, "~ log(MAP) + (1|Plot)"))
  model <- glmer(form, family = binomial, data = PA_mod)
  if (summary(model)[10]$coefficients[8] < .05){
    if (fixef(model)[[2]] < -.8){
      lowMAP <- cbind(lowMAP, PA_mod[asv])
    }
    else if (fixef(model)[[2]] > .8){
      highMAP <- cbind(highMAP, PA_mod[asv])
    }
    else{
      unclassified_MAP <- cbind(unclassified_MAP, PA_mod[asv])
    }
  }
  else{
    unclassified_MAP <- cbind(unclassified_MAP, PA_mod[asv])
  }
}

# DON
for (asv in colnames(PA_mod)[8:ncol(PA_mod)]){
  form <- as.formula(paste0(asv, "~ sqrt(DON_mg.N_kg) + (1|Plot)"))
  model <- glmer(form, family = binomial, data = PA_mod)
  if (summary(model)[10]$coefficients[8] < .05){
    if (fixef(model)[[2]] < -.8){
      lowN <- cbind(lowN, PA_mod[asv])
    }
    else if (fixef(model)[[2]] > .8){
      highN <- cbind(highN, PA_mod[asv])
    }
    else{
      unclassified_N <- cbind(unclassified_N, PA_mod[asv])
    }
  }
  else{
    unclassified_N <- cbind(unclassified_N, PA_mod[asv])
  }
}
```

``` r
# make lists of significant ASVs
lowP_ASVs <- colnames(lowP)[2:length(colnames(lowP))]
highP_ASVs <- colnames(highP)[2:length(colnames(highP))]
unclassified_P_ASVs <- colnames(unclassified_P)[2:ncol(unclassified_P)]

lowMAP_ASVs <- colnames(lowMAP)[2:length(colnames(lowMAP))]
highMAP_ASVs <- colnames(highMAP)[2:length(colnames(highMAP))]
unclassified_MAP_ASVs <- colnames(unclassified_MAP)[2:ncol(unclassified_MAP)]

lowN_ASVs <- colnames(lowN)[2:length(colnames(lowN))]
highN_ASVs <- colnames(highN)[2:length(colnames(highN))]
unclassified_N_ASVs <- colnames(unclassified_N)[2:ncol(unclassified_N)]

# extract significant ASVs with counts
ASV_affin <- as.data.frame(rowid_to_column(ASV_tab_rar, var = "ASV"))
ASV_affin[1] <- paste0("ASV", ASV_affin$ASV)
ASV_affin <- subset(ASV_affin, select = -c(X:Species)) # only counts matrix
names(ASV_affin) <- gsub(x = names(ASV_affin), pattern = "\\.", replacement = "_")  
ASV_affin[,2:ncol(ASV_affin)][ASV_affin[,2:ncol(ASV_affin)] >= 1] <- 1 # presence/absence

lowP_counts <- ASV_affin[ASV_affin$ASV %in% lowP_ASVs,]
highP_counts <- ASV_affin[ASV_affin$ASV %in% highP_ASVs,]
unclass_P_counts <- ASV_affin[ASV_affin$ASV %in% unclassified_P_ASVs,]

lowMAP_counts <- ASV_affin[ASV_affin$ASV %in% lowMAP_ASVs,]
highMAP_counts <- ASV_affin[ASV_affin$ASV %in% highMAP_ASVs,]
unclass_MAP_counts <- ASV_affin[ASV_affin$ASV %in% unclassified_MAP_ASVs,]

lowN_counts <- ASV_affin[ASV_affin$ASV %in% lowN_ASVs,]
highN_counts <- ASV_affin[ASV_affin$ASV %in% highN_ASVs,]
unclass_N_counts <- ASV_affin[ASV_affin$ASV %in% unclassified_N_ASVs,]

# get sums of samples for calculating proportions
sum_df <- data.frame(colSums(ASV_tab_pa[3:ncol(ASV_tab_pa)]))
sum_df$Sample_ID <- rownames(sum_df)
sum_df$Sample_ID <- gsub("\\.", "_", sum_df$Sample_ID)
sum_df <- rename(sum_df, total_ASVs = colSums.ASV_tab_pa.3.ncol.ASV_tab_pa...)
```

### P affinity threshold

``` r
############## Resin P ##############
P_affin_by_sample <- data.frame(colnames(ASV_affin)[2:ncol(ASV_affin)],
                                colSums(highP_counts[2:ncol(highP_counts)]),
                                row.names = NULL)
P_affin_by_sample <- cbind(P_affin_by_sample,
                           colSums(lowP_counts[2:ncol(lowP_counts)])                         )
rownames(P_affin_by_sample) <- NULL
colnames(P_affin_by_sample) <- c("Sample_ID", "High", "Low")

P_meta_affin <- subset(meta_sub, select = c("Sample_ID", "ResinP.mg_kg"))

P_affin_by_sample <- left_join(P_affin_by_sample, P_meta_affin, by = "Sample_ID")

# reshape dataframe with pivot
P_affinities <- as_tibble(P_affin_by_sample)
P_affinities <- pivot_longer(P_affinities, High:Low, names_to = "Affinity",
                             values_to = "Count")

#merge with affinities
P_affinities <- left_join(P_affinities, sum_df, by = "Sample_ID")

# P affinity plot
P_plot <- ggplot(P_affinities, aes(x = ResinP.mg_kg, y = Count/total_ASVs, color = Affinity)) +
  geom_point(alpha = .7) +
  geom_smooth(se = FALSE) +
  geom_vline(xintercept = 1.36, color = "black", linetype = "dashed") +
  annotate("text", 
           label = expression(atop("Resin phosphate = ", paste("1.36 mg P ", kg^-1))),
           x = 0.3, y = .33) +
  labs(x = expression(paste("Resin phosphate (mg P ", kg^-1, ")")),
       y = "Proportion of ASVs in sample") +
  theme(legend.position = "top") +
  scale_x_log10() +
  scale_y_sqrt(limits = c(0, 0.4)) +
  #ylim(0, .4) +
  scale_color_manual(values = c("blue", "red"),
                     name = NULL,
                     breaks = c("Low", "High"),
                     labels = c("Low precipitation affinity",
                                "High precipitation affinity"))
P_plot
```

![](Plots/affinity/P%20affinity%20plot-1.png)

Threshold for strong effect of resin phosphate. Low affinity defined as
effect size \< -0.8 (blue points and lines). High affinity defined as
effect size \> 0.8 (red points and lines). The proportion of ASVs with
low resin P affinity and that of ASVs with high resin P affinity is
equal at 0.92 mg/kg. Axes are scaled: x on log scale, y on square root
scale.

### MAP affinity threshold

``` r
############## MAP ##############
MAP_affin_by_sample <- data.frame(colnames(ASV_affin)[2:ncol(ASV_affin)],
                                colSums(highMAP_counts[2:ncol(highMAP_counts)]),
                                row.names = NULL)
MAP_affin_by_sample <- cbind(MAP_affin_by_sample,
                           colSums(lowMAP_counts[2:ncol(lowMAP_counts)])                         )
rownames(MAP_affin_by_sample) <- NULL
colnames(MAP_affin_by_sample) <- c("Sample_ID", "High", "Low")

MAP_meta_affin <- subset(meta_sub, select = c("Sample_ID", "MAP"))

MAP_affin_by_sample <- left_join(MAP_affin_by_sample, MAP_meta_affin, by = "Sample_ID")

# reshape dataframe with pivot
MAP_affinities <- as_tibble(MAP_affin_by_sample)
MAP_affinities <- pivot_longer(MAP_affinities, High:Low, names_to = "Affinity",
                             values_to = "Count")

# merge with total ASVs present per sample
MAP_affinities <- left_join(MAP_affinities, sum_df, by = "Sample_ID")

# affinity plot for MAP
MAP_plot <- ggplot(MAP_affinities, 
                   aes(x = MAP, y = Count/total_ASVs, color = Affinity)) +
  geom_jitter(width = .8) +
  geom_smooth(se = FALSE) + ylim(0,.45) +
  geom_vline(xintercept = 2567, color = "black", linetype = "dashed") +
  annotate("text", label = "Mean annual precipitation =\n2567 cm",
           x = 2000, y = .25) +
  labs(x = "Mean annual precipitation (cm)",
       y = "Proportion of ASVs in sample") +
  theme(legend.position = "top") +
  scale_y_sqrt() +
  scale_color_manual(values = c("blue", "red"),
                  name = NULL,
                  breaks = c("Low", "High"),
                  labels = c("Low precipitation affinity",
                             "High precipitation affinity"))
MAP_plot
```

![](Plots/affinity/MAP%20affinity%20plot-1.png)

Threshold for strong effect of mean annual precipitation (MAP). Low
affinity defined as effect size \< -0.8 (blue points and lines). High
affinity defined as effect size \> 0.8 (red points and lines). The
proportion of ASVs with low MAP affinity and that of ASVs with high MAP
affinity is equal at 2567 cm. Y axis is square root scaled.

### DON affinity threshold

``` r
############## DON ##############
N_affin_by_sample <- data.frame(colnames(ASV_affin)[2:ncol(ASV_affin)],
                                colSums(highN_counts[2:ncol(highN_counts)]),
                                row.names = NULL)
N_affin_by_sample <- cbind(N_affin_by_sample,
                           colSums(lowN_counts[2:ncol(lowN_counts)])                         )
rownames(N_affin_by_sample) <- NULL
colnames(N_affin_by_sample) <- c("Sample_ID", "High", "Low")

N_meta_affin <- subset(meta_sub, select = c("Sample_ID", "DON_mg.N_kg"))

N_affin_by_sample <- left_join(N_affin_by_sample, N_meta_affin, by = "Sample_ID")

# reshape dataframe with pivot
N_affinities <- as_tibble(N_affin_by_sample)
N_affinities <- pivot_longer(N_affinities, High:Low, names_to = "Affinity",
                             values_to = "Count")

#merge with affinities
N_affinities <- left_join(N_affinities, sum_df, by = "Sample_ID")

# N affinity plot
N_plot <- ggplot(N_affinities,
          aes(x = DON_mg.N_kg, y = Count/total_ASVs, color = Affinity)) +
  geom_point(alpha = .7) +
  geom_smooth(se = FALSE) +
  ylim(0, 0.9) +
  geom_vline(xintercept = 14.35, color = "black", linetype = "dashed") +
  annotate("text", label = expression(paste("DON = 14.35 mg DON ", kg^-1)),
          x = 6, y = .45) +
  labs(x = expression(paste("Dissolved organic nitrogen (mg P ", kg^-1, ")")),
       y = "Proportion of ASVs in sample") +
  theme(legend.position = "top") +
  scale_x_sqrt() +
  scale_y_sqrt() +
  scale_color_manual(values = c("blue", "red"),
                     name = NULL,
                     breaks = c("Low", "High"),
                     labels = c("Low nitrogen affinity",
                                "High nitrogen affinity"))
N_plot
```

![](Plots/affinity/DON%20affinity%20plot-1.png)

Threshold for strong effect of dissolved organic nitrogen (DON). Low
affinity defined as effect size \< -0.8 (blue points and lines). High
affinity defined as effect size \> 0.8 (red points and lines). The
proportion of ASVs with low DON affinity and that of ASVs with high DON
affinity is equal at 15.76 mg/kg. Axes are square root scaled.

### Phosphorus (Resin P) affinity group composition

``` r
ASV_taxa <- ASV_tab_rar[1:8]

ASV_taxa <- rowid_to_column(ASV_taxa, var = "ASV_index")
ASV_taxa <- rename(ASV_taxa, ASV_seq = X)
ASV_taxa$ASV_index <- c(paste0("ASV", ASV_taxa$ASV_index))

highP_taxa <- as.data.frame(names(highP[2:ncol(highP)]))
highP_taxa <- rename(highP_taxa, ASV_index = `names(highP[2:ncol(highP)])`)
highP_taxa <- left_join(highP_taxa, ASV_taxa, by = "ASV_index")

lowP_taxa <- as.data.frame(names(lowP[2:ncol(lowP)]))
lowP_taxa <- rename(lowP_taxa, ASV_index = `names(lowP[2:ncol(lowP)])`)
lowP_taxa <- left_join(lowP_taxa, ASV_taxa, by = "ASV_index")

unclassP_taxa <- as.data.frame(names(unclassified_P[2:ncol(unclassified_P)]))
unclassP_taxa <- rename(unclassP_taxa, ASV_index = `names(unclassified_P[2:ncol(unclassified_P)])`)
unclassP_taxa <- left_join(unclassP_taxa, ASV_taxa, by = "ASV_index")

# Phylum
highP_phylum <- as.data.frame(table(highP_taxa$Phylum))
colnames(highP_phylum) <- c("Phylum", "Count")
highP_phylum$Affinity <- "High"

lowP_phylum <- as.data.frame(table(lowP_taxa$Phylum))
colnames(lowP_phylum) <- c("Phylum", "Count")
lowP_phylum$Affinity <- "Low"

unclassP_phylum <- as.data.frame(table(unclassP_taxa$Phylum))
colnames(unclassP_phylum) <- c("Phylum", "Count")
unclassP_phylum$Affinity <- "Unclassified"

P_phylum <- rbind(highP_phylum, lowP_phylum, unclassP_phylum)

ggplot(P_phylum, aes(fill = Phylum, y = Count, x = Affinity)) + 
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Phosphorus affinity",
       y = "Proportion of classified ASVs") +
  scale_y_continuous(expand = c(0,0)) 
```

![](Plots/affinity/P%20affinity%20phylum-1.png)

``` r
# Class
highP_class <- as.data.frame(table(highP_taxa$Class))
colnames(highP_class) <- c("Class", "Count")
highP_class$Affinity <- "High"

lowP_class <- as.data.frame(table(lowP_taxa$Class))
colnames(lowP_class) <- c("Class", "Count")
lowP_class$Affinity <- "Low"

unclassP_class <- as.data.frame(table(unclassP_taxa$Class))
colnames(unclassP_class) <- c("Class", "Count")
unclassP_class$Affinity <- "Unclassified"

P_class <- rbind(highP_class, lowP_class, unclassP_class)

ggplot(P_class, aes(fill = Class, y = Count, x = Affinity)) +
  geom_bar(position = "Fill", stat = "identity") +
  labs(x = "Phosphorus affinity",
       y = "Proportion of classified ASVs") +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill=guide_legend(ncol = 2))
```

![](Plots/affinity/P%20affinity%20class-1.png)

### Precipitation (MAP) affinity group composition

``` r
highMAP_taxa <- as.data.frame(names(highMAP[2:ncol(highMAP)]))
highMAP_taxa <- rename(highMAP_taxa, ASV_index = `names(highMAP[2:ncol(highMAP)])`)
highMAP_taxa <- left_join(highMAP_taxa, ASV_taxa, by = "ASV_index")

lowMAP_taxa <- as.data.frame(names(lowMAP[2:ncol(lowMAP)]))
lowMAP_taxa <- rename(lowMAP_taxa, ASV_index = `names(lowMAP[2:ncol(lowMAP)])`)
lowMAP_taxa <- left_join(lowMAP_taxa, ASV_taxa, by = "ASV_index")

unclassMAP_taxa <- as.data.frame(names(unclassified_MAP[2:ncol(unclassified_MAP)]))
unclassMAP_taxa <- rename(unclassMAP_taxa,
                        ASV_index = `names(unclassified_MAP[2:ncol(unclassified_MAP)])`)
unclassMAP_taxa <- left_join(unclassMAP_taxa, ASV_taxa, by = "ASV_index")

# Phylum
highMAP_phylum <- as.data.frame(table(highMAP_taxa$Phylum))
colnames(highMAP_phylum) <- c("Phylum", "Count")
highMAP_phylum$Affinity <- "High"

lowMAP_phylum <- as.data.frame(table(lowMAP_taxa$Phylum))
colnames(lowMAP_phylum) <- c("Phylum", "Count")
lowMAP_phylum$Affinity <- "Low"

unclassMAP_phylum <- as.data.frame(table(unclassMAP_taxa$Phylum))
colnames(unclassMAP_phylum) <- c("Phylum", "Count")
unclassMAP_phylum$Affinity <- "Unclassified"

MAP_phylum <- rbind(highMAP_phylum, lowMAP_phylum, unclassMAP_phylum)

ggplot(MAP_phylum, aes(fill = Phylum, y = Count, x = Affinity)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Precipitation affinity",
       y = "Proportion of classified ASVs") +
  scale_y_continuous(expand = c(0,0)) 
```

![](Plots/affinity/MAP%20affinity%20phylum-1.png)

``` r
# Class
highMAP_class <- as.data.frame(table(highMAP_taxa$Class))
colnames(highMAP_class) <- c("Class", "Count")
highMAP_class$Affinity <- "High"

lowMAP_class <- as.data.frame(table(lowMAP_taxa$Class))
colnames(lowMAP_class) <- c("Class", "Count")
lowMAP_class$Affinity <- "Low"

unclassMAP_class <- as.data.frame(table(unclassMAP_taxa$Class))
colnames(unclassMAP_class) <- c("Class", "Count")
unclassMAP_class$Affinity <- "Unclassified"

MAP_class <- rbind(highMAP_class, lowMAP_class, unclassMAP_class)

ggplot(MAP_class, aes(fill = Class, y = Count, x = Affinity)) +
  geom_bar(position = "Fill", stat = "identity") +
  labs(x = "Precipitation affinity",
       y = "Proportion of classified ASVs") +
  guides(fill=guide_legend(ncol = 2)) + 
  scale_y_continuous(expand = c(0,0))
```

![](Plots/affinity/MAP%20affinity%20class-1.png)

### Nitrogen (DON) affinity group composition

``` r
highN_taxa <- as.data.frame(names(highN[2:ncol(highN)]))
highN_taxa <- rename(highN_taxa, ASV_index = `names(highN[2:ncol(highN)])`)
highN_taxa <- left_join(highN_taxa, ASV_taxa, by = "ASV_index")

lowN_taxa <- as.data.frame(names(lowN[2:ncol(lowN)]))
lowN_taxa <- rename(lowN_taxa, ASV_index = `names(lowN[2:ncol(lowN)])`)
lowN_taxa <- left_join(lowN_taxa, ASV_taxa, by = "ASV_index")

unclassN_taxa <- as.data.frame(names(unclassified_N[2:ncol(unclassified_N)]))
unclassN_taxa <- rename(unclassN_taxa,
                        ASV_index = `names(unclassified_N[2:ncol(unclassified_N)])`)
unclassN_taxa <- left_join(unclassN_taxa, ASV_taxa, by = "ASV_index")

# Phylum
highN_phylum <- as.data.frame(table(highN_taxa$Phylum))
colnames(highN_phylum) <- c("Phylum", "Count")
highN_phylum$Affinity <- "High"

lowN_phylum <- as.data.frame(table(lowN_taxa$Phylum))
colnames(lowN_phylum) <- c("Phylum", "Count")
lowN_phylum$Affinity <- "Low"

unclassN_phylum <- as.data.frame(table(unclassN_taxa$Phylum))
colnames(unclassN_phylum) <- c("Phylum", "Count")
unclassN_phylum$Affinity <- "Unclassified"

N_phylum <- rbind(highN_phylum, lowN_phylum, unclassN_phylum)

ggplot(N_phylum, aes(fill = Phylum, y = Count, x = Affinity)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Nitrogen affinity",
       y = "Proportion of classified ASVs") +
  scale_y_continuous(expand = c(0,0))
```

![](Plots/affinity/N%20affinity%20phylum-1.png)

``` r
# Class
highN_class <- as.data.frame(table(highN_taxa$Class))
colnames(highN_class) <- c("Class", "Count")
highN_class$Affinity <- "High"

lowN_class <- as.data.frame(table(lowN_taxa$Class))
colnames(lowN_class) <- c("Class", "Count")
lowN_class$Affinity <- "Low"

unclassN_class <- as.data.frame(table(unclassN_taxa$Class))
colnames(unclassN_class) <- c("Class", "Count")
unclassN_class$Affinity <- "Unclassified"

N_class <- rbind(highN_class, lowN_class, unclassN_class)

ggplot(N_class, aes(fill = Class, y = Count, x = Affinity)) +
  geom_bar(position = "Fill", stat = "identity") +
  labs(x = "Nitrogen affinity",
       y = "Proportion of classified ASVs") +
  guides(fill=guide_legend(ncol = 2)) +
  scale_y_continuous(expand = c(0,0))
```

![](Plots/affinity/N%20affinity%20class-1.png)
