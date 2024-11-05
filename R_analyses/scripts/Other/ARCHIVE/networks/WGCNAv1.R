install.packages("WGCNA")
library(WGCNA)

load(file = "data/DEG.RData")
counts <- counts(agg_cell_types) %>% t()
rownames(counts) <- do.call(paste, list(agg_cell_types$sample, agg_cell_types$treatment, agg_cell_types$sctype_labels, sep = "_"))
                            
input_mat = counts
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)


par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 20
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)


#######################################################################

cor <- temp_cor     # Return cor function to original namespace
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )






module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

splits <- str_split(mME$treatment, "_", simplify = TRUE)

mME$t2 <- paste(splits[,3], splits[,1], sep = "_")

mME %>% ggplot(., aes(x=t2, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")









# pick out a few modules of interest here

modules_of_interest = c("grey","red", "blue", "brown")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id
expr_normalized <- t(counts)
# Get normalized expression for those genes
expr_normalized[1:5,1:10]
#>                       B-3      B-4      B-5      L-3      L-4      L-5      S-3
#> AC149818.2_FG001 7.600901 7.077399 7.803434 7.220840 7.410408 8.028223 7.160846
#> AC149829.2_FG003 8.782014 8.179876 7.900062 8.299778 7.529891 8.631731 8.055118
#> AC182617.3_FG001 8.047244 7.120668 6.885533 7.501391 7.279413 7.809565 7.184253
#> AC186512.3_FG001 6.901539 7.389644 6.975945 6.859593 7.370816 6.633722 7.798843
#> AC186512.3_FG007 7.919688 7.754506 7.670946 7.417760 7.988427 7.904850 7.484542
#>                       S-4      S-5   B_L1.1
#> AC149818.2_FG001 7.401382 7.345322 6.524435
#> AC149829.2_FG003 8.744502 8.142909 8.240407
#> AC182617.3_FG001 8.140134 6.972400 7.777347
#> AC186512.3_FG001 6.949501 6.952659 6.059033
#> AC186512.3_FG007 8.375664 7.762799 6.335663
subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")




