library(dplyr); library(plotly); library(png); library(reticulate)

use_python(Sys.which("python3"), required = TRUE)
work_dir <- "."
data_dir <- "../data"

corr_file <- file.path(work_dir, "all_correlations.csv")
annot_file <- file.path(data_dir, "nutrient_annotation.csv")

pval_thresh <- 0.01
n_top <- 10

pw <- 1600; ph <- 1000; pscale <- 2
pmargin <- list(l = 10, r = 10, t = 200, b = 80)

cat_colors <- c(
  "Minerals" = "#4682B4", "Fat" = "#FFD700",
  "Water-Soluble Vitamins" = "#00CED1", "Carb" = "#32CD32",
  "Protein" = "#FF69B4", "Fat-Soluble Vitamins" = "#8A2BE2",
  "Fiber" = "#00FA9A", "Caffeine" = "#8B4513"
)

sp_colors <- c(
  "Amino Acid" = "#FF6347", "Lipid" = "#4169E1",
  "Peptide" = "#228B22", "Xenobiotics" = "#FF8C00",
  "Nucleotide" = "#DAA520", "Cofactors and Vitamins" = "#20B2AA",
  "Carbohydrate" = "#E377C2", "Energy" = "#7f7f7f",
  "Partially Characterized Molecules" = "#BC8F8F"
)

hex_to_rgba <- function(hex, alpha = 0.4) {
  r <- strtoi(substr(hex, 2, 3), 16)
  g <- strtoi(substr(hex, 4, 5), 16)
  b <- strtoi(substr(hex, 6, 7), 16)
  sprintf("rgba(%d,%d,%d,%.1f)", r, g, b, alpha)
}

#### significant correlations ####

results_df <- read.csv(corr_file, stringsAsFactors = FALSE)
annot <- read.csv(annot_file, stringsAsFactors = FALSE)
annot$nutrient <- gsub('^"|"$', "", annot$nutrient)

sig <- results_df[results_df$pvalue < pval_thresh, ]

sig$category <- annot$category[match(sig$dietary_component, annot$nutrient)]
sig <- sig[!is.na(sig$category), ]
sig <- sig[sig$category != "No use", ]

sig$sub_pathway <- gsub('^"|"$', "", sig$sub_pathway)
sig <- sig[!is.na(sig$sub_pathway) & sig$sub_pathway != "", ]

n_sig <- nrow(sig)

#######################
#### top pathways ####
connections <- sig %>%
  group_by(category, sub_pathway, super_pathway) %>%
  summarise(count = n(), .groups = "drop")

pw_counts <- connections %>%
  group_by(sub_pathway) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total))

top_pw <- pw_counts$sub_pathway[1:n_top]
conn_top <- connections %>% filter(sub_pathway %in% top_pw)

pw_sp_map <- conn_top %>%
  distinct(sub_pathway, super_pathway) %>%
  { setNames(.$super_pathway, .$sub_pathway) }

for (i in seq_along(top_pw)) {
  p <- top_pw[i]
  cat(sprintf("  %2d. %s (%d) -> %s\n", i, p,
      pw_counts$total[pw_counts$sub_pathway == p], pw_sp_map[p]))
}

####################
## bar plot data
sig_top <- sig[sig$sub_pathway %in% top_pw, ]
uniq_metabs <- sig_top[!duplicated(sig_top$metabolite),
               c("metabolite", "super_pathway", "sub_pathway")]
n_uniq <- nrow(uniq_metabs)

class_df <- uniq_metabs %>%
  count(super_pathway, name = "n") %>%
  arrange(desc(n)) %>%
  mutate(pct = n / n_uniq * 100)

bar_levels <- rev(class_df$super_pathway)
class_df$super_pathway <- factor(class_df$super_pathway, levels = bar_levels)
bar_classes <- as.character(bar_levels)

#####################
###### sankey ######
right_order <- c()
for (sp in rev(bar_classes)) {
  pw_in_sp <- top_pw[pw_sp_map[top_pw] == sp]
  pw_in_sp <- pw_in_sp[order(pw_counts$total[match(pw_in_sp, pw_counts$sub_pathway)],
                              decreasing = TRUE)]
  right_order <- c(right_order, pw_in_sp)
}

cat_order <- conn_top %>%
  group_by(category) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(category)

all_nodes <- c(cat_order, right_order)
node_idx <- setNames(seq_along(all_nodes) - 1, all_nodes)

get_node_y <- function(nodes, flows, total) {
  n <- length(nodes)
  pad <- 0.025
  avail <- 0.90 - (n - 1) * pad
  pos <- numeric(n); names(pos) <- nodes
  y <- 0.05
  for (i in seq_along(nodes)) {
    h <- (flows[nodes[i]] / total) * avail
    pos[i] <- y + h / 2
    y <- y + h + pad
  }
  pos
}

l_flows <- sapply(cat_order, function(x) sum(conn_top$count[conn_top$category == x]))
names(l_flows) <- cat_order
r_flows <- pw_counts$total[match(right_order, pw_counts$sub_pathway)]
names(r_flows) <- right_order

l_y <- get_node_y(cat_order, l_flows, sum(l_flows))
r_y <- get_node_y(right_order, r_flows, sum(r_flows))

node_x <- ifelse(all_nodes %in% cat_order, 0.001, 0.999)
node_y <- ifelse(all_nodes %in% cat_order, l_y[all_nodes], r_y[all_nodes])

sources <- node_idx[conn_top$category]
targets <- node_idx[conn_top$sub_pathway]

node_col <- sapply(all_nodes, function(nd) {
  if (nd %in% cat_order) cat_colors[nd] else sp_colors[pw_sp_map[nd]]
})
link_col <- sapply(conn_top$category, function(x) hex_to_rgba(cat_colors[x], 0.4))

fig <- plot_ly(
  type = "sankey", arrangement = "snap",
  node = list(pad = 35, thickness = 20,
              line = list(color = "black", width = 0.5),
              label = all_nodes, color = unname(node_col),
              x = node_x, y = node_y),
  link = list(source = unname(sources), target = unname(targets),
              value = conn_top$count, color = unname(link_col))
) %>% layout(
  title = list(
    text = paste0("Nutrient Categories \u2192 Top 10 Metabolic Pathways<br>",
                   "(p < 0.01, ", n_sig, " significant correlations)"),
    font = list(size = 44, color = "black"), x = 0.5, xanchor = "center"),
  font = list(size = 30, color = "black"),
  margin = pmargin,
  annotations = list(
    list(x = 0.01, y = 1.28, text = "Nutrients", showarrow = FALSE,
         font = list(size = 40, color = "black"), xref = "paper", yref = "paper"),
    list(x = 0.99, y = 1.28, text = "Metabolic Pathways", showarrow = FALSE,
         font = list(size = 40, color = "black"), xref = "paper", yref = "paper"))
)

tmp <- tempfile(fileext = ".png")
save_image(fig, file = tmp, width = pw, height = ph, scale = pscale)

#############################
###### sankey + barplot ######

img <- readPNG(tmp)
ih <- nrow(img); iw <- ncol(img)

out_png <- file.path(work_dir, "combined_final_figure_R.png")
out_pdf <- file.path(work_dir, "combined_final_figure_R.pdf")

draw_fig <- function() {
  layout(matrix(c(1, 2), nrow = 1), widths = c(4.8, 4.8))

  #sankey 
  par(mar = c(0.5, 0.5, 0.5, 0))
  plot(0, 0, type = "n", xlim = c(0, iw), ylim = c(ih, 0),
       xlab = "", ylab = "", axes = FALSE, asp = NA)
  rasterImage(img, 0, ih, iw, 0)

  # barplot 
  par(mar = c(4, 1, 5, 28))
  bpct <- class_df$pct[match(bar_classes, class_df$super_pathway)]
  bn <- class_df$n[match(bar_classes, class_df$super_pathway)]

  bp <- barplot(bpct, horiz = TRUE, col = sp_colors[bar_classes],
                border = "white", xlim = c(0, max(bpct) * 1.4),
                axes = TRUE, names.arg = rep("", length(bpct)),
                space = 0.5, cex.axis = 1.56)

  text(bpct + 0.8, bp, labels = sprintf("%.1f%% (n=%d)", bpct, bn),
       adj = 0, cex = 1.69, col = "black")
  mtext("% of Unique Metabolites", side = 1, line = 2.5, cex = 1.69, font = 2)
  mtext("Chemical Class Distribution", side = 3, line = 3, cex = 1.95, font = 2)
  mtext(sprintf("(Unique metabolites in Top %d pathways, n=%d)", n_top, n_uniq),
        side = 3, line = 1.5, cex = 1.43)

  legend("right", legend = rev(bar_classes),
         fill = sp_colors[rev(bar_classes)], border = "white",
         title = "Chemical Class", cex = 1.3, bty = "o",
         bg = adjustcolor("white", alpha.f = 0.9), box.col = "#cccccc",
         inset = c(-0.55, 0), xpd = TRUE)
}

png(out_png, width = 32, height = 12, units = "in", res = 300, bg = "white")
draw_fig(); dev.off()
pdf(out_pdf, width = 32, height = 12, bg = "white")
draw_fig(); dev.off()

unlink(tmp)

