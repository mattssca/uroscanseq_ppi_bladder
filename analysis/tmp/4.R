library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(scales)

# =============================================================================
# FIXED VISUALIZATION 1: SIGNIFICANCE OVERVIEW HEATMAP
# =============================================================================

create_significance_overview_fixed <- function(significant_hubs) {
  
  # Create summary matrix for all significant genes
  all_sig_genes <- unique(unlist(lapply(significant_hubs, function(x) x$gene)))
  all_sig_genes <- all_sig_genes[!is.na(all_sig_genes)]
  
  if(length(all_sig_genes) == 0) {
    cat("No significant genes found for visualization.\n")
    return(NULL)
  }
  
  # Create matrices for different significance types
  sig_matrix <- matrix(0, nrow = length(all_sig_genes), ncol = length(significant_hubs))
  expr_sig_matrix <- matrix(NA, nrow = length(all_sig_genes), ncol = length(significant_hubs))
  
  rownames(sig_matrix) <- all_sig_genes
  colnames(sig_matrix) <- names(significant_hubs)
  rownames(expr_sig_matrix) <- all_sig_genes
  colnames(expr_sig_matrix) <- names(significant_hubs)
  
  # Fill matrices
  for(subtype in names(significant_hubs)) {
    results <- significant_hubs[[subtype]]
    
    if(nrow(results) > 0) {
      for(i in 1:nrow(results)) {
        gene <- results$gene[i]
        if(!is.na(gene) && gene != "") {
          
          # Significance type based on number of methods
          sig_type <- results$n_methods_significant[i]
          sig_matrix[gene, subtype] <- sig_type
          
          # Store fold changes
          if(!is.na(results$expr_fold_change[i])) {
            expr_sig_matrix[gene, subtype] <- log2(results$expr_fold_change[i])
          }
        }
      }
    }
  }
  
  # Remove rows with all zeros
  sig_matrix_clean <- sig_matrix[rowSums(sig_matrix) > 0, , drop = FALSE]
  expr_sig_matrix_clean <- expr_sig_matrix[rowSums(!is.na(expr_sig_matrix)) > 0, , drop = FALSE]
  
  # Create significance overview plot
  if(nrow(sig_matrix_clean) > 0) {
    max_methods <- max(sig_matrix_clean, na.rm = TRUE)
    sig_colors <- colorRampPalette(c("white", "lightblue", "blue", "darkblue", "red"))(max_methods + 1)
    
    pheatmap(sig_matrix_clean,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = sig_colors,
             main = "Number of Significant Methods per Gene",
             filename = "figs/significance_overview_fixed.pdf",
             width = 8, height = max(6, nrow(sig_matrix_clean) * 0.2))
  }
  
  return(list(
    sig_matrix = sig_matrix_clean,
    expr_sig_matrix = expr_sig_matrix_clean
  ))
}

# =============================================================================
# FIXED VISUALIZATION 2: VOLCANO PLOTS
# =============================================================================

create_volcano_plots_fixed <- function(significant_hubs) {
  
  volcano_plots <- list()
  
  for(subtype in names(significant_hubs)) {
    results <- significant_hubs[[subtype]]
    
    if(nrow(results) > 0) {
      # Prepare data for volcano plot
      volcano_data <- results %>%
        filter(!is.na(expr_fold_change)) %>%
        mutate(
          log2_fc = log2(expr_fold_change),
          # Create a pseudo p-value based on significance
          pseudo_pvalue = case_when(
            n_methods_significant >= 3 ~ 0.001,
            n_methods_significant == 2 ~ 0.01,
            n_methods_significant == 1 ~ 0.05,
            TRUE ~ 0.1
          ),
          neg_log10_p = -log10(pseudo_pvalue),
          significance = case_when(
            n_methods_significant >= 3 ~ "3+ Methods",
            n_methods_significant == 2 ~ "2 Methods",
            n_methods_significant == 1 ~ "1 Method",
            TRUE ~ "Not Sig"
          ),
          label = ifelse(n_methods_significant >= 2 & abs(log2_fc) > 0.5, gene, "")
        )
      
      # Create volcano plot
      p <- ggplot(volcano_data, aes(x = log2_fc, y = neg_log10_p, color = significance)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_text_repel(aes(label = label), 
                        max.overlaps = 15,
                        box.padding = 0.3,
                        size = 3) +
        scale_color_manual(values = c(
          "3+ Methods" = "red",
          "2 Methods" = "orange",
          "1 Method" = "blue", 
          "Not Sig" = "grey"
        )) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "red") +
        theme_minimal() +
        labs(
          title = paste("Significance Plot -", toupper(subtype), "Subtype"),
          x = "Log2 Expression Fold Change",
          y = "-Log10 Pseudo P-value",
          color = "Significance"
        ) +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "bottom"
        )
      
      volcano_plots[[subtype]] <- p
      
      # Save individual plots
      ggsave(paste0("figs/significance_plot_", subtype, ".pdf"), p, width = 10, height = 8)
    }
  }
  
  return(volcano_plots)
}

# =============================================================================
# FIXED VISUALIZATION 3: FOLD CHANGE HEATMAPS
# =============================================================================

create_fold_change_heatmaps_fixed <- function(sig_matrices) {
  
  if(is.null(sig_matrices) || is.null(sig_matrices$expr_sig_matrix)) {
    cat("No expression data available for heatmap.\n")
    return(NULL)
  }
  
  # Expression fold change heatmap
  expr_matrix <- sig_matrices$expr_sig_matrix
  
  # Only keep rows with data in at least 2 subtypes
  complete_rows <- rowSums(!is.na(expr_matrix)) >= 2
  expr_matrix_clean <- expr_matrix[complete_rows, , drop = FALSE]
  
  if(nrow(expr_matrix_clean) > 2) {
    # Replace NAs with 0 for visualization
    expr_matrix_viz <- expr_matrix_clean
    expr_matrix_viz[is.na(expr_matrix_viz)] <- 0
    
    pheatmap(expr_matrix_viz,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             scale = "none",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(-2, 2, length.out = 101),
             main = "Expression Log2 Fold Changes (Significant Genes)",
             filename = "figs/expression_fold_changes_fixed.pdf",
             width = 8, height = max(6, nrow(expr_matrix_viz) * 0.3))
  }
}

# =============================================================================
# FIXED VISUALIZATION 4: TOP SIGNIFICANT HUBS BAR PLOTS
# =============================================================================

create_top_hubs_barplots_fixed <- function(significant_hubs, top_n = 10) {
  
  for(subtype in names(significant_hubs)) {
    results <- significant_hubs[[subtype]]
    
    if(nrow(results) > 0) {
      # Get top significant genes
      top_genes <- results %>%
        filter(n_methods_significant >= 1) %>%
        arrange(desc(n_methods_significant), expr_fdr) %>%
        head(top_n) %>%
        filter(!is.na(expr_fold_change))
      
      if(nrow(top_genes) > 0) {
        # Expression fold change plot
        expr_plot <- ggplot(top_genes, aes(x = reorder(gene, expr_fold_change), 
                                           y = expr_fold_change, 
                                           fill = factor(n_methods_significant))) +
          geom_bar(stat = "identity") +
          geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
          scale_fill_manual(values = c("1" = "lightblue", "2" = "blue", "3" = "darkblue"),
                            name = "# Methods\nSignificant") +
          coord_flip() +
          theme_minimal() +
          labs(
            title = paste("Top Significant Genes -", toupper(subtype)),
            x = "Gene",
            y = "Expression Fold Change"
          )
        
        # Save plot
        ggsave(paste0("figs/top_significant_genes_", subtype, ".pdf"), expr_plot, 
               width = 10, height = max(6, nrow(top_genes) * 0.4))
      }
    }
  }
}

# =============================================================================
# FIXED VISUALIZATION 5: EXPRESSION DISTRIBUTION COMPARISON
# =============================================================================

create_expression_comparison_fixed <- function(significant_hubs) {
  
  # Combine all data
  all_data <- bind_rows(significant_hubs, .id = "subtype") %>%
    filter(!is.na(expr_fold_change))
  
  if(nrow(all_data) > 0) {
    # Overall distribution plot
    dist_plot <- ggplot(all_data, aes(x = log2(expr_fold_change), fill = subtype)) +
      geom_density(alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      facet_wrap(~subtype, scales = "free_y") +
      theme_minimal() +
      labs(
        title = "Distribution of Expression Fold Changes",
        x = "Log2 Expression Fold Change",
        y = "Density"
      ) +
      theme(legend.position = "none")
    
    ggsave("figs/expression_distribution_comparison_fixed.pdf", dist_plot, width = 12, height = 8)
    
    # Methods comparison
    methods_plot <- ggplot(all_data, aes(x = factor(n_methods_significant), y = abs(log2(expr_fold_change)))) +
      geom_boxplot(aes(fill = factor(n_methods_significant))) +
      scale_fill_manual(values = c("1" = "lightblue", "2" = "orange", "3" = "red"),
                        name = "# Methods") +
      theme_minimal() +
      labs(
        title = "Effect Size vs Number of Significant Methods",
        x = "Number of Significant Methods",
        y = "Absolute Log2 Fold Change"
      )
    
    ggsave("figs/methods_vs_effect_size.pdf", methods_plot, width = 8, height = 6)
    
    return(list(dist_plot = dist_plot, methods_plot = methods_plot))
  }
  
  return(NULL)
}

# =============================================================================
# FIXED VISUALIZATION 6: SIGNIFICANCE SUMMARY
# =============================================================================

create_significance_summary_fixed <- function(significant_hubs) {
  
  # Count significant genes per subtype
  summary_data <- data.frame()
  
  for(subtype in names(significant_hubs)) {
    results <- significant_hubs[[subtype]]
    
    if(nrow(results) > 0) {
      counts <- data.frame(
        subtype = subtype,
        one_method = sum(results$n_methods_significant == 1, na.rm = TRUE),
        two_methods = sum(results$n_methods_significant == 2, na.rm = TRUE),
        three_methods = sum(results$n_methods_significant >= 3, na.rm = TRUE),
        total_significant = nrow(results)
      )
    } else {
      counts <- data.frame(
        subtype = subtype,
        one_method = 0,
        two_methods = 0,
        three_methods = 0,
        total_significant = 0
      )
    }
    
    summary_data <- rbind(summary_data, counts)
  }
  
  # Reshape for plotting
  summary_long <- summary_data %>%
    select(-total_significant) %>%
    pivot_longer(-subtype, names_to = "significance_level", values_to = "count")
  
  # Create stacked bar plot
  summary_plot <- ggplot(summary_long, aes(x = subtype, y = count, fill = significance_level)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c(
      "one_method" = "lightblue",
      "two_methods" = "orange", 
      "three_methods" = "red"
    ), labels = c("1 Method", "2 Methods", "3+ Methods")) +
    theme_minimal() +
    labs(
      title = "Number of Significant Genes by Evidence Level",
      x = "Subtype",
      y = "Number of Significant Genes",
      fill = "Evidence Level"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave("figs/significance_summary_fixed.pdf", summary_plot, width = 10, height = 6)
  
  return(summary_plot)
}

# =============================================================================
# FIXED MASTER VISUALIZATION FUNCTION
# =============================================================================

create_all_significance_visualizations_fixed <- function(significant_hubs) {
  
  cat("Creating FIXED significance visualizations...\n")
  
  # Create output directory
  if(!dir.exists("figs")) dir.create("figs")
  
  # Check if we have any significant results
  total_genes <- sum(sapply(significant_hubs, nrow))
  if(total_genes == 0) {
    cat("No significant genes found. Creating summary plots only.\n")
    summary_plot <- create_significance_summary_fixed(significant_hubs)
    return(list(summary_plot = summary_plot))
  }
  
  # 1. Significance overview
  cat("1. Creating significance overview heatmap...\n")
  sig_matrices <- create_significance_overview_fixed(significant_hubs)
  
  # 2. Volcano-like plots
  cat("2. Creating significance plots...\n")
  volcano_plots <- create_volcano_plots_fixed(significant_hubs)
  
  # 3. Fold change heatmaps
  cat("3. Creating fold change heatmaps...\n")
  create_fold_change_heatmaps_fixed(sig_matrices)
  
  # 4. Top hubs bar plots
  cat("4. Creating top hubs bar plots...\n")
  create_top_hubs_barplots_fixed(significant_hubs, top_n = 15)
  
  # 5. Expression comparison plots
  cat("5. Creating expression comparison plots...\n")
  comparison_plots <- create_expression_comparison_fixed(significant_hubs)
  
  # 6. Significance summary
  cat("6. Creating significance summary plots...\n")
  summary_plot <- create_significance_summary_fixed(significant_hubs)
  
  cat("All FIXED visualizations complete! Check the figs/ directory.\n")
  
  return(list(
    sig_matrices = sig_matrices,
    volcano_plots = volcano_plots,
    comparison_plots = comparison_plots,
    summary_plot = summary_plot
  ))
}

# =============================================================================
# RUN FIXED VISUALIZATIONS
# =============================================================================

# Create all fixed visualizations
all_plots_fixed <- create_all_significance_visualizations_fixed(significant_hubs)

# Print summary of what was created
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("FIXED VISUALIZATION SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n")
cat("Created the following visualizations:\n")
cat("1. significance_overview_fixed.pdf - Methods-based significance patterns\n")
cat("2. significance_plot_[subtype].pdf - Significance plots for each subtype\n")
cat("3. expression_fold_changes_fixed.pdf - Expression FC heatmap\n")
cat("4. top_significant_genes_[subtype].pdf - Top genes per subtype\n")
cat("5. expression_distribution_comparison_fixed.pdf - FC distributions\n")
cat("6. methods_vs_effect_size.pdf - Evidence quality vs effect size\n")
cat("7. significance_summary_fixed.pdf - Count of significant genes\n")
cat(paste(rep("=", 60), collapse=""), "\n")
