rm(list = ls())
gc()
graphics.off()

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(deldir)
  library(cowplot)
  library(grid)
})

# =========================================================
# 01 paths
# =========================================================
find_project_root <- function(start = getwd()) {
  current <- tryCatch(
    normalizePath(start, winslash = "/", mustWork = TRUE),
    error = function(e) NA_character_
  )
  if (is.na(current)) stop("Cannot resolve the current working directory.")
  
  repeat {
    if (file.exists(file.path(current, "TNBC_TMA_6k_CosMx.Rproj"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  
  stop(
    "Cannot find TNBC_TMA_6k_CosMx.Rproj. ",
    "Open the R project or run the script from the repository root/subdirectory."
  )
}

proj <- find_project_root()

raw_dir <- file.path(proj, "raw")
out_dir <- file.path(proj, "results", "05_0_1")
data_dir <- file.path(out_dir, "data")
fig_dir  <- file.path(out_dir, "figures")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

save_csv_overwrite <- function(x, filename) {
  out <- file.path(data_dir, paste0(filename, ".csv"))
  if (file.exists(out)) file.remove(out)
  fwrite(x, out)
  invisible(out)
}

# =========================================================
# 02 user settings
# =========================================================
celltype_col <- "BC6k_refined"

tumor_label <- "tumor"

# Multi-page overview layout
overview_ncol <- 4
overview_nrow <- 4
plots_per_page <- overview_ncol * overview_nrow

# PDF size
overview_width  <- 15
overview_height <- 15

# Voronoi settings
expand_ratio <- 0.01

# Local density filter
nn_keep_quantile <- 0.90
nn_keep_mult <- 1.15

polygon_linewidth <- 0.08

# Save one PDF per FOV
save_individual_pdf <- FALSE

individual_width  <- 5
individual_height <- 5

# Background
plot_bg  <- "black"
panel_bg <- "black"

axis_text_col <- "grey75"
axis_line_col <- "grey55"
border_col <- "grey70"

# =========================================================
# 03 cell-type palette
# =========================================================
celltype_cols <- c(
  "B.cell"                      = "#F2A65A",
  "Plasma"                      = "#D98C2B",
  "Plasmablast"                 = "#B87333",
  "Plasmacytoid.dendritic.cell" = "#8C6D31",
  "Conventional.dendritic.cell" = "#8BC34A",
  "Macrophage"                  = "#6DAA2C",
  "Mast.cell"                   = "grey50",
  "Monocyte"                    = "#35B56A",
  "Neutrophil"                  = "#2EB7B0",
  "NK.cell"                     = "#23A6D5",
  "T.cell.CD4"                  = "#4A97C9",
  "T.cell.regulatory"           = "#7B8CC7",
  "T.cell.CD8"                  = "#D864A6",
  "Fibroblast"                  = "#B1AA34",
  "Endothelial"                 = "#C99A2A",
  "Pericyte"                    = "#8E7D22",
  "tumor"                       = "#EE7068"
)

# =========================================================
# 04 helper functions
# =========================================================
harmonize_celltype_names <- function(x) {
  
  x <- as.character(x)
  
  name_map <- c(
    "Dendritic.cell" = "Conventional.dendritic.cell"
  )
  
  idx <- x %in% names(name_map)
  x[idx] <- name_map[x[idx]]
  
  x
}

get_nn1_dist <- function(x, y) {
  
  mat <- as.matrix(data.frame(
    x = x,
    y = y
  ))
  
  d <- as.matrix(dist(mat))
  
  diag(d) <- Inf
  
  apply(
    d,
    1,
    min,
    na.rm = TRUE
  )
}

make_voronoi_dt <- function(
    df,
    x_col = "x",
    y_col = "y",
    id_col = "cell_id",
    expand_ratio = 0.01
) {
  
  stopifnot(
    all(
      c(x_col, y_col, id_col) %in% colnames(df)
    )
  )
  
  dt0 <- copy(
    df[, .(
      cell_id = get(id_col),
      x = get(x_col),
      y = get(y_col)
    )]
  )
  
  dt0 <- dt0[
    !is.na(x) &
      !is.na(y)
  ]
  
  if (nrow(dt0) < 3) {
    stop(
      "Too few cells for Voronoi tessellation."
    )
  }
  
  dup_idx <-
    duplicated(dt0[, .(x, y)]) |
    duplicated(
      dt0[, .(x, y)],
      fromLast = TRUE
    )
  
  if (any(dup_idx)) {
    
    set.seed(1)
    
    dt0[
      dup_idx,
      x := x + runif(.N, -1e-8, 1e-8)
    ]
    
    dt0[
      dup_idx,
      y := y + runif(.N, -1e-8, 1e-8)
    ]
  }
  
  xr <- range(
    dt0$x,
    na.rm = TRUE
  )
  
  yr <- range(
    dt0$y,
    na.rm = TRUE
  )
  
  xpad <- diff(xr) * expand_ratio
  ypad <- diff(yr) * expand_ratio
  
  rw <- c(
    xr[1] - xpad,
    xr[2] + xpad,
    yr[1] - ypad,
    yr[2] + ypad
  )
  
  dd <- deldir::deldir(
    dt0$x,
    dt0$y,
    rw = rw
  )
  
  tt <- deldir::tile.list(dd)
  
  poly_list <- vector(
    "list",
    length(tt)
  )
  
  for (i in seq_along(tt)) {
    
    tile_i <- tt[[i]]
    
    if (length(tile_i$x) < 3) {
      next
    }
    
    poly_list[[i]] <- data.table(
      cell_id = dt0$cell_id[i],
      x = tile_i$x,
      y = tile_i$y,
      polygon_id = i,
      vertex_order = seq_along(tile_i$x)
    )
  }
  
  rbindlist(
    poly_list,
    use.names = TRUE,
    fill = TRUE
  )
}

build_fov_polygon_dt <- function(dt_fov) {
  
  dt_fov <- copy(dt_fov)
  
  if (nrow(dt_fov) < 3) {
    return(NULL)
  }
  
  dt_fov[, nn1_dist := get_nn1_dist(x, y)]
  
  nn_keep_cut <- as.numeric(
    quantile(
      dt_fov$nn1_dist,
      probs = nn_keep_quantile,
      na.rm = TRUE
    )
  ) * nn_keep_mult
  
  dt_fov[
    ,
    keep_polygon := nn1_dist <= nn_keep_cut
  ]
  
  poly_dt <- make_voronoi_dt(
    df = dt_fov[, .(
      cell_id,
      x,
      y
    )],
    x_col = "x",
    y_col = "y",
    id_col = "cell_id",
    expand_ratio = expand_ratio
  )
  
  annotation_dt <- copy(
    dt_fov[
      ,
      !c("x", "y"),
      with = FALSE
    ]
  )
  
  poly_dt <- merge(
    poly_dt,
    annotation_dt,
    by = "cell_id",
    all.x = TRUE
  )
  
  setorder(
    poly_dt,
    polygon_id,
    vertex_order
  )
  
  poly_dt <- poly_dt[
    keep_polygon == TRUE
  ]
  
  list(
    polygons = poly_dt,
    nn_keep_cut = nn_keep_cut
  )
}

theme_overview <- function(show_title = TRUE) {
  
  theme_bw(base_size = 7) +
    theme(
      panel.grid = element_blank(),
      
      axis.title = element_blank(),
      
      axis.text = element_text(
        size = 5.5,
        color = axis_text_col
      ),
      
      axis.ticks = element_line(
        color = axis_line_col,
        linewidth = 0.15
      ),
      
      plot.title = if (show_title) {
        element_text(
          size = 10,
          face = "plain",
          color = "white",
          hjust = 0.5,
          margin = margin(
            b = 2
          )
        )
      } else {
        element_blank()
      },
      
      panel.border = element_rect(
        color = border_col,
        linewidth = 0.35,
        fill = NA
      ),
      
      panel.background = element_rect(
        fill = panel_bg,
        color = panel_bg
      ),
      
      plot.background = element_rect(
        fill = plot_bg,
        color = plot_bg
      ),
      
      legend.position = "none",
      
      plot.margin = margin(
        2,
        2,
        2,
        2,
        unit = "mm"
      )
    )
}

make_fov_plot <- function(
    fov_id,
    dt_all
) {
  
  dt_fov <- copy(
    dt_all[
      FOV == fov_id
    ]
  )
  
  obj <- build_fov_polygon_dt(
    dt_fov
  )
  
  if (is.null(obj)) {
    
    p_empty <- ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = paste0(
          "FOV ",
          fov_id,
          "\nInsufficient cells"
        ),
        color = "white",
        size = 4
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void() +
      theme(
        plot.background = element_rect(
          fill = plot_bg,
          color = plot_bg
        ),
        panel.background = element_rect(
          fill = panel_bg,
          color = panel_bg
        )
      )
    
    return(
      list(
        plot = p_empty,
        nn_keep_cut = NA_real_
      )
    )
  }
  
  poly_dt <- obj$polygons
  
  levels_present <-
    names(celltype_cols)[
      names(celltype_cols) %in%
        unique(poly_dt$celltype_plot)
    ]
  
  poly_dt[
    ,
    celltype_plot := factor(
      celltype_plot,
      levels = levels_present
    )
  ]
  
  p <- ggplot(
    poly_dt,
    aes(
      x = x,
      y = y,
      group = polygon_id,
      fill = celltype_plot
    )
  ) +
    geom_polygon(
      color = "black",
      linewidth = polygon_linewidth
    ) +
    scale_fill_manual(
      values = celltype_cols,
      na.value = "grey50",
      drop = FALSE
    ) +
    coord_fixed() +
    labs(
      title = paste0(
        "FOV ",
        fov_id
      )
    ) +
    theme_overview()
  
  list(
    plot = p,
    nn_keep_cut = obj$nn_keep_cut
  )
}

make_legend_plot <- function() {
  
  legend_dt <- data.table(
    x = seq_along(celltype_cols),
    y = 1,
    celltype = factor(
      names(celltype_cols),
      levels = names(celltype_cols)
    )
  )
  
  ggplot(
    legend_dt,
    aes(
      x = x,
      y = y,
      fill = celltype
    )
  ) +
    geom_tile() +
    scale_fill_manual(
      values = celltype_cols,
      drop = FALSE,
      guide = guide_legend(
        title = "BC6k_refined",
        ncol = 2,
        byrow = TRUE,
        override.aes = list(
          color = NA
        )
      )
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      
      legend.title = element_text(
        size = 9,
        color = "white"
      ),
      
      legend.text = element_text(
        size = 8,
        color = "white"
      ),
      
      legend.key.height = unit(
        4,
        "mm"
      ),
      
      legend.key.width = unit(
        4,
        "mm"
      ),
      
      legend.background = element_rect(
        fill = plot_bg,
        color = plot_bg
      ),
      
      legend.box.background = element_rect(
        fill = plot_bg,
        color = plot_bg
      ),
      
      plot.background = element_rect(
        fill = plot_bg,
        color = plot_bg
      )
    )
}

# =========================================================
# 05 load data
# =========================================================
meta <- as.data.table(
  readRDS(
    file.path(
      raw_dir,
      "TNBC_6k_metadata_Bogang_selected.RDS"
    )
  )
)

xy <- readRDS(
  file.path(
    raw_dir,
    "TNBC_6k_xy_Bogang.RDS"
  )
)

# =========================================================
# 06 prepare coordinates
# =========================================================
xy_dt <- as.data.table(xy)

xy_dt[
  ,
  cell_id := rownames(xy)
]

setcolorder(
  xy_dt,
  c(
    "cell_id",
    setdiff(
      names(xy_dt),
      "cell_id"
    )
  )
)

need_xy <- c(
  "x_mm",
  "y_mm"
)

missing_xy <- setdiff(
  need_xy,
  names(xy_dt)
)

if (length(missing_xy) > 0) {
  stop(
    "Missing coordinate columns in xy object: ",
    paste(
      missing_xy,
      collapse = ", "
    )
  )
}

setnames(
  xy_dt,
  old = c(
    "x_mm",
    "y_mm"
  ),
  new = c(
    "x",
    "y"
  )
)

# =========================================================
# 07 merge metadata and coordinates
# =========================================================
need_meta <- c(
  "cell_id",
  "FOV",
  celltype_col
)

missing_meta <- setdiff(
  need_meta,
  names(meta)
)

if (length(missing_meta) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(
      missing_meta,
      collapse = ", "
    )
  )
}

dt <- merge(
  meta,
  xy_dt[, .(
    cell_id,
    x,
    y
  )],
  by = "cell_id",
  all.x = TRUE
)

dt <- dt[
  !is.na(FOV) &
    !is.na(x) &
    !is.na(y)
]

dt[
  ,
  FOV := as.character(FOV)
]

dt[
  ,
  celltype_plot := harmonize_celltype_names(
    get(celltype_col)
  )
]

# =========================================================
# 08 validate cell-type labels
# =========================================================
unknown_celltypes <- setdiff(
  unique(dt$celltype_plot),
  names(celltype_cols)
)

unknown_celltypes <- unknown_celltypes[
  !is.na(unknown_celltypes)
]

if (length(unknown_celltypes) > 0) {
  
  warning(
    "Cell types without predefined colors: ",
    paste(
      unknown_celltypes,
      collapse = ", "
    )
  )
}

# =========================================================
# 09 FOV summary
# =========================================================
fov_summary <- dt[
  ,
  .(
    n_cells = .N,
    n_tumor = sum(
      get(celltype_col) == tumor_label,
      na.rm = TRUE
    ),
    tumor_fraction = mean(
      get(celltype_col) == tumor_label,
      na.rm = TRUE
    ),
    n_celltypes = uniqueN(
      celltype_plot
    )
  ),
  by = FOV
]

setorder(
  fov_summary,
  FOV
)

save_csv_overwrite(
  fov_summary,
  "05_0_1_FOV_summary"
)

# =========================================================
# 10 FOV cell-type composition
# =========================================================
fov_celltype_composition <- dt[
  ,
  .(
    n = .N
  ),
  by = .(
    FOV,
    celltype_plot
  )
]

fov_celltype_composition[
  ,
  fraction := n / sum(n),
  by = FOV
]

setorder(
  fov_celltype_composition,
  FOV,
  -fraction
)

save_csv_overwrite(
  fov_celltype_composition,
  "05_0_1_FOV_celltype_composition"
)

# =========================================================
# 11 FOV order
# =========================================================
all_fovs <- sort(
  unique(
    dt$FOV
  )
)

fov_order_dt <- data.table(
  plot_order = seq_along(all_fovs),
  FOV = all_fovs
)

save_csv_overwrite(
  fov_order_dt,
  "05_0_1_FOV_plot_order"
)

message(
  "Total FOVs: ",
  length(all_fovs)
)

# =========================================================
# 12 build all FOV plots
# =========================================================
fov_plot_list <- vector(
  "list",
  length(all_fovs)
)

nn_summary_list <- vector(
  "list",
  length(all_fovs)
)

names(fov_plot_list) <- all_fovs

for (i in seq_along(all_fovs)) {
  
  fov_id <- all_fovs[i]
  
  message(
    sprintf(
      "[%d/%d] Processing %s",
      i,
      length(all_fovs),
      fov_id
    )
  )
  
  obj <- tryCatch(
    make_fov_plot(
      fov_id = fov_id,
      dt_all = dt
    ),
    error = function(e) {
      
      warning(
        "Failed FOV ",
        fov_id,
        ": ",
        conditionMessage(e)
      )
      
      p_fail <- ggplot() +
        annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = paste0(
            "FOV ",
            fov_id,
            "\nFailed"
          ),
          color = "white",
          size = 4
        ) +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void() +
        theme(
          plot.background = element_rect(
            fill = plot_bg,
            color = plot_bg
          ),
          panel.background = element_rect(
            fill = plot_bg,
            color = plot_bg
          )
        )
      
      list(
        plot = p_fail,
        nn_keep_cut = NA_real_
      )
    }
  )
  
  fov_plot_list[[i]] <- obj$plot
  
  nn_summary_list[[i]] <- data.table(
    FOV = fov_id,
    nn_keep_cut = obj$nn_keep_cut
  )
  
  if (save_individual_pdf) {
    
    individual_file <- file.path(
      fig_dir,
      paste0(
        "05_0_1_",
        fov_id,
        "_all_cell_types.pdf"
      )
    )
    
    if (file.exists(individual_file)) {
      file.remove(individual_file)
    }
    
    ggsave(
      filename = individual_file,
      plot = obj$plot,
      width = individual_width,
      height = individual_height,
      bg = plot_bg
    )
  }
}

nn_summary_dt <- rbindlist(
  nn_summary_list,
  fill = TRUE
)

save_csv_overwrite(
  nn_summary_dt,
  "05_0_1_FOV_local_density_threshold"
)

# =========================================================
# 13 save separate cell-type legend
# =========================================================
legend_source <- make_legend_plot()

legend_grob <- cowplot::get_legend(
  legend_source
)

legend_pdf <- file.path(
  fig_dir,
  "05_0_1_all_cell_types_legend.pdf"
)

if (file.exists(legend_pdf)) {
  file.remove(legend_pdf)
}

pdf(
  legend_pdf,
  width = 5.0,
  height = 3.2,
  bg = plot_bg
)

grid.newpage()

grid.draw(
  legend_grob
)

dev.off()

# =========================================================
# 14 multi-page overview PDF
# =========================================================
overview_file <- file.path(
  fig_dir,
  "05_0_1_all_FOV_all_cell_types_overview.pdf"
)

if (file.exists(overview_file)) {
  file.remove(overview_file)
}

n_pages <- ceiling(
  length(fov_plot_list) /
    plots_per_page
)

pdf(
  overview_file,
  width = overview_width,
  height = overview_height,
  bg = plot_bg,
  onefile = TRUE
)

for (page_i in seq_len(n_pages)) {
  
  start_i <-
    (page_i - 1) *
    plots_per_page + 1
  
  end_i <- min(
    page_i * plots_per_page,
    length(fov_plot_list)
  )
  
  page_plots <- fov_plot_list[
    start_i:end_i
  ]
  
  n_missing <-
    plots_per_page -
    length(page_plots)
  
  if (n_missing > 0) {
    
    blank_plot <- ggplot() +
      theme_void() +
      theme(
        plot.background = element_rect(
          fill = plot_bg,
          color = plot_bg
        ),
        panel.background = element_rect(
          fill = plot_bg,
          color = plot_bg
        )
      )
    
    page_plots <- c(
      page_plots,
      rep(
        list(blank_plot),
        n_missing
      )
    )
  }
  
  page_grid <- cowplot::plot_grid(
    plotlist = page_plots,
    ncol = overview_ncol,
    nrow = overview_nrow,
    align = "hv",
    axis = "tblr"
  )
  
  page_grid <- cowplot::ggdraw(
    page_grid
  ) +
    theme(
      plot.background = element_rect(
        fill = plot_bg,
        color = plot_bg
      )
    )
  
  print(
    page_grid
  )
}

dev.off()

# =========================================================
# 15 optional overview with FOV summary order
# =========================================================
selection_table <- merge(
  fov_order_dt,
  fov_summary,
  by = "FOV",
  all.x = TRUE
)

setorder(
  selection_table,
  plot_order
)

save_csv_overwrite(
  selection_table,
  "05_0_1_FOV_selection_reference"
)

message("05_0_1 finished successfully.")
message("Total FOVs: ", length(all_fovs))
message("Total overview pages: ", n_pages)
message("Data directory: ", data_dir)
message("Figure directory: ", fig_dir)