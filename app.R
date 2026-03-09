# =========================================================
# app.R
# Anubis — Ancient Egyptian Collagen DB
# Update requested:
# - Round mean_peptide_score to 2 decimals (بعد العلامه) everywhere it appears
# - Make Hagar slightly unique (very subtle): gold frame + tiny soft glow (not too much)
# - Keep Nawal clean with white frame
# - Team tab stays LAST
# - Footer fixed on every page
# =========================================================

# -----------------------------
# 0) Libraries
# -----------------------------
library(shiny)
library(bslib)
library(dplyr)
library(readr)
library(DT)
library(plotly)
library(ggplot2)
library(stringr)
library(purrr)
library(scales)
library(tidyr)
library(htmltools)
library(tibble)

# -----------------------------
# 1) Operators and Helpers
# -----------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_read_csv <- function(path) {
  if (is.null(path) || length(path) == 0 || !file.exists(path)) return(NULL)
  tryCatch(
    suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE)),
    error = function(e) {
      message("Failed to read CSV: ", path, " | ", e$message)
      NULL
    }
  )
}

clean_names_local <- function(df) {
  if (is.null(df)) return(NULL)
  names(df) <- names(df) |>
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") |>
    stringr::str_replace_all("(^_+|_+$)", "") |>
    tolower()
  df
}

normalize_tissue_value <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(
    str_detect(x, "bone|osse|femur|tibia|humerus|skull|tooth|dentin|enamel|mandible|rib") ~ "Bone",
    str_detect(x, "skin|hide|derm|dermis|epiderm|leather") ~ "Skin",
    TRUE ~ NA_character_
  )
}

infer_tissue_from_columns <- function(df) {
  if (is.null(df)) return(rep(NA_character_, 0))
  
  candidates <- c(
    "tissue", "source_tissue", "sample_tissue", "material", "sample_type",
    "dataset", "group", "origin", "source", "source_file", "input_file",
    "file", "filename", "file_name", "db_name", "database", "fasta_path",
    "folder", "path", "project", "subset"
  )
  
  n <- nrow(df)
  inferred <- rep(NA_character_, n)
  
  for (col in intersect(candidates, names(df))) {
    vals <- normalize_tissue_value(df[[col]])
    fill_idx <- is.na(inferred) & !is.na(vals)
    inferred[fill_idx] <- vals[fill_idx]
  }
  
  if ("protein_id" %in% names(df)) {
    vals <- normalize_tissue_value(df$protein_id)
    fill_idx <- is.na(inferred) & !is.na(vals)
    inferred[fill_idx] <- vals[fill_idx]
  }
  
  if ("species" %in% names(df)) {
    vals <- normalize_tissue_value(df$species)
    fill_idx <- is.na(inferred) & !is.na(vals)
    inferred[fill_idx] <- vals[fill_idx]
  }
  
  inferred
}

ensure_tissue_column <- function(df, default_label = "Combined") {
  if (is.null(df)) return(NULL)
  df <- clean_names_local(df)
  
  if (!"tissue" %in% names(df)) {
    inferred <- infer_tissue_from_columns(df)
    if (length(inferred) == nrow(df) && any(!is.na(inferred))) {
      df$tissue <- inferred
    } else {
      df$tissue <- default_label
    }
  } else {
    vals <- normalize_tissue_value(df$tissue)
    df$tissue <- ifelse(is.na(vals), as.character(df$tissue), vals)
    df$tissue[df$tissue %in% c("", "NA", "na")] <- default_label
  }
  
  df$tissue <- dplyr::case_when(
    tolower(df$tissue) == "bone" ~ "Bone",
    tolower(df$tissue) == "skin" ~ "Skin",
    TRUE ~ default_label
  )
  
  df
}

safe_n_distinct <- function(df, col) {
  if (is.null(df) || !col %in% names(df)) return(0)
  dplyr::n_distinct(df[[col]], na.rm = TRUE)
}

format_big <- function(x) scales::comma(x %||% 0)

# Round to 2 decimals (بعد العلامه)
round2 <- function(x) {
  if (is.null(x)) return(x)
  suppressWarnings({
    y <- as.numeric(x)
    ifelse(is.na(y), x, round(y, 2))
  })
}

metric_box <- function(title, value, subtitle = NULL, accent = "#0f172a") {
  div(
    class = "metric-card",
    style = paste0("--accent:", accent, ";"),
    div(class = "metric-label", title),
    div(class = "metric-value", value),
    if (!is.null(subtitle)) div(class = "metric-sub", subtitle)
  )
}

section_title <- function(title, subtitle = NULL) {
  div(
    class = "section-head",
    h2(title),
    if (!is.null(subtitle)) p(subtitle)
  )
}

empty_dt <- function(msg = "No data available.") {
  datatable(
    data.frame(Message = msg),
    rownames = FALSE,
    options = list(dom = "t")
  )
}

# -----------------------------
# 2) Plot helpers
# -----------------------------
plot_palette_main <- c(
  "#0f172a", "#175cd3", "#0f766e", "#c0841a",
  "#7c3aed", "#dc2626", "#0891b2", "#16a34a",
  "#d97706", "#db2777", "#2563eb", "#475569"
)

plot_palette_tissue <- c(
  "Skin" = "#16a34a",
  "Bone" = "#c0841a",
  "Combined" = "#2563eb"
)

base_gg_theme <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#e5e7eb"),
      axis.title.x = element_text(color = "#344054"),
      axis.title.y = element_text(color = "#344054"),
      axis.text.x = element_text(color = "#344054"),
      axis.text.y = element_text(color = "#344054"),
      plot.title = element_text(face = "bold", color = "#101828")
    )
}

ggplotly_clean <- function(p, tooltip = "text") {
  ggplotly(p, tooltip = tooltip) %>%
    layout(
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)",
      legend = list(orientation = "h", y = -0.18)
    )
}

# -----------------------------
# 3) File discovery helpers
# -----------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg_name <- "--file="
  script_path <- sub(file_arg_name, "", args[grep(file_arg_name, args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)))
  }
  
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)))
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

build_search_dirs <- function(user_dir = NULL) {
  dirs <- c(
    user_dir,
    Sys.getenv("DATA_DIR", unset = NA),
    "/content",
    getwd(),
    get_script_dir(),
    file.path(getwd(), "data"),
    file.path(get_script_dir(), "data"),
    path.expand("~/Downloads"),
    path.expand("~")
  )
  
  dirs <- dirs[!is.na(dirs) & nzchar(dirs)]
  unique(normalizePath(dirs, winslash = "/", mustWork = FALSE))
}

find_existing_file <- function(filename, search_dirs) {
  candidates <- unique(c(filename, file.path(search_dirs, filename)))
  exists <- candidates[file.exists(candidates)]
  if (length(exists) == 0) return(NULL)
  normalizePath(exists[1], winslash = "/", mustWork = FALSE)
}

load_csv_with_meta <- function(filename, search_dirs, required = FALSE) {
  path <- find_existing_file(filename, search_dirs)
  df <- safe_read_csv(path)
  list(
    filename = filename,
    found = !is.null(path),
    path = path %||% NA_character_,
    rows = if (!is.null(df)) nrow(df) else 0,
    cols = if (!is.null(df)) ncol(df) else 0,
    data = df,
    required = required
  )
}

bind_tissue_file <- function(filename, tissue_label, search_dirs) {
  out <- load_csv_with_meta(filename, search_dirs)
  df <- out$data
  if (is.null(df)) return(NULL)
  df <- clean_names_local(df)
  df$tissue <- tissue_label
  df
}

# -----------------------------
# 4) Data loaders
# -----------------------------
coalesce_step2_data <- function(search_dirs) {
  skin_split <- bind_tissue_file("step2_skin_peptides_with_species.csv", "Skin", search_dirs)
  bone_split <- bind_tissue_file("step2_bone_peptides_with_species.csv", "Bone", search_dirs)
  
  if (!is.null(skin_split) || !is.null(bone_split)) {
    return(bind_rows(skin_split, bone_split))
  }
  
  combined_meta <- load_csv_with_meta("step2_peptides_with_species.csv", search_dirs, required = TRUE)
  ensure_tissue_column(combined_meta$data)
}

coalesce_step3_species_summary <- function(search_dirs) {
  skin_split <- bind_tissue_file("step3_skin_species_summary.csv", "Skin", search_dirs)
  bone_split <- bind_tissue_file("step3_bone_species_summary.csv", "Bone", search_dirs)
  
  if (!is.null(skin_split) || !is.null(bone_split)) {
    return(clean_names_local(bind_rows(skin_split, bone_split)))
  }
  
  combined_meta <- load_csv_with_meta("step3_species_summary.csv", search_dirs, required = TRUE)
  ensure_tissue_column(combined_meta$data)
}

coalesce_step3_unique <- function(search_dirs) {
  skin_split <- bind_tissue_file("step3_skin_unique_peptides_scored.csv", "Skin", search_dirs)
  bone_split <- bind_tissue_file("step3_bone_unique_peptides_scored.csv", "Bone", search_dirs)
  
  if (!is.null(skin_split) || !is.null(bone_split)) {
    return(clean_names_local(bind_rows(skin_split, bone_split)))
  }
  
  combined_meta <- load_csv_with_meta("step3_unique_peptides_scored.csv", search_dirs, required = TRUE)
  ensure_tissue_column(combined_meta$data)
}

build_diagnostics_table <- function(search_dirs) {
  files <- c(
    "step2_peptides_with_species.csv",
    "step3_species_summary.csv",
    "step3_unique_peptides_scored.csv",
    "step2_skin_peptides_with_species.csv",
    "step2_bone_peptides_with_species.csv",
    "step3_skin_species_summary.csv",
    "step3_bone_species_summary.csv",
    "step3_skin_unique_peptides_scored.csv",
    "step3_bone_unique_peptides_scored.csv"
  )
  
  purrr::map_dfr(files, function(f) {
    meta <- load_csv_with_meta(f, search_dirs)
    tibble(
      filename = meta$filename,
      found = meta$found,
      path = meta$path,
      rows = meta$rows,
      cols = meta$cols
    )
  })
}

# -----------------------------
# 5) Initial data load
# -----------------------------
DEFAULT_DATA_DIR <- "/content"
SEARCH_DIRS <- build_search_dirs(DEFAULT_DATA_DIR)

step2_data_raw <- coalesce_step2_data(SEARCH_DIRS)
step3_species_raw <- coalesce_step3_species_summary(SEARCH_DIRS)
step3_unique_raw <- coalesce_step3_unique(SEARCH_DIRS)
diagnostics_raw <- build_diagnostics_table(SEARCH_DIRS)

# -----------------------------
# 6) Images
# -----------------------------
HOME_IMG_1 <- "https://raw.githubusercontent.com/hagarelazab/ANUBUIS/main/anubis-egyptian-god-Google-Search-03-08-2026_04_02_PM.png"
HOME_IMG_2 <- "https://raw.githubusercontent.com/hagarelazab/ANUBUIS/main/Ancient-Civilizations-Ancient-Egypt-03-08-2026_03_36_PM%20(1).png"

TEAM_IMG_1 <- "https://raw.githubusercontent.com/hagarelazab/ANUBUIS/main/WhatsApp%20Image%202026-03-08%20at%2022.00.33.jpeg"
TEAM_IMG_2 <- "https://raw.githubusercontent.com/hagarelazab/ANUBUIS/main/WhatsApp%20Image%202026-03-08%20at%2021.59.47.jpeg"

# -----------------------------
# 7) Theme
# -----------------------------
premium_theme <- bs_theme(
  version = 5,
  bg = "#f7f8fb",
  fg = "#101828",
  primary = "#0f172a",
  secondary = "#475467",
  success = "#12715b",
  warning = "#b54708",
  base_font = font_google("Inter"),
  heading_font = font_google("Cormorant Garamond")
)

# -----------------------------
# 8) Footer (fixed)
# -----------------------------
app_footer <- tags$footer(
  class = "app-footer",
  div(
    class = "footer-inner",
    div(class = "footer-accent"),
    div(class = "footer-text", HTML("&copy; 2026 Hagar Elazab. All Rights Reserved."))
  )
)

# -----------------------------
# 9) UI
# -----------------------------
ui <- page_navbar(
  theme = premium_theme,
  title = div(class = "brand-wrap", span(class = "brand-title", "Anubis")),
  id = "main_nav",
  selected = "Home",
  
  header = tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      :root {
        --surface: #ffffff;
        --line: #e4e7ec;
        --text: #101828;
        --muted: #667085;
        --shadow: 0 10px 30px rgba(16,24,40,.06);
        --radius: 22px;
        --footer-h: 64px;
      }

      html, body { height: 100%; }

      body {
        background: linear-gradient(180deg, #f8fafc 0%, #f7f8fb 100%);
        padding-bottom: var(--footer-h) !important;
      }

      .navbar {
        background: linear-gradient(90deg, #0f172a 0%, #162033 100%) !important;
        border: 0 !important;
        box-shadow: 0 8px 22px rgba(2, 6, 23, .16);
      }
      .navbar-brand, .nav-link, .navbar-nav .nav-link { color: rgba(255,255,255,.92) !important; }
      .navbar-nav .nav-link.active { color: #ffffff !important; position: relative; font-weight: 600; }
      .navbar-nav .nav-link.active::after {
        content: ''; position: absolute; left: 12px; right: 12px; bottom: 5px;
        height: 2px; border-radius: 999px; background: linear-gradient(90deg, #d4af37, #f4d58d);
      }

      .brand-wrap { display: flex; align-items: center; gap: 10px; }
      .brand-title { font-family: 'Cormorant Garamond', serif; font-size: 2rem; line-height: 1; letter-spacing: .2px; color: #fff; }

      .app-shell { max-width: 1440px; margin: 0 auto; padding: 28px 22px 42px 22px; }

      .hero-minimal {
        background:
          radial-gradient(circle at top left, rgba(200,161,77,.18), transparent 35%),
          linear-gradient(135deg, #0f172a 0%, #111827 50%, #0b1220 100%);
        color: white; padding: 28px 30px; border-radius: 26px;
        box-shadow: 0 18px 40px rgba(15,23,42,.16);
        margin-bottom: 22px; border: 1px solid rgba(255,255,255,.05);
      }
      .hero-title { margin: 0; font-family: 'Cormorant Garamond', serif; font-size: 3.3rem; line-height: 1; letter-spacing: .4px; }
      .hero-sub { margin-top: 10px; font-size: 1.06rem; color: rgba(255,255,255,.80); max-width: 980px; }

      .tissue-switch { display: flex; gap: 10px; flex-wrap: wrap; margin: 18px 0 6px; }
      .tissue-pill {
        display: inline-flex; align-items: center; gap: 8px; padding: 9px 14px;
        border-radius: 999px; background: rgba(255,255,255,.08);
        border: 1px solid rgba(255,255,255,.1); color: rgba(255,255,255,.92);
        font-size: .92rem; font-weight: 600;
      }
      .tissue-dot { width: 8px; height: 8px; border-radius: 50%; background: #d4af37; }

      .metric-grid { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 18px; margin-top: 18px; }
      .metric-card {
        background: var(--surface); border: 1px solid var(--line); border-radius: 22px;
        padding: 18px 18px 16px 18px; box-shadow: var(--shadow);
        position: relative; overflow: hidden;
      }
      .metric-card::before {
        content: ''; position: absolute; left: 0; top: 0; width: 100%; height: 4px;
        background: linear-gradient(90deg, var(--accent), rgba(255,255,255,0));
      }
      .metric-label { color: var(--muted); font-size: .9rem; font-weight: 600; margin-bottom: 10px; }
      .metric-value { color: var(--text); font-size: 2rem; line-height: 1; font-weight: 800; letter-spacing: -.03em; }
      .metric-sub { margin-top: 9px; color: var(--muted); font-size: .88rem; }

      .section-head { margin: 4px 0 16px; }
      .section-head h2 { margin: 0; font-size: 1.55rem; color: var(--text); font-weight: 800; }
      .section-head p { margin: 6px 0 0; color: var(--muted); }

      .panel-card {
        background: var(--surface); border: 1px solid var(--line); border-radius: var(--radius);
        padding: 18px; box-shadow: var(--shadow); margin-bottom: 18px;
      }
      .split-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 18px; }

      .home-story-grid { display: grid; grid-template-columns: minmax(360px, 520px) 1fr; gap: 22px; margin-bottom: 22px; align-items: stretch; }
      .collage-card { background: var(--surface); border: 1px solid var(--line); border-radius: 24px; overflow: hidden; box-shadow: var(--shadow); position: relative; }
      .collage-wrap {
        display: grid; grid-template-columns: 1fr; grid-template-rows: 1.25fr .75fr;
        gap: 12px; padding: 12px;
        background:
          radial-gradient(circle at top left, rgba(200,161,77,.12), transparent 42%),
          linear-gradient(180deg, #0b1220 0%, #0f172a 100%);
      }
      .collage-hero { border-radius: 18px; overflow: hidden; position: relative; min-height: 320px; box-shadow: 0 14px 30px rgba(2, 6, 23, .22); }
      .collage-hero img { width: 100%; height: 100%; object-fit: cover; display: block; }
      .collage-row { display: grid; grid-template-columns: 1fr 1fr; gap: 12px; }
      .collage-small { border-radius: 18px; overflow: hidden; position: relative; min-height: 170px; border: 1px solid rgba(255,255,255,.08); }
      .collage-small img { width: 100%; height: 100%; object-fit: cover; display: block; }
      .collage-label {
        position: absolute; left: 14px; bottom: 14px; padding: 8px 12px; border-radius: 999px;
        background: rgba(15, 23, 42, .70); color: rgba(255,255,255,.92);
        font-size: .85rem; font-weight: 700; backdrop-filter: blur(8px); border: 1px solid rgba(255,255,255,.12);
      }

      .anubis-info-card { background: var(--surface); border: 1px solid var(--line); border-radius: 24px; padding: 24px; box-shadow: var(--shadow); }
      .anubis-kicker {
        display: inline-flex; align-items: center; gap: 8px; padding: 7px 12px;
        border-radius: 999px; background: #f8f3e6; color: #8a6116;
        font-size: .84rem; font-weight: 700; margin-bottom: 14px;
      }
      .anubis-kicker::before { content: ''; width: 8px; height: 8px; border-radius: 50%; background: #c8a14d; }
      .anubis-title { margin: 0 0 12px 0; font-family: 'Cormorant Garamond', serif; font-size: 2.2rem; line-height: 1; color: #101828; }
      .anubis-lead { color: #475467; font-size: 1.02rem; line-height: 1.8; margin-bottom: 16px; }
      .anubis-grid { display: grid; grid-template-columns: 1fr 1fr; gap: 14px; margin-top: 12px; }
      .anubis-mini { background: #fbfcfe; border: 1px solid #e8edf3; border-radius: 18px; padding: 16px; }
      .anubis-mini h4 { margin: 0 0 10px 0; font-size: 1rem; color: #101828; }
      .anubis-mini li { color: #475467; line-height: 1.7; }
      .subtle-badge-row { display: flex; gap: 10px; flex-wrap: wrap; margin-top: 16px; }
      .subtle-badge { background: #fff; border: 1px solid #e6e8ef; color: #344054; border-radius: 999px; padding: 8px 12px; font-size: .88rem; font-weight: 600; }
      .subtle-badge strong { color: #101828; }

      .zoms-card { background: linear-gradient(180deg, #ffffff 0%, #fbfcfe 100%); border: 1px solid #e6e8ef; border-radius: 24px; box-shadow: var(--shadow); padding: 22px 22px; margin: 18px 0 22px 0; position: relative; overflow: hidden; }
      .zoms-card::before{
        content:''; position:absolute; inset:-40px -40px auto auto; width:220px; height:220px; border-radius:999px;
        background: radial-gradient(circle at center, rgba(200,161,77,.20), rgba(200,161,77,0)); transform: rotate(18deg);
      }
      .zoms-kicker{
        display:inline-flex; align-items:center; gap:8px; padding:7px 12px; border-radius:999px;
        background:#eef2ff; color:#1e3a8a; font-size:.84rem; font-weight:800; margin-bottom:12px; position:relative;
      }
      .zoms-kicker::before{ content:''; width:8px;height:8px;border-radius:50%; background:#1e3a8a; }
      .zoms-title{ margin:0 0 10px 0; font-size:1.6rem; font-weight:900; color:#0f172a; }
      .zoms-text{ color:#475467; font-size:1.02rem; line-height:1.9; margin:0; position:relative; }
      .zoms-steps{ display:grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap:12px; margin-top:16px; position:relative; }
      .zoms-step{ background:#ffffff; border:1px solid #e6e8ef; border-radius:18px; padding:14px 14px; }
      .zoms-step strong{ color:#101828; }
      .zoms-step p{ margin:6px 0 0 0; color:#475467; line-height:1.7; }

      /* TEAM TAB */
      .team-shell{ max-width: 1100px; margin: 0 auto; padding: 28px 22px 42px 22px; }
      .team-hero{
        background:
          radial-gradient(circle at top left, rgba(200,161,77,.14), transparent 40%),
          linear-gradient(135deg, #0f172a 0%, #111827 55%, #0b1220 100%);
        color:#fff; border-radius: 26px; padding: 26px 28px;
        box-shadow: 0 18px 40px rgba(15,23,42,.14);
        border: 1px solid rgba(255,255,255,.06);
        margin-bottom: 18px;
      }
      .team-hero h1{ margin:0; font-family:'Cormorant Garamond', serif; font-size: 2.6rem; line-height:1; }
      .team-hero p{ margin: 10px 0 0 0; color: rgba(255,255,255,.82); max-width: 900px; line-height: 1.8; }

      .team-card{
        background: linear-gradient(180deg, #ffffff 0%, #fbfcfe 100%);
        border: 1px solid #e6e8ef; border-radius: 24px; box-shadow: var(--shadow);
        padding: 22px 22px; margin: 0 0 18px 0; position: relative; overflow: hidden;
      }
      .team-card::before{
        content:''; position:absolute; inset:-60px auto auto -60px; width:260px; height:260px; border-radius:999px;
        background: radial-gradient(circle at center, rgba(15,23,42,.08), rgba(15,23,42,0));
      }
      .team-grid{ display:grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap:18px; position:relative; }
      .team-member{
        background:#ffffff; border:1px solid #e6e8ef; border-radius:22px; padding:18px;
        display:flex; flex-direction:column; align-items:center; text-align:center;
        box-shadow: 0 10px 22px rgba(16,24,40,.06);
      }

      .avatar{
        width:160px; height:160px; border-radius:999px; overflow:hidden;
        background:#fff;
      }

      /* Gold frame for Hagar + tiny soft glow (VERY subtle) */
      .avatar.avatar-hagar{
        border: 3px solid rgba(200,161,77,.82);
        box-shadow:
          0 14px 28px rgba(2,6,23,.12),
          0 0 0 1px rgba(200,161,77,.18),
          0 0 24px rgba(200,161,77,.12);
      }

      /* White frame for Nawal (clean, not special) */
      .avatar.avatar-nawal{
        border: 3px solid #ffffff;
        box-shadow:
          0 14px 28px rgba(2,6,23,.12),
          0 0 0 1px rgba(228,231,236,1);
      }

      .avatar img{ width:100%; height:100%; object-fit:cover; display:block; }
      .avatar img.avatar-hagar-img{ object-position: 50% 18%; }
      .avatar img.avatar-nawal-img{ object-position: 50% 20%; }

      .member-name{ margin-top:12px; font-weight:900; color:#101828; letter-spacing:.2px; font-size:1.05rem; }
      .member-affil{ margin-top: 8px; color:#475467; font-size: .94rem; line-height: 1.6; max-width: 520px; }

      .dataTables_wrapper .dataTables_filter input { border-radius: 12px !important; border: 1px solid #d0d5dd !important; padding: 6px 10px !important; }
      table.dataTable { border-collapse: separate !important; border-spacing: 0 !important; }
      table.dataTable thead th { background: #f8fafc !important; color: #344054 !important; border-bottom: 1px solid #eaecf0 !important; font-weight: 700 !important; }

      /* FIXED footer */
      .app-footer{
        position: fixed; left: 0; right: 0; bottom: 0; z-index: 9999;
        height: var(--footer-h);
        display: flex; align-items: center; justify-content: center;
        pointer-events: none;
        background: linear-gradient(180deg, rgba(247,248,251,0) 0%, rgba(247,248,251,.92) 35%, rgba(247,248,251,1) 100%);
        border-top: 1px solid rgba(230,232,239,.85);
        backdrop-filter: blur(10px);
      }
      .footer-inner{ width: 100%; max-width: 1440px; padding: 0 22px; text-align: center; pointer-events: none; }
      .footer-accent{
        height: 2px; width: min(520px, 92%); margin: 0 auto 10px auto; border-radius: 999px;
        background: linear-gradient(90deg, rgba(15,23,42,0), rgba(200,161,77,.95), rgba(15,23,42,0));
      }
      .footer-text{ color: #475467; font-weight: 700; letter-spacing: .2px; line-height: 1; }

      @media (max-width: 1200px) { .home-story-grid { grid-template-columns: 1fr; } .collage-hero { min-height: 280px; } }
      @media (max-width: 1100px) { .metric-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); } .split-grid { grid-template-columns: 1fr; } .zoms-steps{ grid-template-columns: 1fr; } .team-grid{ grid-template-columns: 1fr; } }
      @media (max-width: 700px) { .metric-grid { grid-template-columns: 1fr; } .hero-title { font-size: 2.3rem; } .anubis-grid { grid-template-columns: 1fr; } .collage-row { grid-template-columns: 1fr; } .collage-small { min-height: 200px; } :root { --footer-h: 72px; } .avatar{ width:148px; height:148px; } }
    "))
  ),
  
  footer = app_footer,
  
  # -----------------------------
  # HOME
  # -----------------------------
  nav_panel(
    title = "Home",
    value = "Home",
    div(
      class = "app-shell",
      div(
        class = "hero-minimal",
        h1(class = "hero-title", "Anubis"),
        div(class = "hero-sub",
            "Ancient Egyptian Collagen DB — a dashboard named after the guardian of tombs, mummification, and the passage between worlds."),
        uiOutput("home_tissue_pills")
      ),
      
      div(
        class = "home-story-grid",
        div(class = "collage-card", uiOutput("home_collage_ui")),
        div(
          class = "anubis-info-card",
          div(class = "anubis-kicker", "Mythology context"),
          h2(class = "anubis-title", "Why Anubis?"),
          p(
            class = "anubis-lead",
            "Anubis (Greek name) is associated with death, mummification, tomb protection, and the journey through the Duat. His Egyptian name was Anpu (Inpu)."
          ),
          div(
            class = "anubis-grid",
            div(
              class = "anubis-mini",
              h4("The meaning of black"),
              tags$ul(
                tags$li("Rebirth and regeneration"),
                tags$li("The fertile soil of the Nile"),
                tags$li("The darkened color of a body after mummification")
              )
            ),
            div(
              class = "anubis-mini",
              h4("Main roles of Anubis"),
              tags$ul(
                tags$li("Protector of tombs"),
                tags$li("God of mummification"),
                tags$li("The weighing of the heart")
              )
            )
          ),
          div(
            class = "subtle-badge-row",
            div(class = "subtle-badge", HTML("<strong>Greek name:</strong> Anubis")),
            div(class = "subtle-badge", HTML("<strong>Egyptian name:</strong> Anpu / Inpu")),
            div(class = "subtle-badge", HTML("<strong>Realm:</strong> Duat"))
          )
        )
      ),
      
      div(
        class = "zoms-card",
        div(class = "zoms-kicker", "Method overview"),
        h3(class = "zoms-title", "ZooMS (Zooarchaeology by Mass Spectrometry)"),
        p(
          class = "zoms-text",
          "ZooMS (Zooarchaeology by Mass Spectrometry) is a biomolecular analytical technique that uses collagen type I peptide mass fingerprinting to identify animal species from archaeological and paleontological bone fragments. Collagen type I is highly abundant and structurally stable, allowing it to survive for thousands of years in buried remains. In ZooMS analysis, collagen is extracted from bone, enzymatically digested (commonly with trypsin), and analyzed using matrix-assisted laser desorption/ionization time-of-flight (MALDI-TOF) mass spectrometry. The resulting peptide mass spectra contain species-specific marker peptides that enable taxonomic identification, even when bones are highly fragmented and lack diagnostic morphology. ZooMS has become a powerful tool in archaeology, paleontology, and conservation science for identifying fragmented faunal remains, tracking ancient biodiversity, and studying human–animal interactions in the past."
        ),
        div(
          class = "zoms-steps",
          div(class = "zoms-step", strong("1) Extract collagen"),
              p("Collagen type I is isolated from archaeological bone (often fragmented remains).")),
          div(class = "zoms-step", strong("2) Digest & measure"),
              p("Enzymatic digestion (trypsin) + MALDI-TOF yields a peptide mass fingerprint.")),
          div(class = "zoms-step", strong("3) Identify species"),
              p("Species-specific marker peptides enable taxonomic identification without morphology."))
        )
      ),
      
      uiOutput("overview_metrics")
    )
  ),
  
  # -----------------------------
  # OVERVIEW
  # -----------------------------
  nav_panel(
    title = "Overview",
    value = "Overview",
    div(
      class = "app-shell",
      div(
        class = "split-grid",
        div(
          class = "panel-card",
          section_title("Species distribution", "How many unique species are represented in each tissue subset."),
          plotlyOutput("home_species_bar", height = "340px")
        ),
        div(
          class = "panel-card",
          section_title("Top tissues by unique peptides", "A quick look at total unique peptide signal by section."),
          plotlyOutput("home_unique_tissue_bar", height = "340px")
        )
      )
    )
  ),
  
  nav_panel(
    title = "Data Snapshot",
    value = "Data Snapshot",
    div(
      class = "app-shell",
      div(class = "panel-card",
          section_title("Combined snapshot", "The main combined Step 2 table, moved out of Home for a cleaner front page."),
          DTOutput("home_snapshot"))
    )
  ),
  
  nav_panel(
    title = "Proteins per Species",
    value = "Proteins per Species",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("Proteins per species", "Counts of distinct proteins by species, separated by tissue."),
        fluidRow(
          column(4, selectInput("pps_tissue", "Tissue", choices = c("All", "Skin", "Bone", "Combined"), selected = "All")),
          column(8, uiOutput("pps_species_hint"))
        ),
        plotlyOutput("proteins_species_plot", height = "500px"),
        br(),
        DTOutput("proteins_species_table")
      )
    )
  ),
  
  nav_panel(
    title = "Peptides per Species",
    value = "Peptides per Species",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("Peptides per species", "Counts of total peptide occurrences and unique peptides by species."),
        fluidRow(
          column(4, selectInput("peps_tissue", "Tissue", choices = c("All", "Skin", "Bone", "Combined"), selected = "All")),
          column(4, selectInput("peps_metric", "Metric", choices = c("Total peptide occurrences", "Unique peptides"), selected = "Unique peptides")),
          column(4, sliderInput("peps_top_n", "Top N species", min = 5, max = 30, value = 12, step = 1))
        ),
        plotlyOutput("peptides_species_plot", height = "500px"),
        br(),
        DTOutput("peptides_species_table")
      )
    )
  ),
  
  nav_panel(
    title = "Unique Peptides per Species",
    value = "Unique Peptides per Species",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("Species-specific peptides", "Peptides found in exactly one species, scored and explorable by tissue."),
        fluidRow(
          column(4, selectInput("ups_tissue", "Tissue", choices = c("All", "Skin", "Bone", "Combined"), selected = "All")),
          column(4, selectInput("ups_species", "Species", choices = "All", selected = "All")),
          column(4, sliderInput("ups_top_n", "Top N rows for chart", min = 5, max = 50, value = 15, step = 1))
        ),
        plotlyOutput("unique_peptides_species_plot", height = "500px"),
        br(),
        DTOutput("unique_peptides_species_table")
      )
    )
  ),
  
  nav_panel(
    title = "Peptide Existence",
    value = "Peptide Existence",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("Peptide confidence landscape", "Distribution of peptide evidence scores across species and tissues."),
        fluidRow(
          column(4, selectInput("pex_tissue", "Tissue", choices = c("All", "Skin", "Bone", "Combined"), selected = "All")),
          column(4, selectInput("pex_species", "Species", choices = "All", selected = "All")),
          column(4, selectInput("pex_kind", "View", choices = c("Histogram", "Box plot"), selected = "Histogram"))
        ),
        plotlyOutput("peptide_existence_plot", height = "500px"),
        br(),
        DTOutput("peptide_existence_table")
      )
    )
  ),
  
  nav_panel(
    title = "Existence per Species",
    value = "Existence per Species",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("Existence per species", "Average and total evidence scores at the species level."),
        fluidRow(
          column(4, selectInput("eps_tissue", "Tissue", choices = c("All", "Skin", "Bone", "Combined"), selected = "All")),
          column(4, selectInput("eps_metric", "Metric",
                                choices = c("sum_peptide_score", "mean_peptide_score", "unique_peptide_count"),
                                selected = "sum_peptide_score")),
          column(4, sliderInput("eps_top_n", "Top N species", min = 5, max = 30, value = 12, step = 1))
        ),
        plotlyOutput("existence_species_plot", height = "500px"),
        br(),
        DTOutput("existence_species_table")
      )
    )
  ),
  
  nav_panel(
    title = "Species Summary",
    value = "Species Summary",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("Species summary", "Full species-level summary table with tissue-aware filtering."),
        fluidRow(
          column(4, selectInput("s3sum_tissue", "Tissue", choices = c("All", "Skin", "Bone", "Combined"), selected = "All")),
          column(4, selectInput("s3sum_sort", "Sort by",
                                choices = c("unique_peptide_count", "sum_peptide_score", "mean_peptide_score"),
                                selected = "unique_peptide_count")),
          column(4, selectInput("s3sum_dir", "Direction", choices = c("Descending", "Ascending"), selected = "Descending"))
        ),
        DTOutput("step3_species_summary_table")
      )
    )
  ),
  
  nav_panel(
    title = "Unique Peptides",
    value = "Unique Peptides",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("Unique peptides", "Full peptide-level scored table with species and tissue filters."),
        fluidRow(
          column(4, selectInput("s3up_tissue", "Tissue", choices = c("All", "Skin", "Bone", "Combined"), selected = "All")),
          column(4, selectInput("s3up_species", "Species", choices = "All", selected = "All")),
          column(4, sliderInput("s3up_min_score", "Minimum peptide score", min = 1, max = 5, value = 1, step = 1))
        ),
        DTOutput("step3_unique_peptides_table")
      )
    )
  ),
  
  nav_panel(
    title = "Diagnostics",
    value = "Diagnostics",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title("File diagnostics", "This tab tells you exactly where the app searched and what it found."),
        h4("Search directories"),
        verbatimTextOutput("search_dirs_text"),
        br(),
        DTOutput("diagnostics_table")
      )
    )
  ),
  
  # -----------------------------
  # MEET OUR TEAM (LAST TAB)
  # -----------------------------
  nav_panel(
    title = "Meet our team",
    value = "Team",
    div(
      class = "team-shell",
      div(
        class = "team-hero",
        h1("Meet our team"),
        p("The people behind the Anubis dashboard and the Ancient Egyptian Collagen DB.")
      ),
      div(
        class = "team-card",
        div(
          class = "team-grid",
          div(
            class = "team-member",
            div(class = "avatar avatar-hagar",
                tags$img(src = TEAM_IMG_1, class = "avatar-hagar-img", alt = "Hagar El Azab")),
            div(class = "member-name", "Hagar El Azab"),
            div(class = "member-affil",
                HTML("Department of Biotechnology,<br>Faculty of Agriculture,<br>Cairo University, Giza, Egypt."))
          ),
          div(
            class = "team-member",
            div(class = "avatar avatar-nawal",
                tags$img(src = TEAM_IMG_2, class = "avatar-nawal-img", alt = "Nawal Hassan")),
            div(class = "member-name", "Nawal Hassan"),
            div(class = "member-affil",
                HTML("Undergraduate Studies at the Department of Biotechnology,<br>Faculty of Agriculture,<br>Ain Shams University, Cairo, Egypt."))
          )
        )
      )
    )
  )
)

# -----------------------------
# 10) Server
# -----------------------------
server <- function(input, output, session) {
  
  step2_data <- reactiveVal(step2_data_raw)
  step3_species_summary <- reactiveVal(step3_species_raw)
  step3_unique_peptides <- reactiveVal(step3_unique_raw)
  diagnostics_tbl <- reactiveVal(diagnostics_raw)
  
  step2_tbl <- reactive(step2_data())
  s3_species_tbl <- reactive(step3_species_summary())
  s3_unique_tbl <- reactive(step3_unique_peptides())
  
  available_tissues <- reactive({
    vals <- c(
      if (!is.null(step2_tbl()) && "tissue" %in% names(step2_tbl())) unique(step2_tbl()$tissue) else NULL,
      if (!is.null(s3_species_tbl()) && "tissue" %in% names(s3_species_tbl())) unique(s3_species_tbl()$tissue) else NULL,
      if (!is.null(s3_unique_tbl()) && "tissue" %in% names(s3_unique_tbl())) unique(s3_unique_tbl()$tissue) else NULL
    )
    vals <- unique(na.omit(vals))
    vals <- vals[vals != ""]
    vals <- intersect(c("Skin", "Bone", "Combined"), vals)
    if (length(vals) == 0) vals <- "Combined"
    vals
  })
  
  filter_by_tissue <- function(df, tissue_choice) {
    if (is.null(df)) return(NULL)
    if (!"tissue" %in% names(df)) return(df)
    if (is.null(tissue_choice) || tissue_choice == "All") return(df)
    df %>% filter(tissue == tissue_choice)
  }
  
  observe({
    tissues <- c("All", available_tissues())
    updateSelectInput(session, "pps_tissue", choices = tissues, selected = "All")
    updateSelectInput(session, "peps_tissue", choices = tissues, selected = "All")
    updateSelectInput(session, "ups_tissue", choices = tissues, selected = "All")
    updateSelectInput(session, "pex_tissue", choices = tissues, selected = "All")
    updateSelectInput(session, "eps_tissue", choices = tissues, selected = "All")
    updateSelectInput(session, "s3sum_tissue", choices = tissues, selected = "All")
    updateSelectInput(session, "s3up_tissue", choices = tissues, selected = "All")
  })
  
  observe({
    df <- s3_unique_tbl()
    df <- filter_by_tissue(df, input$ups_tissue)
    species_choices <- if (is.null(df) || !"species" %in% names(df)) "All" else c("All", sort(unique(df$species)))
    updateSelectInput(session, "ups_species", choices = species_choices, selected = "All")
    updateSelectInput(session, "pex_species", choices = species_choices, selected = "All")
    updateSelectInput(session, "s3up_species", choices = species_choices, selected = "All")
  })
  
  output$home_tissue_pills <- renderUI({
    pills <- available_tissues()
    div(class = "tissue-switch", lapply(pills, function(x) div(class = "tissue-pill", span(class = "tissue-dot"), x)))
  })
  
  output$overview_metrics <- renderUI({
    df2 <- step2_tbl()
    proteins <- safe_n_distinct(df2, "protein_id")
    peptides_total <- if (!is.null(df2) && "peptide" %in% names(df2)) nrow(df2) else 0
    unique_peptides <- safe_n_distinct(df2, "peptide")
    species_total <- safe_n_distinct(df2, "species")
    div(
      class = "metric-grid",
      metric_box("Proteins", format_big(proteins), "Distinct protein IDs loaded", "#0f172a"),
      metric_box("Peptide occurrences", format_big(peptides_total), "All rows in Step 2", "#2563eb"),
      metric_box("Unique peptides", format_big(unique_peptides), "Deduplicated peptide sequences", "#0f766e"),
      metric_box("Species", format_big(species_total), paste0("Tissues available: ", paste(available_tissues(), collapse = ", ")), "#c0841a")
    )
  })
  
  output$home_collage_ui <- renderUI({
    div(
      class = "collage-wrap",
      div(class = "collage-hero", tags$img(src = HOME_IMG_1, alt = "Anubis"), div(class = "collage-label", "Anubis")),
      div(
        class = "collage-row",
        div(class = "collage-small", tags$img(src = HOME_IMG_2, alt = "Ancient Egypt"), div(class = "collage-label", "Ancient Egypt")),
        div(
          class = "collage-small",
          div(
            style = paste0(
              "width:100%;height:100%;display:flex;flex-direction:column;justify-content:center;",
              "padding:18px;background:rgba(255,255,255,.06);color:rgba(255,255,255,.92);",
              "backdrop-filter: blur(8px);"
            ),
            div(style = "font-weight:800;font-size:1.05rem;margin-bottom:8px;", "Explore the DB"),
            div(style = "color:rgba(255,255,255,.78);line-height:1.7;",
                "Use the tabs to explore proteins, peptides, species summaries, and tissue subsets.")
          )
        )
      )
    )
  })
  
  output$home_species_bar <- renderPlotly({
    df <- step2_tbl()
    validate(need(!is.null(df), "No Step 2 data found."))
    validate(need(all(c("species", "tissue") %in% names(df)), "Step 2 data must contain species and tissue."))
    plot_df <- df %>%
      distinct(species, tissue) %>%
      count(tissue, name = "n_species") %>%
      mutate(tissue = factor(tissue, levels = c("Skin", "Bone", "Combined")),
             tooltip = paste0("Tissue: ", tissue, "<br>Species: ", scales::comma(n_species)))
    g <- ggplot(plot_df, aes(x = tissue, y = n_species, fill = tissue, text = tooltip)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = scales::comma(n_species)), vjust = -0.4, color = "#344054") +
      scale_fill_manual(values = plot_palette_tissue, drop = FALSE) +
      labs(x = NULL, y = "Distinct species") +
      expand_limits(y = max(plot_df$n_species, na.rm = TRUE) * 1.12) +
      base_gg_theme()
    ggplotly_clean(g)
  })
  
  output$home_unique_tissue_bar <- renderPlotly({
    df <- s3_species_tbl()
    validate(need(!is.null(df), "No species summary found."))
    validate(need(all(c("tissue", "unique_peptide_count") %in% names(df)), "Species summary must contain tissue and unique_peptide_count."))
    plot_df <- df %>%
      group_by(tissue) %>%
      summarise(unique_peptides = sum(unique_peptide_count, na.rm = TRUE), .groups = "drop") %>%
      mutate(tissue = factor(tissue, levels = c("Skin", "Bone", "Combined")),
             tooltip = paste0("Tissue: ", tissue, "<br>Unique peptides: ", scales::comma(unique_peptides)))
    g <- ggplot(plot_df, aes(x = tissue, y = unique_peptides, fill = tissue, text = tooltip)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = scales::comma(unique_peptides)), vjust = -0.4, color = "#344054") +
      scale_fill_manual(values = plot_palette_tissue, drop = FALSE) +
      labs(x = NULL, y = "Unique peptides") +
      expand_limits(y = max(plot_df$unique_peptides, na.rm = TRUE) * 1.12) +
      base_gg_theme()
    ggplotly_clean(g)
  })
  
  output$home_snapshot <- renderDT({
    df <- step2_tbl()
    if (is.null(df)) return(empty_dt("No Step 2 data available."))
    keep <- intersect(c("protein_id", "peptide", "pep_len", "species", "taxid", "tissue"), names(df))
    datatable(df[, keep, drop = FALSE], rownames = FALSE, filter = "top",
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  proteins_by_species <- reactive({
    df <- step2_tbl()
    validate(need(!is.null(df), "No Step 2 data."))
    validate(need(all(c("species", "protein_id") %in% names(df)), "Need species and protein_id in Step 2 data."))
    df <- filter_by_tissue(df, input$pps_tissue)
    grp <- if ("tissue" %in% names(df)) c("tissue", "species") else c("species")
    df %>%
      group_by(across(all_of(grp))) %>%
      summarise(n_proteins = n_distinct(protein_id), .groups = "drop") %>%
      arrange(desc(n_proteins), species)
  })
  
  output$pps_species_hint <- renderUI({
    df <- proteins_by_species()
    div(
      class = "subtle-badge-row",
      div(class = "subtle-badge", HTML(paste0("<strong>", nrow(df), "</strong> species rows"))),
      div(class = "subtle-badge", HTML(paste0("<strong>", length(unique(df$species)), "</strong> distinct species")))
    )
  })
  
  output$proteins_species_plot <- renderPlotly({
    df <- proteins_by_species()
    validate(need(nrow(df) > 0, "No rows to plot."))
    if (!"tissue" %in% names(df)) df$tissue <- "Combined"
    plot_df <- df %>% mutate(
      species = reorder(species, n_proteins),
      tooltip = paste0("Species: ", species, "<br>Proteins: ", scales::comma(n_proteins), "<br>Tissue: ", tissue)
    )
    g <- ggplot(plot_df, aes(x = species, y = n_proteins, fill = species, text = tooltip)) +
      geom_col(width = 0.75) +
      scale_fill_manual(values = rep(plot_palette_main, length.out = nrow(plot_df))) +
      labs(x = NULL, y = "Distinct proteins") +
      base_gg_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ggplotly_clean(g)
  })
  
  output$proteins_species_table <- renderDT({
    datatable(proteins_by_species(), rownames = FALSE, filter = "top",
              options = list(pageLength = 12, scrollX = TRUE))
  })
  
  peptides_by_species <- reactive({
    df <- step2_tbl()
    validate(need(!is.null(df), "No Step 2 data."))
    validate(need(all(c("species", "peptide") %in% names(df)), "Need species and peptide in Step 2 data."))
    df <- filter_by_tissue(df, input$peps_tissue)
    grp <- if ("tissue" %in% names(df)) c("tissue", "species") else c("species")
    df %>%
      group_by(across(all_of(grp))) %>%
      summarise(total_peptide_occurrences = n(), unique_peptides = n_distinct(peptide), .groups = "drop") %>%
      arrange(desc(unique_peptides), species)
  })
  
  output$peptides_species_plot <- renderPlotly({
    df <- peptides_by_species()
    validate(need(nrow(df) > 0, "No rows to plot."))
    metric_col <- if (identical(input$peps_metric, "Total peptide occurrences")) "total_peptide_occurrences" else "unique_peptides"
    plot_df <- df %>%
      arrange(desc(.data[[metric_col]])) %>%
      slice_head(n = input$peps_top_n) %>%
      mutate(
        species = reorder(species, .data[[metric_col]]),
        tooltip = paste0("Species: ", species, "<br>", metric_col, ": ", scales::comma(.data[[metric_col]]), "<br>Tissue: ", tissue)
      )
    g <- ggplot(plot_df, aes(x = species, y = .data[[metric_col]], fill = species, text = tooltip)) +
      geom_col(width = 0.75) +
      scale_fill_manual(values = rep(plot_palette_main, length.out = nrow(plot_df))) +
      labs(x = NULL, y = metric_col) +
      base_gg_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ggplotly_clean(g)
  })
  
  output$peptides_species_table <- renderDT({
    datatable(peptides_by_species(), rownames = FALSE, filter = "top",
              options = list(pageLength = 12, scrollX = TRUE))
  })
  
  unique_species_filtered <- reactive({
    df <- s3_unique_tbl()
    validate(need(!is.null(df), "No unique peptide data."))
    df <- filter_by_tissue(df, input$ups_tissue)
    if (!is.null(input$ups_species) && input$ups_species != "All" && "species" %in% names(df)) df <- df %>% filter(species == input$ups_species)
    df
  })
  
  output$unique_peptides_species_plot <- renderPlotly({
    df <- unique_species_filtered()
    validate(need(nrow(df) > 0, "No rows to plot."))
    validate(need(all(c("species", "peptide") %in% names(df)), "Need species and peptide."))
    plot_df <- df %>%
      group_by(tissue, species) %>%
      summarise(n_unique_peptides = n_distinct(peptide), .groups = "drop") %>%
      arrange(desc(n_unique_peptides)) %>%
      slice_head(n = input$ups_top_n) %>%
      mutate(
        species = reorder(species, n_unique_peptides),
        tooltip = paste0("Species: ", species, "<br>Unique peptides: ", scales::comma(n_unique_peptides), "<br>Tissue: ", tissue)
      )
    g <- ggplot(plot_df, aes(x = species, y = n_unique_peptides, fill = species, text = tooltip)) +
      geom_col(width = 0.75) +
      scale_fill_manual(values = rep(plot_palette_main, length.out = nrow(plot_df))) +
      labs(x = NULL, y = "Species-specific peptides") +
      base_gg_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    ggplotly_clean(g)
  })
  
  output$unique_peptides_species_table <- renderDT({
    df <- unique_species_filtered()
    if (is.null(df) || nrow(df) == 0) return(empty_dt("No matching unique peptides."))
    keep <- intersect(c("peptide", "species", "tissue", "peptide_score", "n_proteins"), names(df))
    datatable(df[, keep, drop = FALSE], rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })
  
  peptide_existence_filtered <- reactive({
    df <- s3_unique_tbl()
    validate(need(!is.null(df), "No unique peptide data."))
    validate(need("peptide_score" %in% names(df), "peptide_score column is missing."))
    df <- filter_by_tissue(df, input$pex_tissue)
    if (!is.null(input$pex_species) && input$pex_species != "All" && "species" %in% names(df)) df <- df %>% filter(species == input$pex_species)
    df
  })
  
  output$peptide_existence_plot <- renderPlotly({
    df <- peptide_existence_filtered()
    validate(need(nrow(df) > 0, "No rows to plot."))
    if (!"tissue" %in% names(df)) df$tissue <- "Combined"
    
    if (identical(input$pex_kind, "Histogram")) {
      plot_df <- df %>% mutate(
        tissue = factor(tissue, levels = c("Skin", "Bone", "Combined")),
        tooltip = paste0("Score: ", peptide_score, "<br>Tissue: ", tissue)
      )
      g <- ggplot(plot_df, aes(x = peptide_score, fill = tissue, text = tooltip)) +
        geom_histogram(position = "identity", alpha = 0.7, bins = 20) +
        scale_fill_manual(values = plot_palette_tissue, drop = FALSE) +
        labs(x = "Peptide score", y = "Count") +
        base_gg_theme()
      ggplotly_clean(g)
    } else {
      plot_df <- df %>% mutate(tooltip = paste0("Species: ", species, "<br>Score: ", peptide_score, "<br>Tissue: ", tissue))
      g <- ggplot(plot_df, aes(x = species, y = peptide_score, fill = species, text = tooltip)) +
        geom_boxplot(outlier_alpha = 0.7) +
        scale_fill_manual(values = rep(plot_palette_main, length.out = dplyr::n_distinct(plot_df$species))) +
        labs(x = NULL, y = "Peptide score") +
        base_gg_theme() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
      ggplotly_clean(g)
    }
  })
  
  output$peptide_existence_table <- renderDT({
    df <- peptide_existence_filtered()
    keep <- intersect(c("peptide", "species", "tissue", "peptide_score", "n_proteins"), names(df))
    datatable(df[, keep, drop = FALSE], rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })
  
  # ---- IMPORTANT CHANGE: mean_peptide_score rounded to 2 decimals ----
  existence_species_filtered <- reactive({
    df <- s3_species_tbl()
    validate(need(!is.null(df), "No species summary data."))
    validate(need(all(c("species", "unique_peptide_count", "sum_peptide_score", "mean_peptide_score") %in% names(df)),
                  "Need species, unique_peptide_count, sum_peptide_score, mean_peptide_score."))
    df <- filter_by_tissue(df, input$eps_tissue)
    if (!"tissue" %in% names(df)) df$tissue <- "Combined"
    
    # round mean_peptide_score to 2 decimals (and keep numeric)
    df$mean_peptide_score <- round(as.numeric(df$mean_peptide_score), 2)
    
    df
  })
  
  output$existence_species_plot <- renderPlotly({
    df <- existence_species_filtered()
    validate(need(nrow(df) > 0, "No rows to plot."))
    metric_col <- input$eps_metric
    
    plot_df <- df %>%
      arrange(desc(.data[[metric_col]])) %>%
      slice_head(n = input$eps_top_n) %>%
      mutate(
        species = reorder(species, .data[[metric_col]]),
        tooltip = paste0(
          "Species: ", species,
          "<br>", metric_col, ": ",
          if (metric_col == "mean_peptide_score") sprintf("%.2f", .data[[metric_col]]) else scales::comma(.data[[metric_col]]),
          "<br>Tissue: ", tissue
        )
      )
    
    g <- ggplot(plot_df, aes(x = species, y = .data[[metric_col]], fill = species, text = tooltip)) +
      geom_col(width = 0.75) +
      scale_fill_manual(values = rep(plot_palette_main, length.out = nrow(plot_df))) +
      labs(x = NULL, y = metric_col) +
      base_gg_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
    
    ggplotly_clean(g)
  })
  
  output$existence_species_table <- renderDT({
    df <- existence_species_filtered()
    if ("mean_peptide_score" %in% names(df)) df$mean_peptide_score <- sprintf("%.2f", df$mean_peptide_score)
    datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 12, scrollX = TRUE))
  })
  
  output$step3_species_summary_table <- renderDT({
    df <- s3_species_tbl()
    if (is.null(df)) return(empty_dt("No species summary available."))
    df <- filter_by_tissue(df, input$s3sum_tissue)
    
    if ("mean_peptide_score" %in% names(df)) df$mean_peptide_score <- round(as.numeric(df$mean_peptide_score), 2)
    
    dir_desc <- identical(input$s3sum_dir, "Descending")
    sort_col <- input$s3sum_sort
    if (sort_col %in% names(df)) df <- df %>% arrange(if (dir_desc) desc(.data[[sort_col]]) else .data[[sort_col]])
    
    if ("mean_peptide_score" %in% names(df)) df$mean_peptide_score <- sprintf("%.2f", df$mean_peptide_score)
    
    datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })
  
  output$step3_unique_peptides_table <- renderDT({
    df <- s3_unique_tbl()
    if (is.null(df)) return(empty_dt("No unique peptides available."))
    df <- filter_by_tissue(df, input$s3up_tissue)
    if (!is.null(input$s3up_species) && input$s3up_species != "All" && "species" %in% names(df)) df <- df %>% filter(species == input$s3up_species)
    if ("peptide_score" %in% names(df)) df <- df %>% filter(peptide_score >= input$s3up_min_score)
    datatable(df, rownames = FALSE, filter = "top", options = list(pageLength = 15, scrollX = TRUE))
  })
  
  output$search_dirs_text <- renderText({ paste(SEARCH_DIRS, collapse = "\n") })
  
  output$diagnostics_table <- renderDT({
    datatable(diagnostics_tbl(), rownames = FALSE, filter = "top", options = list(pageLength = 20, scrollX = TRUE))
  })
}

# -----------------------------
# 11) Run app
# -----------------------------
shinyApp(ui, server)