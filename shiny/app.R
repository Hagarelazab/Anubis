# Anubis - Ancient Egyptian Collagen DB
# Updated runnable version generated from the full uploaded app.R reference.
# Place this app.R in the same folder as your CSV/XML data files, or put the data files in a subfolder named data/.

required_packages <- c(
  "shiny", "bslib", "dplyr", "readr", "DT", "plotly", "ggplot2",
  "stringr", "purrr", "scales", "tidyr", "htmltools", "tibble", "xml2"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Please install the missing R packages before running the app: ",
    paste(missing_packages, collapse = ", "),
    "\nRun this in R: install.packages(c(",
    paste(sprintf('\"%s\"', missing_packages), collapse = ", "),
    "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
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
  library(xml2)
})

options(shiny.maxRequestSize = 1024 * 1024^2)

# Operators and Helpers
`%||%` <- function(a, b) if (!is.null(a)) a else b

first_non_empty <- function(...) {
  vals <- list(...)
  for (val in vals) {
    if (is.null(val)) next
    x <- trimws(as.character(val)[1] %||% "")
    if (nzchar(x)) return(x)
  }
  ""
}

read_taxonomy_label_from_xml <- function(path) {
  if (is.null(path) || !file.exists(path)) return("")
  doc <- tryCatch(xml2::read_xml(path), error = function(e) NULL)
  if (is.null(doc)) return("")
  labels <- tryCatch(xml2::xml_attr(xml2::xml_find_all(doc, ".//taxon"), "label"), error = function(e) character(0))
  labels <- trimws(as.character(labels))
  labels <- labels[nzchar(labels)]
  if (length(labels) == 0) "" else labels[1]
}


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
  # Use nested calls instead of the base pipe |> so the app works on older R versions.
  names(df) <- tolower(
    stringr::str_replace_all(
      stringr::str_replace_all(names(df), "[^A-Za-z0-9]+", "_"),
      "(^_+|_+$)",
      ""
    )
  )
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

round2 <- function(x) {
  if (is.null(x)) return(x)
  suppressWarnings({
    y <- as.numeric(x)
    ifelse(is.na(y), x, round(y, 2))
  })
}

# Robust E-value parser. Keeps scientific notation correctly (e.g. 4.4e-05),
# removes harmless text symbols such as <, >, commas, and whitespace,
# and returns numeric values for correct sorting/ranking.
parse_evalue <- function(x) {
  if (is.null(x)) return(NA_real_)
  x_chr <- trimws(as.character(x))
  x_chr <- gsub(",", "", x_chr, fixed = TRUE)
  x_chr <- gsub("<", "", x_chr, fixed = TRUE)
  x_chr <- gsub(">", "", x_chr, fixed = TRUE)
  x_chr <- gsub("^[Ee]=", "", x_chr)
  x_chr[x_chr %in% c("", "NA", "NaN", "NULL", "Inf", "-Inf")] <- NA_character_
  suppressWarnings(as.numeric(x_chr))
}

format_evalue <- function(x, digits = 3) {
  y <- parse_evalue(x)
  ifelse(is.na(y), NA_character_, formatC(y, format = "e", digits = digits))
}

# DT renderer for E-values.
# Important: keeps the underlying value numeric for correct sorting,
# but displays it in scientific notation so 4.4 becomes 4.40e+00
# and 0.046 becomes 4.60e-02.
evalue_column_defs <- function(target_index, digits = 2) {
  list(list(
    targets = target_index,
    render = DT::JS(sprintf(
      "function(data, type, row, meta) { if (type === 'display' || type === 'filter') { if (data === null || data === '' || isNaN(Number(data))) { return ''; } return Number(data).toExponential(%d); } return data; }",
      as.integer(digits)
    ))
  ))
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



table_mini_card <- function(label, value, subtitle = NULL) {
  div(
    class = "table-min i-card",
    div(class = "table-mini-label", label),
    div(class = "table-mini-value", value),
    if (!is.null(subtitle)) div(class = "table-mini-sub", subtitle)
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

# Plot helpers
plot_palette_main <- c(
  "#0f172a", "#175cd3", "#0f766e", "#c0841a",
  "#7c3aed", "#dc2626", "#0891b2", "#16a34a",
  "#d97706", "#db2777", "#B68D40", "#475569"
)

plot_palette_tissue <- c(
  "Skin" = "#15803d",
  "Bone" = "#b7791f",
  "Combined" = "#B68D40"
)

# Premium chart helpers
# These functions only change the visual design of the charts.
# They do not change any biological/statistical calculations.
make_plot_ymax <- function(x, pad = 1.22) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0 || max(x, na.rm = TRUE) <= 0) return(1)
  max(x, na.rm = TRUE) * pad
}

base_gg_theme <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.margin = margin(t = 4, r = 0, b = 0, l = 0),
      legend.text = element_text(color = "#344054", size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#eef2f7", linewidth = 0.45),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(color = "#344054", size = 12, face = "bold", margin = margin(t = 12)),
      axis.title.y = element_text(color = "#344054", size = 12, face = "bold", margin = margin(r = 12)),
      axis.text.x = element_text(color = "#344054", size = 11, face = "bold"),
      axis.text.y = element_text(color = "#475467", size = 10),
      plot.title = element_text(face = "bold", color = "#101828"),
      plot.margin = margin(t = 16, r = 18, b = 12, l = 12)
    )
}

ggplotly_clean <- function(p, tooltip = "text") {
  ggplotly(p, tooltip = tooltip) %>%
    layout(
      paper_bgcolor = "rgba(0,0,0,0)",
      plot_bgcolor = "rgba(0,0,0,0)",
      hoverlabel = list(
        bgcolor = "#0f172a",
        bordercolor = "#0f172a",
        font = list(color = "#ffffff", size = 12)
      ),
      legend = list(
        orientation = "h",
        y = -0.18,
        x = 0.5,
        xanchor = "center",
        font = list(color = "#344054", size = 12)
      ),
      margin = list(l = 58, r = 28, t = 22, b = 72)
    ) %>%
    config(displayModeBar = FALSE, responsive = TRUE)
}

# File discovery helpers
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg_name <- "--file="
  script_path <- sub(file_arg_name, "", args[grep(file_arg_name, args)])
  if (length(script_path) > 0 && nzchar(script_path[1])) {
    return(dirname(normalizePath(script_path[1], winslash = "/", mustWork = FALSE)))
  }

  frame_files <- vapply(sys.frames(), function(x) {
    ofile <- x$ofile
    if (is.null(ofile)) "" else as.character(ofile)[1]
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[length(frame_files)], winslash = "/", mustWork = FALSE)))
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

build_search_dirs <- function(user_dir = NULL) {
  app_dir <- get_script_dir()
  wd <- getwd()
  
  dirs <- c(
    user_dir,
    Sys.getenv("DATA_DIR", unset = NA),
    Sys.getenv("APP_DATA_DIR", unset = NA),
    wd,
    app_dir,
    file.path(wd, "data"),
    file.path(app_dir, "data"),
    file.path(dirname(app_dir), "data"),
    path.expand("~/xtandem_test"),
    file.path(path.expand("~/xtandem_test"), "data"),
    path.expand("~/Downloads"),
    file.path(path.expand("~/Downloads"), "Shiny"),
    "/content",
    "/content/data"
  )
  
  dirs <- dirs[!is.na(dirs) & nzchar(dirs)]
  unique(normalizePath(dirs, winslash = "/", mustWork = FALSE))
}

find_existing_file <- function(filename, search_dirs) {
  candidates <- unique(c(filename, file.path(search_dirs, filename)))
  exists <- candidates[file.exists(candidates)]
  if (length(exists) > 0) {
    return(normalizePath(exists[1], winslash = "/", mustWork = FALSE))
  }
  
  recursive_dirs <- search_dirs[dir.exists(search_dirs)]
  found <- unlist(lapply(recursive_dirs, function(d) {
    tryCatch(
      list.files(
        d,
        pattern = paste0("^", stringr::str_replace_all(filename, "([\\^\\$\\.\\|\\?\\*\\+\\(\\)\\[\\]\\{\\}])", "\\\\\\1"), "$"),
        recursive = TRUE,
        full.names = TRUE
      ),
      error = function(e) character(0)
    )
  }), use.names = FALSE)
  
  found <- found[file.exists(found)]
  if (length(found) == 0) return(NULL)
  normalizePath(found[1], winslash = "/", mustWork = FALSE)
}

find_existing_files_by_patterns <- function(patterns, search_dirs) {
  recursive_dirs <- unique(search_dirs[dir.exists(search_dirs)])
  if (length(recursive_dirs) == 0) return(character(0))
  
  found <- unlist(lapply(recursive_dirs, function(d) {
    unlist(lapply(patterns, function(p) {
      tryCatch(
        list.files(
          d,
          pattern = p,
          recursive = TRUE,
          full.names = TRUE,
          ignore.case = TRUE
        ),
        error = function(e) character(0)
      )
    }), use.names = FALSE)
  }), use.names = FALSE)
  
  found <- unique(found[file.exists(found)])
  if (length(found) == 0) return(character(0))
  normalizePath(found, winslash = "/", mustWork = FALSE)
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

# Data loaders
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

# Spectra identification helpers
# FIX 2026-05-16:
# X! Tandem labels do not always include UniProt OS=... species text.
# Some labels look like CO4A1_MOUSE Collagen... or sp|P02452|CO1A1_HUMAN ...
# The old function returned NA for those labels, so species ranking became empty.
# This version supports both OS=... and UniProt mnemonic suffixes such as _HUMAN/_MOUSE.
uniprot_species_map <- c(
  HUMAN = "Homo sapiens",
  MOUSE = "Mus musculus",
  RAT   = "Rattus norvegicus",
  BOVIN = "Bos taurus",
  SHEEP = "Ovis aries",
  PIG   = "Sus scrofa",
  HORSE = "Equus caballus",
  CHICK = "Gallus gallus",
  CANLF = "Canis lupus familiaris",
  FELCA = "Felis catus",
  RABIT = "Oryctolagus cuniculus",
  CAPHI = "Capra hircus",
  CAMEL = "Camelus",
  DROME = "Camelus dromedarius"
)

extract_uniprot_species_code <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(NA_character_)
  y <- trimws(as.character(x)[1])

  # Keep the UniProt mnemonic token, whether the label is:
  # sp|P02452|CO1A1_HUMAN Collagen... OR CO1A1_HUMAN Collagen...
  if (grepl("^(sp|tr)\\|", y)) {
    parts <- strsplit(y, "\\|")[[1]]
    token <- if (length(parts) >= 3) parts[3] else y
  } else {
    token <- strsplit(y, "\\s+")[[1]][1]
  }

  token <- sub("\\s.*$", "", token)
  token <- sub(";.*$", "", token)
  token <- sub(",.*$", "", token)
  token <- toupper(token)

  m <- stringr::str_match(token, "_([A-Z0-9]{3,10})$")
  code <- m[, 2]
  ifelse(is.na(code) | code == "", NA_character_, code)
}

extract_species_from_label <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(NA_character_)
  y <- trimws(as.character(x)[1])

  # 1) Standard UniProt text, e.g. OS=Homo sapiens OX=9606
  m <- stringr::str_match(y, "OS=([^=]+?)(?:\\sOX=|\\sGN=|\\sPE=|\\sSV=|$)")
  out <- trimws(m[, 2])
  if (!is.na(out) && nzchar(out)) return(out)

  # 2) Alternative organism/taxon formats sometimes found in FASTA headers
  m <- stringr::str_match(y, "(?:Organism|organism|Taxon|taxon)=([^;|]+)")
  out <- trimws(m[, 2])
  if (!is.na(out) && nzchar(out)) return(out)

  # 3) UniProt mnemonic suffix, e.g. CO1A1_HUMAN or CO4A1_MOUSE
  code <- extract_uniprot_species_code(y)
  if (!is.na(code) && nzchar(code)) {
    mapped <- uniprot_species_map[code]
    if (length(mapped) > 0 && !is.na(mapped[1]) && nzchar(mapped[1])) return(unname(mapped[1]))
    return(code)  # fallback: still rank by the species code instead of returning NA
  }

  NA_character_
}

extract_protein_from_label <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(NA_character_)
  y <- trimws(as.character(x)[1])
  if (grepl("^(sp|tr)\\|", y)) {
    parts <- strsplit(y, "\\|")[[1]]
    y <- if (length(parts) >= 3) parts[3] else y
  }
  y <- stringr::str_replace(y, "\\sOS=.*$", "")
  y <- stringr::str_replace(y, "\\sOX=.*$", "")
  y <- trimws(y)
  ifelse(is.na(y) | y == "", NA_character_, y)
}

standardize_spectra_columns <- function(df) {
  if (is.null(df)) return(NULL)
  df <- clean_names_local(df)
  
  pick_first <- function(candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit) == 0) return(NULL)
    hit[1]
  }
  
  protein_col <- pick_first(c("protein", "protein_id", "protein_accession", "accession", "hit_protein", "label", "description"))
  peptide_col <- pick_first(c("peptide", "peptide_sequence", "sequence", "peptide_seq", "seq"))
  species_col <- pick_first(c("species", "organism", "species_name", "taxon_name"))
  evalue_col <- pick_first(c("evalue", "e_value", "expect", "eval", "expectation_value"))
  
  out <- tibble(
    protein = if (!is.null(protein_col)) as.character(df[[protein_col]]) else NA_character_,
    peptide = if (!is.null(peptide_col)) as.character(df[[peptide_col]]) else NA_character_,
    species = if (!is.null(species_col)) as.character(df[[species_col]]) else NA_character_,
    evalue = if (!is.null(evalue_col)) parse_evalue(df[[evalue_col]]) else NA_real_
  ) %>%
    mutate(
      species = ifelse(is.na(species) & !is.na(protein), vapply(protein, extract_species_from_label, character(1)), species),
      protein = ifelse(!is.na(protein), vapply(protein, extract_protein_from_label, character(1)), protein)
    ) %>%
    filter(!(is.na(protein) & is.na(peptide) & is.na(species) & is.na(evalue))) %>%
    distinct()
  
  if (nrow(out) == 0) return(NULL)
  out
}

safe_read_spectra_table <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  
  tryCatch(
    {
      if (ext == "csv") {
        suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
      } else if (ext %in% c("tsv", "txt")) {
        suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )
}

parse_xtandem_xml <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  
  doc <- tryCatch(xml2::read_xml(path), error = function(e) NULL)
  if (is.null(doc)) return(NULL)
  
  protein_nodes <- xml2::xml_find_all(doc, ".//protein")
  domain_nodes <- xml2::xml_find_all(doc, ".//domain")
  
  if (length(protein_nodes) == 0 || length(domain_nodes) == 0) return(NULL)
  
  protein_tbl <- tibble(
    protein_id = xml2::xml_attr(protein_nodes, "id"),
    protein_label = xml2::xml_attr(protein_nodes, "label")
  ) %>%
    mutate(
      protein_label = as.character(protein_label),
      protein = vapply(protein_label, extract_protein_from_label, character(1)),
      species = vapply(protein_label, extract_species_from_label, character(1))
    ) %>%
    distinct(protein_id, .keep_all = TRUE) %>%
    select(protein_id, protein, species)
  
  domain_tbl <- tibble(
    domain_id = xml2::xml_attr(domain_nodes, "id"),
    peptide = xml2::xml_attr(domain_nodes, "seq"),
    evalue = parse_evalue(xml2::xml_attr(domain_nodes, "expect"))
  ) %>%
    mutate(
      protein_id = stringr::str_replace(domain_id, "\\.[0-9]+$", "")
    ) %>%
    select(protein_id, peptide, evalue) %>%
    distinct()
  
  out <- domain_tbl %>%
    left_join(protein_tbl, by = "protein_id") %>%
    transmute(
      protein = protein,
      peptide = peptide,
      species = species,
      evalue = evalue
    ) %>%
    filter(!(is.na(protein) & is.na(peptide) & is.na(species) & is.na(evalue))) %>%
    distinct()
  
  if (nrow(out) == 0) return(NULL)
  out
}

build_species_rank_table <- function(df, selected_peptide = NULL) {

  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  df <- df %>%
    mutate(
      protein = trimws(as.character(protein)),
      peptide = trimws(as.character(peptide)),
      species = trimws(as.character(species)),
      evalue = parse_evalue(evalue)
    )

  if (!is.null(selected_peptide) &&
      nzchar(selected_peptide) &&
      "peptide" %in% names(df)) {

    df <- df %>%
      filter(peptide == selected_peptide)
  }

  df <- df %>%
    filter(
      !is.na(species),
      species != "",
      !is.na(protein),
      protein != ""
    )

  if (nrow(df) == 0) {
    return(NULL)
  }

  species_rank <- df %>%
    group_by(species) %>%
    summarise(

      proteins_count =
        dplyr::n_distinct(protein),

      peptides_count =
        dplyr::n_distinct(peptide),

      best_evalue =
        if (all(is.na(evalue))) {
          NA_real_
        } else {
          min(evalue, na.rm = TRUE)
        },

      .groups = "drop"
    ) %>%
    arrange(
      desc(proteins_count),
      desc(peptides_count),
      is.na(best_evalue),
      best_evalue
    ) %>%
    mutate(
      species_rank = row_number()
    ) %>%
    select(
      species_rank,
      species,
      proteins_count,
      peptides_count,
      best_evalue
    ) %>%
    slice_head(n = 3)

  return(species_rank)
}


build_overall_species_rank_table <- function(df) {

  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  df <- df %>%
    mutate(
      protein = trimws(as.character(protein)),
      peptide = trimws(as.character(peptide)),
      species = trimws(as.character(species)),
      evalue = parse_evalue(evalue)
    )

  df <- df %>%
    filter(
      !is.na(species),
      species != "",
      !is.na(protein),
      protein != ""
    )

  if (nrow(df) == 0) {
    return(NULL)
  }

  overall_rank <- df %>%
    group_by(species) %>%
    summarise(

      proteins_count =
        dplyr::n_distinct(protein),

      peptides_count =
        dplyr::n_distinct(peptide),

      best_evalue =
        if (all(is.na(evalue))) {
          NA_real_
        } else {
          min(evalue, na.rm = TRUE)
        },

      .groups = "drop"
    ) %>%
    arrange(
      desc(proteins_count),
      desc(peptides_count),
      is.na(best_evalue),
      best_evalue
    ) %>%
    mutate(
      species_rank = row_number()
    ) %>%
    select(
      species_rank,
      species,
      proteins_count,
      peptides_count,
      best_evalue
    ) %>%
    slice_head(n = 3)

  return(overall_rank)
}


load_builtin_spectra_results <- function(search_dirs) {
  # FIX 2026-05-09:
  # The old app stopped at the first valid spectra file.
  # Your folder has a tiny spectra_identification_results.csv (3 rows) AND a real output.xml (~24,695 rows).
  # This function now scans all candidates and strongly prefers output.xml when it is valid.

  exact_candidates <- c(
    "output.xml",
    "spectra_identification_results.csv",
    "xtandem_results.csv",
    "tandem_results.csv",
    "spectra_results.csv"
  )

  exact_paths <- unique(unlist(lapply(exact_candidates, function(f) {
    path <- find_existing_file(f, search_dirs)
    if (is.null(path)) character(0) else path
  }), use.names = FALSE))

  pattern_hits <- find_existing_files_by_patterns(
    patterns = c("output\\.xml$", "spectra", "xtandem", "tandem"),
    search_dirs = search_dirs
  )

  candidate_paths <- unique(c(exact_paths, pattern_hits))
  candidate_paths <- candidate_paths[file.exists(candidate_paths)]
  if (length(candidate_paths) == 0) return(NULL)

  parse_candidate <- function(path) {
    ext <- tolower(tools::file_ext(path))
    df <- if (ext == "xml") {
      parse_xtandem_xml(path)
    } else {
      standardize_spectra_columns(safe_read_spectra_table(path))
    }
    if (is.null(df) || nrow(df) == 0) return(NULL)
    list(
      path = normalizePath(path, winslash = "/", mustWork = FALSE),
      basename = basename(path),
      ext = ext,
      rows = nrow(df),
      data = df
    )
  }

  parsed_candidates <- purrr::compact(purrr::map(candidate_paths, parse_candidate))
  if (length(parsed_candidates) == 0) return(NULL)

  # Strongly prefer valid output.xml, because this is the real X! Tandem result file in your Diagnostics tab.
  output_xml_candidates <- parsed_candidates[vapply(parsed_candidates, function(x) identical(tolower(x$basename), "output.xml"), logical(1))]
  if (length(output_xml_candidates) > 0) {
    best_idx <- which.max(vapply(output_xml_candidates, function(x) x$rows, numeric(1)))
    best <- output_xml_candidates[[best_idx]]
  } else {
    # If no output.xml is found, fall back to the valid spectra file with the largest row count.
    best_idx <- which.max(vapply(parsed_candidates, function(x) x$rows, numeric(1)))
    best <- parsed_candidates[[best_idx]]
  }

  message("Using built-in spectra file: ", best$path, " (", best$rows, " rows)")
  attr(best$data, "source_path") <- best$path
  attr(best$data, "source_rows") <- best$rows
  attr(best$data, "source_file") <- best$basename
  best$data
}


parse_mgf_file <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)

  lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character(0))
  if (length(lines) == 0) return(NULL)

  rows <- list()
  in_block <- FALSE
  current <- list(title = NA_character_, pepmass = NA_real_, charge = NA_character_, n_peaks = 0L)
  spectrum_index <- 0L

  for (line in lines) {
    x <- trimws(line)
    if (!nzchar(x)) next

    if (toupper(x) == "BEGIN IONS") {
      in_block <- TRUE
      current <- list(title = NA_character_, pepmass = NA_real_, charge = NA_character_, n_peaks = 0L)
      next
    }

    if (toupper(x) == "END IONS" && in_block) {
      spectrum_index <- spectrum_index + 1L
      rows[[length(rows) + 1L]] <- tibble(
        spectrum_id = spectrum_index,
        spectrum_title = current$title,
        precursor_mz = current$pepmass,
        charge = current$charge,
        n_peaks = current$n_peaks
      )
      in_block <- FALSE
      next
    }

    if (!in_block) next

    if (startsWith(toupper(x), "TITLE=")) {
      current$title <- sub("^TITLE=", "", x, ignore.case = TRUE)
    } else if (startsWith(toupper(x), "PEPMASS=")) {
      val <- sub("^PEPMASS=", "", x, ignore.case = TRUE)
      current$pepmass <- suppressWarnings(as.numeric(strsplit(val, "\\s+")[[1]][1]))
    } else if (startsWith(toupper(x), "CHARGE=")) {
      current$charge <- sub("^CHARGE=", "", x, ignore.case = TRUE)
    } else if (grepl("^[0-9]+(\\.[0-9]+)?\\s+[0-9]+(\\.[0-9]+)?", x)) {
      current$n_peaks <- current$n_peaks + 1L
    }
  }

  if (length(rows) == 0) return(NULL)
  out <- bind_rows(rows)
  attr(out, "upload_kind") <- "mgf"
  attr(out, "source_file") <- basename(path)
  out
}

looks_like_mgf_file <- function(path) {
  if (is.null(path) || !file.exists(path)) return(FALSE)
  lines <- tryCatch(readLines(path, n = 80, warn = FALSE), error = function(e) character(0))
  any(trimws(toupper(lines)) == "BEGIN IONS")
}

read_uploaded_spectra_file <- function(fileinfo) {
  if (is.null(fileinfo) || nrow(fileinfo) == 0) return(NULL)
  
  ext <- tolower(tools::file_ext(fileinfo$name[1]))
  path <- fileinfo$datapath[1]
  if (!file.exists(path)) return(NULL)

  # Accept MGF even when the browser/RStudio truncates or hides the extension.
  # MGF is raw spectra, so this previews spectra metadata instead of protein IDs.
  if (ext == "mgf" || looks_like_mgf_file(path)) {
    return(parse_mgf_file(path))
  }
  
  if (ext == "xml") {
    out <- parse_xtandem_xml(path)
    if (!is.null(out)) attr(out, "upload_kind") <- "identification"
    return(out)
  }
  
  out <- standardize_spectra_columns(safe_read_spectra_table(path))
  if (!is.null(out)) attr(out, "upload_kind") <- "identification"
  out
}


xml_escape_local <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x <- gsub("'", "&apos;", x, fixed = TRUE)
  x
}

find_xtandem_executable <- function(search_dirs = character(0)) {
  env_bin <- Sys.getenv("XTANDEM_BIN", unset = "")
  if (nzchar(env_bin) && file.exists(env_bin)) {
    return(normalizePath(env_bin, winslash = "/", mustWork = FALSE))
  }

  if (.Platform$OS.type == "windows") {
    command_hits <- c(Sys.which("xtandem.exe"), Sys.which("tandem.exe"), Sys.which("xtandem"), Sys.which("tandem"))
    extra_names <- c("xtandem.exe", "tandem.exe", "xtandem", "tandem")
  } else {
    command_hits <- c(Sys.which("xtandem"), Sys.which("tandem"), Sys.which("xtandem.exe"), Sys.which("tandem.exe"))
    extra_names <- c("xtandem", "tandem", "xtandem.exe", "tandem.exe")
  }

  command_hits <- command_hits[nzchar(command_hits) & file.exists(command_hits)]
  if (length(command_hits) > 0) {
    return(normalizePath(command_hits[1], winslash = "/", mustWork = FALSE))
  }

  extra_hits <- unlist(lapply(search_dirs[dir.exists(search_dirs)], function(d) {
    unlist(lapply(extra_names, function(nm) {
      p <- file.path(d, nm)
      if (file.exists(p)) p else character(0)
    }), use.names = FALSE)
  }), use.names = FALSE)

  extra_hits <- unique(extra_hits[file.exists(extra_hits)])
  if (length(extra_hits) == 0) return(NULL)
  normalizePath(extra_hits[1], winslash = "/", mustWork = FALSE)
}

get_xtandem_runtime <- function(search_dirs) {
  taxonomy_xml_path <- find_existing_file(
    Sys.getenv("XTANDEM_TAXONOMY_XML", unset = "taxonomy.xml"),
    search_dirs
  )
  default_input_xml_path <- find_existing_file(
    Sys.getenv("XTANDEM_DEFAULT_INPUT", unset = "default_input.xml"),
    search_dirs
  )
  taxonomy_label_value <- first_non_empty(
    Sys.getenv("XTANDEM_TAXON", unset = ""),
    read_taxonomy_label_from_xml(taxonomy_xml_path)
  )

  list(
    executable = find_xtandem_executable(search_dirs),
    taxonomy_xml = taxonomy_xml_path,
    default_input_xml = default_input_xml_path,
    taxonomy_label = taxonomy_label_value,
    threads = suppressWarnings(as.integer(Sys.getenv("XTANDEM_THREADS", unset = "1"))),
    spectrum_path_type = Sys.getenv("XTANDEM_SPECTRUM_PATH_TYPE", unset = "mgf"),
    parent_error_minus = suppressWarnings(as.numeric(Sys.getenv("XTANDEM_PARENT_ERROR_MINUS", unset = "10"))),
    parent_error_plus = suppressWarnings(as.numeric(Sys.getenv("XTANDEM_PARENT_ERROR_PLUS", unset = "10"))),
    parent_error_units = Sys.getenv("XTANDEM_PARENT_ERROR_UNITS", unset = "ppm"),
    fragment_error = suppressWarnings(as.numeric(Sys.getenv("XTANDEM_FRAGMENT_ERROR", unset = "0.4"))),
    fragment_error_units = Sys.getenv("XTANDEM_FRAGMENT_ERROR_UNITS", unset = "Daltons"),
    cleavage_site = Sys.getenv("XTANDEM_CLEAVAGE_SITE", unset = "[RK]|{P}"),
    maximum_missed_cleavages = suppressWarnings(as.integer(Sys.getenv("XTANDEM_MAX_MISSED_CLEAVAGES", unset = "1"))),
    xsl_path = find_existing_file(
      Sys.getenv("XTANDEM_XSL_PATH", unset = "tandem-input-style.xsl"),
      search_dirs
    )
  )
}

validate_xtandem_runtime <- function(runtime, taxon_override = NULL) {
  problems <- character(0)

  if (is.null(runtime$executable) || !file.exists(runtime$executable)) {
    problems <- c(problems, "XTANDEM_BIN is missing or the X! Tandem executable was not found.")
  }
  if (is.null(runtime$taxonomy_xml) || !file.exists(runtime$taxonomy_xml)) {
    problems <- c(problems, "XTANDEM_TAXONOMY_XML is missing or taxonomy.xml was not found.")
  }
  if (is.null(runtime$default_input_xml) || !file.exists(runtime$default_input_xml)) {
    problems <- c(problems, "XTANDEM_DEFAULT_INPUT is missing or default_input.xml was not found.")
  }

  taxon_value <- first_non_empty(taxon_override, runtime$taxonomy_label, read_taxonomy_label_from_xml(runtime$taxonomy_xml))
  if (!nzchar(taxon_value)) {
    problems <- c(problems, "No X! Tandem taxon/database label was provided.")
  }

  problems
}


write_xtandem_default_xml <- function(
  default_input_xml,
  patched_default_input_xml,
  taxonomy_xml,
  output_xml_path,
  taxon_label
) {
  lines <- tryCatch(
    readLines(default_input_xml, warn = FALSE),
    error = function(e) character(0)
  )

  if (length(lines) == 0) {
    stop("Failed to read default_input.xml.", call. = FALSE)
  }

  replace_or_append_note <- function(lines, label, value) {
    note_prefix <- sprintf('<note type="input" label="%s">', label)
    replacement <- sprintf('<note type="input" label="%s">%s</note>', label, xml_escape_local(value))
    hit <- which(grepl(note_prefix, lines, fixed = TRUE))
    if (length(hit) > 0) {
      lines[hit[1]] <- replacement
    } else {
      bioml_end <- tail(which(grepl("</bioml>", lines, fixed = TRUE)), 1)
      if (length(bioml_end) == 0 || is.na(bioml_end)) {
        lines <- c(lines, replacement)
      } else {
        lines <- append(lines, replacement, after = bioml_end - 1)
      }
    }
    lines
  }

  lines <- replace_or_append_note(
    lines,
    "list path, taxonomy information",
    normalizePath(taxonomy_xml, winslash = "/", mustWork = FALSE)
  )
  lines <- replace_or_append_note(lines, "protein, taxon", taxon_label)
  lines <- replace_or_append_note(
    lines,
    "output, path",
    normalizePath(output_xml_path, winslash = "/", mustWork = FALSE)
  )

  writeLines(lines, con = patched_default_input_xml, useBytes = TRUE)
  normalizePath(patched_default_input_xml, winslash = "/", mustWork = FALSE)
}

write_xtandem_input_xml <- function(
  input_xml_path,
  default_input_xml,
  taxonomy_xml,
  spectrum_path,
  output_xml_path,
  taxon_label,
  xsl_path = NULL,
  threads = 1L,
  spectrum_path_type = "mgf",
  parent_error_minus = 10,
  parent_error_plus = 10,
  parent_error_units = "ppm",
  fragment_error = 0.4,
  fragment_error_units = "Daltons",
  cleavage_site = "[RK]|{P}",
  maximum_missed_cleavages = 1L
) {
  notes <- c(
    sprintf('<note type="input" label="list path, default parameters">%s</note>', xml_escape_local(normalizePath(default_input_xml, winslash = "/", mustWork = FALSE))),
    sprintf('<note type="input" label="list path, taxonomy information">%s</note>', xml_escape_local(normalizePath(taxonomy_xml, winslash = "/", mustWork = FALSE))),
    sprintf('<note type="input" label="protein, taxon">%s</note>', xml_escape_local(taxon_label)),
    sprintf('<note type="input" label="spectrum, path">%s</note>', xml_escape_local(normalizePath(spectrum_path, winslash = "/", mustWork = FALSE))),
    sprintf('<note type="input" label="spectrum, path type">%s</note>', xml_escape_local(spectrum_path_type)),
    sprintf('<note type="input" label="output, path">%s</note>', xml_escape_local(normalizePath(output_xml_path, winslash = "/", mustWork = FALSE))),
    sprintf('<note type="input" label="spectrum, threads">%s</note>', as.integer(threads %||% 1L)),
    sprintf('<note type="input" label="protein, cleavage site">%s</note>', xml_escape_local(cleavage_site)),
    sprintf('<note type="input" label="scoring, maximum missed cleavage sites">%s</note>', as.integer(maximum_missed_cleavages %||% 1L)),
    sprintf('<note type="input" label="spectrum, parent monoisotopic mass error minus">%s</note>', as.character(parent_error_minus %||% 10)),
    sprintf('<note type="input" label="spectrum, parent monoisotopic mass error plus">%s</note>', as.character(parent_error_plus %||% 10)),
    sprintf('<note type="input" label="spectrum, parent monoisotopic mass error units">%s</note>', xml_escape_local(parent_error_units %||% "ppm")),
    sprintf('<note type="input" label="spectrum, fragment monoisotopic mass error">%s</note>', as.character(fragment_error %||% 0.4)),
    sprintf('<note type="input" label="spectrum, fragment monoisotopic mass error units">%s</note>', xml_escape_local(fragment_error_units %||% "Daltons"))
  )

  if (!is.null(xsl_path) && file.exists(xsl_path)) {
    notes <- c(
      notes,
      sprintf('<note type="input" label="output, xsl path">%s</note>', xml_escape_local(normalizePath(xsl_path, winslash = "/", mustWork = FALSE)))
    )
  }

  xml_lines <- c(
    '<?xml version="1.0"?>',
    '<bioml>',
    notes,
    '</bioml>'
  )

  writeLines(xml_lines, con = input_xml_path, useBytes = TRUE)
  normalizePath(input_xml_path, winslash = "/", mustWork = FALSE)
}

run_xtandem_search <- function(fileinfo, runtime, taxon_override = NULL) {
  if (is.null(fileinfo) || nrow(fileinfo) == 0) {
    stop("Upload an MGF file first.", call. = FALSE)
  }

  original_name <- as.character(fileinfo$name[1] %||% "")
  ext <- tolower(tools::file_ext(original_name))
  source_path <- fileinfo$datapath[1]

  if (!file.exists(source_path)) {
    stop("Uploaded file is no longer available on disk.", call. = FALSE)
  }
  if (!(ext == "mgf" || looks_like_mgf_file(source_path))) {
    stop("Phase 1 supports only MGF upload.", call. = FALSE)
  }

  issues <- validate_xtandem_runtime(runtime, taxon_override = taxon_override)
  if (length(issues) > 0) {
    stop(paste(issues, collapse = "\n"), call. = FALSE)
  }

  taxon_label <- first_non_empty(taxon_override, runtime$taxonomy_label, read_taxonomy_label_from_xml(runtime$taxonomy_xml))
  job_dir <- file.path(tempdir(), paste0("xtandem_job_", as.integer(as.numeric(Sys.time())), "_", sample(100000:999999, 1)))
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

  mgf_name <- basename(original_name)
  if (!nzchar(mgf_name) || tolower(tools::file_ext(mgf_name)) != "mgf") {
    mgf_name <- paste0(tools::file_path_sans_ext(ifelse(nzchar(mgf_name), mgf_name, "upload")), ".mgf")
  }
  job_mgf <- file.path(job_dir, mgf_name)
  ok_copy <- file.copy(source_path, job_mgf, overwrite = TRUE)
  if (!isTRUE(ok_copy)) {
    stop("Failed to stage the uploaded MGF file for X! Tandem.", call. = FALSE)
  }

  output_xml <- file.path(job_dir, "output.xml")
  input_xml <- file.path(job_dir, "input.xml")
  patched_default_input_xml <- file.path(job_dir, "default_input.xml")
  log_txt <- file.path(job_dir, "xtandem_stdout.log")

  write_xtandem_default_xml(
    default_input_xml = runtime$default_input_xml,
    patched_default_input_xml = patched_default_input_xml,
    taxonomy_xml = runtime$taxonomy_xml,
    output_xml_path = output_xml,
    taxon_label = taxon_label
  )

  write_xtandem_input_xml(
    input_xml_path = input_xml,
    default_input_xml = patched_default_input_xml,
    taxonomy_xml = runtime$taxonomy_xml,
    spectrum_path = job_mgf,
    output_xml_path = output_xml,
    taxon_label = taxon_label,
    xsl_path = runtime$xsl_path,
    threads = runtime$threads,
    spectrum_path_type = runtime$spectrum_path_type,
    parent_error_minus = runtime$parent_error_minus,
    parent_error_plus = runtime$parent_error_plus,
    parent_error_units = runtime$parent_error_units,
    fragment_error = runtime$fragment_error,
    fragment_error_units = runtime$fragment_error_units,
    cleavage_site = runtime$cleavage_site,
    maximum_missed_cleavages = runtime$maximum_missed_cleavages
  )

  exit_status <- tryCatch(
    suppressWarnings(
      system2(
        command = runtime$executable,
        args = c(input_xml),
        stdout = log_txt,
        stderr = log_txt,
        wait = TRUE
      )
    ),
    error = function(e) {
      stop(paste("Failed to start X! Tandem:", e$message), call. = FALSE)
    }
  )

  if (!file.exists(log_txt)) {
    writeLines("X! Tandem finished with no console output.", con = log_txt, useBytes = TRUE)
  }

  if (!file.exists(output_xml)) {
    log_excerpt <- tryCatch(
      paste(utils::tail(readLines(log_txt, warn = FALSE), 20), collapse = "\n"),
      error = function(e) ""
    )
    stop(
      paste(
        "X! Tandem finished but output.xml was not created.",
        paste0("Exit status: ", exit_status %||% NA_integer_),
        if (nzchar(log_excerpt)) paste("Last log lines:", log_excerpt, sep = "\n") else "No X! Tandem log output was captured.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  parsed <- parse_xtandem_xml(output_xml)
  if (is.null(parsed) || nrow(parsed) == 0) {
    stop(
      paste(
        "X! Tandem created output.xml, but no protein/peptide/species rows were parsed.",
        "Check whether the database taxon matches your taxonomy.xml entry and whether the FASTA database contains the target species.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  attr(parsed, "upload_kind") <- "xtandem_identification"
  attr(parsed, "source_file") <- basename(output_xml)
  attr(parsed, "job_dir") <- normalizePath(job_dir, winslash = "/", mustWork = FALSE)
  attr(parsed, "input_xml_path") <- normalizePath(input_xml, winslash = "/", mustWork = FALSE)
  attr(parsed, "patched_default_input_xml") <- normalizePath(patched_default_input_xml, winslash = "/", mustWork = FALSE)
  attr(parsed, "output_xml_path") <- normalizePath(output_xml, winslash = "/", mustWork = FALSE)
  attr(parsed, "stdout_log_path") <- normalizePath(log_txt, winslash = "/", mustWork = FALSE)
  attr(parsed, "mgf_path") <- normalizePath(job_mgf, winslash = "/", mustWork = FALSE)
  attr(parsed, "taxon_label") <- taxon_label
  parsed
}

build_xtandem_runtime_table <- function(runtime) {
  tibble(
    setting = c(
      "X! Tandem executable",
      "taxonomy.xml",
      "default_input.xml",
      "taxonomy label",
      "threads",
      "parent mass error minus",
      "parent mass error plus",
      "parent mass error units",
      "fragment mass error",
      "fragment mass error units",
      "cleavage site",
      "max missed cleavages"
    ),
    value = c(
      runtime$executable %||% NA_character_,
      runtime$taxonomy_xml %||% NA_character_,
      runtime$default_input_xml %||% NA_character_,
      runtime$taxonomy_label %||% NA_character_,
      runtime$threads %||% NA_integer_,
      runtime$parent_error_minus %||% NA_real_,
      runtime$parent_error_plus %||% NA_real_,
      runtime$parent_error_units %||% NA_character_,
      runtime$fragment_error %||% NA_real_,
      runtime$fragment_error_units %||% NA_character_,
      runtime$cleavage_site %||% NA_character_,
      runtime$maximum_missed_cleavages %||% NA_integer_
    )
  )
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
    "step3_bone_unique_peptides_scored.csv",
    "spectra_identification_results.csv",
    "xtandem_results.csv",
    "tandem_results.csv",
    "spectra_results.csv",
    "output.xml",
    "taxonomy.xml",
    "default_input.xml",
    "tandem-input-style.xsl"
  )
  
  purrr::map_dfr(files, function(f) {
    path <- find_existing_file(f, search_dirs)
    ext <- if (!is.null(path)) tolower(tools::file_ext(path)) else ""
    df <- if (is.null(path)) {
      NULL
    } else if (ext == "xml") {
      parse_xtandem_xml(path)
    } else {
      safe_read_spectra_table(path) %||% safe_read_csv(path)
    }
    
    tibble(
      filename = f,
      found = !is.null(path),
      path = path %||% NA_character_,
      rows = if (!is.null(df)) nrow(df) else 0,
      cols = if (!is.null(df)) ncol(df) else 0
    )
  })
}

# Initial data load
DEFAULT_DATA_DIR <- Sys.getenv("DATA_DIR", unset = getwd())
SEARCH_DIRS <- build_search_dirs(DEFAULT_DATA_DIR)

step2_data_raw <- coalesce_step2_data(SEARCH_DIRS)
step3_species_raw <- coalesce_step3_species_summary(SEARCH_DIRS)
step3_unique_raw <- coalesce_step3_unique(SEARCH_DIRS)
diagnostics_raw <- build_diagnostics_table(SEARCH_DIRS)

# Images
HOME_IMG_1 <- "https://raw.githubusercontent.com/hagarelazab/ANUBUIS/main/anubis-egyptian-god-Google-Search-03-08-2026_04_02_PM.png"
HOME_IMG_2 <- "https://raw.githubusercontent.com/hagarelazab/ANUBUIS/main/Ancient-Civilizations-Ancient-Egypt-03-08-2026_03_36_PM%20(1).png"

# Team images
# IMPORTANT:
# Shiny cannot reliably display absolute local image paths directly in the browser.
# This exposes /Users/hagarhaitham/Downloads as a Shiny static resource called "hagar_downloads".
HAGAR_DOWNLOADS_DIR <- "/Users/hagarhaitham/Downloads"

if (dir.exists(HAGAR_DOWNLOADS_DIR)) {
  try(addResourcePath("hagar_downloads", HAGAR_DOWNLOADS_DIR), silent = TRUE)
}

find_hagar_image <- function() {
  candidates <- c(
    "hagar_team.png",
    "hagar_team.jpg",
    "hagar_team.jpeg",
    "hagar picture.png",
    "hagar picture.jpg",
    "hagar picture.jpeg",
    "Hagar picture.png",
    "Hagar picture.jpg",
    "Hagar picture.jpeg",
    "Hagar.png",
    "Hagar.jpg",
    "Hagar.jpeg",
    "hagar.png",
    "hagar.jpg",
    "hagar.jpeg"
  )

  if (dir.exists(HAGAR_DOWNLOADS_DIR)) {
    hit <- candidates[file.exists(file.path(HAGAR_DOWNLOADS_DIR, candidates))]
    if (length(hit) > 0) {
      return(paste0("hagar_downloads/", utils::URLencode(hit[1], reserved = TRUE)))
    }
  }

  # fallback if no matching file is found
  "hagar_team.png"
}

TEAM_IMG_1 <- find_hagar_image()
TEAM_IMG_2 <- "https://raw.githubusercontent.com/hagarelazab/ANUBUIS/main/WhatsApp%20Image%202026-03-08%20at%2021.59.47.jpeg"

# Theme
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

app_footer <- tags$footer(
  class = "app-footer",
  div(
    class = "footer-inner",
    div(class = "footer-accent"),
    div(class = "footer-text", HTML("&copy; 2026 Hagar Elazab. All Rights Reserved."))
  )
)

# UI
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
        background: linear-gradient(180deg, rgba(255,255,255,.98) 0%, rgba(251,252,254,.96) 100%);
        border: 1px solid rgba(228,231,236,.92);
        border-radius: 24px;
        padding: 22px 22px 18px 22px;
        box-shadow: 0 18px 45px rgba(16,24,40,.08);
        margin-bottom: 18px;
        position: relative;
        overflow: hidden;
      }
      .panel-card::before {
        content: '';
        position: absolute;
        inset: 0 0 auto 0;
        height: 4px;
        background: linear-gradient(90deg, rgba(15,23,42,0), rgba(212,175,55,.75), rgba(37,99,235,.55), rgba(15,23,42,0));
        opacity: .75;
      }
      .panel-card .section-head { position: relative; z-index: 1; }
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
        width:125px;
        height:125px; border-radius:999px; overflow:hidden;
        background:#fff;
      }

      .avatar.avatar-hagar{
        border: 3px solid rgba(200,161,77,.82);
        box-shadow:
          0 14px 28px rgba(2,6,23,.12),
          0 0 0 1px rgba(200,161,77,.18),
          0 0 24px rgba(200,161,77,.12);
      }

      .avatar.avatar-nawal{
        border: 3px solid #ffffff;
        box-shadow:
          0 14px 28px rgba(2,6,23,.12),
          0 0 0 1px rgba(228,231,236,1);
      }

      .avatar img{ width:100%; height:100%; display:block; }
      .avatar img.avatar-hagar-img{
        object-fit: cover !important;
        object-position: center center !important;
        transform: none !important;
        background: #ffffff;
      }
      .avatar img.avatar-nawal-img{
        object-fit: cover !important;
        object-position: center center !important;
        transform: none !important;
      }

      .member-name{ margin-top:12px; font-weight:900; color:#101828; letter-spacing:.2px; font-size:1.05rem; }
      .member-affil{ margin-top: 8px; color:#475467; font-size: .94rem; line-height: 1.6; max-width: 520px; }


      .snapshot-note {
        display: flex; align-items: flex-start; gap: 12px;
        background: linear-gradient(180deg, #fbfcfe 0%, #ffffff 100%);
        border: 1px solid #e6e8ef; border-radius: 18px;
        padding: 14px 16px; color: #475467; line-height: 1.7;
        margin: 4px 0 18px 0;
      }
      .snapshot-note-icon {
        width: 28px; height: 28px; border-radius: 999px;
        display: inline-flex; align-items: center; justify-content: center;
        background: #f8f3e6; color: #8a6116; font-weight: 900;
        flex: 0 0 auto;
      }
      .snapshot-table-wrap .dataTables_wrapper { margin-top: 6px; }
      .snapshot-table-wrap table.dataTable tbody td { padding-top: 12px !important; padding-bottom: 12px !important; }



      .table-metric-grid {
        display: grid;
        grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: 14px;
        margin: 18px 0 18px 0;
      }
      .table-mini-card {
        background: linear-gradient(180deg, #ffffff 0%, #fbfcfe 100%);
        border: 1px solid #e6e8ef;
        border-radius: 18px;
        padding: 14px 15px;
        box-shadow: 0 10px 24px rgba(16,24,40,.055);
        position: relative;
        overflow: hidden;
      }
      .table-mini-card::before {
        content: '';
        position: absolute;
        left: 0; top: 0; bottom: 0;
        width: 4px;
        background: linear-gradient(180deg, #d4af37, #0f172a);
        opacity: .85;
      }
      .table-mini-label {
        color: #667085;
        font-size: .78rem;
        font-weight: 800;
        text-transform: uppercase;
        letter-spacing: .04em;
        margin-bottom: 8px;
      }
      .table-mini-value {
        color: #101828;
        font-size: 1.45rem;
        line-height: 1;
        font-weight: 900;
        letter-spacing: -.02em;
      }
      .table-mini-sub {
        margin-top: 7px;
        color: #667085;
        font-size: .82rem;
        font-weight: 600;
      }
      .professional-table-wrap {
        border: 1px solid #e6e8ef;
        border-radius: 20px;
        overflow: hidden;
        box-shadow: 0 12px 28px rgba(16,24,40,.055);
        background: #ffffff;
      }
      .professional-table-wrap .dataTables_wrapper {
        padding: 14px 14px 10px 14px;
      }
      .professional-table-wrap table.dataTable tbody tr:hover {
        background: #fff8e6 !important;
      }
      .professional-table-wrap table.dataTable tbody td {
        vertical-align: middle;
        border-top: 1px solid #f1f3f6 !important;
      }
      .professional-table-wrap table.dataTable thead th {
        text-transform: uppercase;
        letter-spacing: .03em;
        font-size: .78rem;
      }

      .dataTables_wrapper .dataTables_filter input { border-radius: 12px !important; border: 1px solid #d0d5dd !important; padding: 6px 10px !important; }
      table.dataTable { border-collapse: separate !important; border-spacing: 0 !important; }
      table.dataTable thead th { background: #f8fafc !important; color: #344054 !important; border-bottom: 1px solid #eaecf0 !important; font-weight: 700 !important; }

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
      @media (max-width: 700px) { .metric-grid { grid-template-columns: 1fr; } .hero-title { font-size: 2.3rem; } .anubis-grid { grid-template-columns: 1fr; } .collage-row { grid-template-columns: 1fr; } .collage-small { min-height: 200px; } :root { --footer-h: 72px; } .avatar{ width:118px; height:118px; } }
    "))
  ),
  
  footer = app_footer,
  
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
      div(
        class = "panel-card snapshot-table-wrap",
        section_title("Combined snapshot", "A clean preview of the Step 2 peptide table with quick dataset indicators."),
        uiOutput("snapshot_metrics"),
        div(
          class = "snapshot-note",
          div(class = "snapshot-note-icon", "i"),
          div(HTML("This table is the main combined peptide reference. Use the filters under each column to quickly isolate a protein, peptide, species, taxid, or tissue subset."))
        ),
        DTOutput("home_snapshot")
      )
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
        uiOutput("species_summary_metrics"),
        div(class = "professional-table-wrap", DTOutput("step3_species_summary_table"))
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
        uiOutput("unique_peptides_metrics"),
        div(class = "professional-table-wrap", DTOutput("step3_unique_peptides_table"))
      )
    )
  ),
  
  nav_panel(
    title = "Spectra Identification",
    value = "Spectra Identification",
    div(
      class = "app-shell",
      div(
        class = "panel-card",
        section_title(
          "Spectra identification",
          "Table 1 shows protein, peptide, and E-value. Click any row to see Table 2 with species ranking for that peptide."
        ),
        radioButtons(
          "spectra_mode",
          "Choose source",
          choices = c(
            "Explore our data" = "builtin",
            "Upload your own" = "upload"
          ),
          selected = "builtin",
          inline = TRUE
        ),
        conditionalPanel(
          condition = "input.spectra_mode === 'builtin'",
          p(
            class = "text-muted",
            "Built-in mode reads local spectra results from the app search folders, including XTandem output.xml. The app prefers valid output.xml; if no output.xml is found, it selects the largest valid spectra file."
          ),
          uiOutput("spectra_builtin_source_info"),
          h4("Top 3 overall ranked species"),
          DTOutput("spectra_builtin_overall_species_table"),
          br(),
          h4("Table 1: Protein / Peptide / E-value"),
          DTOutput("spectra_builtin_table"),
          br(),
          h4("Table 2: Species ranking for selected peptide"),
          uiOutput("spectra_builtin_selected_peptide"),
          DTOutput("spectra_builtin_species_table")
        ),
        conditionalPanel(
          condition = "input.spectra_mode === 'upload'",
          fileInput(
            "spectra_upload",
            "Upload MGF for X! Tandem search",
            accept = c(".mgf")
          ),
          fluidRow(
            column(
              4,
              textInput(
                "spectra_upload_taxon",
                "X! Tandem taxon / database label",
                value = first_non_empty(Sys.getenv("XTANDEM_TAXON", unset = ""), read_taxonomy_label_from_xml(find_existing_file(Sys.getenv("XTANDEM_TAXONOMY_XML", unset = "taxonomy.xml"), SEARCH_DIRS)))
              )
            ),
            column(
              3,
              numericInput(
                "spectra_upload_top_n",
                "Top N peptides",
                min = 5,
                max = 100,
                value = 25,
                step = 1
              )
            ),
            column(
              5,
              br(),
              actionButton("run_xtandem_btn", "Run X! Tandem on uploaded MGF", class = "btn btn-primary")
            )
          ),
          p(
            class = "text-muted",
            "Upload an MGF file, run X! Tandem, then review identified peptides and ranked species from the generated output."
          ),
          uiOutput("spectra_upload_status"),
          br(),
          h4("Top 3 overall ranked species"),
          DTOutput("spectra_upload_overall_species_table"),
          br(),
          h4("Table 1: Protein / Peptide / E-value"),
          DTOutput("spectra_upload_table"),
          br(),
          h4("Table 2: Species ranking for selected peptide"),
          uiOutput("spectra_upload_selected_peptide"),
          DTOutput("spectra_upload_species_table")
        )
      )
    )
  ),
  
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

# Server
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
      metric_box("Peptide occurrences", format_big(peptides_total), "All rows in Step 2", "#B68D40"),
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
    y_max <- make_plot_ymax(plot_df$n_species, 1.26)
    g <- ggplot(plot_df, aes(x = tissue, y = n_species, fill = tissue, text = tooltip)) +
      geom_col(width = 0.55, alpha = 0.96) +
      geom_text(
        aes(y = n_species + y_max * 0.035, label = scales::comma(n_species)),
        color = "#0f172a",
        fontface = "bold",
        size = 4.2
      ) +
      scale_fill_manual(values = plot_palette_tissue, drop = FALSE) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0))) +
      coord_cartesian(ylim = c(0, y_max), clip = "off") +
      labs(x = NULL, y = "Distinct species") +
      base_gg_theme()
    ggplotly_clean(g)
  })
  
  output$home_unique_tissue_bar <- renderPlotly({

    df <- step2_tbl()
    validate(need(!is.null(df), "No Step 2 data found."))
    validate(need(all(c("tissue", "peptide") %in% names(df)),
                  "Step 2 data must contain tissue and peptide."))
    plot_df <- df %>%
      filter(!is.na(peptide), !is.na(tissue)) %>%
      group_by(tissue) %>%
      summarise(unique_peptides = n_distinct(peptide), .groups = "drop") %>%
      mutate(tissue = factor(tissue, levels = c("Skin", "Bone", "Combined")),
             tooltip = paste0("Tissue: ", tissue, "<br>Unique peptides: ", scales::comma(unique_peptides)))
    y_max <- make_plot_ymax(plot_df$unique_peptides, 1.26)
    g <- ggplot(plot_df, aes(x = tissue, y = unique_peptides, fill = tissue, text = tooltip)) +
      geom_col(width = 0.55, alpha = 0.96) +
      geom_text(
        aes(y = unique_peptides + y_max * 0.035, label = scales::comma(unique_peptides)),
        color = "#0f172a",
        fontface = "bold",
        size = 4.2
      ) +
      scale_fill_manual(values = plot_palette_tissue, drop = FALSE) +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0))) +
      coord_cartesian(ylim = c(0, y_max), clip = "off") +
      labs(x = NULL, y = "Unique peptides") +
      base_gg_theme()
    ggplotly_clean(g)
  })
  
  output$snapshot_metrics <- renderUI({
    df <- step2_tbl()
    if (is.null(df)) return(NULL)
    protein_total <- safe_n_distinct(df, "protein_id")
    peptide_total <- safe_n_distinct(df, "peptide")
    species_total <- safe_n_distinct(df, "species")
    tissue_total <- safe_n_distinct(df, "tissue")
    div(
      class = "metric-grid",
      metric_box("Rows", format_big(nrow(df)), "Total table records", "#0f172a"),
      metric_box("Proteins", format_big(protein_total), "Distinct protein IDs", "#175cd3"),
      metric_box("Peptides", format_big(peptide_total), "Distinct peptide sequences", "#0f766e"),
      metric_box("Species", format_big(species_total), paste0("Tissue groups: ", tissue_total), "#c0841a")
    )
  })
  
  output$home_snapshot <- renderDT({
    df <- step2_tbl()
    if (is.null(df)) return(empty_dt("No Step 2 data available."))
    keep <- intersect(c("protein_id", "peptide", "pep_len", "species", "taxid", "tissue"), names(df))
    datatable(
      df[, keep, drop = FALSE],
      rownames = FALSE,
      filter = "top",
      class = "stripe hover compact",
      options = list(
        pageLength = 10,
        lengthMenu = c(10, 15, 25, 50),
        scrollX = TRUE,
        autoWidth = TRUE
      )
    )
  })
  
  proteins_by_species <- reactive({
    df <- step2_tbl()
    validate(need(!is.null(df), "No Step 2 data."))
    validate(need(all(c("species", "protein_id") %in% names(df)), "Need species and protein_id in Step 2 data."))
    df <- filter_by_tissue(df, input$pps_tissue)
    grp <- if ("tissue" %in% names(df)) c("tissue", "species") else c("species")
    df %>%
      group_by(across(all_of(grp))) %>%
      summarise(n_proteins = n_distinct(protein_id, na.rm = TRUE), .groups = "drop") %>%
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

    # Professional horizontal chart
    # Fixes overlapping species labels and removes the legend text that was
    # sitting on top of the x-axis labels in the previous vertical version.
    plot_df <- df %>%
      arrange(desc(n_proteins), species) %>%
      mutate(
        bar_group = ifelse(row_number() <= 3, "Top 3", "Other"),
        fill_color = ifelse(bar_group == "Top 3", "#c0841a", "#162033"),
        species_label = stringr::str_wrap(species, width = 24),
        species_label = factor(species_label, levels = rev(species_label)),
        tooltip = paste0(
          "Species: ", species,
          "<br>Proteins: ", scales::comma(n_proteins),
          "<br>Tissue: ", tissue,
          "<br>Group: ", bar_group
        )
      )

    x_max <- make_plot_ymax(plot_df$n_proteins, 1.18)

    g <- ggplot(plot_df, aes(x = species_label, y = n_proteins, text = tooltip)) +
      geom_col(aes(fill = fill_color), width = 0.64, alpha = 0.98, show.legend = FALSE) +
      geom_text(
        aes(y = n_proteins + x_max * 0.025, label = scales::comma(n_proteins)),
        color = "#0f172a",
        fontface = "bold",
        size = 3.8,
        hjust = 0
      ) +
      scale_fill_identity() +
      scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0))) +
      coord_flip(ylim = c(0, x_max), clip = "off") +
      labs(x = NULL, y = "Distinct proteins") +
      base_gg_theme() +
      theme(
        axis.text.y = element_text(size = 10.5, color = "#344054", face = "bold"),
        axis.text.x = element_text(size = 10, color = "#475467"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "#edf1f5"),
        plot.margin = margin(t = 18, r = 42, b = 18, l = 18)
      )

    ggplotly_clean(g) %>%
      layout(
        showlegend = FALSE,
        margin = list(l = 145, r = 55, t = 22, b = 55)
      )
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
      summarise(total_peptide_occurrences = n(),
                unique_peptides = n_distinct(peptide, na.rm = TRUE),
                .groups = "drop") %>%
      arrange(desc(unique_peptides), species)
  })
  
  output$peptides_species_plot <- renderPlotly({
    df <- peptides_by_species()
    validate(need(nrow(df) > 0, "No rows to plot."))

    metric_col <- if (identical(input$peps_metric, "Total peptide occurrences")) {
      "total_peptide_occurrences"
    } else {
      "unique_peptides"
    }

    metric_label <- if (identical(metric_col, "total_peptide_occurrences")) {
      "Total peptide occurrences"
    } else {
      "Unique peptides"
    }

    plot_df <- df %>%
      arrange(desc(.data[[metric_col]])) %>%
      slice_head(n = input$peps_top_n) %>%
      mutate(
        value = as.numeric(.data[[metric_col]])
      ) %>%
      arrange(value) %>%
      mutate(
        rank_desc = rank(-value, ties.method = "first"),
        bar_color = ifelse(rank_desc <= 3, "#c0841a", "#0f172a"),
        species_label = as.character(species),
        label_text = scales::comma(value),
        tooltip = paste0(
          "Species: ", species_label,
          "<br>", metric_label, ": ", scales::comma(value),
          "<br>Tissue: ", tissue
        )
      )

    y_max <- make_plot_ymax(plot_df$value, 1.22)

    plot_ly(
      data = plot_df,
      x = ~species_label,
      y = ~value,
      type = "bar",
      text = ~label_text,
      textposition = "outside",
      hoverinfo = "text",
      hovertext = ~tooltip,
      marker = list(
        color = plot_df$bar_color,
        line = list(color = "rgba(15,23,42,0.12)", width = 1)
      )
    ) %>%
      layout(
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor = "rgba(0,0,0,0)",
        showlegend = FALSE,
        margin = list(l = 70, r = 30, t = 26, b = 115),
        bargap = 0.32,
        hoverlabel = list(
          bgcolor = "#0f172a",
          bordercolor = "#0f172a",
          font = list(color = "#ffffff", size = 12)
        ),
        xaxis = list(
          title = "",
          tickangle = -45,
          tickfont = list(color = "#344054", size = 11),
          showgrid = FALSE,
          zeroline = FALSE
        ),
        yaxis = list(
          title = list(text = metric_label, font = list(color = "#344054", size = 13)),
          range = c(0, y_max),
          tickformat = ",",
          gridcolor = "#eef2f7",
          zeroline = FALSE,
          tickfont = list(color = "#475467", size = 11)
        )
      ) %>%
      config(displayModeBar = FALSE, responsive = TRUE)
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
      summarise(n_unique_peptides = n_distinct(peptide, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(n_unique_peptides)) %>%
      slice_head(n = input$ups_top_n) %>%
      mutate(value = as.numeric(n_unique_peptides)) %>%
      arrange(value) %>%
      mutate(
        rank_desc = rank(-value, ties.method = "first"),
        bar_color = ifelse(rank_desc <= 3, "#c0841a", "#0f172a"),
        species_label = as.character(species),
        label_text = scales::comma(value),
        tooltip = paste0(
          "Species: ", species_label,
          "<br>Species-specific peptides: ", scales::comma(value),
          "<br>Tissue: ", tissue
        )
      )

    y_max <- make_plot_ymax(plot_df$value, 1.22)

    plot_ly(
      data = plot_df,
      x = ~species_label,
      y = ~value,
      type = "bar",
      text = ~label_text,
      textposition = "outside",
      hoverinfo = "text",
      hovertext = ~tooltip,
      marker = list(
        color = plot_df$bar_color,
        line = list(color = "rgba(15,23,42,0.12)", width = 1)
      )
    ) %>%
      layout(
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor = "rgba(0,0,0,0)",
        showlegend = FALSE,
        margin = list(l = 70, r = 30, t = 26, b = 115),
        bargap = 0.32,
        hoverlabel = list(
          bgcolor = "#0f172a",
          bordercolor = "#0f172a",
          font = list(color = "#ffffff", size = 12)
        ),
        xaxis = list(
          title = "",
          tickangle = -45,
          tickfont = list(color = "#344054", size = 11),
          showgrid = FALSE,
          zeroline = FALSE
        ),
        yaxis = list(
          title = list(text = "Species-specific peptides", font = list(color = "#344054", size = 13)),
          range = c(0, y_max),
          tickformat = ",",
          gridcolor = "#eef2f7",
          zeroline = FALSE,
          tickfont = list(color = "#475467", size = 11)
        )
      ) %>%
      config(displayModeBar = FALSE, responsive = TRUE)
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
    if ("peptide_score" %in% names(df)) df$peptide_score <- suppressWarnings(as.numeric(df$peptide_score))
    df
  })
  
  output$peptide_existence_plot <- renderPlotly({
    df <- peptide_existence_filtered()
    validate(need(nrow(df) > 0, "No rows to plot."))
    if (!"tissue" %in% names(df)) df$tissue <- "Combined"

    df <- df %>%
      mutate(
        peptide_score = suppressWarnings(as.numeric(peptide_score)),
        tissue = as.character(tissue)
      ) %>%
      filter(!is.na(peptide_score))

    validate(need(nrow(df) > 0, "No numeric peptide scores to plot."))

    if (identical(input$pex_kind, "Histogram")) {
      plot_df <- df %>%
        mutate(
          tooltip = paste0(
            "Peptide score: ", peptide_score,
            "<br>Tissue: ", tissue
          )
        )

      plot_ly(
        plot_df,
        x = ~peptide_score,
        type = "histogram",
        nbinsx = 20,
        marker = list(
          color = "#0f172a",
          line = list(color = "#ffffff", width = 1)
        ),
        opacity = 0.92,
        hovertemplate = paste(
          "Peptide score bin: %{x}<br>",
          "Count: %{y}<extra></extra>"
        )
      ) %>%
        layout(
          paper_bgcolor = "rgba(0,0,0,0)",
          plot_bgcolor = "rgba(0,0,0,0)",
          showlegend = FALSE,
          bargap = 0.18,
          margin = list(l = 70, r = 35, t = 20, b = 80),
          xaxis = list(
            title = list(text = "Peptide score", font = list(color = "#344054", size = 13)),
            gridcolor = "rgba(0,0,0,0)",
            zeroline = FALSE,
            tickfont = list(color = "#475467", size = 11)
          ),
          yaxis = list(
            title = list(text = "Count", font = list(color = "#344054", size = 13)),
            gridcolor = "#eef2f7",
            zeroline = FALSE,
            tickfont = list(color = "#475467", size = 11)
          ),
          hoverlabel = list(bgcolor = "#0f172a", bordercolor = "#0f172a", font = list(color = "#ffffff"))
        ) %>%
        config(displayModeBar = FALSE, responsive = TRUE)
    } else {
      plot_df <- df %>%
        filter(!is.na(species), species != "") %>%
        group_by(species) %>%
        mutate(median_score = median(peptide_score, na.rm = TRUE)) %>%
        ungroup() %>%
        arrange(desc(median_score)) %>%
        mutate(
          species = factor(species, levels = unique(species)),
          tooltip = paste0(
            "Species: ", species,
            "<br>Peptide score: ", peptide_score,
            "<br>Tissue: ", tissue
          )
        )

      validate(need(nrow(plot_df) > 0, "No species rows to plot."))

      plot_ly(
        plot_df,
        x = ~peptide_score,
        y = ~species,
        type = "box",
        orientation = "h",
        boxpoints = "outliers",
        marker = list(color = "#0f172a", opacity = 0.55),
        line = list(color = "#c0841a", width = 2),
        fillcolor = "rgba(15,23,42,0.12)",
        hoverinfo = "text",
        text = ~tooltip
      ) %>%
        layout(
          paper_bgcolor = "rgba(0,0,0,0)",
          plot_bgcolor = "rgba(0,0,0,0)",
          showlegend = FALSE,
          margin = list(l = 170, r = 35, t = 20, b = 70),
          xaxis = list(
            title = list(text = "Peptide score", font = list(color = "#344054", size = 13)),
            gridcolor = "#eef2f7",
            zeroline = FALSE,
            tickfont = list(color = "#475467", size = 11)
          ),
          yaxis = list(
            title = "",
            autorange = "reversed",
            gridcolor = "rgba(0,0,0,0)",
            tickfont = list(color = "#344054", size = 11)
          ),
          hoverlabel = list(bgcolor = "#0f172a", bordercolor = "#0f172a", font = list(color = "#ffffff"))
        ) %>%
        config(displayModeBar = FALSE, responsive = TRUE)
    }
  })
  
  output$peptide_existence_table <- renderDT({
    df <- peptide_existence_filtered()
    keep <- intersect(c("peptide", "species", "tissue", "peptide_score", "n_proteins"), names(df))
    datatable(df[, keep, drop = FALSE], rownames = FALSE, filter = "top",
              options = list(pageLength = 15, scrollX = TRUE))
  })
  
  existence_species_filtered <- reactive({
    df <- s3_species_tbl()
    validate(need(!is.null(df), "No species summary data."))
    validate(need(all(c("species", "unique_peptide_count", "sum_peptide_score", "mean_peptide_score") %in% names(df)),
                  "Need species, unique_peptide_count, sum_peptide_score, mean_peptide_score."))
    df <- filter_by_tissue(df, input$eps_tissue)
    if (!"tissue" %in% names(df)) df$tissue <- "Combined"
    
    df$mean_peptide_score <- round(as.numeric(df$mean_peptide_score), 2)
    
    df
  })
  
  output$existence_species_plot <- renderPlotly({
    df <- existence_species_filtered()
    validate(need(nrow(df) > 0, "No rows to plot."))
    metric_col <- input$eps_metric

    plot_df <- df %>%
      mutate(value = suppressWarnings(as.numeric(.data[[metric_col]]))) %>%
      filter(!is.na(value), !is.na(species), species != "") %>%
      arrange(desc(value)) %>%
      slice_head(n = input$eps_top_n) %>%
      mutate(
        rank = row_number(),
        group = ifelse(rank <= 3, "Top 3", "Other"),
        color = ifelse(rank <= 3, "#c0841a", "#0f172a"),
        label = ifelse(metric_col == "mean_peptide_score", sprintf("%.2f", value), scales::comma(value)),
        species = factor(species, levels = rev(species)),
        tooltip = paste0(
          "Species: ", species,
          "<br>", metric_col, ": ", label,
          "<br>Tissue: ", tissue
        )
      )

    validate(need(nrow(plot_df) > 0, "No numeric species values to plot."))

    x_max <- max(plot_df$value, na.rm = TRUE) * 1.18

    plot_ly(
      plot_df,
      x = ~value,
      y = ~species,
      type = "bar",
      orientation = "h",
      marker = list(color = ~color, line = list(color = "rgba(255,255,255,0.9)", width = 1)),
      text = ~label,
      textposition = "outside",
      cliponaxis = FALSE,
      hoverinfo = "text",
      hovertext = ~tooltip
    ) %>%
      layout(
        paper_bgcolor = "rgba(0,0,0,0)",
        plot_bgcolor = "rgba(0,0,0,0)",
        showlegend = FALSE,
        margin = list(l = 175, r = 70, t = 20, b = 70),
        xaxis = list(
          title = list(text = metric_col, font = list(color = "#344054", size = 13)),
          range = c(0, x_max),
          tickformat = ",",
          gridcolor = "#eef2f7",
          zeroline = FALSE,
          tickfont = list(color = "#475467", size = 11)
        ),
        yaxis = list(
          title = "",
          gridcolor = "rgba(0,0,0,0)",
          tickfont = list(color = "#344054", size = 11)
        ),
        hoverlabel = list(bgcolor = "#0f172a", bordercolor = "#0f172a", font = list(color = "#ffffff"))
      ) %>%
      config(displayModeBar = FALSE, responsive = TRUE)
  })
  
  output$existence_species_table <- renderDT({
    df <- existence_species_filtered()

    dt <- datatable(df, rownames = FALSE, filter = "top",
                    options = list(pageLength = 12, scrollX = TRUE))
    if ("mean_peptide_score" %in% names(df)) dt <- formatRound(dt, "mean_peptide_score", 2)
    dt
  })
  

  output$species_summary_metrics <- renderUI({
    df <- s3_species_tbl()
    if (is.null(df)) return(NULL)
    df <- filter_by_tissue(df, input$s3sum_tissue)

    species_n <- if ("species" %in% names(df)) dplyr::n_distinct(df$species, na.rm = TRUE) else nrow(df)
    unique_total <- if ("unique_peptide_count" %in% names(df)) sum(suppressWarnings(as.numeric(df$unique_peptide_count)), na.rm = TRUE) else NA_real_
    score_total <- if ("sum_peptide_score" %in% names(df)) sum(suppressWarnings(as.numeric(df$sum_peptide_score)), na.rm = TRUE) else NA_real_
    best_species <- if ("species" %in% names(df) && "unique_peptide_count" %in% names(df) && nrow(df) > 0) {
      df %>% mutate(unique_peptide_count = suppressWarnings(as.numeric(unique_peptide_count))) %>% arrange(desc(unique_peptide_count)) %>% slice(1) %>% pull(species)
    } else "—"

    div(
      class = "table-metric-grid",
      table_mini_card("Species", format_big(species_n), "after filters"),
      table_mini_card("Unique peptides", ifelse(is.na(unique_total), "—", format_big(unique_total)), "total signal"),
      table_mini_card("Peptide score", ifelse(is.na(score_total), "—", format_big(score_total)), "summed evidence"),
      table_mini_card("Top species", best_species, "by unique peptides")
    )
  })

  output$unique_peptides_metrics <- renderUI({
    df <- s3_unique_tbl()
    if (is.null(df)) return(NULL)
    df <- filter_by_tissue(df, input$s3up_tissue)
    if (!is.null(input$s3up_species) && input$s3up_species != "All" && "species" %in% names(df)) df <- df %>% filter(species == input$s3up_species)
    if ("peptide_score" %in% names(df)) {
      df$peptide_score <- suppressWarnings(as.numeric(df$peptide_score))
      df <- df %>% filter(is.na(peptide_score) | peptide_score >= input$s3up_min_score)
    }

    peptide_n <- if ("peptide" %in% names(df)) dplyr::n_distinct(df$peptide, na.rm = TRUE) else nrow(df)
    species_n <- if ("species" %in% names(df)) dplyr::n_distinct(df$species, na.rm = TRUE) else NA_integer_
    protein_total <- if ("n_proteins" %in% names(df)) sum(suppressWarnings(as.numeric(df$n_proteins)), na.rm = TRUE) else NA_real_
    avg_score <- if ("peptide_score" %in% names(df)) mean(suppressWarnings(as.numeric(df$peptide_score)), na.rm = TRUE) else NA_real_

    div(
      class = "table-metric-grid",
      table_mini_card("Peptides", format_big(peptide_n), "visible rows"),
      table_mini_card("Species", ifelse(is.na(species_n), "—", format_big(species_n)), "represented"),
      table_mini_card("Protein links", ifelse(is.na(protein_total), "—", format_big(protein_total)), "total matches"),
      table_mini_card("Avg score", ifelse(is.nan(avg_score) || is.na(avg_score), "—", round(avg_score, 2)), "filtered peptides")
    )
  })

  output$step3_species_summary_table <- renderDT({
    df <- s3_species_tbl()
    if (is.null(df)) return(empty_dt("No species summary available."))
    df <- filter_by_tissue(df, input$s3sum_tissue)

    numeric_cols <- intersect(c("unique_peptide_count", "sum_peptide_score", "mean_peptide_score"), names(df))
    for (cc in numeric_cols) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))

    dir_desc <- identical(input$s3sum_dir, "Descending")
    sort_col <- input$s3sum_sort
    if (sort_col %in% names(df)) df <- df %>% arrange(if (dir_desc) desc(.data[[sort_col]]) else .data[[sort_col]])

    dt <- datatable(
      df,
      rownames = FALSE,
      filter = "top",
      class = "display compact nowrap premium-table",
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50),
        scrollX = TRUE,
        dom = "<'row align-items-center mb-2'<'col-sm-6'l><'col-sm-6'f>>rt<'row align-items-center mt-2'<'col-sm-5'i><'col-sm-7'p>>",
        autoWidth = TRUE
      )
    )
    if ("mean_peptide_score" %in% names(df)) dt <- formatRound(dt, "mean_peptide_score", 2)
    if ("species" %in% names(df)) dt <- formatStyle(dt, "species", fontWeight = "700", color = "#101828")
    if ("unique_peptide_count" %in% names(df)) dt <- formatStyle(dt, "unique_peptide_count", fontWeight = "700", color = "#0f172a")
    if ("tissue" %in% names(df)) dt <- formatStyle(dt, "tissue", color = "#475467", fontWeight = "600")
    dt
  })
  
  output$step3_unique_peptides_table <- renderDT({
    df <- s3_unique_tbl()
    if (is.null(df)) return(empty_dt("No unique peptides available."))
    df <- filter_by_tissue(df, input$s3up_tissue)
    if (!is.null(input$s3up_species) && input$s3up_species != "All" && "species" %in% names(df)) df <- df %>% filter(species == input$s3up_species)
    if ("peptide_score" %in% names(df)) {
      df$peptide_score <- suppressWarnings(as.numeric(df$peptide_score))
      df <- df %>% filter(is.na(peptide_score) | peptide_score >= input$s3up_min_score)
    }
    if ("n_proteins" %in% names(df)) df$n_proteins <- suppressWarnings(as.numeric(df$n_proteins))

    dt <- datatable(
      df,
      rownames = FALSE,
      filter = "top",
      class = "display compact nowrap premium-table",
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50),
        scrollX = TRUE,
        dom = "<'row align-items-center mb-2'<'col-sm-6'l><'col-sm-6'f>>rt<'row align-items-center mt-2'<'col-sm-5'i><'col-sm-7'p>>",
        autoWidth = TRUE
      )
    )
    if ("peptide" %in% names(df)) dt <- formatStyle(dt, "peptide", fontWeight = "700", color = "#101828")
    if ("species" %in% names(df)) dt <- formatStyle(dt, "species", color = "#0f172a", fontWeight = "600")
    if ("peptide_score" %in% names(df)) dt <- formatStyle(dt, "peptide_score", fontWeight = "800", color = styleInterval(c(2, 4), c("#475467", "#b7791f", "#0f172a")))
    if ("tissue" %in% names(df)) dt <- formatStyle(dt, "tissue", color = "#475467", fontWeight = "600")
    dt
  })
  
  spectra_builtin_tbl <- reactive({
    load_builtin_spectra_results(SEARCH_DIRS)
  })
  

  output$spectra_builtin_source_info <- renderUI({
    df <- spectra_builtin_tbl()
    if (is.null(df) || nrow(df) == 0) {
      return(div(class = "text-muted", "No built-in spectra source loaded yet."))
    }
    src <- attr(df, "source_path") %||% "Unknown source"
    src_rows <- attr(df, "source_rows") %||% nrow(df)
    div(
      class = "subtle-badge-row",
      div(class = "subtle-badge", HTML(paste0("<strong>Using spectra file:</strong> ", htmltools::htmlEscape(basename(src))))),
      div(class = "subtle-badge", HTML(paste0("<strong>Rows loaded:</strong> ", scales::comma(src_rows)))),
      div(class = "subtle-badge", HTML(paste0("<strong>Full path:</strong> ", htmltools::htmlEscape(src))))
    )
  })

  xtandem_runtime <- reactive({
    get_xtandem_runtime(SEARCH_DIRS)
  })

  xtandem_upload_result <- reactiveVal(NULL)
  xtandem_upload_error <- reactiveVal(NULL)
  xtandem_upload_last_run <- reactiveVal(NULL)
  observe({
    runtime <- xtandem_runtime()
    current_taxon <- trimws(as.character(input$spectra_upload_taxon %||% ""))
    fallback_taxon <- first_non_empty(runtime$taxonomy_label, read_taxonomy_label_from_xml(runtime$taxonomy_xml))
    if (!nzchar(current_taxon) && nzchar(fallback_taxon)) {
      updateTextInput(session, "spectra_upload_taxon", value = fallback_taxon)
    }
  })


  observeEvent(input$spectra_upload, {
    xtandem_upload_result(NULL)
    xtandem_upload_error(NULL)
    xtandem_upload_last_run(NULL)
  })

  observeEvent(input$run_xtandem_btn, {
    xtandem_upload_result(NULL)
    xtandem_upload_error(NULL)

    if (is.null(input$spectra_upload) || nrow(input$spectra_upload) == 0) {
      xtandem_upload_error("Upload an MGF file first.")
      return()
    }

    run_name <- tolower(tools::file_ext(input$spectra_upload$name[1] %||% ""))
    if (!(run_name == "mgf" || looks_like_mgf_file(input$spectra_upload$datapath[1]))) {
      xtandem_upload_error("Phase 1 supports only MGF upload.")
      return()
    }

    result <- tryCatch(
      withProgress(message = "Running X! Tandem", value = 0, {
        incProgress(0.15, detail = "Checking server configuration")
        runtime <- xtandem_runtime()

        incProgress(0.35, detail = "Preparing X! Tandem input XML")
        out <- run_xtandem_search(
          fileinfo = input$spectra_upload,
          runtime = runtime,
          taxon_override = input$spectra_upload_taxon
        )

        incProgress(0.85, detail = "Parsing output.xml")
        out
      }),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      xtandem_upload_error(conditionMessage(result))
      xtandem_upload_last_run(Sys.time())
      return()
    }

    xtandem_upload_result(result)
    xtandem_upload_last_run(Sys.time())
  })

  spectra_uploaded_tbl <- reactive({
    xtandem_upload_result()
  })

  spectra_builtin_table1 <- reactive({
    df <- spectra_builtin_tbl()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df %>%
      transmute(
        protein = protein,
        peptide = peptide,
        e_value = evalue
      ) %>%
      arrange(peptide, is.na(e_value), e_value, protein) %>%
      distinct()
  })

  spectra_upload_table1 <- reactive({
    df <- spectra_uploaded_tbl()
    if (is.null(df) || nrow(df) == 0) return(NULL)

    df %>%
      transmute(
        protein = protein,
        peptide = peptide,
        e_value = evalue
      ) %>%
      arrange(peptide, is.na(e_value), e_value, protein) %>%
      distinct() %>%
      slice_head(n = input$spectra_upload_top_n %||% 25)
  })

  selected_builtin_peptide <- reactive({
    df <- spectra_builtin_table1()
    idx <- input$spectra_builtin_table_rows_selected
    if (is.null(df) || nrow(df) == 0 || length(idx) == 0) return(NULL)
    if (is.na(df$peptide[idx[1]]) || !nzchar(as.character(df$peptide[idx[1]]))) NULL else as.character(df$peptide[idx[1]])
  })

  selected_upload_peptide <- reactive({
    df <- spectra_upload_table1()
    idx <- input$spectra_upload_table_rows_selected
    if (is.null(df) || nrow(df) == 0 || length(idx) == 0) return(NULL)
    if (!"peptide" %in% names(df)) return(NULL)
    if (is.na(df$peptide[idx[1]]) || !nzchar(as.character(df$peptide[idx[1]]))) NULL else as.character(df$peptide[idx[1]])
  })


  output$spectra_builtin_overall_species_table <- renderDT({
    df <- build_overall_species_rank_table(spectra_builtin_tbl())

    if (is.null(df) || nrow(df) == 0) {
      return(empty_dt("No built-in spectra ranking available."))
    }

    dt <- datatable(
      df,
      rownames = FALSE,
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        order = list(list(0, "asc")),
        columnDefs = evalue_column_defs(4, digits = 2)
      )
    )
    dt
  })

  output$spectra_upload_overall_species_table <- renderDT({
    uploaded <- spectra_uploaded_tbl()

    if (is.null(uploaded)) {
      err <- xtandem_upload_error()
      if (!is.null(err)) {
        return(empty_dt(err))
      }
      return(empty_dt("Run X! Tandem on an uploaded MGF first."))
    }

    df <- build_overall_species_rank_table(uploaded)

    if (is.null(df) || nrow(df) == 0) {
      return(empty_dt("No ranked species found in the uploaded X! Tandem results."))
    }

    dt <- datatable(
      df,
      rownames = FALSE,
      options = list(
        dom = "t",
        paging = FALSE,
        searching = FALSE,
        info = FALSE,
        scrollX = TRUE,
        order = list(list(0, "asc")),
        columnDefs = evalue_column_defs(4, digits = 2)
      )
    )
    dt
  })

  output$spectra_upload_status <- renderUI({
    err <- xtandem_upload_error()
    res <- xtandem_upload_result()
    last_run <- xtandem_upload_last_run()

    if (!is.null(err)) {
      return(
        div(
          class = "subtle-badge-row",
          div(class = "subtle-badge", HTML(paste0("<strong>Status:</strong> Failed"))),
          div(class = "subtle-badge", style = "border-color:#fecdca;color:#b42318;", htmltools::htmlEscape(err))
        )
      )
    }

    if (is.null(res)) {
      return(
        div(
          class = "subtle-badge-row",
          div(class = "subtle-badge", HTML("<strong>Status:</strong> Waiting for uploaded MGF")),
          if (!is.null(last_run)) div(class = "subtle-badge", HTML(paste0("<strong>Last attempt:</strong> ", format(last_run, "%Y-%m-%d %H:%M:%S"))))
        )
      )
    }

    div(
      class = "subtle-badge-row",
      div(class = "subtle-badge", HTML("<strong>Status:</strong> Search complete")),
      div(class = "subtle-badge", HTML(paste0("<strong>Taxon:</strong> ", htmltools::htmlEscape(attr(res, "taxon_label") %||% "")))),
      div(class = "subtle-badge", HTML(paste0("<strong>Rows parsed:</strong> ", scales::comma(nrow(res))))),
      div(class = "subtle-badge", HTML(paste0("<strong>Job folder:</strong> ", htmltools::htmlEscape(attr(res, "job_dir") %||% ""))))
    )
  })

  output$spectra_builtin_selected_peptide <- renderUI({
    pep <- selected_builtin_peptide()
    if (is.null(pep)) {
      return(div(class = "text-muted", "Click a row in Table 1 to see species details."))
    }
    div(class = "subtle-badge-row",
        div(class = "subtle-badge", HTML(paste0("<strong>Selected peptide:</strong> ", pep))))
  })

  output$spectra_upload_selected_peptide <- renderUI({
    pep <- selected_upload_peptide()
    res <- spectra_uploaded_tbl()

    if (is.null(res)) {
      err <- xtandem_upload_error()
      if (!is.null(err)) {
        return(div(class = "text-muted", "Fix the upload/server settings above, then run X! Tandem again."))
      }
      return(div(class = "text-muted", "Run X! Tandem on an uploaded MGF, then click a row in Table 1 to see species details."))
    }

    if (is.null(pep)) {
      return(div(class = "text-muted", "Click a row in Table 1 to see species details."))
    }

    div(class = "subtle-badge-row",
        div(class = "subtle-badge", HTML(paste0("<strong>Selected peptide:</strong> ", pep))))
  })

  output$spectra_builtin_table <- renderDT({
    df <- spectra_builtin_table1()
    if (is.null(df) || nrow(df) == 0) {
      return(empty_dt(
        "No built-in spectra identification file found. Add spectra_identification_results.csv, xtandem_results.csv, tandem_results.csv, spectra_results.csv, or output.xml to one of the search folders."
      ))
    }

    dt <- datatable(
      df,
      rownames = FALSE,
      filter = "top",
      selection = "single",
      options = list(
        pageLength = 100,
        scrollX = TRUE,
        order = list(list(1, "asc"), list(2, "asc")),
        columnDefs = evalue_column_defs(2, digits = 2)
      )
    )
    dt
  })

  output$spectra_upload_table <- renderDT({
    if (is.null(input$spectra_upload)) {
      return(empty_dt("Upload an MGF file, then run X! Tandem."))
    }

    df <- spectra_upload_table1()
    if (is.null(df) || nrow(df) == 0) {
      err <- xtandem_upload_error()
      if (!is.null(err)) {
        return(empty_dt(err))
      }
      return(empty_dt("No uploaded X! Tandem results yet. Click 'Run X! Tandem on uploaded MGF'."))
    }

    dt <- datatable(
      df,
      rownames = FALSE,
      filter = "top",
      selection = "single",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        order = list(list(1, "asc"), list(2, "asc")),
        columnDefs = evalue_column_defs(2, digits = 2)
      )
    )
    dt
  })

  output$spectra_builtin_species_table <- renderDT({
    pep <- selected_builtin_peptide()
    df <- build_species_rank_table(spectra_builtin_tbl(), pep)

    if (is.null(pep)) {
      return(empty_dt("Select a row from Table 1 first."))
    }

    if (is.null(df) || nrow(df) == 0) {
      return(empty_dt("No species found for the selected peptide."))
    }

    dt <- datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(0, "asc")),
        columnDefs = evalue_column_defs(4, digits = 2)
      )
    )
    dt
  })

  output$spectra_upload_species_table <- renderDT({
    uploaded <- spectra_uploaded_tbl()
    pep <- selected_upload_peptide()

    if (is.null(uploaded)) {
      err <- xtandem_upload_error()
      if (!is.null(err)) {
        return(empty_dt(err))
      }
      return(empty_dt("Run X! Tandem on an uploaded MGF first."))
    }

    df <- build_species_rank_table(uploaded, pep)

    if (is.null(pep)) {
      return(empty_dt("Select a row from Table 1 first."))
    }

    if (is.null(df) || nrow(df) == 0) {
      return(empty_dt("No species found for the selected peptide."))
    }

    dt <- datatable(
      df,
      rownames = FALSE,
      filter = "top",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(0, "asc")),
        columnDefs = evalue_column_defs(4, digits = 2)
      )
    )
    dt
  })
}

shinyApp(ui, server)
