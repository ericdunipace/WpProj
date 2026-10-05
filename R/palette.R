# JAMA color palette, same values as ggsci::pal_jama("default").
# Kept internal so plotting does not depend on ggsci.
jama_colors <- c("#374E55", "#DF8F44", "#00A1D5", "#B24745",
                 "#79AF97", "#6A6599", "#80796B")

pal_jama <- function(n) {
  if (n > length(jama_colors)) {
    warning("The JAMA palette has ", length(jama_colors), " colors but ", n,
            " are needed; extra groups will be shown as NA.")
  }
  jama_colors[seq_len(n)]
}

jama_scale <- function(aesthetics, ...) {
  # ggplot2 >= 3.5.0 deprecated `scale_name`, older versions require it
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::discrete_scale(aesthetics, palette = pal_jama, ...)
  } else {
    ggplot2::discrete_scale(aesthetics, scale_name = "jama", palette = pal_jama, ...)
  }
}

scale_color_jama <- function(...) jama_scale("colour", ...)

scale_fill_jama <- function(...) jama_scale("fill", ...)
