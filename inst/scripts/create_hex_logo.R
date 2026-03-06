# Create Hexagonal Logo for EpigenicR
# Install hexSticker if needed: install.packages("hexSticker")

library(hexSticker)
library(ggplot2)

# Option 1: Simple text-based hex sticker
# ========================================
sticker(
  # No subplot - just text
  subplot = ~ plot.new(),
  s_x = 1,
  s_y = 1,
  s_width = 1,
  s_height = 1,
  # Package name
  package = "EpigenicR",
  p_size = 20,
  p_color = "#FFFFFF",
  p_family = "sans",
  p_y = 1,
  # Hex colors
  h_fill = "#2C3E50",      # Dark blue-grey
  h_color = "#3498DB",      # Bright blue
  h_size = 1.5,
  # Output
  filename = "man/figures/logo.png",
  dpi = 300
)

message("✓ Simple logo created: man/figures/logo.png")


# Option 2: Logo with DNA/Epigenome visualization
# ================================================
# Create a simple epigenetic visualization
# Use x positions that stay within hexagon bounds
x_pos <- seq(1.8, 7.2, length.out = 8)

epigenetic_plot <- ggplot() +
  # Simulate histone modifications as colored bars
  geom_segment(aes(x = x_pos, xend = x_pos, 
                   y = c(0.3, rep(0, 6), 0.3), 
                   yend = c(2.5, 1.5, 3, 2, 2.8, 1.8, 2.3, 2.6)),
               color = c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", 
                        "#E74C3C", "#3498DB", "#2ECC71", "#F39C12"),
               size = 3, lineend = "round") +
  # Add methylation marks as points
  geom_point(aes(x = c(x_pos[2], x_pos[4], x_pos[6], x_pos[7]), 
                 y = c(3.5, 3.3, 3.4, 3.6)),
             color = "#9B59B6", size = 4, shape = 19) +
  # DNA backbone
  geom_line(aes(x = x_pos, y = rep(0.2, 8)), 
            color = "#34495E", size = 1.5, alpha = 0.7) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  coord_cartesian(clip = "off")

sticker(
  epigenetic_plot,
  s_x = 1,
  s_y = 0.9,
  s_width = 1.6,
  s_height = 1.2,
  # Package name
  package = "EpigenicR",
  p_size = 18,
  p_color = "#FFFFFF",
  p_family = "sans",
  p_y = 1.5,
  # Hex colors
  h_fill = "#2C3E50",
  h_color = "#3498DB",
  h_size = 1.5,
  # Output
  filename = "man/figures/logo_with_plot.png",
  dpi = 300
)

message("✓ Plot-based logo created: man/figures/logo_with_plot.png")


# Option 3: Logo with custom icon/image (if you have one)
# ========================================================
# If you have a custom icon or image file, you can use it like this:
# 
# sticker(
#   "path/to/your/icon.png",
#   s_x = 1,
#   s_y = 0.9,
#   s_width = 0.6,
#   s_height = 0.6,
#   package = "EpigenicR",
#   p_size = 18,
#   p_color = "#FFFFFF",
#   p_y = 1.45,
#   h_fill = "#2C3E50",
#   h_color = "#3498DB",
#   h_size = 1.5,
#   filename = "man/figures/logo_custom.png",
#   dpi = 300
# )


# Color Palette Suggestions for Epigenomics
# ==========================================
# Dark backgrounds: #2C3E50, #34495E, #1A1A2E
# Accent colors:
#   - Blue (DNA): #3498DB, #2980B9
#   - Red (H3K27ac): #E74C3C, #C0392B
#   - Green (H3K4me3): #2ECC71, #27AE60
#   - Orange (H3K36me3): #F39C12, #E67E22
#   - Purple (Methylation): #9B59B6, #8E44AD

message("\n=== Hex Logo Creation Complete ===")
message("Two versions created:")
message("1. man/figures/logo.png - Simple text version")
message("2. man/figures/logo_with_plot.png - With epigenetic visualization")
message("\nTo use in README, add at the top:")
message('# EpigenicR <img src="man/figures/logo.png" align="right" height="139" />')
