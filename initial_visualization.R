#just for use in initial data visualization
data_full <- read_csv("COMPLETE_COMBINED_MUSCLE_FIBER_DATA.csv")
ggplot(df_phalloidin, aes(x = X, y = Y)) +
  geom_point(size = 0.1, alpha = 0.3, colour = "magenta") +
  facet_wrap(~Sample, scales = "free") +
  theme_dark() +
  labs(
    title = "Raw fiber positions — all samples (Phalloidin C1)",
    x = "X (pixels)", y = "Y (pixels)"
  ) +
  theme(
    strip.text    = element_text(size = 7, colour = "white"),
    axis.text     = element_text(size = 5),
    panel.spacing = unit(0.3, "lines")
  )

ggsave("raw_samples_overview.png", width = 16, height = 20, dpi = 150)

# Detach sf to avoid the ggsave conflict
detach("package:sf", unload = TRUE)

# Rebuild the plot
p <- ggplot(df_phalloidin, aes(x = X, y = Y)) +
  geom_point(size = 0.1, alpha = 0.3, colour = "magenta") +
  facet_wrap(~Sample, scales = "free") +
  theme_dark() +
  labs(
    title = "Raw fiber positions — all samples (Phalloidin C1)",
    x = "X (pixels)", y = "Y (pixels)"
  ) +
  theme(
    strip.text    = element_text(size = 7, colour = "white"),
    axis.text     = element_text(size = 5),
    panel.spacing = unit(0.3, "lines")
  )

# Save using png() instead of ggsave
png("raw_samples_overview.png", width = 16, height = 20, units = "in", res = 150)
print(p)
dev.off()