# Explanatory plot --------------------------------------------------------
rm(list=ls())
library(tibble)
library(ggplot2)
library(dplyr)
theme_set(theme_light())

# Defining 11 unique survival patterns ------------------------------------
# Underlying: 0--3
start <- 0; end <- 3; GP <- 1.5
surv.patterns <- tribble(
  ~id, ~pre.surgery.time, ~had.surgery, ~post.surgery.time, ~died,
  # A: Has surgery in GP, dies after GP.
  "A",                 1,            T,                2.5,     T, 
  # B: Has surgery in GP, censored at end
  "B",              0.74,            T,                end,     F,  
  # C: Has surgery and dies in GP
  "C",              0.25,            T,                1.1,     T,
  # D: Has surgery and is censored in GP
  "D",              0.50,            T,               0.80,     F,
  # E: Has surgery in GP, censored after GP (but not at end)
  "E",              0.33,            T,               2.00,     T,
  # F: Has surgery and dies after GP,
  "F",              1.75,            T,               2.75,     T,
  # G: Has surgery and is censored after GP (but not at end)
  "G",              1.95,            T,               2.25,     F,
  # H: Never has surgery, and censored at end
  "H",               end,            F,               end,      F,
  # I: Never has surgery, dies after GP
  "I",               2.0,            F,               2.0,      T,
  # J: Never has surgery, dies in GP
  "J",               0.4,            F,               0.4,      T,
  # K: Never has surgery, censored in GP
  "K",               1.2,            F,               1.2,      F
)
surv.patterns$id <- forcats::fct_inorder(surv.patterns$id)

ggplot(surv.patterns, aes(x = forcats::fct_rev(id), group = id), lwd = 4) + 
  geom_hline(yintercept = c(start, GP, end), lty = 3, col = 'grey') + 
  geom_pointrange(aes(ymin = 0, ymax = pre.surgery.time, y = -100), 
                  colour = 'red2', lwd = 3) + 
  geom_pointrange(data = surv.patterns[surv.patterns$had.surgery,],
                  aes(y = -100, ymin = pre.surgery.time, ymax = post.surgery.time),
                  colour = 'steelblue', lwd = 3) +
  geom_point(data = surv.patterns[surv.patterns$died,],
             aes(y = post.surgery.time),
             colour = 'black', size = 4, pch = 15) + 
  labs(x = 'Survival pattern\n', y = '\nTime from enrollment', title = 'Observed') + 
  scale_y_continuous(breaks = c(start, GP, end),
                     labels = c('Start', "GP", "End")) + 
  coord_flip(ylim = c(start, end)) + 
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9)
  ) -> P1

# Censoring version -------------------------------------------------------
ggplot(surv.patterns, aes(x = forcats::fct_rev(id), group = id), lwd = 4) + 
  geom_hline(yintercept = c(start, GP, end), lty = 3, col = 'grey') + 
  geom_pointrange(aes(ymin = 0, ymax = pmin(pre.surgery.time, GP), y = -100), 
                  colour = 'red2', lwd = 3) +
  geom_pointrange(
    data = surv.patterns[surv.patterns$id%in%LETTERS[6:9],],
    aes(ymin = pmin(pre.surgery.time, GP), y = -100,
                      ymax = post.surgery.time), lwd = 3.2, colour = 'black') +
  geom_pointrange(
    data = surv.patterns[surv.patterns$id%in%LETTERS[6:9],],
    aes(ymin = pmin(pre.surgery.time, GP), y = -100,
        ymax = post.surgery.time), lwd = 2, colour = 'red2') +  
  geom_point(data = surv.patterns[surv.patterns$id%in%LETTERS[6:9]&surv.patterns$died,],
             aes(y = post.surgery.time),
             colour = 'black', size = 4, pch = 15) + 
  labs(x = 'Artificial survival pattern\n', y = '\nTime from enrollment',
       title = 'Control clone arm') + 
  scale_y_continuous(breaks = c(start, GP, end),
                     labels = c('Start', "GP", "End")) + 
  coord_flip(ylim = c(start, end)) + 
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9)
  ) -> P2

# Surgery version here
legend.set <- tribble(
  ~id, ~label, ~y,
  "A", "Post-surgery", 10,
  "A", "Post-surgery", 20,
  "A", "Pre-surgery", 10,
  "A", "Pre-surgery", 20
)
ggplot(surv.patterns, aes(x = forcats::fct_rev(id), group = id), lwd = 4) + 
  geom_hline(yintercept = c(start, GP, end), lty = 3, col = 'grey') + 
  geom_pointrange(aes(ymin = 0, ymax = pmin(pre.surgery.time, GP), y = -100), 
                  colour = 'red2', lwd = 3) + 
  geom_pointrange(data = surv.patterns[surv.patterns$had.surgery & surv.patterns$pre.surgery.time < GP, ],
                  aes(ymin = pre.surgery.time, ymax = pmin(GP, post.surgery.time), y = -100),
                  colour = 'steelblue', lwd = 3) + 
  geom_pointrange(
    data = surv.patterns[surv.patterns$id%in%c('A', 'B', 'E'),],
    aes(ymin = GP, ymax = post.surgery.time, y = -100), lwd = 3.2, colour = 'black') +
  geom_pointrange(data = surv.patterns[surv.patterns$id%in%c('A', 'B', 'E'), ],
                  aes(ymin = GP, ymax = post.surgery.time, y = -100),
                  colour = 'steelblue', lwd = 2.5) + 
  geom_point(data = surv.patterns[surv.patterns$id%in%c('A', 'B', 'E', 'C')&surv.patterns$died,],
             aes(y = post.surgery.time),
             colour = 'black', size = 4, pch = 15) +
  geom_line(data = legend.set, aes(y = y, colour = label), lwd = 1.15) +
  scale_color_manual(values = c('steelblue', 'red2'))+
  labs(x = 'Artificial survival pattern\n', y = '\nTime from enrollment',
       title = 'Surgery clone arm', colour = NULL) + 
  scale_y_continuous(breaks = c(start, GP, end),
                     labels = c('Start', "GP", "End")) + 
  coord_flip(ylim = c(start, end)) + 
  theme(
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.key = element_blank(),
    legend.background = element_rect(fill = NA),
    legend.position = 'inside',
    legend.position.inside = c(0.75, 0.5)
  ) -> P3


# Arrange -----------------------------------------------------------------
rescaler <- function(x){
  x + theme(
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    title = element_text(size = 16.5),
    legend.text = element_text(size = 13)
  )
}

library(ggpubr)
ggarrange(rescaler(P1), rescaler(P2), rescaler(P3), nrow = 1)
ggsave('./docs-final/Figure1.png')
ggsave('./docs-final/Figure1b.tiff', device = 'tiff', dpi = 300L)
