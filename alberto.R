
library(tidyverse)

d_anto <- read_csv("data/anthocyanins_longFormat.csv")

d_anto %>%
  ggplot(aes(x = tiempo, y = concentration, colour = processing)) +
  # geom_line() +
  geom_smooth(method = "lm", linetype = 2) +
  geom_point() +
  scale_y_log10() +
  facet_wrap(compound ~ Temp + sweetener,
             scales = "free") +
  theme(legend.position = "top")




anthocyanins$processing %>% unique()

anthocyanins %>%
  # group_by(compound, Temp, sweetener, tiempo, processing) %>%
  # summarize()
  filter(Temp == "4", sweetener == "ST",
         compound == "Delphinidin.3.O.sambubioside.5.O.glucoside",
         tiempo < 0.002)

