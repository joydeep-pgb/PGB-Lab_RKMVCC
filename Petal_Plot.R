library(dplyr); library(ggplot2)
petals = 5
petal_angle = 360/petals

# Data in R

data <- data.frame(
  
  category = c("Metabolism", "Lipid metabolism", "AMPK signaling pathway", 
               
               "PPAR signaling pathway", "Lipid biosynthesis proteins"),
  
  value = c(37, 17, 8, 6, 6)
  
)
  
  
data %>% mutate(petal = row_number(),
         theta0 = petal * petal_angle) |>
  reframe(theta = theta0 + c(0, -petal_angle/2,  0, 
                             petal_angle/2, 0),
          r = value^0.5 * c(0, 0.6, 1, 0.6, 0), .by = c(category, value, petal, theta0)) |> 
  ggplot(aes(theta, r + 0.1, group = petal, fill = petal |> as.character())) +
  ggforce::stat_bspline(geom = "area", n = 1000) +
  guides(fill = "none") +
  coord_radial() +
  theme_void()
