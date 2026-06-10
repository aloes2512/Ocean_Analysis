# from R4DS page 400 nested data
library(dslabs)
library (tidyverse)
library(modelr)
gapminder %>%  head() %>% names()
nz <- filter(gapminder, country == "New Zealand")
nz %>%
  ggplot(aes(year, life_expectancy)) +
  geom_line() +
  ggtitle("Full data = ")
germ <- filter(gapminder, country == "Germany")
germ %>%
  ggplot(aes(year, life_expectancy)) +
  geom_line() +
  ggtitle("Full data = ")
nz_mod <- lm(life_expectancy ~ year, data = nz)
nz %>%
  add_predictions(nz_mod) %>%
  ggplot(aes(year, pred)) +
  geom_line() +
  ggtitle("Linear trend + ")
germ_mod <- lm(life_expectancy ~ year, data = germ)
germ %>%
  add_predictions(germ_mod) %>%
  ggplot(aes(year, pred)) +
  geom_line() +
  ggtitle("Linear trend + ")

nz %>%
  add_residuals(nz_mod) %>%
  ggplot(aes(year, resid)) +
  geom_hline(yintercept = 0, color = "white", size = 3) +
  geom_line() +
  ggtitle("Remaining pattern")
germ %>%
  add_residuals(germ_mod) %>%
  ggplot(aes(year, resid)) +
  geom_hline(yintercept = 0, color = "white", size = 3) +
  geom_line() +
  ggtitle("Remaining pattern")
# Erzeugung einer liste vo Daten
# jede Zeile entspricht dem tibble eines gruppenelements
by_country <- gapminder %>%
  group_by(country, continent) %>%
  nest()
head(by_country,2)
by_country$data[[1]] # a tibble 57 x 7
# model fitting function
country_model <- function(df) {
  lm(life_expectancy ~ year, data = df)
}
#Anwendung der Funktion auf Länder_Liste
models <- map(by_country$data, country_model)
# füge die Modellparameter der Länderliste hinzu
by_country <- by_country %>%
  mutate(model = map(data, country_model))
head(by_country,1)
