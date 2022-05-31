library(magrittr)
library(tidyverse)
library(bayesplot)

setwd("~/Desktop/exam")

df = read_tsv("alignmaf_3.nexus.run1.p", skip=1)
burnin = df$Gen %>% 
  max() %>% 
  multiply_by(0.25) %>% 
  floor()


df2 = df %>% 
  filter(Gen > burnin) 

mcmc_trace(df2, pars=c("pi(A)",'pi(T)'))
mcmc_trace(df2,pars = 'pi(T)')
mcmc_trace(df2,pars='pi(G)')
mcmc_trace(df2,pars='pi(C)')
color_scheme_set("mix-blue-red")
mcmc_trace(df2, pars = c("pi(A)", "pi(T)","pi(G)","pi(C)"), 
           facet_args = list(ncol = 1, strip.position = "left"))

#plots of frequency posterior probability
mcmc_areas(
  df, 
  pars = c("pi(A)", "pi(T)"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

mcmc_areas(
  df, 
  pars = c("pi(G)","pi(C)"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

mcmc_dens(
  df2, 
  pars = c("r(A<->C)",'r(G<->T)','r(C<->T)',"r(A<->G)",'r(A<->T)','r(C<->G)'),
  facet_args = list(ncol = 2, strip.position = "left"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

mcmc_areas(
  df2, 
  pars = c("TL"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)
  
mcmc_areas(
  df2, 
  pars = c("pi(A)","pi(T)"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

mcmc_areas(
  df2, 
  pars = c("pi(A)","pi(T)"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)

mcmc_areas(
  df2, 
  pars = c("pi(A)","pi(T)"),
  prob = 0.8, # 80% intervals
  prob_outer = 0.99, # 99%
  point_est = "mean"
)
  
  
  
burnin = df$Gen %>% 
  max() %>% 
  multiply_by(0.25) %>% 
  floor()


df2 = df %>% 
  filter(Gen > burnin) 
  
mcmc_intervals(df2, prob_outer = 1)
mcmc_areas(df2, prob_outer = 1)

df2long = pivot_longer(df2, cols = c("CG", "AC"))

ggplot(df2long) + 
  geom_density(mapping=aes(x=value, fill=name), alpha=0.3) + 
  labs(x="Substitution rate")

ggplot(df2, aes(x=CG, y=AC)) + 
  geom_point(col="blue") + 
  geom_abline(intercept=0, slope=1, lty=2, col="red") +
  xlim(0,0.25) + 
  ylim(0,0.25) + 
  labs(x="CG rate", y ="AC rate") 

ggplot(df2, aes(x=CG, y= AC)) + 
  geom_hex(col="blue") + 
  geom_abline(intercept=0, slope=1, lty=2, col="red") +
  xlim(0,0.25) + 
  ylim(0, 0.25) + 
  labs(x="CG rate", y ="AC rate") 

df2 %>%
  nrow()
df2 %>%
  filter(AC>CG) %>%
  nrow()

df3 = df %>% 
  filter(Gen > 75000) %>%
  select(Codon_1st = `m{1}`,
         ) %>%
  pivot_longer(cols=c("Codon_1st"))

ggplot(df3) + 
  geom_density(mapping=aes(x=value, fill=name), alpha=0.3) + 
  labs(x="Relative substitution rate")

