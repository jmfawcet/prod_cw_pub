# Custom Family -----------------------------------------------------------

# Define the custom family
circle_mix <- custom_family(
  "circle_mix",
  dpars = c('mu', 'kappa', 'theta'),
  links = c('identity', 'log', 'logit'),
  ub = c(NA, NA, NA),
  type = "real"
)

stan_funs = "
  /* von Mises log-PDF of a single response
     * for kappa > 100 the normal approximation is used
     * for reasons of numerial stability
     * Args: 
     *   y: the response vector between -pi and pi 
     *   mu: location parameter vector
     *   kappa: precision parameter
     * Returns:  
     *   a scalar to be added to the log posterior 
     */ 
     real von_mises_real_lpdf(real y, real mu, real kappa) {
       if (kappa < 100) {
         return von_mises_lpdf(y | mu, kappa);
       } else {
         return normal_lpdf(y | mu, sqrt(1 / kappa));
       }
     }
     
    /* von Mises log-PDF of a response vector
     * for kappa > 100 the normal approximation is used
     * for reasons of numerial stability
     * Args: 
     *   y: the response vector between -pi and pi 
     *   mu: location parameter vector
     *   kappa: precision parameter
     * Returns:  
     *   a scalar to be added to the log posterior 
     */ 
     real von_mises_vector_lpdf(vector y, vector mu, real kappa) {
       if (kappa < 100) {
         return von_mises_lpdf(y | mu, kappa);
       } else {
         return normal_lpdf(y | mu, sqrt(1 / kappa));
       }
     }
  
  real uniform_circ_lpdf(real y, real mu) {
    return -log(2.0 * pi());
  }

  real circle_mix_lpdf(real y, real mu, real kappa, real theta) {
    real log_prob;
    
    log_prob = log_sum_exp(
                  bernoulli_lpmf(1 | theta) + 
                    von_mises_real_lpdf(y | mu, kappa),
                  bernoulli_lpmf(0 | theta) + 
                    uniform_circ_lpdf(y | mu));

    return(log_prob);
  }
"

stanvars = stanvar(scode = stan_funs, block = "functions")

# Load data ------------------------------------------------------------

comb_test_dat = readRDS('data/comb_test_dat.rds')

# Note that Label and WP effectively contrasts E1 (no label) to 
# E2 and E3 (labels before or concurrent). Importantly, it would seem
# that WP and LABEL are inverted such that WP == 1 is e1 and WP == 2
# are e2 and e3 and LABEL == Y is no label and LABEL == N is label;
# this can be verified in the data frame itself. We invert it here, but
# keep this note for transparency (results should be unchanged)
comb_test_dat = comb_test_dat %>%
  mutate(WP = 1 - WP, LABEL = ifelse(LABEL=='Y', 'N', 'Y'))

# ANOVAs ------------------------------------------------------------------

# Analyzing colour and position separately
comb_test_dat %>%
  pivot_longer(cols=c('raddiff_col', 'raddiff_pos')) -> comb_test_dat_comb_dv

sum_stat = comb_test_dat %>% 
  mutate(WP = factor(WP)) %>%
  group_by(WP, SID, COND) %>% 
  summarize(abs_diff_col=mean(abs_diff_col), abs_diff_pos=mean(abs_diff_pos))

ezANOVA(sum_stat, .(abs_diff_col), .(SID), .(COND), between=.(WP))
ezStats(sum_stat, .(abs_diff_col), .(SID), .(COND), between=.(WP))

ezANOVA(sum_stat, .(abs_diff_pos), .(SID), .(COND), between=.(WP))
ezStats(sum_stat, .(abs_diff_pos), .(SID), .(COND), between=.(WP))

# Comparing colour and position as an IV
comb_test_dat %>%
  pivot_longer(cols=c('abs_diff_col', 'abs_diff_pos')) -> comb_test_dat_comb_dv_deg

combdv_sum_stat = comb_test_dat_comb_dv_deg %>% 
  group_by(WP, SID, COND, name) %>% 
  summarize(abs_diff=mean(abs(value)))

ezANOVA(combdv_sum_stat, .(abs_diff), .(SID), .(COND, name), between=.(WP))
ezStats(combdv_sum_stat, .(abs_diff), .(SID), .(COND, name), between=.(WP))

# Combining colour and position
combdv_sum_stat = comb_test_dat_comb_dv_deg %>% 
  mutate(WP = factor(WP)) %>%
  group_by(WP, SID, COND) %>% 
  summarize(abs_diff=mean(abs(value)))

ezANOVA(combdv_sum_stat, .(abs_diff), .(SID), .(COND), between=.(WP))
ezStats(combdv_sum_stat, .(abs_diff), .(SID), .(COND), between=.(WP))

# Fit Von Mises Model -----------------------------------------------------

# Slightly modified data set for us in mixture models, etc. tid is for
# use as variable precision variable
comb_test_dat2 = comb_test_dat %>%
  mutate(tid=1:n()) %>%
  pivot_longer(cols=c(raddiff_col, raddiff_pos), names_to='DV',values_to='diff') %>%
  mutate(DV = ifelse(DV=='raddiff_col', 'col', 'pos')) %>%
  mutate(dnumber = 1:n())

# Priors used, assuming intercept is removed
vm1_cond.prior = c(prior(normal(0, 1.5), class = b, dpar = "kappa" ),
                   prior(normal(0, 1.5), class = sd, dpar='kappa'),
                   set_prior('lkj(4)', class = "cor"))

# Note that tid makes this variable precision, with
# similar encoding precision for colour and location in a trial,
# it is also possible to refit allowing differential attention to
# location and colour by using dnumber instead (which assigns
# a different value to each response in a trial).
comb_vm1_cond = brm(bf(diff~1,
                       kappa ~ LABEL:COND - 1 + (COND - 1 |SID) + (1|tid),
                       mu~0,
                       family=von_mises),
                    data=comb_test_dat2,
                    backend='cmdstanr',
                    chains = 8, iter = 6000, cores=8,
                    prior = vm1_cond.prior,
                    control=list(adapt_delta=.95))

# Fit the mixture model ---------------------------------------------------

# Priors used, assuming intercept is removed
m1_cond.prior = c(prior(normal(0, 1.5), class = b, dpar = "kappa" ),
                  prior(logistic(0, .5), class = b, dpar='theta'),
                  prior(normal(0, 1.5), class = sd, dpar='kappa'),
                  prior(normal(0, 1), class = sd, dpar='theta'),
                  set_prior('lkj(4)', class = "cor"))

# Note that tid makes this variable precision, with
# similar encoding precision for colour and location in a trial,
# it is also possible to refit allowing differential attention to
# location and colour by using dnumber instead (which assigns
# a different value to each response in a trial).
comb_m1 = brm(bf(diff~1,
                 kappa ~ LABEL:COND - 1 + (COND - 1 |P|SID) + (1|tid),
                 theta ~ LABEL:COND - 1 + (COND - 1 |P|SID),
                 mu~0, family=circle_mix),
              data=comb_test_dat2,
              backend='cmdstanr',
              chains = 8, iter = 6000, cores=8,
              stanvars = stanvars,
              prior = m1_cond.prior)

# To fit models inclusive of position/colour as a predictor, simply modify
# the above equation to LABEL:COND:name - 1. Also, the above models can take
# a long time to run; when originally fit, with a decent computer they could take 
# several days each. The model fits were too large to push to github, unfortunately.
# 
# These models exclude image as a random effect, because to include it made them
# run far too long, but including them that term is as easy as + (COND - 1 | I | IMAGE)