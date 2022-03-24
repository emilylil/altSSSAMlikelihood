# altSSSAMlikelihood

# A simple simulation/estimation procedure to check if a stock assessment model can estimate the catch and proportions-at-age, or catch-at-age when the underlying distribution is "correct" or "incorrect."

Options for assumed distribution:
1. Multivariate Log Normal (Catch-at-age)
2. Univariate Log normal (Catch) and Multinomial (Proportions-at-age)
3. Univariate Log normal (Catch) and Dirichlet Multinomial (Proportions-at-age)

Wish list for additional distributions:
4. Univariate Log normal (Catch) and Logistic Normal (Proportions-at-age)
5. Univariate Log normal (Catch) and multivariate-Tweedie (Proportions-at-age)

Resources/Inspiration:
Thorson, James T., et al. "Model-based estimates of effective sample size in stock assessment models using the Dirichlet-multinomial distribution." Fisheries Research 192 (2017): 84-93.

Fisch, Nicholas, et al. "Assessing likelihoods for fitting composition data within stock assessments, with emphasis on different degrees of process and observation error." Fisheries Research 243 (2021): 106069.

Albertsen, Christoffer Moesgaard, Anders Nielsen, and Uffe Høgsbro Thygesen. "Choosing the observational likelihood in state-space stock assessment models." Canadian Journal of Fisheries and Aquatic Sciences 74.5 (2017): 779-789.
