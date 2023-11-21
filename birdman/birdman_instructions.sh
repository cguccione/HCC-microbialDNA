To run Birdman the following steps are needed:

1. Install birdman conda env: /panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/conda_env/birdman_env.yml

2. Compile Stan model
- Ensure you have complied this model by activating birdman environment, opening python and running the following
```python
import cmdstanpy
cmdstanpy.CmdStanModel(stan_file='/panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/birdman/negative_binomial_stan/negative_binomial_single.stan')
``` 
- 
