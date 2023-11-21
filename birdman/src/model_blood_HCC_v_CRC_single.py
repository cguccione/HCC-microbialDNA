from pkg_resources import resource_filename

import biom
from birdman import SingleFeatureModel
import numpy as np
import pandas as pd

'''
Written by Hazel Dilmore (https://github.com/ahdilmore/MARS_Birdman/blob/main/) 
Modified by Caitlin edit
Goal: Run birdman on HCC data
'''

#Run specific variables
group_name = 'blood_HCC_v_CRC'
metric = 'tumor_type'
#Need to edit class below

#Note: See birdman_instructions.sh to ensure model is compiled
MODEL_PATH = "/panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/birdman/negative_binomial_stan/negative_binomial_single.stan"

MD = pd.read_table("/panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/processed_data/metadata/metadata_" + group_name + ".tsv", sep="\t", index_col='sample_name')

# NAME CLASS SOMETHING RELEVANT TO YOUR MODEL 
class blood_HCC_v_CRC_ModelSingle(SingleFeatureModel):
    def __init__(
        self,
        table: biom.Table,
        feature_id: str,
        # OPTIONAL: CHANGE PARAMETERS  
        beta_prior: float = 5.0,
        inv_disp_sd: float = 0.5,
        num_iter: int = 500,
        num_warmup: int = 500,
        **kwargs
    ):
        super().__init__(
            table=table,
            feature_id=feature_id,
            model_path=MODEL_PATH,
            num_iter=num_iter,
            num_warmup=num_warmup,
            **kwargs
        )


        D = table.shape[0]
        A = np.log(1 / D) 
	    # REPLACE WITH YOUR PATSY STYLE FORMULA 
        self.create_regression(formula=metric, metadata=MD)

        param_dict = {
            "depth": np.log(table.sum(axis="sample")),
            "B_p": beta_prior,
            "inv_disp_sd": inv_disp_sd,
	    "A": A
        }
        self.add_parameters(param_dict)

        self.specify_model(
            params=["beta_var", "inv_disp"],
            dims={
                "beta_var": ["covariate"],
                "log_lhood": ["tbl_sample"],
                "y_predict": ["tbl_sample"]
            },
            coords={
                "covariate": self.colnames,
                "tbl_sample": self.sample_names,
            },
            include_observed_data=True,
            posterior_predictive="y_predict",
            log_likelihood="log_lhood"

        )
