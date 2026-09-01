# Latent Space Rating Relation Model (LSRRM)

This repository provides the source code for the Latent Space Rating Relation Model (LSRRM) developed in:

> Leng, C. H., Böckenholt, U., Lee, H. W., & Yao, G. (2025). Item response models for rating relational data. *Psychometrika, 90*(3), 1067–1096. https://doi.org/10.1017/psy.2025.10016

## Overview

The code in this repository implements the LSRRM for analyzing relational rating data proposed by Leng et al. (2025). The model jointly represents individual response tendencies and latent relational structures, allowing relational responses to be analyzed within an item response modeling framework.

## Getting Started

If you are using this repository for the first time, we recommend starting with `demo.R`.

The example script provides a basic workflow for fitting the model. To apply the model to your own data, modify the **Setting** section of `demo.R`.

In particular:

1. Import your relational response data and assign it to the object `y`.
2. Specify the model and estimation settings in `Args`.
3. Run the remaining sections of the script to estimate the model and obtain the corresponding outputs.

Additional modifications can be made to the model settings and estimation options according to the structure of the data and the purpose of the analysis.

## Citation

If you use this software, source code, or any part of the implementation in your research, please cite:

> Leng, C. H., Böckenholt, U., Lee, H. W., & Yao, G. (2025). Item response models for rating relational data. *Psychometrika, 90*(3), 1067–1096. https://doi.org/10.1017/psy.2025.10016

## License and Copyright

Copyright and intellectual property rights associated with this software belong to the Licensor and are protected under applicable intellectual property laws.

The software is provided for research and academic use subject to the terms specified by the Licensor. Redistribution, modification, or incorporation of the software into other projects should comply with the applicable license terms.

Use of this software in academic or research work requires appropriate citation of the associated publication:

> Leng, C. H., Böckenholt, U., Lee, H. W., & Yao, G. (2025). Item response models for rating relational data. *Psychometrika, 90*(3), 1067–1096. https://doi.org/10.1017/psy.2025.10016
