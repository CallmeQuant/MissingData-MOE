# MissingData_MOE

## Overview

This repository contains the `SelvarMix` package, as described in [1]. Since the package is on archive and contains some deprecated parts (but not serious), some modifications are made in `SelvarMix_fix` file. 

## Installation

To install the `SelvarMix_fix` package in R, you need to have the `.tar.gz` file available on your local machine. You can install the package using the following command:

```R
install.packages(path/to/.tar.gzfile, repos = NULL, type = "source")
```

### Important Notes:

- **Correct File Path**: Make sure to adjust `"path/to.tar.gzfile"` to the actual location of the `.tar.gz` file on your local system. The file path you provide should point directly to where you saved the package locally.
  
  For example:
  - If the file is saved in `C:/Users/username/Documents/`, update the command as follows:
    ```R
    install.packages("C:/Users/username/Documents/SelvarMix_fix.tar.gz", repos = NULL, type = "source")
    ```

- **Dependencies**: If the package requires any dependencies that are not installed on your system, R will prompt you to install them. Ensure all required packages are installed for the `SelvarMix_fix` package to work correctly.

## Quarto File (`selvarclustmv3.qmd`)

The Quarto file (`selvarclustmv3.qmd`) is included in this repository to support experimenting with an extended version of `SelvarMix`. This file contains code, analysis, and notes related to ongoing experiments.

## Overview Structure

- **deprecated/**: Contains files that are no longer used in the current project.
- **SelvarMix_fix.tar.gz**: The tarball file used for package installation.
- **amputation.R**: Contains code for producing missing values mechanism.
- **selvarclustmv3.qmd**: Quarto file for experimenting with SelvarMix extensions.

## Reference

[1] Celeux, G., Maugis-Rabusseau, C. & Sedki, M. Variable selection in model-based clustering and discriminant analysis with a regularization approach. Adv Data Anal Classif 13, 259–278 (2019). https://doi.org/10.1007/s11634-018-0322-5
