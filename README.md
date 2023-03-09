# ProteoformQuant

ProteoformQuant is a Python tool for quantitative analysis of proteoforms from mass spectrometry data.

## Setup/Installation

To use ProteoformQuant, you need to create a new Conda environment using Mamba package manager.

### 1. Install Mamba

Mamba is a fast package manager for Conda environments. You can install it using Conda itself by running the following command:

# ProteoformQuant

ProteoformQuant is a Python tool for quantitative analysis of proteoforms from mass spectrometry data.

## Setup/Installation

To use ProteoformQuant, you can create a Conda environment using the provided 'environment.yml' file. This can be done using Conda directly, however we recommend using the Mamba package manager.

### 1. Install Conda (and Mamba)

If not already done, install Conda (https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html#regular-installation)

If you wish to use Manda to create the Conda environment (faster) you can install mamba by running the following command in a terminal:

```bash
conda install mamba -n base -c conda-forge
```

### 2. Install Conda (and Mamba)

Next, create the environment using either Conda or Mamba by running the following command in the folder

```bash
# With Conda
conda env create --file environment.yml
#With Manba
mamba env create --file environment.yml
```

## Usage

```bash
import foobar

# returns 'words'
foobar.pluralize('word')

# returns 'geese'
foobar.pluralize('goose')

# returns 'phenomenon'
foobar.singularize('phenomena')
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
