# LOggBe
LOggBe is a Python package for evolutionary analysis of bacterial core genes.

## Installation
To install LOggBe, follow these steps:

1. Clone the LOggBe repository:
   
```sh
git clone https://github.com/rolesucsd/LOggBe.git
```

2. Install the required dependencies:

```sh
pip install -r requirements.txt
conda create -n iqtree -c bioconda iqtree
```

## Usage
To run LOggBe, use the following command:

```sh
logbee /path/to/alignment_directory /path/to/output_directory --jobs 10
```

## Output 
LOggBe generates the following outputs:

* Filtered Alignment: A processed alignment file containing the filtered sequences.
* Phylogenetic Tree: A tree file representing the evolutionary relationships among the core genes.
* Rates Summary: A text file named "rates.txt" that provides the relative rate of each core gene.

The information from the rates summary can be combined with gene annotation data to draw conclusions about the speed of evolution for different types of genes.

![Example Output](LOggBe.png)

We appreciate your interest in LOggBe. If you encounter any issues or have suggestions for improvement, please don't hesitate to contact us.
