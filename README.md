# Protein Sequence Descriptor Converter

This Python script converts protein sequences into numerical descriptors based on amino acid properties. It supports two descriptor types:

1. **VHSE (Vector of Hydrophobic, Steric, and Electronic properties)**: Based on the paper: https://doi.org/10.1002/bip.20296
2. **Z-scales**: Based on the paper: https://pubs.acs.org/doi/full/10.1021/jm9700575

## Requirements

- Python 3.x
- pandas
- scikit-learn

To install the required packages:

```bash
pip install pandas scikit-learn
```

## Files

- `descriptor_converter.py`: Main Python script
- `example.csv`: Input file containing protein sequences (one sequence per line, no header)
- `data_vhse.csv`: Output file with calculated descriptors

## Usage

1. Prepare your input file (example.csv) with one protein sequence per line
2. Run the script:

```bash
python descriptor_converter.py
```

3. The script will generate `data_vhse.csv` containing the numerical descriptors

## Customization

You can modify the descriptor type by changing the parameter in the script:

```python
# For VHSE descriptors (8 values per amino acid)
data_descriptors = get_descriptors(data_list, descriptor='vhse')

# For Z-scales descriptors (5 values per amino acid)
data_descriptors = get_descriptors(data_list, descriptor='zscales')
```

## Functions

- `assign_descriptor(sequence, descriptor)`: Converts a single protein sequence to descriptors
- `get_descriptors(data_list, descriptor)`: Processes a list of sequences
- `norm_features(data)`: Normalizes the descriptors using MinMaxScaler

## Output

The output file contains numerical descriptors where:
- Each row represents one protein sequence
- Each column represents a descriptor value
- For VHSE: 8 values per amino acid
- For Z-scales: 5 values per amino acid

## Example

Input protein sequence: "ACDE"
- With VHSE, this would generate 4×8 = 32 descriptor values
- With Z-scales, this would generate 4×5 = 20 descriptor values

## Citation

If using this tool for research, please cite the relevant descriptor papers:

- VHSE: Mei H, Liao ZH, Zhou Y, Li SZ. A new set of amino acid descriptors and its application in peptide QSARs. Biopolymers. 2005.
- Z-scales: Sandberg M, Eriksson L, Jonsson J, Sjöström M, Wold S. New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids. Journal of Medicinal Chemistry. 1998.
