from padelpy import from_smiles

# calculate molecular descriptors for propane
descriptors = from_smiles('CCC')

# calculate molecular descriptors for propane and butane
# descriptors = from_smiles(['CCC', 'CCCC'])

# in addition to descriptors, calculate PubChem fingerprints
desc_fp = from_smiles('CCC', fingerprints=True)

# only calculate fingerprints
fingerprints = from_smiles('CCC', fingerprints=True, descriptors=False)

# setting the number of threads, this uses one cpu thread to compute descriptors
# descriptors = from_smiles(['CCC', 'CCCC'], threads=1)

# save descriptors to a CSV file
_ = from_smiles('CCC', output_csv='descriptors.csv')