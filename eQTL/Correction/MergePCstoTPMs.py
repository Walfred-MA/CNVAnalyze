# Read the eigenvec file
with open('pca_results.eigenvec', 'r') as f:
    lines = f.readlines()

# Extract the headers from the second column
headers = [line.strip().split()[1] for line in lines]

# Extract the data excluding the first column
data = [line.strip().split()[2:] for line in lines]

# Transpose the data
transposed_data = list(zip(*data))

# Write the transposed data to a new CSV file
with open('transformed_pca_results.csv', 'w') as out:
    # Write headers
    out.write('\t'.join(headers) + '\n')
    
    # Write data rows
    for row in transposed_data:
        out.write('\t'.join(row) + '\n')
