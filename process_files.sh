#!/bin/bash

# Loop through all .txt files in the current directory
for file in *.txt; do
    # Skip if no .txt files are found
    [ -e "$file" ] || continue  

    echo "Processing: $file"

    # Step 1: Remove lines that have no allele information
    trimmed_file="trimmed_$file"
    grep -v '\s-\s-' "$file" > "$trimmed_file"

    # Step 2: Remove non-autosomal chromosome data (keep only 1-22)
    pruned_file="pruned_$file"
    awk '$2 >= 1 && $2 <= 22' "$trimmed_file" > "$pruned_file"

    echo "Processed file saved as: $pruned_file"
done

echo "Processing completed for all files."

