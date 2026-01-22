#!/bin/bash

# Check if a filename was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <filename>"
    exit 1
fi

# File to process
file="$1"

# Use sed to remove ~~...~~ patterns
# This will overwrite the file in-place
sed -i.bak 's/~~[^~]*~~//g' "$file"

echo "All ~~...~~ patterns removed. Original file saved as $file.bak"

