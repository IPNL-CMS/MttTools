#! /bin/bash

# Start a new analysis

read -e -p "Enter the analysis name: " ANALYSIS_NAME
date=`LANG=C date +%d%b%y`

cat > analysis.json <<EOL
{
  "name": "${ANALYSIS_NAME}",
  "date": "${date}"
}
EOL
