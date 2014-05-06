# Quick How-to

1. Extract data using 'extractDataFromSkim.py'
2. Extract MC using 'extractMCFromSkim.py'
3. Extract MC systematics using 'extractMCSystematicsFromSkim.py' (can take a bit longer)
4. Merge samples in a single dataset, using 'mergeFiles.py -p <path to folder where theta files are located>'
5. Generate generator systematics (matching & scale) using 'computeGeneratorSystematics.py -p <path to folder where theta files are located>'
6. Finally, create the final root file used by theta for limit computation using the script 'createThetaFile.py -p <path to folder where theta files are located>' -o filename.root
