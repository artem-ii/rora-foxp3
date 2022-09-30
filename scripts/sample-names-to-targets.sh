#! /bin/zsh


RAW_DATA_PATH="/Volumes/artemii-4tb/data/microarray/treg-microarray/raw-data/"
ANALYSIS_PATH="/Users/artemii/Documents/RORaFoxp3-microarray-analysis/treg-microarray/"
CONFIG_PATH="${ANALYSIS_PATH}config/"
TARGETS_PATH="${CONFIG_PATH}targets.txt"

# Create targets file with microarray files names
ls $RAW_DATA_FOLDER > $TARGETS_PATH