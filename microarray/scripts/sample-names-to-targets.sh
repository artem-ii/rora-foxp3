#! /bin/bash


RAW_DATA_PATH="/Volumes/artemii-4tb/data/microarray/treg-microarray/raw-data/"
CONFIG_PATH="/Volumes/artemii-4tb/data/microarray/treg-microarray/raw-data/"
TARGETS_PATH="${CONFIG_PATH}targets.txt"

# Create targets file with microarray files names
ls $RAW_DATA_FOLDER > $TARGETS_PATH