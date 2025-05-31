#!/bin/bash

# --- Function to build Apptainer images ---
build_apptainer_image() {
  container_sif="$1"
  container_def="$2"
  algo_name="$3"
  image_name="${container_sif##*/}"

  # Check if the .sif image already exists
  if [ -f "$container_sif" ]; then
    # Prompt for force rebuild
    read -p "A .sif image for $algo_name already exists. Force rebuild? (y/N) " -n 1 -r
    echo    # Move to a new line after input
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      # Force rebuild
      apptainer build --force "$container_sif" "$container_def"
    else
      echo "Skipping rebuild for $algo_name image."
    fi
  else
    # Build the image (no existing image found)
    echo "Building $algo_name image..."
    apptainer build "$container_sif" "$container_def"
  fi
}

# --- Build algorithm images ---
# Loop through each directory (algorithm) in the "algorithms" folder
for algo_name in algorithms/*; do
  # echo "Building $algo_name image..."
  # Check if it's a directory and not the "base" folder
  if [ -d "$algo_name" ] && [ "${algo_name##*/}" != "base" ]; then 
    # Construct the full paths for the command
    container_sif="algorithms/${algo_name##*/}/container.sif"
    container_def="algorithms/${algo_name##*/}/container.def"
    
    # Build the algorithm image
    build_apptainer_image "$container_sif" "$container_def" "${algo_name##*/}"
  fi
done