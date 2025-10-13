#!/bin/bash
# Download small example dataset for testing MetaMAG Explorer
# This downloads 3 samples from each environment (9 total) for quick testing

echo "==============================================="
echo "MetaMAG Explorer - Example Dataset Downloader"
echo "==============================================="
echo ""
echo "This will download a small subset of data for testing:"
echo "- 3 human gut samples (~6 GB)"
echo "- 3 plant rhizosphere samples (~6 GB)"  
echo "- 3 rumen samples (~6 GB)"
echo "Total size: ~18 GB"
echo ""

# Create directories
mkdir -p example_data/human_gut
mkdir -p example_data/plant_rhizosphere
mkdir -p example_data/rumen

# Check if SRA toolkit is installed
if ! command -v prefetch &> /dev/null; then
    echo "ERROR: SRA toolkit not found!"
    echo "Please install from: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit"
    exit 1
fi

echo "Starting downloads..."
echo ""

# Download 3 human gut samples
echo "[1/3] Downloading human gut samples..."
for SRA in SRR9654993 SRR9654994 SRR9654995; do
    echo "  Downloading $SRA..."
    prefetch $SRA --output-directory example_data/human_gut
    fasterq-dump $SRA --split-files -O example_data/human_gut
    rm -rf example_data/human_gut/$SRA
done

# Download 3 plant samples  
echo "[2/3] Downloading plant rhizosphere samples..."
for SRA in SRR15883089 SRR15883091 SRR15883095; do
    echo "  Downloading $SRA..."
    prefetch $SRA --output-directory example_data/plant_rhizosphere
    fasterq-dump $SRA --split-files -O example_data/plant_rhizosphere
    rm -rf example_data/plant_rhizosphere/$SRA
done

# Download 3 rumen samples
echo "[3/3] Downloading rumen samples..."
for SRA in SRR14129954 SRR14129955 SRR14129957; do
    echo "  Downloading $SRA..."
    prefetch $SRA --output-directory example_data/rumen
    fasterq-dump $SRA --split-files -O example_data/rumen
    rm -rf example_data/rumen/$SRA
done

echo ""
echo "==============================================="
echo "Download complete!"
echo "Data saved in: example_data/"
echo "==============================================="
