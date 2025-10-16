#!/bin/bash
# install_metamag.sh - User-friendly MetaMAG installation script

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║     MetaMAG Explorer Installation          ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════╝${NC}"
echo ""

# Function to check if conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        echo -e "${RED} Conda is not installed!${NC}"
        echo "Please install Miniconda first from: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    echo -e "${GREEN} Conda is installed${NC}"
}

# Function to create conda environment
create_environment() {
    echo -e "\n${BLUE}Creating MetaMAG conda environment...${NC}"
    
    if conda env list | grep -q "^metamag "; then
        echo -e "${BLUE}Environment 'metamag' already exists. Would you like to:${NC}"
        echo "1) Update existing environment"
        echo "2) Remove and recreate"
        echo "3) Skip environment creation"
        read -p "Enter choice (1-3): " choice
        
        case $choice in
            1)
                echo "Updating existing environment..."
                ;;
            2)
                echo "Removing existing environment..."
                conda env remove -n metamag -y
                conda create -n metamag python=3.9 -y
                ;;
            3)
                echo "Skipping environment creation..."
                return
                ;;
            *)
                echo "Invalid choice. Exiting."
                exit 1
                ;;
        esac
    else
        conda create -n metamag python=3.9 -y
        echo -e "${GREEN}✓ Created 'metamag' environment${NC}"
    fi
}

# Function to display installation menu
show_menu() {
    echo -e "\n${BLUE}Select installation option:${NC}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "1) Quick Install    - Essential tools only (FastQC, Fastp, Trimming)"
    echo "2) Standard Install - Core pipeline tools"
    echo "3) Full Install     - All tools including annotation"
    echo "4) Custom Install   - Choose specific modules"
    echo "5) Use Existing     - Detect tools already on system"
    echo "6) Verify Install   - Check installation status"
    echo "0) Exit"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
}

# Function to install tools
install_tools() {
    local option=$1
    
    echo -e "\n${BLUE}Activating metamag environment...${NC}"
    eval "$(conda shell.bash hook)"
    conda activate metamag
    
    case $option in
        1)  # Quick Install
            echo -e "${BLUE}Installing essential tools...${NC}"
            python setup_tools.py --tools fastqc fastp trimmomatic
            ;;
        2)  # Standard Install
            echo -e "${BLUE}Installing core pipeline tools...${NC}"
            python setup_tools.py --steps qc trimming host_removal assembly binning evaluation
            ;;
        3)  # Full Install
            echo -e "${BLUE}Installing all tools...${NC}"
            python setup_tools.py --all
            ;;
        4)  # Custom Install
            show_custom_menu
            ;;
        5)  # Use Existing
            echo -e "${BLUE}Scanning for existing tools...${NC}"
            python setup_tools.py --use-existing
            ;;
        6)  # Verify
            echo -e "${BLUE}Verifying installation...${NC}"
            python setup_tools.py --verify
            ;;
    esac
}

# Function to show custom installation menu
show_custom_menu() {
    echo -e "\n${BLUE}Select modules to install:${NC}"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "Available modules:"
    echo "  qc           - Quality control (FastQC, MultiQC)"
    echo "  preprocessing - Trimming and filtering (Fastp, BWA, Samtools)"
    echo "  assembly     - Assembly tools (IDBA, MEGAHIT, MetaQUAST)"
    echo "  binning      - Binning tools (MetaBAT2, MaxBin2, DAS Tool)"
    echo "  evaluation   - Quality assessment (CheckM2, dRep)"
    echo "  taxonomy     - Classification (GTDB-Tk, Kraken2, Bracken)"
    echo "  annotation   - Functional annotation (EggNOG, dbCAN, Prodigal)"
    echo ""
    read -p "Enter modules to install (space-separated): " modules
    
    if [ -n "$modules" ]; then
        echo -e "${BLUE}Installing selected modules: $modules${NC}"
        python setup_tools.py --update --steps $modules
    fi
}

# Main installation flow
main() {
    # Step 1: Check prerequisites
    check_conda
    
    # Step 2: Create/update environment
    create_environment
    
    # Step 3: Show menu and get user choice
    while true; do
        show_menu
        read -p "Enter your choice (0-6): " choice
        
        if [ "$choice" = "0" ]; then
            echo -e "${GREEN}Installation complete. Exiting.${NC}"
            break
        elif [ "$choice" -ge 1 ] && [ "$choice" -le 6 ]; then
            install_tools $choice
            
            # Post-installation message
            echo -e "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
            echo -e "${GREEN}Step completed!${NC}"
            echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
            
            read -p "Press Enter to continue..."
        else
            echo -e "${RED}Invalid choice. Please try again.${NC}"
        fi
    done
    
    # Final instructions
    echo -e "\n${BLUE}════════════════════════════════════════${NC}"
    echo -e "${GREEN}Installation Complete!${NC}"
    echo -e "${BLUE}════════════════════════════════════════${NC}"
    echo ""
    echo "To use MetaMAG Explorer:"
    echo "1. Activate the environment: conda activate metamag"
    echo "2. Run the pipeline: metamag --help"
    echo ""
    echo "Quick start example:"
    echo "  metamag run --samples samples.txt --steps qc trimming"
    echo ""
}

# Run main function
main
