#!/bin/bash

# Color definitions
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

download_if_needed() {
    local URL="$1"

    local SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    local FILENAME="$(basename "$URL")"
    local LOCAL_FILE="$SCRIPT_DIR/${FILENAME}"

    if [ ! -f "$LOCAL_FILE" ]; then
        printf "${BLUE}File '$FILENAME' does not exist locally. Downloading...${NC}\n"
        wget -qO "$LOCAL_FILE" "$URL"

        # Check if file is empty or very small (potential captcha/download failure)
        if [ ! -s "$LOCAL_FILE" ] || [ $(stat -f%z "$LOCAL_FILE" 2>/dev/null || stat -c%s "$LOCAL_FILE" 2>/dev/null) -lt 1024 ]; then
            printf "${YELLOW}WARNING: Download may have failed for ${RED}'$FILENAME'!${NC}\n"
            printf "${YELLOW}The file is empty or suspiciously small. This might be due to:${NC}\n"
            printf "${YELLOW}  - A captcha requirement${NC}\n"
            printf "${YELLOW}  - Network issues${NC}\n"
            printf "${YELLOW}  - Invalid URL or access restrictions${NC}\n"
            printf "${YELLOW}Please check the file and try downloading manually if needed.${NC}\n"
        else
            printf "${GREEN}Successfully downloaded '$FILENAME'${NC}\n"
        fi
    else
        printf "${BLUE}File '$FILENAME' exists locally. Not downloading${NC}\n"
    fi
}

download_if_needed "https://3d-asteroids.space/data/asteroids/models/i/25143_Itokawa_200k.obj"
download_if_needed "https://3d-asteroids.space/data/comets/models/67P_ESA_NAVCAM_Jul2015data_256k.obj"
download_if_needed "https://pdssbn.astro.umd.edu/holdings/nh-a-lorri_mvic-5-geophys-v1.0/data/shape_models/mu69_merged.obj"
download_if_needed "https://sbn.psi.edu/pds/shape-models/files/a8567/a8567.tab.obj"
download_if_needed "https://data.darts.isas.jaxa.jp/pub/hayabusa2/paper/Watanabe_2019/SHAPE_SFM_3M_v20180804.obj"
download_if_needed "https://sbn.psi.edu/pds/shape-models/files/RADAR/4179toutatis.tab.obj"
download_if_needed "https://sbn.psi.edu/pds/shape-models/files/hartley2/hartley2_2012_cart.obj"
