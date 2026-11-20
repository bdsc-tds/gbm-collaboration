#!/bin/bash

# ==============================================================================
# CONFIGURATION & HELP
# ==============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG_FILE="${SCRIPT_DIR}/symlink_config.yml"
CSV_DIR_NAME="data"
MAX_FILENAME_LEN=250

RED='\033[0;31m'
BRED='\033[1;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'
CHECK="✔"
CROSS="✘"

# Default flags
CLEAN_LINKS=false
LINK_CAPTIONS=true
SELECTED_KEY=""

# Argument parsing
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      echo -e "${BLUE}${BOLD}Figure Organizer Help${NC}"
      echo -e "Config File: ${YELLOW}$CONFIG_FILE${NC}"
      echo -e "\nUsage:"
      echo -e "  $0 [options]              Organize all panels."
      echo -e "  $0 [options] <panel_key>  Organize only a specific panel (e.g., figure1)."
      echo -e "\nOptions:"
      echo -e "  -c, --clean               Remove all symlinks from target folders before relinking."
      echo -e "  --nocaption               Do not link <stem>.caption.json sidecars (see README.md)."
      echo -e "  -h, --help                Show this help message."
      echo -e "\nEnvironment:"
      echo -e "  BASE_DIR=<dir>            Panel folder instead of the config's BASE_DIR."
      echo -e "\nYAML Syntax:"
      echo -e "  - path/to/file.pdf        (Links PDF and CSV if found)"
      echo -e "  - path/to/file.pdf --nocsv (Links PDF only)"
      exit 0
      ;;
    --nocaption)
      LINK_CAPTIONS=false
      shift
      ;;
    -c|--clean)
      CLEAN_LINKS=true
      shift
      ;;
    -*)
      echo -e "${RED}Unknown option: $1${NC}"
      echo -e "Use -h or --help for usage information."
      exit 1
      ;;
    *)
      if [[ -z "$SELECTED_KEY" ]]; then
        SELECTED_KEY="$1"
      fi
      shift
      ;;
  esac
done

if [ ! -f "$CONFIG_FILE" ]; then
    echo -e "${RED}Error: '$CONFIG_FILE' not found.${NC}"
    exit 1
fi

# Extract BASE_DIR
RAW_BASE=$(grep "^BASE_DIR:" "$CONFIG_FILE" | sed 's/^BASE_DIR://' | tr -d '"' | tr -d "'" | tr -d '\r' | xargs)
RAW_BASE="${BASE_DIR:-$RAW_BASE}"  # BASE_DIR=<dir> ./symlink_figures.sh overrides the config
CLEAN_BASE="${RAW_BASE%/}"

# ==============================================================================
# PROCESS FILE FUNCTION
# ==============================================================================
process_file() {
    local raw_line="$1"
    local panel_dir="$2"

    # 1. Clean Path String (remove the "- " prefix)
    local content=$(echo "$raw_line" | sed 's/^[[:space:]]*- //' | tr -d '\r' | xargs)
    content=$(echo "${content%%#*}" | xargs)   # drop trailing comment
    [ -z "$content" ] && return

    # 2. Check for --nocsv flag
    local skip_csv=false
    local item_path="$content"
    if [[ "$content" == *"--nocsv"* ]]; then
        skip_csv=true
        # Remove the flag from the path string and trim
        item_path=$(echo "${content//--nocsv/}" | xargs)
    fi

    # 3. Construct Paths & Link Names
    local full_source_path=""
    local link_name=""

    # Remove instances of "../" for clean link name generation
    local clean_item_path="${item_path//..\//}"

    if [[ "$item_path" == /* ]]; then
        # ABSOLUTE PATH
        full_source_path="$item_path"
        link_name=$(basename "$item_path")
    else
        # RELATIVE PATH (readlink -f will handle the actual ../ resolution later)
        full_source_path="${CLEAN_BASE}/${item_path}"
        link_name="${clean_item_path//\//_}" # Replace / with _
        link_name="${link_name#_}"           # Remove leading _
    fi

    # 4. Filename Safety Check (truncate if too long)
    if [ ${#link_name} -gt $MAX_FILENAME_LEN ]; then
        local ext="${link_name##*.}"
        local base_no_ext="${link_name%.*}"
        link_name="${base_no_ext: -200}.${ext}"
    fi

    # 5. Verify Source File Existence
    local filename=$(basename "$full_source_path")
    if [ ! -f "$full_source_path" ]; then
        echo -e "    ${RED}${CROSS} File not found: ${filename}${NC}"
        echo -e "      ${YELLOW}Tried: $full_source_path${NC}"
        return
    fi

    # 6. Direct CSV Handling
    # If the target file is intrinsically a .csv file, put it directly in 'data/' and skip figure logic
    local file_ext="${filename##*.}"
    if [[ "${file_ext,,}" == "csv" ]]; then
        mkdir -p "${panel_dir}/${CSV_DIR_NAME}"
        local csv_link_path="${panel_dir}/${CSV_DIR_NAME}/${link_name}"
        local abs_csv_target=$(readlink -f "$full_source_path")

        if ln -sf "$abs_csv_target" "$csv_link_path"; then
            echo -e "    ${GREEN}${CHECK} Linked CSV (direct): ${link_name}${NC}"
        else
            echo -e "    ${RED}${CROSS} Failed to link CSV: ${link_name}${NC}"
        fi
        return # Exit the function, skip the figure logic below
    fi

    # 7. CSV Logic for Figure Files
    local csv_status=""
    local csv_exists=false
    local stem="${filename%.*}"
    local source_dir=$(dirname "$full_source_path")
    local csv_source="${source_dir}/${stem}.csv"

    if [ "$skip_csv" = true ]; then
        csv_status="${YELLOW} (csv skipped)${NC}"
    elif [ -f "$csv_source" ]; then
        csv_status="${GREEN} (+ csv ${CHECK} )${NC}"
        csv_exists=true
    else
        csv_status="${BRED} (${CROSS} csv missing)${NC}"
    fi

    # 8. Create Links for Figures
    local fig_link_path="${panel_dir}/${link_name}"
    local abs_fig_target=$(readlink -f "$full_source_path")

    if ln -sf "$abs_fig_target" "$fig_link_path"; then
        if [ "$csv_exists" = true ]; then
            mkdir -p "${panel_dir}/${CSV_DIR_NAME}"
            local link_stem="${link_name%.*}"
            local csv_link_path="${panel_dir}/${CSV_DIR_NAME}/${link_stem}.csv"
            local abs_csv_target=$(readlink -f "$csv_source")
            ln -sf "$abs_csv_target" "$csv_link_path"
        fi
        # caption sidecar (README.md), linked next to the figure
        local cap_source="${source_dir}/${stem}.caption.json"
        if [ "$LINK_CAPTIONS" = true ] && [ -f "$cap_source" ]; then
            ln -sf "$(readlink -f "$cap_source")" "${panel_dir}/${link_name%.*}.caption.json"
            csv_status="${csv_status}${GREEN} (+ caption)${NC}"
        fi
        echo -e "    ${GREEN}${CHECK} Linked: ${link_name}${csv_status}"
    else
        echo -e "    ${RED}${CROSS} Failed to link: ${link_name}${NC}"
    fi
}

# ==============================================================================
# MAIN PARSER
# ==============================================================================
echo -e "\n${BLUE}${BOLD}=== Figure Organizer ===${NC}"
[ -n "$CLEAN_BASE" ] && echo -e "${BLUE}Base Dir: $CLEAN_BASE${NC}"

CURRENT_KEY=""
MATCH_FOUND=0
IS_ACTIVE=0

while IFS= read -r line || [ -n "$line" ]; do
    # Skip comments and empty lines
    [[ "$line" =~ ^[[:space:]]*# ]] && continue
    [[ -z "${line// /}" ]] && continue
    [[ "$line" =~ ^BASE_DIR: ]] && continue

    # Detect Figure Sections (e.g., figure1:)
    if [[ "$line" =~ ^([a-zA-Z0-9_-]+): ]]; then
        CURRENT_KEY="${BASH_REMATCH[1]}"
        if [[ -z "$SELECTED_KEY" || "$CURRENT_KEY" == "$SELECTED_KEY" ]]; then
            echo -e "\n${BLUE}${BOLD}Panel: $CURRENT_KEY${NC}"
            PANEL_BASE_DIR="${SCRIPT_DIR}/figures/$CURRENT_KEY"  # beside this script, whatever the current dir
            mkdir -p "$PANEL_BASE_DIR"

            # Clean Symlinks Option
            if [[ "$CLEAN_LINKS" == true ]]; then
                # Deletes all symlinks (-type l) inside the target directory and its subdirectories
                find "$PANEL_BASE_DIR" -type l -delete
                echo -e "    ${YELLOW}Cleaned previous symlinks.${NC}"
            fi

            IS_ACTIVE=1
            MATCH_FOUND=1
        else
            IS_ACTIVE=0
        fi
        continue
    fi

    # Detect List Items (lines starting with -)
    if [[ "$line" =~ ^[[:space:]]*- ]] && [[ "$IS_ACTIVE" -eq 1 ]]; then
        process_file "$line" "$PANEL_BASE_DIR"
    fi
done < "$CONFIG_FILE"

if [[ "$MATCH_FOUND" -eq 0 ]] && [[ -n "$SELECTED_KEY" ]]; then
    echo -e "\n${RED}Error: Key '$SELECTED_KEY' not found.${NC}"
fi

echo -e "\n${GREEN}${BOLD}Done!${NC}"