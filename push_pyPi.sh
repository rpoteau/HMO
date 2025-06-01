#!/bin/bash

# ANSI colors
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
CYAN='\033[1;36m'
RED='\033[1;31m'
WHITE_BG_BLACK_TEXT='\033[30;47m'
RESET='\033[0m'

PYPROJECT="pyproject.toml"
DIST_DIR="dist"

# nice utility
print_padded_line_wbg() {
    # Usage: print_padded_line "your message" width
    local msg="$1"
    local width="$2"
    local msg_len=${#msg}
    local pad_len=0
    if (( msg_len < width )); then
        pad_len=$((width - msg_len))
        pad=$(printf '%*s' "$pad_len" "")
        msg="$msg$pad"
    fi
    echo -e "${WHITE_BG_BLACK_TEXT}${msg}${RESET}"
}
#
# Clear separator line
SEPARATOR_RAW="---------------------------------------------------------------------------------------------"
SEPARATOR_WIDTH=${#SEPARATOR_RAW}
SEPARATOR="${WHITE_BG_BLACK_TEXT}${SEPARATOR_RAW}${RESET}"

# Read current version from pyproject.toml
CURRENT_VERSION=$(grep "^version" $PYPROJECT | head -n1 | cut -d '"' -f2)
echo -e "$SEPARATOR"
print_padded_line_wbg "                    Current version in pyproject.toml: $CURRENT_VERSION" "$SEPARATOR_WIDTH"
echo -e "$SEPARATOR"
echo

# Check if any .tar.gz exists in dist/
if [ -d "$DIST_DIR" ]; then
    echo -e "$SEPARATOR"
    ARCHIVES=$(ls dist/*.tar.gz 2>/dev/null)
    if [ -z "$ARCHIVES" ]; then
        print_padded_line_wbg "${YELLOW}No tar.gz archive found.${RESET}" "$SEPARATOR_WIDTH"
    else
        print_padded_line_wbg "Archives found in dist/: $ARCHIVES" "$SEPARATOR_WIDTH"
    fi
    echo -e "$SEPARATOR"
else
    echo -e "$SEPARATOR"
    echo -e "${WHITE_BG_BLACK_TEXT}No dist/ directory.  ${RESET}"
    echo -e "$SEPARATOR"
fi
echo

# Get the latest published version on PyPI (optional)
PACKAGE_NAME=$(grep "^name" $PYPROJECT | head -n1 | cut -d '"' -f2)
echo -e "$SEPARATOR"
print_padded_line_wbg "Querying PyPI for $PACKAGE_NAME..." "$SEPARATOR_WIDTH"
echo -e "$SEPARATOR"
LATEST_PYPI=$(curl -s https://pypi.org/pypi/$PACKAGE_NAME/json | jq -r '.info.version')

if [ "$LATEST_PYPI" != "null" ]; then
    echo -e "${CYAN}Latest published version on PyPI:${RESET} ${YELLOW}$LATEST_PYPI${RESET}"
else
    echo -e "${RED}Package does not exist on PyPI (or PyPI error).${RESET}"
fi

# Ask whether to increment the version
echo -e "$SEPARATOR"
print_padded_line_wbg "Do you want to increment the $CURRENT_VERSION version? (y/n) "  "$SEPARATOR_WIDTH"
echo -e "$SEPARATOR"
read -r REPLY

if [[ "$REPLY" =~ ^[Yy]$ ]]; then
    echo -e "${CYAN}Which level? ([p]atch / [m]inor / [M]ajor)${RESET}"
    read -r LEVEL

    IFS='.' read -r MAJOR MINOR PATCH <<< "$CURRENT_VERSION"

    case $LEVEL in
        p)
            PATCH=$((PATCH + 1))
            ;;
        m)
            MINOR=$((MINOR + 1))
            PATCH=0
            ;;
        M)
            MAJOR=$((MAJOR + 1))
            MINOR=0
            PATCH=0
            ;;
        *)
            echo -e "${RED}Unknown type, version not modified.${RESET}"
            ;;
    esac

    NEW_VERSION="$MAJOR.$MINOR.$PATCH"
    echo -e "${GREEN}Updating version $CURRENT_VERSION to: $NEW_VERSION${RESET}"
    # Update pyproject.toml in place
    sed -i "s/^version = \".*\"/version = \"$NEW_VERSION\"/" $PYPROJECT
    echo "     - in  pyproject.toml   ... Done"
    # Update  hmo/__init__.py in place
    if grep -q "^__version__ *= *" hmo/__init__.py; then
        sed -i "s/^__version__ *= *.*/__version__ = \"$NEW_VERSION\"/" hmo/__init__.py
    else
        echo "__version__ = \"$NEW_VERSION\"" >> hmo/__init__.py
    fi
    echo "     - in  hmo/__init__.py  ... Done"
    # Update __last_update__ field automatically
    today=$(date +%Y-%m-%d)
    if grep -q "^__last_update__ *= *" hmo/__init__.py; then
        sed -i "s/^__last_update__ *= *.*/__last_update__ = \"$today\"/" hmo/__init__.py
    else
        echo "__last_update__ = \"$today\"" >> hmo/__init__.py
    fi
    echo "last update field in hmo/__init__.py  ... Set to $today"

    # --- GIT SECTION ---
    echo -e "$SEPARATOR"
    print_padded_line_wbg "Git commit and tag...  " "$SEPARATOR_WIDTH"
    echo -e "$SEPARATOR"
    echo

    git add -A
    print_padded_line_wbg "Git status before commit:" "$SEPARATOR_WIDTH"
    git status
    echo "Proceed with commit? (y/n)"
    read -r CONFIRM
    if [[ "$CONFIRM" =~ ^[Yy]$ ]]; then
        git commit -m "Bump version: $CURRENT_VERSION → $NEW_VERSION"
        git tag "v$NEW_VERSION"
        git push
        git push --tags
    else
        echo "Commit cancelled."
        exit 1
    fi
    echo

    echo -e "$SEPARATOR"
    print_padded_line_wbg "Removing old builds: rm -rf build dist *.egg-info" "$SEPARATOR_WIDTH"
    echo -e "$SEPARATOR"
    rm -rf build dist *.egg-info
    echo

    # Build the package
    echo -e "$SEPARATOR"
    print_padded_line_wbg "Building the package: python -m build" "$SEPARATOR_WIDTH"
    echo -e "$SEPARATOR"
    python -m build
    echo

    # Upload to PyPI
    echo -e "$SEPARATOR"
    print_padded_line_wbg "Uploading to PyPI: twine upload dist/*" "$SEPARATOR_WIDTH"
    echo -e "$SEPARATOR"
    twine upload dist/*
    echo

    # Reinstall in editable mode
    echo -e "$SEPARATOR"
    print_padded_line_wbg "Reinstalling in editable mode: pip install -e ." "$SEPARATOR_WIDTH"
    echo -e "$SEPARATOR"
    pip install -e .
    echo

    echo -e "$SEPARATOR"
    echo -e "${GREEN}🎉 Process completed!${RESET}"
    echo -e "$SEPARATOR"
else
    echo -e "${YELLOW}Version $CURRENT_VERSION kept unchanged.${RESET}"
fi

