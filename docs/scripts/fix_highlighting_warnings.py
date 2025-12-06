#!/usr/bin/env python3
"""
Fix syntax highlighting warnings in RST files.
Changes code blocks with mathematical symbols from 'python' to 'text'
"""

import re
from pathlib import Path

def fix_highlighting_in_file(filepath):
    """Fix code blocks that contain mathematical symbols."""
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original_content = content
    fixes_made = 0
    
    # Find all code blocks marked as python
    # Pattern: .. code-block:: python followed by content
    code_block_pattern = r'(\.\. code-block:: python\n(?:   :.*\n)*\n)((?:   .*\n)*)'
    
    def check_and_fix_block(match):
        nonlocal fixes_made
        header = match.group(1)
        block_content = match.group(2)
        
        # Check if block contains mathematical symbols
        math_symbols = ['×', '≈', '₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉', 
                        '→', '├', '└', '│', '✓', '✗', 'Δ', 'Σ', 'Π']
        
        has_math = any(symbol in block_content for symbol in math_symbols)
        
        if has_math:
            # Replace 'python' with 'text'
            new_header = header.replace('python', 'text')
            fixes_made += 1
            return new_header + block_content
        else:
            return match.group(0)
    
    # Apply the fix
    content = re.sub(code_block_pattern, check_and_fix_block, content)
    
    if content != original_content:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return fixes_made
    
    return 0

def main():
    """Fix all RST files in the documentation."""
    
    # Find the documentation directory
    docs_dirs = [
        Path("/mnt/user-data/outputs/HPC-ScaleTest-FINAL-COMPLETE-DOCS/docs/source"),
        Path("docs/source"),  # If run from package root
        Path("source"),       # If run from docs directory
    ]
    
    docs_dir = None
    for d in docs_dirs:
        if d.exists():
            docs_dir = d
            break
    
    if not docs_dir:
        print("❌ Could not find documentation directory")
        return 1
    
    print("="*80)
    print("Fixing Syntax Highlighting Warnings")
    print("="*80)
    print(f"\nScanning: {docs_dir}")
    
    # Find all RST files
    rst_files = list(docs_dir.rglob("*.rst"))
    print(f"Found {len(rst_files)} RST files")
    
    total_fixes = 0
    files_fixed = 0
    
    for rst_file in rst_files:
        fixes = fix_highlighting_in_file(rst_file)
        if fixes > 0:
            print(f"  ✓ {rst_file.name}: Fixed {fixes} code blocks")
            total_fixes += fixes
            files_fixed += 1
    
    print("\n" + "="*80)
    print(f"✅ Complete!")
    print(f"   Files modified: {files_fixed}")
    print(f"   Code blocks fixed: {total_fixes}")
    print("="*80)
    
    if total_fixes > 0:
        print("\nRecommendation: Rebuild documentation with 'make clean html'")
    
    return 0

if __name__ == "__main__":
    exit(main())

