#!/usr/bin/env python3
import sys
import os
import subprocess
import re
from dataclasses import dataclass
from typing import List, Tuple, Optional

# Global settings and paths
debug = False
USE_IUPAC = True
LOCALMAX_PATH ="/Users/ss61/Desktop/JT_code_autoseed_110325/localmax-motif-3" # Path to localmax-motif
EXTENDER_PATH ="/Users/ss61/Desktop/JT_code_autoseed_110325/motif_code_GPT/seedextender_flex_debug_final"  # Path to seedextender_flex
SPACEK_PATH = "/Users/ss61/Desktop/JT_code_autoseed_110325/spacek40"  # Path to spacek40
seed_history = []  # Global list for convergence checking

@dataclass
class LocalMaxResult:
    kmer: str
    gap: int
    position: int
    background: int
    signal: int
    is_localmax: bool
    shift_ratio: float
    ic: float
    avg_var: float
    exp_var: float

@dataclass
class ExtendedSeedResult:
    seed: str
    bg_count: int
    sig_count: int
    parent_count: int
    parent_kmer: str
    extension: str
    shift_ratio: float = 0.0
    ic: float = 0.0
    avg_var: float = 0.0
    exp_var: float = 0.0
    max_type: str = ""
    confidence: str = ""
    is_localmax: bool = True

def parse_localmax_output(output: str) -> List[LocalMaxResult]:
    results = []
    for line in output.splitlines():  # Don't skip first line
        if not line or line.startswith("Total") or line.startswith("kmer"):  # Skip actual headers
            continue
        fields = line.split()
        if len(fields) < 10:
            continue
        try:
            results.append(LocalMaxResult(
                kmer=fields[0],
                gap=int(fields[1]),
                position=int(fields[2]),
                background=int(fields[3]),
                signal=int(fields[4]),
                is_localmax="Localmax" in fields[5],
                shift_ratio=float(fields[6]) if fields[6] != "NA" else 0.0,
                ic=float(fields[7]),
                avg_var=float(fields[8].strip('%'))/100,
                exp_var=float(fields[9].strip('%'))/100
            ))
        except (ValueError, IndexError):
            continue
    return results

def parse_extender_output(output: str) -> List[ExtendedSeedResult]:
    results = []
    grouped_results = {}
    for line in output.splitlines():
        if "Full_kmer" in line or not line:
            continue
        fields = line.split()
        if len(fields) < 9:
            continue
        try:
            # Handle 'NA' values for numeric fields
            shift_ratio = 0.0 if fields[6] == 'NA' else float(fields[6])
            ic = 0.0 if fields[7] == 'NA' else float(fields[7].strip('%'))/100
            avg_var = 0.0 if fields[8] == 'NA' else float(fields[8].strip('%'))/100
            # Fixed the logic for exp_var
            exp_var = 0.0 if (len(fields) <= 9 or fields[9] == 'NA') else float(fields[9].strip('%'))/100

            res = ExtendedSeedResult(
                seed=fields[0],
                bg_count=int(fields[1]),
                sig_count=int(fields[2]),
                parent_count=int(fields[3]),
                parent_kmer=fields[4],
                extension=fields[5],
                shift_ratio=shift_ratio,
                ic=ic,
                avg_var=avg_var,
                exp_var=exp_var,
                max_type=fields[10] if len(fields) > 10 else "",
                confidence=fields[11] if len(fields) > 11 else "",
                is_localmax=True
            )
            if res.parent_kmer not in grouped_results or res.sig_count > grouped_results[res.parent_kmer].sig_count:
                grouped_results[res.parent_kmer] = res
        except (ValueError, IndexError) as e:
            if debug:
                print(f"Failed to parse line: {line}")
                print(f"Error: {e}")
            continue
    return sorted(grouped_results.values(), key=lambda x: x.sig_count, reverse=True)

def parse_motif_output(output: str) -> str:
    if debug:
        print("\nParsing PFM from motif output:")
        print(output)
    IUPAC_MAP = {1: 'A', 2: 'C', 3: 'M', 4: 'G', 5: 'R', 6: 'S', 7: 'V', 8: 'T', 9: 'W', 10: 'Y', 11: 'H', 12: 'K', 13: 'D', 14: 'B', 15: 'N'}
    matrix = []
    skipped = 0
    started = False
    for line in output.splitlines():
        if "Motif from all matches" in line:
            started = True
            if debug: print("\nFound motif start")
            continue
        if "Match statistics:" in line:
            break
        if started and line and not line.startswith("Total"):
            nums = line.split()
            if len(nums) > 1:
                if skipped == 0 and nums[0].isdigit() and int(nums[0]) == 1:
                    skipped = 1
                    continue
                values = [float(x) for x in nums[1:]]
                matrix.append(values)
    if len(matrix) != 4:
        if debug: print(f"Error: Got {len(matrix)} rows instead of 4 ACGT rows")
        return ""
    if debug:
        print("\nFinal matrix:")
        for row in matrix:
            print(row)
    
    # Start with reasonable cutoffs
    primary_cutoff = 0.4    # For single base
    secondary_cutoff = 0.5  # For two base redundant
    tertiary_cutoff = 0.6   # For three base redundant
    
    iter = 0
    while primary_cutoff <= 0.95 and primary_cutoff > 0 and iter < 50:
        if debug: print(f"\nTrying cutoffs: primary={primary_cutoff}, secondary={secondary_cutoff}, tertiary={tertiary_cutoff}")
        consensus = []
        for col in range(len(matrix[0])):
            col_vals = [row[col] for row in matrix]
            total = sum(col_vals)
            if total <= 0:
                consensus.append('N')
                if debug: print("All zeros -> N")
                continue
            
            # Normalize to get frequencies
            freqs = [val/total for val in col_vals]
            max_freq = max(freqs)
            max_idx = freqs.index(max_freq)
            
            # Sort frequencies in descending order
            sorted_freqs = sorted(freqs, reverse=True)
            
            # Calculate IUPAC based on new rules
            if max_freq >= primary_cutoff:
                # Use single base (most abundant > 50%)
                nuc_value = 1 << max_idx
            elif sorted_freqs[1] >= secondary_cutoff * max_freq:
                # Second most abundant base is >50% of max base, use two-base code
                indices = [i for i, freq in enumerate(freqs) if freq >= secondary_cutoff * max_freq]
                nuc_value = sum(1 << i for i in indices)
            elif sorted_freqs[2] >= tertiary_cutoff * max_freq:
                # Third most abundant base is >50% of max base, use three-base code
                indices = [i for i, freq in enumerate(freqs) if freq >= tertiary_cutoff * max_freq]
                nuc_value = sum(1 << i for i in indices)
            else:
                # No clear pattern, use N
                nuc_value = 15
            
            base = IUPAC_MAP.get(nuc_value, 'N')
            consensus.append(base)
        
        consensus = ''.join(consensus)
        # Define the low information characters
        low_info = set('BDHVN')

        # Trim left end until max 2 low info characters remain
        while len(consensus) >= 3:
            if (consensus[0] in low_info and
                consensus[1] in low_info and
                consensus[2] in low_info):
                consensus = consensus[1:]  # Remove just one character at a time
            else:
                break

        # Trim right end until max 2 low info characters remain
        while len(consensus) >= 3:
            if (consensus[-3] in low_info and
                consensus[-2] in low_info and
                consensus[-1] in low_info):
                consensus = consensus[:-1]  # Remove just one character at a time
            else:
                break
        
        ic = sum(2 if c in 'ACGT' else 1 if c in 'RYMKWS' else 0.415 if c in 'BDHV' else 0 for c in consensus)
        if debug: print(f"\nConsensus: {consensus} (IC: {ic:.2f}, Length excluding N: {len(consensus)-consensus.upper().count('N')})")
        
        if 8 <= ic <= 20:
            if debug: print("IC in range 8-20")
            if ic/(len(consensus)-consensus.upper().count('N')+0.001) >= 1:
                if debug: print("IC per base also acceptable, returning consensus")
                return consensus
        
        iter += 1
        if ic > 20:
            # Too high information content, need MORE redundancy
            primary_cutoff += 0.025  # Make it harder to use single bases
            secondary_cutoff -= 0.025  # Make it easier to use two-base codes
            tertiary_cutoff -= 0.025  # Make it easier to use three-base codes
        else:
            # Too low information content, need LESS redundancy
            primary_cutoff -= 0.025  # Make it easier to use single bases
            secondary_cutoff += 0.025  # Make it harder to use two-base codes
            tertiary_cutoff += 0.025  # Make it harder to use three-base codes
    
    if debug: print("Failed to find consensus with IC 8-20")
    return ""
   
def run_localmax_kmer(bg_file: str, signal_file: str, min_len: int, max_len: int,
                      count_threshold: int, length_cutoff: float, top_n: int) -> List[LocalMaxResult]:
    cmd = [LOCALMAX_PATH, bg_file, signal_file, str(min_len), str(max_len), str(count_threshold), str(length_cutoff)]
    full_cmd = f"{' '.join(cmd)}"
    if debug: print(f"\nRunning command: {full_cmd}")
    result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True, check=False)
    with open("short_seeds_filtered.tmp", 'w') as f:
        f.write(result.stdout)
    return parse_localmax_output(result.stdout)
    

def run_seedextender(bg_file: str, signal_file: str, seeds_file: str, max_extension: int,
                    min_count: int, length_cutoff: float, robust: bool = True) -> List[ExtendedSeedResult]:
   base_cmd = [EXTENDER_PATH, bg_file, signal_file, seeds_file, str(max_extension), str(min_count), str(length_cutoff)]
   if robust:
       base_cmd.append("robust")
       
   cmd_str = f"{' '.join(base_cmd)} | grep Localmax"
   
   if debug:
       print(f"\nRunning seedextender command: {cmd_str}")
       print(f"Input seeds file: {seeds_file}")
       os.system(f"cat {seeds_file}")
   print("Running seedextender with command:", cmd_str)

   result = subprocess.run(cmd_str, shell=True, capture_output=True, text=True, check=False)
   
   with open("extended_seeds.tmp", 'w') as f:
       f.write(result.stdout)

   return parse_extender_output(result.stdout)

def run_localmax_motif(bg_file: str, signal_file: str, seed: str, multinomial: int,
                       lambda_val: Optional[float] = None) -> Tuple[str, str]:
    cmd = [LOCALMAX_PATH, bg_file, signal_file, seed, str(multinomial)]
    if lambda_val is not None:
        cmd.append(str(lambda_val))
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    with open("motif.tmp", 'w') as f:
        f.write(result.stdout)
    return result.stdout, parse_motif_output(result.stdout)

def seeds_converged(old_seed: str, new_seed: str) -> bool:
    global seed_history
    if new_seed == old_seed:
        return True
    if new_seed in seed_history:
        return True
    seed_history.append(old_seed)
    return False

def create_merged_svg(final_seeds, seed_counts):
    """Create merged SVG with a table showing seeds, counts and percentages:
                        Seed              Background                Signal
       Original     CAAC (seed1)       100 (0.33%)           2000 (0.33%)
       Refined      TCAACC (seed2)      82 (0.25%)           3200 (0.45%)
    """
    def get_text_width(text: str, font_size: int, font_family: str = "Courier") -> int:
        """Estimate text width in pixels. For monospace fonts only."""
        if font_family == "Courier":
            return len(text) * (font_size * 0.6)  # Courier character width is ~0.6 of font size
        raise ValueError("Width calculation only supported for Courier font")

    def svg_table(data: List[List[str]], column_widths: List[int], row_height: int,
                 start_x: int, start_y: int, font_size: int = 20,
                 font_family: str = "Courier", font_color: str = "black",
                 show_grid: bool = False, grid_color: str = "#CCCCCC",
                 h_padding: int = 15, alignments: List[str] = None) -> str:
        if not data or not data[0]:
            return ""
        
        n_cols = len(data[0])
        n_rows = len(data)
        if len(column_widths) != n_cols:
            raise ValueError("Number of column widths must match number of columns in data")
        
        if alignments is None:
            alignments = ["left"] + ["right"] * (n_cols - 1)
        
        x_positions = []
        curr_x = start_x
        for width in column_widths:
            curr_x += width
            x_positions.append(curr_x)
        
        svg = []
        
        if show_grid:
            curr_x = start_x
            for width in column_widths:
                curr_x += width
                svg.append(
                    f'<line x1="{curr_x}" y1="{start_y-row_height/2}" '
                    f'x2="{curr_x}" y2="{start_y + (n_rows-0.5)*row_height}" '
                    f'stroke="{grid_color}" stroke-width="1"/>'
                )
            
            for row_idx in range(n_rows+1):
                y = start_y + (row_idx * row_height) - row_height/2
                svg.append(
                    f'<line x1="{start_x}" y1="{y}" '
                    f'x2="{x_positions[-1]}" y2="{y}" '
                    f'stroke="{grid_color}" stroke-width="1"/>'
                )
        
        for row_idx, row in enumerate(data):
            y = start_y + (row_idx * row_height)
            for col_idx, cell in enumerate(row):
                if alignments[col_idx] == "left":
                    x = x_positions[col_idx] - column_widths[col_idx] + h_padding
                    anchor = "start"
                elif alignments[col_idx] == "center":
                    x = x_positions[col_idx] - column_widths[col_idx]/2
                    anchor = "middle"
                else:  # right
                    x = x_positions[col_idx] - h_padding
                    anchor = "end"
                
                if row_idx == 0:
                    x = x_positions[col_idx] - column_widths[col_idx]/2
                    anchor = "middle"
                
                svg.append(
                    f'<text x="{x}" y="{y}" '
                    f'font-size="{font_size}" font-family="{font_family}" '
                    f'text-anchor="{anchor}" fill="{font_color}">{cell}</text>'
                )
        
        return "\n".join(svg)

    # Font settings
    title_font_size = 24
    table_font_size = 20
    filename_font_size = 16
    font_family = "Courier"
    title_color = "black"
    table_color = "black"
    filename_color = "gray"

    # Vertical spacing and dimensions
    initial_offset = 20     # Start 20 points lower
    logo_height = 120       # Height for each motif logo
    extra_spacing = 30      # Additional space between motif groups
    group_height = logo_height + extra_spacing
    merged_width = 3000     # Total SVG width

    # Horizontal positioning
    left_margin = 10        # Left edge margin
    table_position = 550    # Where table text begins
    
    # Vertical positioning within each group
    title_offset = 0        # Distance from group top to title
    header_offset = 20      # Distance from group top to table header
    filename_offset = logo_height-15  # Distance from group top to filename

    # Initialize SVG
    merged_svg_filename = "merged_logos.svg"
    current_y = initial_offset

    merged_svg = (
        '<?xml version="1.0" standalone="no"?>\n'
        '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" '
        '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n'
        f'<svg version="1.1" xmlns="http://www.w3.org/2000/svg" '
        f'width="{merged_width}" height="{len(final_seeds) * group_height + initial_offset}">\n'
    )

    # Process each motif and create its group
    for i, (orig_seed, refined_seed, pfm_filename) in enumerate(final_seeds, start=1):
        svg_filename = pfm_filename.replace(".pfm", ".svg")
        
        try:
            with open(svg_filename, "r") as f:
                content = f.read()
        except Exception as e:
            print(f"Error reading {svg_filename}: {e}")
            continue

        # Clean SVG tags
        content = re.sub(r'<\?xml.*?\?>', '', content, flags=re.DOTALL)
        content = re.sub(r'<!DOCTYPE.*?>', '', content, flags=re.DOTALL)
        content = re.sub(r'<svg[^>]*>', '', content, flags=re.DOTALL)
        content = re.sub(r'</svg>', '', content, flags=re.DOTALL)
        content = content.strip()

        # Get statistics
        orig_sig = seed_counts.get(orig_seed, {}).get('signal', 0)
        orig_bg = seed_counts.get(orig_seed, {}).get('background', 0)
        orig_sig_pct = seed_counts.get(orig_seed, {}).get('signal_percent', 0)
        orig_bg_pct = seed_counts.get(orig_seed, {}).get('background_percent', 0)
        
        refined_sig = seed_counts.get(refined_seed, {}).get('signal', 0)
        refined_bg = seed_counts.get(refined_seed, {}).get('background', 0)
        refined_sig_pct = seed_counts.get(refined_seed, {}).get('signal_percent', 0)
        refined_bg_pct = seed_counts.get(refined_seed, {}).get('background_percent', 0)

        # Calculate column widths based on content
        seed_width = max(
            get_text_width("Seed", table_font_size, font_family),
            get_text_width(f"Original: {orig_seed}", table_font_size, font_family),
            get_text_width(f"Refined:  {refined_seed}", table_font_size, font_family)
        )

        bg_width = max(
            get_text_width("Background", table_font_size, font_family),
            get_text_width(f"{orig_bg} ({orig_bg_pct:.2f}%)", table_font_size, font_family),
            get_text_width(f"{refined_bg} ({refined_bg_pct:.2f}%)", table_font_size, font_family)
        )

        sig_width = max(
            get_text_width("Signal", table_font_size, font_family),
            get_text_width(f"{orig_sig} ({orig_sig_pct:.2f}%)", table_font_size, font_family),
            get_text_width(f"{refined_sig} ({refined_sig_pct:.2f}%)", table_font_size, font_family)
        )

        extra_width = 20
        column_widths = [seed_width + extra_width, bg_width + extra_width, sig_width + extra_width]

        # Format title
        title_text = f"{orig_seed} -> {refined_seed}"

        # Format data for table
        data = [
            ["Seed", "Background", "Signal"],
            [f"Original: {orig_seed}",
             f"{orig_bg} ({orig_bg_pct:.2f}%)",
             f"{orig_sig} ({orig_sig_pct:.2f}%)"],
            [f"Refined : {refined_seed}",
             f"{refined_bg} ({refined_bg_pct:.2f}%)",
             f"{refined_sig} ({refined_sig_pct:.2f}%)"]
        ]

        # Add group to SVG
        merged_svg += f'  <g transform="translate({left_margin}, {current_y})">\n'
        merged_svg += f'    <g transform="translate(0, 0)">\n{content}\n    </g>\n'
        
        # Add text elements
        merged_svg += (
            f'    <text x="{left_margin}" y="{title_offset}" '
            f'font-size="{title_font_size}" font-family="{font_family}" '
            f'fill="{title_color}">{title_text}</text>\n'
        )

        # Add the table
        merged_svg += svg_table(
            data=data,
            column_widths=column_widths,
            row_height=25,
            start_x=table_position,
            start_y=header_offset,
            font_size=table_font_size,
            font_family=font_family,
            font_color=table_color,
            show_grid=False,
            h_padding=10,
            alignments=["left", "right", "right"]
        ) + "\n"

        # Add filename
        merged_svg += (
            f'    <text x="{left_margin}" y="{filename_offset}" '
            f'font-size="{filename_font_size}" font-family="{font_family}" '
            f'fill="{filename_color}">{pfm_filename}</text>\n'
        )

        merged_svg += "  </g>\n"
        current_y += group_height

    merged_svg += "</svg>\n"

    # Save the SVG
    with open(merged_svg_filename, "w") as out:
        out.write(merged_svg)

    if debug:
        print("Merged SVG generated:", merged_svg_filename)
        print(f"Total height: {current_y}")

def check_motif_similarity(pfm1: str, pfm2: str, spacing: int = 6) -> float:
    """Check similarity between two PFM files using spacek40.
    Returns similarity score or 0.0 if error."""
    try:
        cmd = f"{SPACEK_PATH} --dist -o=similarity.tmp {pfm1} {pfm2} gapped {spacing}"
        if debug:
            print(f"Running similarity command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if debug:
            print("Spacek output:", result.stdout)
            if result.stderr:
                print("Spacek stderr:", result.stderr)
        for line in result.stdout.splitlines():
            if "gapped" in line:
                # Extract similarity score
                score = line.split('\t')[1]  # Use tab as separator
                if debug:
                    print(f"Found similarity line: {line}")
                    print(f"Extracted score: {score}")
                return 0.0 if 'e' in score.lower() else float(score)
    except Exception as e:
        if debug:
            print(f"Error checking motif similarity: {e}")
    return 0.0
    
    
def run(args: List[str]):
    global debug, seed_history
    if len(args) < 6:
        print("Usage: ./autoseed.py <background_file> <signal_file> <shortest_kmer_length> "
              "<longest_kmer_length> <top_N_seeds> <max_extension> [count_threshold (default: 100)] "
              "[length_diff_cutoff (default: 0.35)] [multinomial (default: 1)] "
              "[similarity_cutoff (default: 0.20)]")
        return
    bg_file = args[0]
    signal_file = args[1]
    min_len = int(args[2])
    max_len = int(args[3])
    top_n = int(args[4])
    max_extension = int(args[5])
    count_threshold = int(args[6]) if len(args) > 6 else 100
    length_cutoff = float(args[7]) if len(args) > 7 else 0.35
    multinomial = int(args[8]) if len(args) > 8 else 1
    similarity_cutoff = float(args[9]) if len(args) > 9 else 0.20
    
    if debug:
        print("\n=== Initial Parameters ===")
        print(f"Background file: {bg_file}")
        print(f"Signal file: {signal_file}")
        print(f"Kmer length range: {min_len}-{max_len}")
        print(f"Top N seeds: {top_n}")
        print(f"Max extension: {max_extension}")
        print(f"Count threshold: {count_threshold}")
        print(f"Length cutoff: {length_cutoff}")
        print(f"Multinomial: {multinomial}")
        print("\n=== Starting Analysis ===")
    
    print("\n1. Running localmax-motif to find short seeds...")
    initial_results = run_localmax_kmer(bg_file, signal_file, min_len, max_len, count_threshold, length_cutoff, top_n+1)
    
    print(f"Found {len(initial_results)} initial seeds")
    if not initial_results:
        print("No initial seeds found. Exiting.")
        return
        
    if debug:
        print("\nTop seeds by signal count:")
        for i, res in enumerate(sorted(initial_results, key=lambda x: x.signal, reverse=True), 1):
            print(f"  {i}. {res.kmer} (sig: {res.signal}, bg: {res.background})")
    
    filtered_results = sorted(initial_results, key=lambda x: x.signal, reverse=True)[:top_n]
    with open("short_seeds_filtered.tmp", 'w') as f:
        for result in filtered_results:
            f.write(f"{result.kmer}\t{result.gap}\t{result.position}\t{result.background}\t{result.signal}\t{result.is_localmax}\n")
    
    print(f"\n2. Running seedextender...")
    extended_results = run_seedextender(bg_file, signal_file, "short_seeds_filtered.tmp", max_extension, count_threshold, length_cutoff, False)
    
    print(f"Found {len(extended_results)} extended seeds")
    if not extended_results:
        print("No extended seeds found. Exiting.")
        return
    
    seed_counts = {}    # Initialize the dictionary here
    
    if debug:
        print("\nStored counts:")
        for seed, counts in seed_counts.items():
            print(f"  {seed}: sig={counts['signal']}, bg={counts['background']}")
        print("\nSeed    Extended to")
        for i, res in enumerate(extended_results, 1):
            print(f"  {i}. {res.parent_kmer} -> {res.seed} (sig: {res.sig_count}, bg: {res.bg_count})")
    
    print("\n3. Starting PWM refinement for each extended seed...")
    
    final_seeds = []
    accepted_pwms = []
    
    for idx, result in enumerate(extended_results, start=1):
        print("\n-----")
        print(f"Processing seed {idx}/{len(extended_results)}: {result.seed}")
        if debug:
            print(f"Initial seed: {result.seed}")
        current_seed = result.seed
        current_seed_history = []
        iter_count = 0
        final_motif_output = ""
        max_iterations = 20
        
        # Capture original seed match counts and percentages
        if iter_count == 0:
            orig_motif_output, _ = run_localmax_motif(bg_file, signal_file, result.seed, multinomial)
            for line in orig_motif_output.splitlines():
                if line.startswith("Background:"):
                    parts = line.split()
                    bg_count = int(parts[1])
                    bg_percent = float(parts[6].strip('()%'))
                    seed_counts[result.seed] = {'signal': 0, 'background': bg_count, 'background_percent': bg_percent}
                elif line.startswith("Signal:"):
                    parts = line.split()
                    sig_count = int(parts[1])
                    sig_percent = float(parts[6].strip('()%'))
                    seed_counts[result.seed].update({'signal': sig_count, 'signal_percent': sig_percent})

        while iter_count <= max_iterations:
            motif_output, new_seed = run_localmax_motif(bg_file, signal_file, current_seed, multinomial)

            # Extract counts and percentages
            for line in motif_output.splitlines():
                if line.startswith("Background:"):
                    parts = line.split()
                    bg_count = int(parts[1])
                    bg_percent = float(parts[6].strip('()%'))
                    # Update counts for the new seed
                    seed_counts[new_seed] = {'signal': 0, 'background': bg_count, 'background_percent': bg_percent}
                elif line.startswith("Signal:"):
                    parts = line.split()
                    sig_count = int(parts[1])
                    sig_percent = float(parts[6].strip('()%'))
                    # Update the existing dict with signal info
                    seed_counts[new_seed].update({'signal': sig_count, 'signal_percent': sig_percent})

                    
            print(f"Iteration {iter_count+1}: Seed change {current_seed} -> {new_seed}", flush=True)
            if iter_count == max_iterations:
                print(f"Seed error: failed to converge in {iter_count+1} iterations")
                final_motif_output = ""
                break
            if not new_seed:
                print("Seed error: Empty seed")
                final_motif_output = ""
                break
            if iter_count > 0 and seeds_converged(current_seed, new_seed):
                print(f"Converged after {iter_count+1} iterations")
                final_motif_output = motif_output
                break
            current_seed_history.append(current_seed)
            current_seed = new_seed
            iter_count += 1
            final_motif_output = motif_output
        if not final_motif_output or not new_seed:
            print(f"Seed error for seed {result.seed}: no valid consensus found.")
            continue
        pwm_matrix = []
        skipped = 0
        started = False
        for line in final_motif_output.splitlines():
            if "Motif from all matches" in line:
                started = True
                continue
            if "Match statistics:" in line:
                break
            if started and line and not line.startswith("Total"):
                nums = line.split()
                if len(nums) > 1:
                    if skipped == 0 and nums[0].isdigit() and int(nums[0]) == 1:
                        skipped = 1
                        continue
                    pwm_matrix.append(nums[1:])
        
        if len(pwm_matrix) == 4:
            original_seed = current_seed_history[0] if current_seed_history else result.seed
            filename = f"refined_model_{idx}_for_seed_{original_seed}.pfm"
            
            # Write current PWM
            with open(filename, "w") as f:
                for row in pwm_matrix:
                    f.write("\t".join(row) + "\n")
                    
            # Initialize similarity check result
            is_similar = False
            
            # Only check similarity if we have previous PWMs
            if accepted_pwms:
                for prev_pwm in accepted_pwms:
                    if prev_pwm == filename:  # Skip self-comparison
                        continue
                    similarity = check_motif_similarity(prev_pwm, filename)
                    if debug:
                        print(f"Comparing {filename} with {prev_pwm}: similarity = {similarity}")
                    if similarity > similarity_cutoff:
                        if debug:
                            print(f"Rejecting motif - similarity {similarity:.3f} to {prev_pwm}")
                        is_similar = True
                        break
            
            # Handle acceptance or rejection
            if not is_similar:
                print(f"Seed {original_seed} -> {new_seed} ACCEPTED as model {len(accepted_pwms)+1}")
                accepted_pwms.append(filename)
                svg_filename = filename.replace(".pfm", ".svg")
                cmd = f"{SPACEK_PATH} -noname --logo {filename} {svg_filename}"
                subprocess.run(cmd, shell=True)
                print(f"Generated SVG logo: {svg_filename}")
                final_seeds.append((original_seed, new_seed, filename))
            else:
                print(f"Seed {original_seed} REJECTED due to similarity with previous model")
                os.remove(filename)  # Clean up rejected PWM file
                
        else:
            print(f"PWM matrix extraction failed for seed {result.seed}")
            
    print("\n=== Final Results ===")
    for original, refined, fname in final_seeds:
        print(f"{original}\t{refined}\t{fname}")
     
    if debug:
        print("\n=== Debug: All Seed Counts ===")
        for seed, counts in seed_counts.items():
            print(f"{seed}: signal={counts['signal']} ({counts.get('signal_percent', 0):.4f}%), "
                  f"background={counts['background']} ({counts.get('background_percent', 0):.4f}%)")

        print("\n=== Debug: Final Seeds List ===")
        for original, refined, fname in final_seeds:
            orig_counts = seed_counts[original]
            refined_counts = seed_counts[refined]
            print(f"{original} -> {refined}: {fname}")
            print(f"  Original: bg={orig_counts['background']} ({orig_counts.get('background_percent', 0):.4f}%), "
                  f"sig={orig_counts['signal']} ({orig_counts.get('signal_percent', 0):.4f}%)")
            print(f"  Refined:  bg={refined_counts['background']} ({refined_counts.get('background_percent', 0):.4f}%), "
                  f"sig={refined_counts['signal']} ({refined_counts.get('signal_percent', 0):.4f}%)")

    print("\n----------------")
    if final_seeds:
        print(f"A TOTAL OF {len(final_seeds)} MODELS FOUND")
        print("\nCreating merged SVG file...")
        create_merged_svg(final_seeds, seed_counts)
    else:
        print("NO MODELS FOUND")
        
    # creates merged svg
    create_merged_svg(final_seeds, seed_counts)

if __name__ == "__main__":
    if "-debug" in sys.argv:
        debug = True
        sys.argv.remove("-debug")
    run(sys.argv[1:])

