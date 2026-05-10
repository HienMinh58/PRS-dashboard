import os
import subprocess

def compute_ld_from_target(bfile_prefix, chromosomes, window_kb=1000, r2_threshold=0.1, out_dir="/app/results/ld_target"):
    """
    Computes LD (Linkage Disequilibrium) statistics from target genotype data using PLINK2.
    
    Args:
        bfile_prefix (str): Prefix of the target PLINK files (.bed, .bim, .fam).
        chromosomes (str or list): Chromosome(s) to process.
        window_kb (int): Max distance in kb between variants to compute LD.
        r2_threshold (float): Minimum r2 to report.
        out_dir (str): Output directory for LD reports.
        
    Returns:
        dict: A summary containing success status, output paths, and parameters used.
    """
    # Use fallback if /app doesn't exist (e.g. running locally outside Docker)
    if not os.path.exists("/app") and out_dir.startswith("/app"):
        out_dir = "." + out_dir[4:]
        
    os.makedirs(out_dir, exist_ok=True)
    out_prefix = os.path.join(out_dir, "target_ld")
    
    # Normalize chromosomes
    if chromosomes in (None, "", "1"):
        chrom_list = ["1"]
    elif chromosomes == "1-22":
        chrom_list = [str(c) for c in range(1, 23)]
    elif isinstance(chromosomes, (list, tuple)):
        chrom_list = [str(c) for c in chromosomes]
    else:
        chrom_list = [str(chromosomes)]
        
    cmd = [
        "plink2",
        "--bfile", bfile_prefix,
        "--r2-unphased",
        "--ld-window-kb", str(window_kb),
        "--ld-window-r2", str(r2_threshold),
        "--out", out_prefix
    ]
    
    if chrom_list == [str(c) for c in range(1, 23)]:
        cmd.append("--autosome")
    else:
        cmd.extend(["--chr", ",".join(chrom_list)])
        
    summary = {
        "success": False,
        "target_prefix": bfile_prefix,
        "chromosomes": chrom_list,
        "window_kb": window_kb,
        "r2_threshold": r2_threshold,
        "out_dir": out_dir,
        "out_files": [],
        "command": cmd,
        "error": None
    }
    
    try:
        # Check if plink2 is available, otherwise try plink
        try:
            subprocess.run(["plink2", "--version"], check=False, capture_output=True)
        except FileNotFoundError:
            cmd[0] = "plink"
            # Fallback to --r2 for PLINK 1.9
            if "--r2-unphased" in cmd:
                cmd[cmd.index("--r2-unphased")] = "--r2"

        process = subprocess.run(cmd, check=True, capture_output=True, text=True)
        summary["success"] = True
        
        # PLINK2 outputs .vcor, PLINK1.9 outputs .ld
        expected_out = f"{out_prefix}.vcor"
        if not os.path.exists(expected_out):
            if os.path.exists(f"{out_prefix}.ld"):
                expected_out = f"{out_prefix}.ld"
                
        summary["out_files"].append(expected_out)
        
    except FileNotFoundError:
        summary["error"] = "PLINK/PLINK2 executable not found in PATH."
    except subprocess.CalledProcessError as e:
        summary["error"] = e.stderr if e.stderr else e.stdout
        
    return summary
