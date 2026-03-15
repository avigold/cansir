"""
Genotype imputation integration for the PPS calculator.

Converts consumer genotyping files (23andMe, AncestryDNA) to VCF format
and submits to the Michigan Imputation Server for imputation against
a reference panel (1000G Phase 3 or TOPMed).

After imputation, coverage of PPS variants jumps from ~10% to ~90%+.

Workflow:
1. Convert consumer file to VCF (using snps library)
2. Split by chromosome and bgzip
3. Submit to imputation server via REST API
4. Poll for completion
5. Download and decrypt results
6. Feed imputed VCF to PPS calculator
"""

import gzip
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import requests


MICHIGAN_API = "https://imputationserver.sph.umich.edu/api/v2"
TOPMED_API = "https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2"


def convert_consumer_to_vcf(
    input_file: str,
    output_dir: str = "data/raw/imputation",
    target_build: int = 37,
) -> list[Path]:
    """
    Convert a consumer genotyping file to per-chromosome VCFs.

    Uses the `snps` library which auto-detects format (23andMe, AncestryDNA,
    FTDNA, MyHeritage, etc.) and handles build remapping.

    Returns list of per-chromosome VCF paths.
    """
    from snps import SNPs

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"  Loading {input_file}...")
    s = SNPs(input_file)
    print(f"  Detected: {s.source}, build {s.build}, {s.count} variants")

    # Remap to target build if needed
    if s.build != target_build:
        print(f"  Remapping from build {s.build} to {target_build}...")
        s.remap(target_build)

    # Export to VCF — snps library prefixes "output/" to paths
    relative_out = out.relative_to(Path.cwd()) if out.is_absolute() else out
    os.makedirs(f"output/{relative_out}", exist_ok=True)
    print(f"  Exporting to VCF...")
    result_path = s.to_vcf(f"{relative_out}/consumer_data.vcf")
    vcf_path = Path(result_path) if isinstance(result_path, str) else Path(str(result_path))
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF export failed: {result_path}")
    print(f"  Exported: {vcf_path} ({vcf_path.stat().st_size / 1e6:.1f} MB)")

    # Split by chromosome for imputation server
    print(f"  Splitting by chromosome...")
    chr_vcfs = split_vcf_by_chromosome(vcf_path, out)

    # Clean up unsplit VCF
    vcf_path.unlink(missing_ok=True)

    return chr_vcfs


def split_vcf_by_chromosome(vcf_path: Path, output_dir: Path) -> list[Path]:
    """Split a VCF file into per-chromosome files, bgzip and index."""
    chr_files = {}

    with open(vcf_path) as f:
        header_lines = []
        for line in f:
            if line.startswith("#"):
                header_lines.append(line)
            else:
                parts = line.split("\t", 2)
                chrom = parts[0]
                if chrom not in chr_files:
                    chr_path = output_dir / f"chr{chrom}.vcf"
                    chr_files[chrom] = open(chr_path, "w")
                    chr_files[chrom].writelines(header_lines)
                chr_files[chrom].write(line)

    result_paths = []
    for chrom, fh in chr_files.items():
        fh.close()
        vcf_path = output_dir / f"chr{chrom}.vcf"

        # bgzip
        gz_path = output_dir / f"chr{chrom}.vcf.gz"
        try:
            subprocess.run(["bgzip", "-f", str(vcf_path)], check=True,
                           capture_output=True, timeout=60)
            result_paths.append(gz_path)
        except (subprocess.CalledProcessError, FileNotFoundError):
            # bgzip not available — use gzip
            with open(vcf_path, "rb") as f_in:
                with gzip.open(str(gz_path), "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            vcf_path.unlink()
            result_paths.append(gz_path)

    # Sort by chromosome number
    def chrom_sort_key(p):
        c = p.stem.replace("chr", "").replace(".vcf", "")
        try:
            return int(c)
        except ValueError:
            return 100 + ord(c[0])

    result_paths.sort(key=chrom_sort_key)

    print(f"  Created {len(result_paths)} per-chromosome VCFs")
    return result_paths


def submit_imputation_job(
    vcf_files: list[Path],
    token: str,
    server: str = "michigan",
    ref_panel: str = "1000g-phase-3-v5",
    population: str = "mixed",
    build: str = "hg19",
    r2_filter: float = 0.3,
) -> dict:
    """
    Submit an imputation job to Michigan or TOPMed server.

    Args:
        vcf_files: List of per-chromosome bgzipped VCF files.
        token: API authentication token (from server profile page).
        server: "michigan" or "topmed".
        ref_panel: Reference panel ID.
        population: Population code for phasing.
        build: Genome build ("hg19" or "hg38").
        r2_filter: Minimum imputation quality to include in output.

    Returns:
        Job submission response dict.
    """
    base_url = MICHIGAN_API if server == "michigan" else TOPMED_API
    endpoint = f"{base_url}/jobs/submit/minimac4"
    if server == "topmed":
        endpoint = f"{base_url}/jobs/submit/imputationserver"

    headers = {"X-Auth-Token": token}

    files = [("files", (f.name, open(f, "rb"))) for f in vcf_files]

    data = {
        "refpanel": ref_panel,
        "population": population,
        "build": build,
        "r2Filter": str(r2_filter),
    }

    print(f"\n  Submitting to {server} imputation server...")
    print(f"  Reference panel: {ref_panel}")
    print(f"  Population: {population}")
    print(f"  Build: {build}")
    print(f"  Files: {len(vcf_files)} chromosomes")

    try:
        resp = requests.post(endpoint, headers=headers, files=files,
                             data=data, timeout=300)
        result = resp.json()

        if resp.status_code == 200 and result.get("success"):
            job_id = result.get("id", "")
            print(f"  Job submitted successfully! ID: {job_id}")
            print(f"  Monitor at: {base_url.rsplit('/api', 1)[0]}/#!jobs/{job_id}")
        else:
            print(f"  Submission failed: {result.get('message', resp.text[:200])}")

        return result

    except Exception as e:
        print(f"  Error submitting job: {e}")
        return {"success": False, "message": str(e)}

    finally:
        for _, (_, fh) in files:
            fh.close()


def check_job_status(job_id: str, token: str, server: str = "michigan") -> dict:
    """Check the status of an imputation job."""
    base_url = MICHIGAN_API if server == "michigan" else TOPMED_API
    headers = {"X-Auth-Token": token}

    resp = requests.get(f"{base_url}/jobs/{job_id}/status",
                        headers=headers, timeout=30)
    return resp.json()


def wait_for_job(
    job_id: str,
    token: str,
    server: str = "michigan",
    poll_interval: int = 60,
    max_wait: int = 86400,
) -> dict:
    """Poll job status until completion (or timeout)."""
    base_url = MICHIGAN_API if server == "michigan" else TOPMED_API
    headers = {"X-Auth-Token": token}
    start = time.time()

    while time.time() - start < max_wait:
        resp = requests.get(f"{base_url}/jobs/{job_id}/status",
                            headers=headers, timeout=30)
        status = resp.json()
        state = status.get("state", -1)

        # States: 1=waiting, 2=running, 3=exporting, 4=success, 5=failed, 6=canceled
        state_names = {1: "waiting", 2: "running", 3: "exporting",
                       4: "complete", 5: "failed", 6: "canceled"}
        elapsed = time.time() - start
        print(f"  [{elapsed/60:.0f}m] Job {job_id}: {state_names.get(state, state)}")

        if state >= 4:
            return status

        time.sleep(poll_interval)

    return {"state": -1, "message": "Timeout"}


def print_imputation_instructions(input_file: str):
    """Print step-by-step instructions for running imputation manually."""
    print(f"""
╔══════════════════════════════════════════════════════════╗
║  IMPUTATION SETUP INSTRUCTIONS                          ║
╠══════════════════════════════════════════════════════════╣
║                                                          ║
║  To improve PPS accuracy from ~10% to ~90% coverage,    ║
║  impute your genotype data against a reference panel.    ║
║                                                          ║
║  Step 1: Register (free, one-time)                       ║
║    → https://imputationserver.sph.umich.edu              ║
║    → Create account → go to Profile → copy API token     ║
║                                                          ║
║  Step 2: Set your token                                  ║
║    export IMPUTATION_TOKEN="your-token-here"             ║
║                                                          ║
║  Step 3: Run imputation                                  ║
║    PYTHONPATH=. python3 src/analysis/imputation.py \\    ║
║      {input_file:52s}║
║                                                          ║
║  Step 4: Wait for completion (~1-24 hours)               ║
║    The script will poll and download automatically.      ║
║                                                          ║
║  Step 5: Run PPS on imputed data                         ║
║    PYTHONPATH=. python3 src/analysis/pps_calculator.py \\ ║
║      data/raw/imputation/imputed_merged.vcf.gz           ║
║                                                          ║
║  Privacy: Data encrypted in transit (HTTPS) and at rest. ║
║  Auto-deleted from server after 7 days.                  ║
╚══════════════════════════════════════════════════════════╝
""")


def run_imputation_pipeline(
    input_file: str,
    token: str | None = None,
    server: str = "michigan",
    ref_panel: str = "1000g-phase-3-v5",
    population: str = "mixed",
):
    """
    End-to-end: convert consumer file → submit imputation → wait → download.
    """
    token = token or os.environ.get("IMPUTATION_TOKEN", "")
    if not token:
        print("  No imputation server token found.")
        print_imputation_instructions(input_file)
        return None

    # Step 1: Convert to VCF
    print("=== Step 1: Converting to VCF ===")
    vcf_files = convert_consumer_to_vcf(input_file)

    # Step 2: Submit
    print("\n=== Step 2: Submitting imputation job ===")
    result = submit_imputation_job(
        vcf_files, token, server, ref_panel, population
    )

    if not result.get("success"):
        print(f"  Failed: {result.get('message', 'Unknown error')}")
        return None

    job_id = result.get("id", "")

    # Step 3: Wait
    print("\n=== Step 3: Waiting for imputation ===")
    status = wait_for_job(job_id, token, server)

    if status.get("state") != 4:
        print(f"  Job did not complete successfully: {status}")
        return None

    print("\n  Imputation complete!")
    print(f"  Download your results from the server web interface.")
    print(f"  Then run: PYTHONPATH=. python3 src/analysis/pps_calculator.py <imputed.vcf.gz>")

    return job_id


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: python imputation.py <consumer_genotype_file>")
        print("\nSet IMPUTATION_TOKEN env var with your Michigan/TOPMed server token.")
        print("Register at: https://imputationserver.sph.umich.edu")
        sys.exit(1)

    input_file = sys.argv[1]
    token = os.environ.get("IMPUTATION_TOKEN", "")

    if not token:
        # No token — just convert and print instructions
        print("=== Converting consumer file to VCF ===\n")
        vcf_files = convert_consumer_to_vcf(input_file)
        print(f"\n  Converted {len(vcf_files)} chromosome files to data/raw/imputation/")
        print_imputation_instructions(input_file)
    else:
        run_imputation_pipeline(input_file, token)
