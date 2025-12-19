**VariantCore: High-Performance Genomic Data Structures**

VariantCore is a lightweight, memory-efficient library for parsing VCF and BED files. It is designed for clinical pipelines where data integrity and memory footprint are critical.

**Why this exists? (Engineering Philosophy)**

Most ad-hoc bioinformatics scripts lack type safety and consume excessive memory. I built this library to demonstrate how Domain-Driven Design can improve pipeline reliability.

**Usage**

I think you might have hit enter too soon!

If you meant **"Show me the whole file,"** here is the complete, polished `README.md` with the new installation instructions and usage examples combined.

You can copy-paste this entire block to overwrite your current `README.md`.

```markdown
# VariantCore: High-Performance Genomic Data Structures

**VariantCore** is a lightweight, strictly-typed library for parsing VCF and BED files. It is designed for clinical bioinformatics pipelines where **data integrity** and **memory efficiency** are critical.

## Why this exists (Engineering Philosophy)

Most ad-hoc bioinformatics scripts lack type safety and consume excessive memory. I built this library to demonstrate how **Domain-Driven Design** can improve pipeline reliability without the overhead of heavy dependencies.

## Installation

You can install `VariantCore` directly from GitHub using pip:

```bash
pip install git+[https://github.com/shivabioinformatics/variant-core.git](https://github.com/shivabioinformatics/variant-core.git)

```

## Usage

### Reading VCFs

```python
from variant_core import VCFReader

# Lazy loading with generators keeps memory usage low
reader = VCFReader("data/sample.vcf")

for variant in reader:
    if variant.is_snp():
        print(f"Found SNP: {variant.chrom}:{variant.pos}")

```

### Reading BED Files

```python
from variant_core import BEDReader

# Automatically handles whitespace and 0-based coordinates
bed = BEDReader("data/targets.bed")

for region in bed:
    print(f"Target Region: {region}")

```

## Testing

To run the test suite:

```bash
pip install -r requirements.txt
pytest

```
