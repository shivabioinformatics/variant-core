**VariantCore: High-Performance Genomic Data Structures**

VariantCore is a lightweight, memory-efficient library for parsing VCF and BED files. It is designed for clinical pipelines where data integrity and memory footprint are critical.

**Why this exists? (Engineering Philosophy)**

Most ad-hoc bioinformatics scripts lack type safety and consume excessive memory. I built this library to demonstrate how Domain-Driven Design can improve pipeline reliability.

**Usage**

## Installation

You can install `VariantCore` directly from GitHub using pip:

```bash
pip install git+[https://github.com/shivabioinformatics/variant-core.git](https://github.com/shivabioinformatics/variant-core.git)```


### Reading VCFs
```python
from variant_core import VCFReader

# Lazy loading with generators keeps memory usage low
reader = VCFReader("data/sample.vcf")

for variant in reader:
    if variant.is_snp():
        print(f"Found SNP: {variant.chrom}:{variant.pos}")



