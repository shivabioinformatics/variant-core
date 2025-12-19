**VariantCore: High-Performance Genomic Data Structures**

VariantCore is a lightweight, memory-efficient library for parsing VCF and BED files. It is designed for clinical pipelines where data integrity and memory footprint are critical.

**Why this exists? (Engineering Philosophy)**

Most ad-hoc bioinformatics scripts lack type safety and consume excessive memory. I built this library to demonstrate how Domain-Driven Design can improve pipeline reliability.


**Usage**

```
from variant_core.vcf import VCFReader
```
