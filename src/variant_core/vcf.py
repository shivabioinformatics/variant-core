from dataclasses import dataclass
from typing import Optional, Generator, Tuple

@dataclass(frozen=True, slots=True)
class Variant:
    """
    Immutable representation of a Genomic Variant.
    """ 
    chrom: str 
    pos: int 
    id: str 
    ref: str 
    alt: str 
    qual: Optional[float]
    filter: str 
    info: tuple
    format_fields: tuple  # Renamed from 'format_file' since it represents columns, not a file
    samples: tuple
    
    # I use __post_init__ to catch bad data immediately at creation time.
    # This aligns with the 'Fail Fast' principle to prevent silent errors downstream.
    def __post_init__(self): 
        # 1. Quality Check
        if self.qual is not None and self.qual < 0: 
            raise ValueError(f"Quality must not be negative, got {self.qual}")
        
        # 2. Position Check
        if self.pos < 1: 
            raise ValueError(f"Position must be >= 1, got {self.pos}")
        
        # 3. Base Character Check
        self._base_validation()
        
    def _base_validation(self):
        # I explicitly added the comma ',' here to handle multi-allelic variants (e.g. ALT="G,T")
        # without this, the parser would crash on standard VCFs.
        valid_bases = set("GTCAN.,")
        
        for base in (self.ref + self.alt).upper():
            if base not in valid_bases:
                raise ValueError(f"Invalid base '{base}' in ref/alt")
    
    def is_snp(self) -> bool:
        """Returns True if the variant is a Single Nucleotide Polymorphism."""
        return len(self.ref) == 1 and len(self.alt) == 1
    
    def __str__(self) -> str:
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"
        
    
class VCFReader:
    """
    I implemented this as a Generator to handle large files lazily.
    This prevents memory spikes when running on limited cloud instances.
    """
    def __init__(self, filepath: str): 
        self.filepath = filepath 
    
    def __iter__(self) -> Generator[Variant, None, None]: 
        with open(self.filepath, 'r') as f: 
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                yield self._parse_line(line)
        
    def _parse_line(self, line: str) -> Variant: 
        fields = line.strip().split('\t')
        if len(fields) < 8:
            raise ValueError("Insufficient columns")
        
        return Variant(
            chrom=fields[0],
            pos=int(fields[1]), 
            id=fields[2],
            ref=fields[3],
            alt=fields[4],
            qual=None if fields[5] == "." else float(fields[5]), 
            filter=fields[6],
            info=tuple(fields[7].split(";")),
            format_fields=tuple(fields[8].split(":")),
            samples=tuple(fields[9:])
        )
