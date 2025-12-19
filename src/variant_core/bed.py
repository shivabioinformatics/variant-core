from dataclasses import dataclass 
from typing import Generator

@dataclass(frozen=True, slots=True) 
class Region: 
    """
    Immutable representation of a genomic region.
    Note: BED format uses 0-based, half-open coordinates (Start is included, End is excluded).
    """
    chrom: str 
    start: int 
    end: int 
    name: str  # I keep this generic as a string to handle standard BED formats safely
    
    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end} ({self.name})"
    
class BEDReader: 
    """
    Memory-efficient BED Parser. 
    I implemented this to handle large files lazily, stripping headers and comments automatically.
    """
    def __init__(self, filepath: str) -> None: 
        self.filepath = filepath
    
    def __iter__(self) -> Generator[Region, None, None]: 
        # I use a generator here to avoid loading millions of regions into memory at once.
        # This is critical for efficient processing in limited cloud environments.
        with open(self.filepath, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines, comments (#), and UCSC track lines
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                yield self._parse_line(line)
    
    def _parse_line(self, line: str) -> Region: 
        fields = line.split('\t')
        
        # Basic validation: BED requires at least chrom, start, end
        if len(fields) < 3:
            raise ValueError(f"Invalid BED line: expected 3+ fields, got {len(fields)}")
        
        # I handle the optional 'name' field (Column 4) safely. 
        # If the file is BED3 (only coords), I default to "."
        name_field = fields[3] if len(fields) > 3 else "."
        
        return Region(
            chrom=fields[0],
            start=int(fields[1]),
            end=int(fields[2]),
            name=name_field
        )