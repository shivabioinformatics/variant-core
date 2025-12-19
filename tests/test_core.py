import pytest
from variant_core import Variant, BEDReader

def test_variant_raises_error_on_negative_quality():
    """
    Test that the __post_init__ validation logic catches invalid data.
    This ensures bad data from upstream pipelines doesn't silently corrupt our database.
    """
    with pytest.raises(ValueError, match="Quality must not be negative"):
        Variant(
            chrom="chr1", pos=100, id=".", ref="A", alt="T",
            qual=-10.0,  # <--- This invalid value should trigger the crash
            filter="PASS", info=(), format_fields=(), samples=()
        )

def test_variant_allows_multiallelic_bases():
    """
    Test that our validator correctly handles comma-separated ALT alleles (e.g., A -> G,T).
    This was a critical edge case I handled to support standard VCF specs.
    """
    v = Variant(
        chrom="chr1", pos=100, id=".", ref="A", 
        alt="G,T",  # <--- This includes a comma
        qual=50.0, filter="PASS", info=(), format_fields=(), samples=()
    )
    # If this crashes, the test fails. If it passes, our logic is sound.
    assert v.alt == "G,T"

def test_bed_parser_handles_basic_region(tmp_path):
    """
    Test that the BEDReader correctly parses a standard 4-column line.
    I use tmp_path (a pytest fixture) to create a real file on the fly, ensuring
    the test is self-contained and doesn't rely on external data files.
    """
    # 1. Create a dummy BED file
    d = tmp_path / "data"
    d.mkdir()
    p = d / "test.bed"
    p.write_text("chr1\t100\t200\tgeneA\n")

    # 2. Read it back using my parser
    reader = BEDReader(str(p))
    regions = list(reader)

    # 3. Assertions (The Proof)
    assert len(regions) == 1
    assert regions[0].chrom == "chr1"
    assert regions[0].start == 100
    assert regions[0].name == "geneA"