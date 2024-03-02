from pathlib import Path

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import LatchAuthor, NextflowFileParameter, NextflowMetadata

NextflowMetadata(
    display_name="COMBINE-LAB/quantaf",
    author=LatchAuthor(),
    parameters={
        "permitlist": NextflowFileParameter(
            type=LatchFile,
            default=LatchFile("latch://1721.account/quantaf/pl_sheet.tsv"),
            path=Path("/root/pl_sheet.tsv"),
        ),
        "reference": NextflowFileParameter(
            type=LatchFile,
            default=LatchFile("latch://1721.account/quantaf/ref_sheet.tsv"),
            path=Path("/root/ref_sheet.tsv"),
        ),
        "sample": NextflowFileParameter(
            type=LatchFile,
            default=LatchFile("latch://1721.account/quantaf/sample_sheet.tsv"),
            path=Path("/root/sample_sheet.tsv"),
        ),
    },
)
