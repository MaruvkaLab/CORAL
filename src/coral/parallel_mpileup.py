"""The exact ``samtools mpileup`` invocation CORAL uses, in one place.

``coral.pileup_manager.Pileup.generate`` runs::

    samtools mpileup -f <ref> -B -d 100 <bam>...

The parallel extractor runs the same command with ``-r <chrom>`` added, and that
identity is what makes a chromosome's pileup its slice of the whole-genome pileup.
Keeping the options here -- in a module that imports nothing from CORAL -- means
the equivalence tests can check the command itself without the rest of the stack.
"""

from __future__ import annotations

#: Options that must match ``Pileup.generate`` exactly.
MPILEUP_OPTS = ["-B", "-d", "100"]


def mpileup_cmd(ref_fasta, bams, region=None):
    """Build the mpileup argv. ``region`` adds ``-r``; everything else is fixed."""
    cmd = ["samtools", "mpileup", "-f", ref_fasta] + MPILEUP_OPTS
    if region is not None:
        cmd += ["-r", region]
    return cmd + list(bams)
