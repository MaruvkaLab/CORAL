import pytest

from coral import alignment_manager as am
from coral import genome_manager as gm


def test_index_is_built_once_and_rebuilt_for_another_command(tmp_path, monkeypatch):
    genome = gm.Genome("ref", None, str(tmp_path), verbose=False)
    (tmp_path / "ref").mkdir()
    open(genome.fasta_path + ".fai", "w").close()
    calls = []
    monkeypatch.setattr(gm, "run_cmd", lambda cmd, verbose=True: calls.append(cmd[-1]))

    genome.index("bwa index {ref}")
    genome.index("bwa index {ref}")
    assert calls == [f"bwa index {genome.fasta_path}"]

    genome.index("minimap2 -x sr -d {ref}.mmi {ref}")
    assert calls[-1] == f"minimap2 -x sr -d {genome.fasta_path}.mmi {genome.fasta_path}"


def test_aligner_commands():
    assert am.aligner_commands("minimap2") == (am.ALIGNERS["minimap2"]["index"], am.ALIGNERS["minimap2"]["align"])
    assert am.aligner_commands("bwa", aligner_cmd="my_aligner {ref} {fq} {cores}") == (None, "my_aligner {ref} {fq} {cores}")
    assert am.aligner_commands("bwa", aligner_cmd="x {ref} {fq} {cores}", index_cmd="x-index {ref}")[0] == "x-index {ref}"
    with pytest.raises(ValueError):
        am.aligner_commands("unknown")
