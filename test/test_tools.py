import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import mdtraj as mdt
import pytest
from crossflow.filehandling import FileHandler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from alphafix import tools

TEST_DIR = Path(__file__).resolve().parent


class DummyResponse:
    def __init__(self, status_code=200, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload or {}
        self.text = text

    def json(self):
        return self._payload


class DummySession:
    def __init__(self, payload=None):
        self.payload = payload or [{"pdbUrl": "https://example.com/model.pdb"}]

    def get(self, url):
        if url.startswith("https://alphafold.com/api/prediction/"):
            return DummyResponse(200, payload=self.payload)
        return DummyResponse(200, text="ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")


class DummySubprocessTask:
    def __init__(self, cmd):
        self.cmd = cmd

    def set_inputs(self, inputs):
        self.inputs = inputs

    def set_outputs(self, outputs):
        self.outputs = outputs

    def __call__(self, *args):
        if "alias" in self.cmd:
            return ""
        if "blastp" in self.cmd:
            out = FileHandler().create("x.csv")
            out.write_text("P00519,95.0\n")
            return out
        if "tleap" in self.cmd:
            prmtop = FileHandler().create("system.prmtop")
            prmtop.write_text("PRMTOP")
            inpcrd = FileHandler().create("system.inpcrd")
            inpcrd.write_text("INPCRD")
            return prmtop, inpcrd, "tleap ok"
        if "pdb4amber" in self.cmd:
            out = FileHandler().create("out.pdb")
            out.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
            renum = FileHandler().create("out_renum.txt")
            renum.write_text("1 1 1\n")
            return out, renum
        if "ambpdb" in self.cmd:
            out = FileHandler().create("x.pdb")
            out.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
            return out
        raise AssertionError(self.cmd)


@pytest.fixture
def sample_pdb():
    return TEST_DIR / "test1.pdb"


def test_aliased_and_check_available(monkeypatch):
    monkeypatch.setattr(tools, "SubprocessTask", DummySubprocessTask)
    assert tools._aliased("blastp") == ""

    monkeypatch.setattr(tools.shutil, "which", lambda cmd: None)
    monkeypatch.setattr(tools, "_aliased", lambda cmd: "")
    with pytest.raises(FileNotFoundError):
        tools._check_available("blastp")


def test_check_exists(tmp_path):
    path = tmp_path / "exists.txt"
    path.write_text("ok")
    tools._check_exists(path)
    with pytest.raises(FileNotFoundError):
        tools._check_exists(tmp_path / "missing.txt")


def test_pdbify_and_trajify(sample_pdb):
    traj = mdt.load_pdb(sample_pdb)
    pdbout = tools._pdbify(traj)
    content = pdbout.read_text()
    assert content.strip()
    assert "ATOM" in content or "HEADER" in content

    traj2 = tools._trajify(sample_pdb)
    assert traj2.topology.n_atoms > 0

    with pytest.raises(TypeError):
        tools._trajify(object())


def test_residue_and_atom_ids(sample_pdb):
    traj = mdt.load_pdb(sample_pdb)
    atom = traj.topology.atom(0)
    expected = f"{atom.residue.name}{atom.residue.resSeq}.{atom.residue.chain.chain_id}"
    assert tools._residue_id(atom.residue) == expected
    assert tools._atom_id(atom) == f"{expected}@{atom.name}"


def test_unique_chain_ids(sample_pdb):
    traj = mdt.load_pdb(sample_pdb)
    out = tools.unique_chain_ids(traj)
    assert out.topology.n_chains == traj.topology.n_chains
    assert [c.chain_id for c in out.topology.chains] == [c.chain_id for c in traj.topology.chains]


def test_bumps_and_hetify(sample_pdb):
    result = tools.bumps(sample_pdb)
    assert isinstance(result, str)

    het = tools.hetify(TEST_DIR / "test_nonprot.pdb")
    assert "HETATM" in het.read_text()


def test_smith_waterman_and_aln_score():
    aln1, aln2 = tools.smith_waterman("ACGT", "AC")
    assert len(aln1) == len(aln2)
    assert tools.aln_score(("ACGT", "ACGT")) == (4, 0, 0)


def test_sequences_and_match_align(sample_pdb):
    traj = mdt.load_pdb(sample_pdb)
    seqs = tools.sequences(traj.topology)
    assert isinstance(seqs, list) and len(seqs) > 0

    pdbout, alignment = tools.match_align(sample_pdb, sample_pdb)
    content = pdbout.read_text()
    assert content.strip()
    assert "ATOM" in content or "HEADER" in content
    assert len(alignment[0]) == len(alignment[2])


def test_sp_search(monkeypatch):
    monkeypatch.setattr(tools, "_check_available", lambda cmd: None)
    monkeypatch.setattr(tools, "SubprocessTask", DummySubprocessTask)
    matches = tools.sp_search("ACGT")
    assert matches[0]["uniprotAccession"] == "P00519"


def test_uni_find_frompdb_and_fromblastp(sample_pdb, monkeypatch):
    out = tools.uni_find_frompdb(sample_pdb)
    assert isinstance(out, dict)

    monkeypatch.setattr(tools, "sp_search", lambda seq: [{"uniprotAccession": "P00519", "percent_identity": 95.0}])
    assert tools.uni_find_fromblastp(sample_pdb, "A") == "P00519"


def test_uni_find(sample_pdb, monkeypatch):
    monkeypatch.setattr(tools, "uni_find_frompdb", lambda pdb_in: {"A": "P00519"})
    monkeypatch.setattr(tools, "uni_find_fromblastp", lambda pdb_in, chain_id: "P00519")
    result = tools.uni_find(sample_pdb)
    assert result["A"][0] == "P00519"


def test_alpha_get(monkeypatch):
    monkeypatch.setattr(tools.requests, "Session", lambda: DummySession())
    out = tools.alpha_get("P00519")
    assert "ATOM" in out.read_text()


def test_uniprot_diff_and_alpha_check(sample_pdb, monkeypatch):
    monkeypatch.setattr(tools, "alpha_get", lambda uniprot_id: sample_pdb)
    monkeypatch.setattr(tools, "match_align", lambda pdb_in, pdb_ref, cutoff=0.02, renumber=False, align=True: (pdb_ref, ("AAAA", "||||", "AAAA")))
    log = tools.uniprot_diff(sample_pdb, "P00519")
    assert isinstance(log, str)

    monkeypatch.setattr(tools, "uni_find", lambda pdb_in: {"A": ("P00519", "from user"), "B": ("P00519", "from user")})
    monkeypatch.setattr(tools, "uniprot_diff", lambda *args, **kwargs: "ok")
    log2 = tools.alpha_check(sample_pdb, ["P00519", "P00519"])
    assert "ALPHACHECK" in log2


def test_alpha_fix(monkeypatch, sample_pdb):
    monkeypatch.setattr(tools, "alpha_get", lambda uniprot_id: sample_pdb)
    monkeypatch.setattr(tools, "uniprot_diff", lambda *args, **kwargs: "ok")
    monkeypatch.setattr(tools, "rest_min", lambda pdbin, pdbref=None, kr=1.0, maxcyc=200: (pdbin, "ok"))

    def fake_match_align(pdb_in, pdb_ref, cutoff=0.02, renumber=False, align=True):
        return sample_pdb, ("AAAA", "||||", "AAAA")

    monkeypatch.setattr(tools, "match_align", fake_match_align)
    out, log = tools.alpha_fix(sample_pdb, ["P00519", "P00519"])
    assert "ALPHAFIX" in log
    assert out.read_text().startswith("HEADER") or "ATOM" in out.read_text()


def test_leap_and_ambpdb_and_renumbered_pdb(tmp_path, monkeypatch):
    monkeypatch.setattr(tools, "_check_available", lambda cmd: None)
    monkeypatch.setattr(tools, "SubprocessTask", DummySubprocessTask)

    pdb = tmp_path / "system.pdb"
    pdb.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
    prmtop, inpcrd, stdout = tools.leap(pdb, ["protein.ff14SB"])
    assert prmtop.read_text() == "PRMTOP"
    assert inpcrd.read_text() == "INPCRD"
    assert "tleap ok" in stdout

    out = tools.ambpdb(pdb, prmtop)
    assert "ATOM" in out.read_text()

    renum_map = tmp_path / "renum.txt"
    renum_map.write_text("1 1 1\n")
    renumbered = tools.renumbered_pdb(pdb, renum_map)
    assert renumbered.read_text().startswith("ATOM")


def test_rest_min_and_rest_min_omm(tmp_path, monkeypatch):
    monkeypatch.setattr(tools, "rest_min_omm", lambda pdbin, pdbref=None, kr=1.0, maxcyc=200: (FileHandler().create("out.pdb"), "ok"))
    monkeypatch.setattr(tools, "renumbered_pdb", lambda pdb_tmp, renumbering: pdb_tmp)
    monkeypatch.setattr(tools, "bumps", lambda pdbin, cutoff=0.2: "")

    pdbin = tmp_path / "in.pdb"
    pdbin.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
    out2, log2 = tools.rest_min(pdbin)
    assert out2.read_text() == "" or out2.read_text().startswith("ATOM")
    assert "ok" in log2

    out3, log3 = tools.rest_min_omm(pdbin)
    assert out3.read_text() == "" or out3.read_text().startswith("ATOM")
    assert log3 == "ok"
