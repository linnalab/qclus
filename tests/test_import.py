import qclus as qc


def test_public_api():
    assert callable(qc.run_qclus)
    assert callable(qc.quickstart_qclus)
    assert callable(qc.utils.fraction_unspliced_from_bam)
    assert callable(qc.utils.fraction_unspliced_from_loom)


def test_version():
    assert qc.__version__ != "unknown"
