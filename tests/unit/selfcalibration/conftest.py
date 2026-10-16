import pytest

from .helpers import FakeImager


@pytest.fixture
def ms(tmp_path) -> str:
    path = tmp_path / "target.ms"
    path.mkdir()
    (path / "table.dat").write_text("visibilities")
    return str(path)


@pytest.fixture
def imager(tmp_path, casa) -> FakeImager:
    return FakeImager(output=str(tmp_path / "out_"), write_model=casa.write_model)
