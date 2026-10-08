# Releasing QClus


Releases are made by hand from a maintainer's machine. Nothing is published automatically.


A release on PyPI cannot be undone, and a version number can never be used twice. If a release turns out to be wrong, the fix is a new version.



## What you need


- Push rights to this repository.
- The owner or maintainer role for the `qclus` project on PyPI, and an account on TestPyPI.
- An API token for each of the two sites in `~/.pypirc`, under `[pypi]` and `[testpypi]`. This file must never be committed.
- `uv` and the GitHub command line tool `gh`.



## Steps


Replace `X.Y.Z` with the new version throughout.


**1. Check that the tests pass.** The latest run on `main` must be green on GitHub. If you have the real sample, also run the regression test:

```bash
QCLUS_TEST_DATA=<data folder> QCLUS_TEST_SAMPLE=<sample name> python -m pytest tests/test_regression.py
```


**2. Set the version.** In a pull request, set `version = "X.Y.Z"` in `pyproject.toml` and give the changelog section its date. Merge it once the tests pass.


**3. Build from a clean copy of the merged commit.** Building from an export keeps stray local files out of the package.

```bash
git switch main && git pull
src=$(mktemp -d)
git archive HEAD | tar -x -C "$src"
uv build "$src" --out-dir dist
uvx twine check --strict dist/*
```


**4. Upload to TestPyPI.** TestPyPI is a practice copy of PyPI, used to test the files before the real upload.

```bash
uvx twine upload -r testpypi dist/*
```


**5. Install from TestPyPI and test.** Take only QClus from TestPyPI and its dependencies from PyPI: TestPyPI holds unrelated packages under the names of real ones.

```bash
env=$(mktemp -d)
uv venv "$env/venv"
py="$env/venv/bin/python"
uv pip install --python "$py" --no-deps --index-url https://test.pypi.org/simple/ "qclus==X.Y.Z"
"$py" -c "import importlib.metadata as m; print('\n'.join(r for r in m.requires('qclus') if 'extra ==' not in r))" > "$env/requirements.txt"
uv pip install --python "$py" -r "$env/requirements.txt" "pytest>=8"
cp -R tests "$env/tests"
(cd "$env" && venv/bin/qclus --version && venv/bin/python -m pytest tests)
```

The tests are copied out of the repository so that they run against the installed package and not the source folder.


**6. Upload to PyPI.** This uploads the same files that were just tested.

```bash
uvx twine upload -r pypi dist/*
```


**7. Tag the commit and create the GitHub release.** Use the changelog section as the release notes.

```bash
git tag -a vX.Y.Z -m "QClus X.Y.Z"
git push origin vX.Y.Z
gh release create vX.Y.Z --title "QClus X.Y.Z" --notes-file <notes file> dist/*
```


**8. Open the next version.** In a pull request, set the version in `pyproject.toml` to the next one with `.dev0` at the end, for example `X.Y.(Z+1).dev0`. An install from GitHub is then never mistaken for the release.
