from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("qclus")
except PackageNotFoundError:
    # Imported from a source checkout that is not installed.
    __version__ = "unknown"
