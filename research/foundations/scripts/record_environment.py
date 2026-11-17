"""Print the pinned verification environment in an archive-friendly form."""

from importlib.metadata import version
import platform
import sys

print("python_version=" + platform.python_version())
print("python_implementation=" + platform.python_implementation())
print("python_executable=" + sys.executable)
print("platform=" + platform.platform())
for package in ("sympy", "mpmath"):
    print(f"{package}_version={version(package)}")
