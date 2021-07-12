import sys
from setuptools import setup
from setuptools_rust import RustExtension

setup(
    name="counts_to_csv",
    version="0.1.0",
    packages=["counts_to_csv"],
    install_requires=["setuptools_rust"],
    rust_extensions=[RustExtension("counts_to_csv.counts_to_csv")],
    include_package_data=True,
    zip_safe=False,
)