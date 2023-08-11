from setuptools import find_packages, setup

setup(
    name="polya",
    version="0.0.1",
    packages=find_packages(exclude=["tests.*", "tests", "figs", "examples", "media"]),
    author="Killian Sheriff",
    author_email="ksheriff@mit.edu",
    description="A python implementation of Polya's enumeration theory.",
    long_description="",
    long_description_content_type="text/markdown",
    license="MIT",
    keywords=[],
    url="https://github.com/killiansheriff/polya",
    install_requires=["numpy", "sympy", "igraph"],
    include_package_data=True,
)
