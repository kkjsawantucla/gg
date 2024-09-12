from setuptools import setup, find_packages

setup(
    name="gg",
    version=0.1,
    description="GG is an open-source code for building graph based grand canonical basin hopping calculators",
    download_url="https://github.com/kkjsawantucla/gg",
    author="Kaustubh Sawant, Geng Sun",
    python_requires=">=3.7",
    packages=find_packages(include=["gg", "gg.*"]),
    license="LGPLv2.1+",
    install_requires=["numpy", "ase", "networkx", "pandas"],
    zip_safe=False,
)
