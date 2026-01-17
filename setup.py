from setuptools import setup, find_packages

setup(
    name="gg",
    version="0.1.0",
    description="GG is an open-source code for building graph based modifiers and perform global optimization",
    download_url="https://github.com/kkjsawantucla/gg",
    author="Kaustubh Sawant, Geng Sun",
    python_requires=">=3.10",
    packages=find_packages(include=["gg", "gg.*"]),
    license="LGPLv2.1+",
    install_requires=["numpy", "ase", "networkx", "openai", "pandas", "PyYAML", "scipy", "rich"],
    scripts=["bin/add_mono", "bin/add_bi"],
    zip_safe=False,
)
