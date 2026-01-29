from setuptools import setup, find_packages

setup(
    name="neosv-trace",
    version="0.1.0",
    description="Structural-variant driven neoantigen prediction pipeline",
    author="Greyson Wintergerst",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "pyensembl",
        "biopython",
    ],
    entry_points={
        "console_scripts": [
            "neosv-trace=neosv_trace.main:main",
        ]
    },
    python_requires=">=3.8",
)
