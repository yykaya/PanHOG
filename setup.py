from setuptools import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="panhog",
    version="0.3.0",
    description="Phylogeny-Aware Pangenome Classification Toolkit",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/yykaya/PanHOG",
    author="Yasin Kaya",
    author_email="yyasinkkaya@gmail.com",

    # These helper modules are imported by PanHOG at runtime and MUST ship with
    # the package, otherwise the installed `panhog` console script cannot run the
    # Ka/Ks (dN/dS), phylogeny-aware, gene-tree validation, or compartment Ka/Ks
    # analyses.
    py_modules=["PanHOG", "PangeneHOG", "panhog_dnds", "panhog_phylo",
                "panhog_genetrees", "panhog_kaks_compartments", "panhog_pansummary"],

    entry_points={
        'console_scripts': [
            'panhog=PanHOG:main',
            'pangenehog=PangeneHOG:main',
        ],
    },

    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn",
        "biopython>=1.80",
        "pyyaml",
    ],
    extras_require={
        # Enables the YN00 and ML dN/dS substitution models.
        "full": ["scipy"],
        # Everything needed to run the test suite.
        "dev": ["pytest", "scipy"],
    },

    project_urls={
        "Source": "https://github.com/yykaya/PanHOG",
        "Bug Tracker": "https://github.com/yykaya/PanHOG/issues",
    },
    keywords=[
        "bioinformatics", "pangenome", "phylogenetics", "orthology",
        "HOG", "OrthoFinder", "dN/dS", "KaKs", "comparative genomics",
    ],

    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.8',
    license="MIT",
)
