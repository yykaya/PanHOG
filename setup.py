from setuptools import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="panhog",
    version="0.2.0",
    description="Phylogeny-Aware Pangenome Classification Toolkit",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/yykaya/PanHOG",
    author="Yasin Kaya",
    author_email="yyasinkkaya@gmail.com",

    py_modules=["PanHOG", "PangeneHOG"],

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
        "biopython",
        "pyyaml",
    ],

    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.8',
    license="MIT",
)
