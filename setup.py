from setuptools import setup, find_packages

# Read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="panhog",
    version="0.1.0",  # Start with a version
    description="A tool for pangenome analysis using HOGs",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/yykaya/PanHOG",
    author="Yasin Kaya",
    author_email="yyasinkkaya@gmail.com",

  
    py_modules=["PanHOG"],

  
    entry_points={
        'console_scripts': [
            'panhog=PanHOG:main',
        ],
    },

    
    install_requires=[], 

    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS :: MacOS X",
    ],
    python_requires='>=3.6',
)
