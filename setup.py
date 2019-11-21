from setuptools import setup, find_packages
from nanovar.version import __version__
import os

current_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(current_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='nanovar',
    version=__version__,
    packages=find_packages(),
    package_data={'nanovar.model': ['ANN.E63B400L2N12-5D0.4-0.3SGDsee31_het_v1.h5'],
                  'nanovar.gaps': ['hg19_filter.bed', 'hg38_filter.bed', 'mm10_filter.bed'],
                  'nanovar.css': ['*.css'],
                  'nanovar.js': ['*.js']},
    include_package_data=True,
    scripts=['nanovar/nanovar'],
    url='https://github.com/cytham/nanovar',
    download_url='https://github.com/cytham/nanovar/releases',
    license='gpl-3.0',
    author='Tham Cheng Yong',
    author_email='chengyong.tham@u.nus.edu',
    description='NanoVar is a long-read structural variant caller.',
    keywords=['nanovar', 'structural variant caller', 'sv', 'nanopore', 'long read', 'low coverage', 'low depth'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    install_requires=['biopython==1.75', 'pybedtools==0.8.0', 'scipy==1.2.1', 'matplotlib==2.2.4', 'numpy==1.16.3',
                      'keras==2.2.4', 'tensorflow==1.13.1', 'natsort==6.2.0'],
    python_requires='>=3.7',
    classifiers=[
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)"
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
