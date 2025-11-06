from setuptools import setup, find_packages

setup(
    name='MGCalibrator',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas>=1.5',
        'numpy>=1.24',
        'matplotlib>=3.7',
        'seaborn>=0.12',
        'pysam>=0.20',
        'biopython',
    ],
    entry_points={
        'console_scripts': [
            'mgcalibrator=MGCalibrator.cli:main'
        ]
    },
    author='Nimrod de Wit',
    author_email='nimrod.de.wit@rivm.nl',
    description='CLI tool for computing absolute abundances in metagenomic BAM files',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/NimroddeWit/MGCalibrator',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
    ],
    python_requires='>=3.8',
)
