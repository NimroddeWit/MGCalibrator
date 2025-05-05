from setuptools import setup, find_packages

setup(
    name='absolutifier',
    version='0.1',
    packages=find_packages(),
    install_requires=['biopython', 'pandas', 'numpy'],
    entry_points={
        'console_scripts': [
            'absolutifier=absolutifier.cli:main'
        ]
    },
)
