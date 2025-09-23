from setuptools import setup, find_packages

setup(
    name='MGCalibrator',
    version='0.1',
    packages=find_packages(),
    install_requires=['biopython', 'pandas', 'numpy', 'matplotlib', 'seaborn'],
    entry_points={
        'console_scripts': [
            'MGCalibrator=MGCalibrator.cli:main'
        ]
    },
)
