from setuptools import setup, find_packages
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

from cgeisolate.version import __version__

setup(
    name='cgeisolate',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=__version__,
    packages=find_packages(),
    data_files=[],
    include_package_data=True,
    url='https://github.com/fredmol/cgeisolate',
    license='Apache-2.0',
    install_requires=[
        'weasyprint>=63.0',
        'Jinja2>=3.1.4',
        'matplotlib>=3.9.2',
    ],
    author='Malte B. Hallgren',
    scripts=['bin/cgeisolate'],
    author_email='malhal@food.dtu.dk',
    description='cgeisolate - K-mer Gene Typer',
    package_data={
        'cgeisolate': [
            'templates/*',
            'assets/*',
        ],
    },
)